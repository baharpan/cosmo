#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <libgen.h>
#include <boost/algorithm/string.hpp>
#include "sdsl/sd_vector.hpp"
#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <bits/stdc++.h>
#include <cstdio>
#include <cstdlib>
using namespace std;
using namespace sdsl;

class reduction {

public:
  int num_color = 0;
  map<size_t, vector<size_t>> color_map;
  vector<size_t> build_backup;



  map<string,size_t> index_maker(){
    cerr<<"Making the pair map"<<endl;
    cerr<<"================================================="<<endl;
    ifstream f ("pair" , ifstream::in);
    map<string,size_t> pairs;
    string s;
    size_t index = 0;
    while( getline(f,s) ){
      locale loc;
      string ss = "";
      for (string::size_type i = 0; i < s.length(); ++i)
        ss += toupper(s[i],loc);
        pairs.insert(pair<string,size_t>(ss,index));
        index++;
      }
      f.close();
      return pairs;
    }


    map<string,size_t> pairs = index_maker();
    vector<int> online_kmers{vector<int>(pairs.size(),-1)};
    map<int,vector<int>> reads_of_kmer;



    vector<string> kmer_counter_read(int k, string read)
    {
      string kmer;
      vector<string> found_kmers_per_read;
      size_t pos = 0;
      if (k > read.size()) cout<<"ERROR: k should be smaller than read length"<<endl;
      while(pos != read.size()-k+1){
        kmer = read.substr(pos,k);
        map<string,size_t>::iterator find = pairs.find(kmer);
        if (find != pairs.end())
          found_kmers_per_read.push_back(kmer);

          pos++;
        }
      return found_kmers_per_read;
    }


    int which_color ( vector<string> &subreadit){
      int max = -1;
      for (size_t i = 0; i < subreadit.size(); i++){
        if (online_kmers[pairs[subreadit[i]]] > max ){
          max = online_kmers[pairs[subreadit[i]]];}
          online_kmers[pairs[subreadit[i]]] += 1;
    }

    while(max != color_map.size() + 1 ){
      for (size_t i = 0; i < subreadit.size(); i++){
		    for (size_t j = 0; j < color_map[max+1].size(); j++){
             if (find (reads_of_kmer[pairs[subreadit[i]]].begin(), reads_of_kmer[pairs[subreadit[i]]].end(),
            color_map[max+1][j]) != reads_of_kmer[pairs[subreadit[i]]].end()) { max++; i = -1; j = -1;  break;}
          }

        }
        return max+1;
      }


    }
    void test(map<size_t, vector<size_t>> color_map, map<int,vector<int>> reads_of_kmer){
      for (size_t i = 0; i < color_map.size(); i++){
        for (size_t k = 0; k < color_map[i].size()-1; k++){
          for (size_t j = 0; j < reads_of_kmer.size(); j++){
            if (find(reads_of_kmer[j].begin(), reads_of_kmer[j].end(), color_map[i][k]) != reads_of_kmer[j].end()
              && find(reads_of_kmer[j].begin(), reads_of_kmer[j].end(), color_map[i][k+1]) != reads_of_kmer[j].end()){
                cerr<<"ERROR: "<<color_map[i][k]<<" and "<< color_map[i][k+1]<<" should not have same label since kmer "<<j<<" is present in both"<<endl;
                exit(0);}

        }
      }
    }
    cerr<<"Test Passed"<<endl;
  }

    void build (vector<string> found_kmers_per_read, size_t read){
      int color = which_color(found_kmers_per_read);
      if (color < num_color ){
        for (vector<string>::iterator it = found_kmers_per_read.begin(); it!= found_kmers_per_read.end(); ++it){
          build_backup.push_back(pairs[*it] + color * pairs.size());
	        reads_of_kmer[pairs[*it]].push_back(read);
	}
          color_map[color].push_back(read);
            }
      else{
        num_color++;

        for (vector<string>::iterator it = found_kmers_per_read.begin(); it!= found_kmers_per_read.end(); ++it){
          build_backup.push_back(pairs[*it] + color * pairs.size());
          vector<int> new_reads;
          new_reads.push_back(read);
	        reads_of_kmer[pairs[*it]].push_back(read);}
        vector<size_t> new_colomn;
        new_colomn.push_back(read);
        color_map[color].push_back(read);
          }
      }


    void parser(int k, ifstream& f , int num_reads)
    {
      size_t count = 0;
      string s;

      if (!getline(f, s) || !(s[0] == '>'||s[0] == '@')) {
          cerr << "Make sure that the input is fastq file" << endl;
          exit(EXIT_FAILURE);
        }
    bool mode = true;
    s = "";
    while(getline(f, s) && s.size() > 0){
      string read;
      if (s[0] == '>' || s[0] == '@') {
        mode = true;
        s = "";
        }

      else if (s[0] == '+') {
        mode = false;
        }

      else if (mode == true) {
        if (count % 10000 == 0) cout<<"Processing read number "<< count <<
        " ("<<(double)count * 100/ num_reads<<" % is processed)\n";
        read += s;
        vector<string> found_kmers_per_read = kmer_counter_read(k, read);
        build (found_kmers_per_read , count);
        count++;
        }
    	}
    }


    void build_recolored_matrix (int k, ifstream& f , int num_reads){
      parser(k, f, num_reads);
      size_t n = (num_color ) * pairs.size();
      size_t m = build_backup.size();
      cerr<<"================================================="<<endl;
      cerr << "Matrix dimensions: n = " << n<< " m = " << m <<endl;
      sdsl::sd_vector_builder b_builder(n, m);
      std::sort(build_backup.begin(), build_backup.end());
      for(vector<size_t>::iterator it = build_backup.begin(); it!= build_backup.end(); ++it){
      	b_builder.set(*it);
      }

      sdsl::sd_vector<> b(b_builder);
      cerr << "Total size of reduced matrix: "<<size_in_mega_bytes(b) <<" Mb"<<endl;
      cerr<<"Writing the matrix to output_matrix_reduced"<<endl;
      string outfilename = "output_matrix_reduced";
      sdsl::store_to_file(b, outfilename);
      ofstream labels;
      labels.open("labels.txt");
      labels<<"Labels\tColors\n";
      for (map<size_t,vector<size_t>>::iterator it = color_map.begin(); it != color_map.end(); ++it){
        labels<<it->first<<"\t";
        for (size_t j = 0; j < it->second.size(); ++j)
          labels<<it->second[j]<<" ";
        labels<<endl;
      }
      labels.close();
      ofstream frequency;
      frequency.open("frequency.txt");
      frequency<<"kmer\treads\tfrequency\n";
      for (auto i = reads_of_kmer.begin(); i != reads_of_kmer.end(); i++){
        frequency<<i->first<<"\t";
        for (size_t j = 0; j < i->second.size(); ++j)
        frequency<<i->second[j]<<" ";
        frequency<<" : "<<i->second.size()<<endl;
        }
      //test(color_map, reads_of_kmer);
      }

};
