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
  map<unsigned long long, vector<unsigned long long>> color_map;
  map<unsigned long long, set<unsigned long long>> kmer_map;
  set<unsigned long long> build_backup;



  map<string,unsigned long long> index_maker(){
    cerr<<"Making the pair map"<<endl;
    cerr<<"================================================="<<endl;
    ifstream f ("pair" , ifstream::in);
    map<string,unsigned long long> pairs;
    string s;
    unsigned long long index = 0;
    while( getline(f,s) ){
      locale loc;
      string ss = "";
      for (size_t i = 0; i < s.length(); ++i)
        ss += toupper(s[i],loc);
        pairs.insert(pair<string,unsigned long long>(ss,index));
        index++;
      }
      f.close();
      return pairs;
    }


    map<string,unsigned long long> pairs = index_maker();
    vector<int> online_kmers{vector<int>(pairs.size(),-1)};


    vector<string> kmer_counter_read(int k, string read)
    {
      string kmer;
      vector<string> found_kmers_per_read;
      size_t pos = 0;
      if (k > read.size()) cout<<"ERROR: k should be smaller than read length"<<endl;
      while(pos != read.size()-k+1){
        kmer = read.substr(pos,k);
        map<string,unsigned long long>::iterator find = pairs.find(kmer);
        if (find != pairs.end())
          found_kmers_per_read.push_back(kmer);

          pos++;
        }
      return found_kmers_per_read;
    }


    int which_color ( vector<string> &subreadit){
      int max = -1;
      vector<size_t> counter;
      for (size_t i = 0; i < subreadit.size(); i++){
        if (find(counter.begin(), counter.end(), pairs[subreadit[i]]) == counter.end()){
          if (online_kmers[pairs[subreadit[i]]] > max ){
          max = online_kmers[pairs[subreadit[i]]];}

          online_kmers[pairs[subreadit[i]]] += 1;
          counter.push_back(pairs[subreadit[i]]);}
        }


        if (max == -1) return max+1;
        while(max != color_map.size() + 1 ){
          for (size_t i = 0; i < subreadit.size(); i++){
            auto it = kmer_map[max+1].find(pairs[subreadit[i]]);
            if (it != kmer_map[max+1].end())
                    { max++; i = -1;  }


                }
        return max+1;
      }
    }



    void build (vector<string> found_kmers_per_read, size_t read){
      int color = which_color(found_kmers_per_read);
      if (color == num_color ) num_color++;
        for (vector<string>::iterator it = found_kmers_per_read.begin(); it!= found_kmers_per_read.end(); ++it){
          build_backup.insert(pairs[*it] + color * pairs.size());
          kmer_map[color].insert(pairs[*it]);
	}

          color_map[color].push_back(read);
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
      cout<<"Finish Processing"<<endl;
    }


    void build_recolored_matrix (int k, ifstream& f , int num_reads){

      parser(k, f, num_reads);
      unsigned long long n = (num_color ) * pairs.size();
      unsigned long long m = build_backup.size();


      cerr<<"================================================="<<endl;
      cerr << "Matrix dimensions: n = " << n<< " m = " << m <<endl;
      sdsl::sd_vector_builder b_builder(n, m);

      //std::sort(build_backup.begin(), build_backup.end());
      for(set<unsigned long long>::iterator it = build_backup.begin(); it!= build_backup.end(); ++it){
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
      for (map<unsigned long long,vector<unsigned long long>>::iterator it = color_map.begin(); it != color_map.end(); ++it){
        labels<<it->first<<"\t";
        for (size_t j = 0; j < it->second.size(); ++j)
          labels<<it->second[j]<<" ";
        labels<<endl;
      }
      labels.close();

      }

};
