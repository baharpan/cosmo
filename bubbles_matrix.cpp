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
#include "tclap/CmdLine.h"
#include <cstdio>
#include <cstdlib>
using namespace std;
using namespace sdsl;


struct parameters {
  string i;
  size_t k = 0;
  size_t n = 0;
};


parameters parse_arguments(int argc, char **argv) {
  parameters params;
  TCLAP::CmdLine cmd(" ");
  TCLAP::ValueArg<string> input("i", "file", "input file", true, " ", "input file", cmd);
  TCLAP::ValueArg<size_t> k_value("k", "kmer_length", "Length of edges (node is k-1)", true, 32, "k value", cmd);
  TCLAP::ValueArg<size_t> num_reads("n", "num_reads", "number of the reads", true, 0, "number of reads", cmd);

  cmd.parse( argc, argv );
  params.k = k_value.getValue();
  params.n = num_reads.getValue();
  params.i = input.getValue();
  return params;
}

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


    vector<vector<string>> kmer_counter_read(size_t k, string read)
    {
      string kmer;
      vector<vector<string>> found_kmers_per_read;
      size_t pos = 0;
      if (k > read.size()) {
        cout<<"ERROR: k should be smaller than minimum read length"<<endl;
        exit(EXIT_FAILURE);
      }
      size_t subread = 0;
      vector<string> new_subread;
      found_kmers_per_read.push_back(new_subread);
      while(pos != read.size()-k+1){
        kmer = read.substr(pos,k);
        map<string,unsigned long long>::iterator find = pairs.find(kmer);
        if (find != pairs.end()){
          if (std::find(found_kmers_per_read[subread].begin(), found_kmers_per_read[subread].end(), kmer) == found_kmers_per_read[subread].end() ){
            found_kmers_per_read[subread].push_back(kmer);


          }
          else{
            found_kmers_per_read.push_back(new_subread);
            subread++;
            found_kmers_per_read[subread].push_back(kmer);
            //cerr<<kmer<<" "<<subread<<endl;


          }
        }
          pos++;
        }
      return found_kmers_per_read;
    }






    void build (vector<string> found_kmers_per_subread, size_t read , size_t num_reads){
        for (vector<string>::iterator it = found_kmers_per_subread.begin(); it!= found_kmers_per_subread.end(); ++it){
          //build_backup.insert(pairs[*it] + read * pairs.size());
          //cout<<*it<<" "<<read<<endl;
          build_backup.insert(read  + num_reads * pairs[*it]);


	}


      }


    size_t parser(size_t k, ifstream& f , size_t num_reads, ofstream& repeats, bool isBuild)
    {
      size_t count = 0;
      size_t numberOfSubreads = 0;
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
        if (count % 10000 == 0 && isBuild) cout<<"Processing read number "<< count <<
        " ("<<(double)count * 100/ num_reads<<" % is processed)\n";
        read += s;

        vector<vector<string>> found_kmers_per_read = kmer_counter_read(k, read);
        if (isBuild){
        if (found_kmers_per_read.size() > 1) repeats<<found_kmers_per_read[1].front().substr(0,k-1)<<" "<<numberOfSubreads+found_kmers_per_read.size() - 1<<endl;}
        for (size_t subRead = 0; subRead < found_kmers_per_read.size(); subRead++){
          if (isBuild){
          build (found_kmers_per_read[subRead] , numberOfSubreads, num_reads);}
          numberOfSubreads++;
        }

        count++;
        }
    	}
      cout<<"Finish Processing"<<endl;
      return (numberOfSubreads);

    }


    void build_recolored_matrix (size_t k, ifstream& f , size_t num_reads){
      ofstream repeats;
      repeats.open("repeats");
      size_t numberOfSubreads = parser(k, f, num_reads, repeats, false);
      f.clear() ;
      f.seekg(0);
      numberOfSubreads = parser(k, f, numberOfSubreads, repeats , true);
      unsigned long long n = (numberOfSubreads ) * pairs.size();
      unsigned long long m = build_backup.size();


      cerr<<"================================================="<<endl;
      cerr << "Matrix dimensions: n = " << n<< " m = " << m <<endl;

      sdsl::sd_vector_builder b_builder(n, m);

      //std::sort(build_backup.begin(), build_backup.end());
      for(set<unsigned long long>::iterator it = build_backup.begin(); it!= build_backup.end(); ++it){
      	b_builder.set(*it);
      }
      sdsl::sd_vector<> b(b_builder);
      cerr << "Total size of matrix: "<<size_in_mega_bytes(b) <<" Mb"<<endl;
      cerr<<"Writing the matrix to output_matrix"<<endl;
      string outfilename = "output_matrix";
      sdsl::store_to_file(b, outfilename);
      cerr<<"Number of subreads: "<<numberOfSubreads<<endl;
      ofstream of;
      of.open("subreadInfo");
      of<<numberOfSubreads<<endl;
    }

      int main(int argc, char* argv[]){
        auto params = parse_arguments(argc, argv);
        size_t k = params.k;
        size_t num_reads = params.n;
        ifstream fastqfile(params.i);
        build_recolored_matrix (k, fastqfile , num_reads);
        fastqfile.close();
        cerr<<"Number of reads:  "<<num_reads <<endl;
        cerr<<"Number of unique kmers:  "<<pairs.size()<<endl;


        cerr<<"DONE!"<<endl;

      }
