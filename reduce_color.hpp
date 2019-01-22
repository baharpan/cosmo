#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <libgen.h> // basename
#include <boost/algorithm/string.hpp>

#include "sdsl/sd_vector.hpp"
#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>
#include "io.hpp"
#include "debruijn_graph.hpp"
#include "algorithm.hpp"
#include <sys/timeb.h>
#include <future>
using namespace std;
using namespace sdsl;

#include <sys/timeb.h>

#include <utility>
#include <ctime>

// TCLAP
#include "tclap/CmdLine.h"

#include <cstdio>
#include <cstdlib>
#include "kmer.hpp"
class reduction {

public:

  vector<string> load_input(){
    vector<string> files;
    ifstream infile("list.txt");
    string line;
    while (std::getline(infile, line)) {
      files.push_back(line);
    }
  return (files);
}

  std::map<string,size_t> index_maker(){
    cerr<<"Making the pair map"<<endl;
    ifstream f ("pair" , std::ifstream::in);
     std::map<string,size_t> pairs;
    string s;
    size_t index=0;
    while(std::getline(f,s)){
      std::locale loc;
      string ss="";
      for (std::string::size_type i=0; i<s.length(); ++i)
        ss += std::toupper(s[i],loc);
        pairs.insert(std::pair<string,size_t>(ss,index++));
      }
      std::map<string,size_t>::iterator it = pairs.begin();
      return pairs;
    }


    inline bool parser(std::ifstream f, std::vector<string>& fastq)
    {
    std::string s;
    if (!std::getline(f, s) || !(s[0] == '>'||s[0] == '@')) {
           std::cerr << "ERROR: File doesn't seem like a FASTQ file." << std::endl;
            exit(EXIT_FAILURE);
        }
    bool mode=true;
    s="";
    fastq.push_back(s);
    while(std::getline(f, s) && s.size() > 0){
      if (s[0] == '>' || s[0] == '@') {
        mode = true;
        s = "";
        fastq.push_back(s);
        }

      else if (s[0] == '+') {
        mode = false;
        }

      else if (mode == true) {
        *(fastq.end()-1) += s;
        }
    	}
    	return fastq.size();
    }

    bool read_okay(const std::string& readstr)
    {
        for (const char& c : readstr) {
            if (c == 'A' || c == 'G' || c == 'C' || c == 'T') {
                continue;
            } else {
                return false;
            }
        }

        return true;
    }

    void  okay_fastq(vector<string> fastq , vector<string>& filter_fastq){
    	for (std::vector<std::string >::iterator rf = fastq.begin(); rf != fastq.end(); ++rf) {
    		if (read_okay(*rf))
    			filter_fastq.push_back(*rf);
    	}
    }

    void kmer_counter(size_t k, std::vector<std::string>& fastq, std::vector<std::vector<std::string>>& found_kmers, std::map<string,size_t>& pairs)
    {
    	string kmer;
    	for (std::vector<string>::iterator it = fastq.begin(); it != fastq.end(); ++it){
    		size_t pos = 0;
    		vector <string> newset;
    		found_kmers.push_back(newset);
    		while(pos != it->length()-k+1){
    			kmer = it->substr(pos,k);
                map<string,size_t>::iterator find = pairs.find(kmer);
      		if (find != pairs.end()) {
    			     found_kmers.back().push_back(kmer);
             }
    		  pos++;
    		}
    	}

    }

    bool is_compatible( std::vector<string> &subreadit , size_t color, std::map<size_t,std::vector<size_t>> &sort , std::map<string,size_t>& pairs){
      if (sort.empty()) return false;
      for (size_t i = 0; i < subreadit.size(); i++){
        if (std::find(sort[color].begin(), sort[color].end(), pairs[subreadit[i]] + color * pairs.size()) != sort[color].end()){
          return false;
        }
    }
    return true;
    }

    void build_recolored_matrix (){
      vector<string> files = load_input();
      std::map<string,size_t> pairs = index_maker();
      vector<string> fast[files.size()];
      cerr<<"Making the found kmers and fastq vectors"<<endl;
      vector<vector<string>> found_kmer[files.size()];
      vector<string> filter_fastq[files.size()];
     // size_t found_kmer_size [files.size()];
      vector<vector<string>> found_kmers;
      vector<string> fastq;
      for(size_t i = 0; i < files.size(); i++){
        cerr<<"currently in file " << i << endl;
        parser(ifstream(files[i]), fast[i]);
        kmer_counter (32, fast[i], found_kmer[i],pairs);
       // found_kmer_size[i] = found_kmer[i].size();
        found_kmers.insert( found_kmers.end(), found_kmer[i].begin(), found_kmer[i].end() );
        fastq.insert( fastq.end(), fast[i].begin(), fast[i].end() );
      }
      std::ofstream subreads_info;
      subreads_info.open("subreads_info");
      unsigned long long number_of_subreads = found_kmers.size();
      subreads_info<<number_of_subreads;
      cerr<<"Number of subreads 		" <<found_kmers.size()<<endl;

      cerr<<"Non filtered fastq file size: 	"<<fastq.size()<<endl;

      int i = 0;
      std::cerr << "=== Determining size of color matrix===" << std::endl;
      unsigned long long tot_kmers = 0;
      for (std::vector<std::string >::iterator rf = fastq.begin(); rf != fastq.end(); ++i, ++rf) {
        tot_kmers += rf->size() - 32+1;
      }

      std::cerr<<"Unique kmers: 		"<<pairs.size()<<endl;
      std::cerr << "Total kmers: 		" << tot_kmers << std::endl;
      cerr<<"constructing b"<<endl;
      map<size_t,std::vector<size_t>> sort;
      vector<size_t> sort_backup;
      unsigned subread = 0;
      size_t count = 0;
      vector<vector<string> >::iterator begin = found_kmers.begin();
      size_t num_color = 0;
      std::map<size_t,vector<size_t>> colors_map;
      size_t index = 0;
      for(std::vector<vector<string>>::iterator subreadit = found_kmers.begin(); subreadit != found_kmers.end(); ++subreadit, ++index){
          count++;
          if (subread % 10000 == 0) cout<<"on read: "<<subread<<endl;
          size_t color = 0;
          bool need_color = true;
          while ( color != num_color + 1){
            if (is_compatible(*subreadit , color, sort, pairs)){
              for (vector<string>::iterator it = subreadit->begin(); it!= subreadit->end(); ++it){
                sort[color].push_back(pairs[*it] + color * pairs.size());
                sort_backup.push_back(pairs[*it] + color * pairs.size());
                }
              colors_map[color].push_back(index);
              need_color = false;
              break;
            }
            else color++;
          }

          if (need_color){
            num_color++;
            std::vector<size_t> new_column;
            for (vector<string>::iterator it = subreadit->begin(); it!= subreadit->end(); ++it){

            new_column.push_back(pairs[*it] + num_color * pairs.size());
            sort_backup.push_back(pairs[*it] + num_color * pairs.size());}
            sort.insert(make_pair(num_color, new_column));
            std::vector<size_t> new_color_vector;
            new_color_vector.push_back(index);
            colors_map.insert(make_pair(num_color, new_color_vector));
              }
      subread++;}

      cerr<<"Size of color map is: "<<colors_map.size()<<" Should be equal to new number of colors: "<<num_color + 1 <<" and should be less than old number of colors which was: "<<number_of_subreads<<endl;
      cerr<<"Size of the sortd set:		"<<sort_backup.size()<<" should be equal to value of m"<<endl;
      size_t n = (num_color +1) * pairs.size();
      size_t m = sort_backup.size();
      std::cerr << "stack allocating sdsl::sd_vector_builder base object with n=" << n<< " m=" << m << std::endl;
      sdsl::sd_vector_builder b_builder(n, m);
      for(std::vector<size_t>::iterator it = sort_backup.begin(); it!= sort_backup.end(); ++it){
      	b_builder.set(*it);
      }


      std::cerr << "sd_vector items: 		" << b_builder.items()  <<" should be equal to sd_vector capacity"<< endl;
      cerr<< "sd_vector capacity: 		" << b_builder.capacity() << std::endl;
      sdsl::sd_vector<> b(b_builder);
      cerr << "Total size: 			"<<size_in_mega_bytes(b) <<" Mb"<<endl;
      std::string outfilename = "output_matrix_reduced";
      sdsl::store_to_file(b, outfilename);
      std::cerr << "DONE!" << std::endl;
      }

};
