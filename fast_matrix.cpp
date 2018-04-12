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
#include <iostream>
using std::cout;
using std::endl;
#include <fstream>
using std::ifstream;
#include <cstring>
unsigned long long global_t;



int getMilliCount(){
    timeb tb;
    ftime(&tb);
    int nCount = tb.millitm + (tb.time & 0xfffff) * 1000;
    return nCount;
}


int getMilliSpan(int nTimeStart){
    int nSpan = getMilliCount() - nTimeStart;
    if(nSpan < 0)
        nSpan += 0x100000 * 1000;
    return nSpan;
}



int main(){
    ifstream file ("kmers.txt");
    string s;
    std::string delimiter = "\t";
    string delimiter_2 = ","  ;
    unsigned long long tot_kmers=0;
    unsigned long long unique_kmers=0;
    ifstream sub("subreads_info");
    string line;
    unsigned long long num_reads=0;
    while (std::getline(sub, line)) {
        num_reads = stol(line)/4;
    }

    vector<unsigned long long> temp;


    while (std::getline(file, s)){
        unique_kmers++;
        unsigned long long pos = 0;
        std::string index;
        std::string read;
        while ((pos = s.find(delimiter)) != std::string::npos) {
            index = s.substr(0, pos);
            s.erase(0, pos + delimiter.length());
            unsigned long long pos2=0;
            while ((pos2 = s.find(delimiter_2)) != std::string::npos) {
              tot_kmers++;
              read = s.substr(0, pos2);
              temp.push_back((stol(read)) + num_reads* (stol(index))); //cout<<(stol(read)-1) + num_reads* (stol(index))<<endl;
              s.erase(0, pos2 + delimiter_2.length());
              }
        temp.push_back((stol(s)) + num_reads* (stol(index)));

            tot_kmers++;

        }
      }

    std::cerr<<"Unique kmers:               "<<unique_kmers<<endl;
    std::cerr<<"temp:               "<<temp.size()<<endl;
    std::sort(temp.begin(),temp.end());
    std::cerr << "Total kmers:              " << tot_kmers << " should be equal to temp size"<<std::endl;

    cerr<<"constructing b"<<endl;
   unsigned long long n = unique_kmers * num_reads;
    unsigned long long m = temp.size();
    std::cerr << "stack allocating sdsl::sd_vector_builder base object with n=" << n<< " m=" << m << std::endl;
    sdsl::sd_vector_builder b_builder(n, m);

   for (std::vector<unsigned long long>::iterator it = temp.begin(); it!= temp.end(); ++it){

         b_builder.set(*it);}
    std::cerr << "sd_vector items:          " << b_builder.items()  <<" should be equal to sd_vector capacity"<< endl;
    cerr<< "sd_vector capacity:             " << b_builder.capacity() << std::endl;
    std::cerr << "Total wall time:          " << (getMilliCount() - global_t) / 1000.0 << " s" << std::endl;
    sdsl::sd_vector<> b(b_builder);
    cerr << "Total size:                    "<<size_in_mega_bytes(b) <<" Mb"<<endl;

    std::string outfilename = "output_matrix";
    sdsl::store_to_file(b, outfilename);
    std::cerr << "DONE!" << std::endl;
}
