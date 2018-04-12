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
//
unsigned long long global_t;
bool trace = false;
string file_extension = ".dbg";

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


void index_maker(ifstream f , std::map<string,size_t>& pairs){

string s;
size_t index=1;
while(std::getline(f,s)){
  std::locale loc;
  string ss="";
  for (std::string::size_type i=0; i<s.length(); ++i)
    ss+=std::toupper(s[i],loc);
    pairs.insert(std::pair<string,size_t>(ss,index++));
  }
std::map<string,size_t>::iterator it = pairs.begin();

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

void kmer_counter(size_t k,std::vector<std::string>& fastq,std::vector<std::vector<std::string>>& found_kmers, std::map<string,size_t> pairs)
{
	std::string kmer;
	for (std::vector<string>::iterator it = fastq.begin(); it != fastq.end(); ++it){
		size_t pos=0;
		vector <string> newset;
		found_kmers.push_back(newset);

    while(pos!= it->length()-k+1){
		  	kmer=it->substr(pos,k);
        map<string,size_t>::iterator findd = pairs.find(kmer);

  		  if (findd!=pairs.end()) {

			      if(std::find((found_kmers.back()).begin(),(found_kmers.back()).end(),kmer)!=(found_kmers.back()).end()){

				       vector <string> newset;
				       found_kmers.push_back(newset);
			        }

			     found_kmers.back().push_back(kmer);
        }

    pos++;
		}
	}

}



int main(){

vector<string> files;
ifstream infile("list.txt");
string line;
while (std::getline(infile, line)) {
    files.push_back(line);
  }



cerr<<"===Constructing the color matrix==="<<endl;
std::map<string,size_t> pairs;
cerr<<"Making the pair map"<<endl;
index_maker(ifstream ("pair") ,  pairs);

vector<string> fast[files.size()];
//cerr<<"the number of the reads is: "<<fastq.size()<<endl;
cerr<<"Making the found kmers and fastq vectors"<<endl;
vector<vector<string>> found_kmer[files.size()];
vector<string> filter_fastq[files.size()];
//size_t found_kmer_size [files.size()];
//okay_fastq(fastq ,filter_fastq);
//cerr<<"Filtered fastq file size: 	"<<filter_fastq.size()<<endl;
//cerr<<"Non filtered fastq file size: 	"<<fastq.size()<<endl;
vector<vector<string>> found_kmers;
vector<string> fastq;
cerr<<"Making final found kmers"<<endl;
for(size_t i=0; i< files.size(); i++){
cerr<<"in file "<<i<<endl;
parser(ifstream(files[i]), fast[i]);

kmer_counter (32, fast[i], found_kmer[i],pairs);
//found_kmer_size[i] = found_kmer[i].size();
//kmer_counter (32, fast[1], found_kmer[1]);
//found_kmers.reserve( found_kmer[0].size() + found_kmer[1].size() );
found_kmers.insert( found_kmers.end(), found_kmer[i].begin(), found_kmer[i].end() );
//found_kmers.insert( found_kmers.end(), found_kmer[1].begin(), found_kmer[1].end() );


//fastq.reserve( fast[0].size() + fast[1].size() ); // preallocate memory
fastq.insert( fastq.end(), fast[i].begin(), fast[i].end() );
//fastq.insert( fastq.end(), fast[1].begin(), fast[1].end() );
}
std::ofstream subreads_info;
subreads_info.open("subreads_info");
size_t number_of_subreads=found_kmers.size();
subreads_info<<number_of_subreads;
cerr<<"Number of subreads 		" <<found_kmers.size()<<endl;
map<int,int> check_point;
//cerr<<"the number of subreads is: "<<found_kmers.size()<<endl;
//okay_fastq(fastq ,filter_fastq);
//cerr<<"Filtered fastq file size: 	"<<filter_fastq.size()<<endl;
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
size_t n = number_of_subreads * pairs.size();
std::vector<size_t> sort;
unsigned subread=0;
size_t count=0;
vector<vector<string> >::iterator begin = found_kmers.begin();

//size_t pointer = found_kmer_size[0];
cerr<<"in while loop"<<endl;


while (!found_kmers.empty()){
    count++;
    vector<vector<string> >::iterator  subreadit = found_kmers.begin();
		/* for(int i=0; i< files.size(); i++){  //Multiple input samples
		if(std::find(found_kmer[i].begin(), found_kmer[i].end(), *subreadit) != found_kmer[i].end()) {

     while (j != files.size()){
      if (count / pointer ==0){
            //cerr<<"here"<<endl;
			     check_point[count]=j;
          // cerr<<"j "<<j<<endl;
           break;
         }
      else{ j++;
            pointer+=found_kmer_size[j];
          }
      }

}
}*/
		if (subread % 10000 == 0) cout<<"on read: "<<subread<<endl;

		  for (vector<string>::iterator it = subreadit->begin(); it!= subreadit->end(); ++it){
        sort.push_back(pairs[*it]+subread*pairs.size());//}
			//b_builder.set(pairs[*it]+subread*pairs.size());
			   }
      found_kmers.erase (found_kmers.begin());
	subread++;}



cerr<<"Size of the sortd set:		"<<sort.size()<<" should be equal to value of m"<<endl;


size_t m = sort.size();

std::cerr << "stack allocating sdsl::sd_vector_builder base object with n=" << n<< " m=" << m << std::endl;
sdsl::sd_vector_builder b_builder(n, m);
std::sort (sort.begin(), sort.end());
for(std::vector<size_t>::iterator it = sort.begin(); it!= sort.end(); ++it){
//	cout<<*it<<endl;
	b_builder.set(*it);}



std::cerr << "sd_vector items: 		" << b_builder.items()  <<" should be equal to sd_vector capacity"<< endl;
cerr<< "sd_vector capacity: 		" << b_builder.capacity() << std::endl;
std::cerr << "Total wall time: 		" << (getMilliCount() - global_t) / 1000.0 << " s" << std::endl;
sdsl::sd_vector<> b(b_builder);
cerr << "Total size: 			"<<size_in_mega_bytes(b) <<" Mb"<<endl;

std::string outfilename = "output_matrix";
sdsl::store_to_file(b, outfilename);
std::cerr << "DONE!" << std::endl;


}
