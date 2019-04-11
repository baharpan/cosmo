#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <boost/algorithm/string.hpp>
#include "sdsl/sd_vector.hpp"
#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sys/timeb.h>
#include <future>
#include <cstdio>
#include <cstdlib>
#include "reduce_color.hpp"

using namespace std;
using namespace sdsl;



int main(int argc, char* argv[]){
  if (argc < 4) {
         cerr << "Error: More arguments needed" << endl;
         return -1;
     }

  ifstream fastqfile(argv[1]);
  int num_reads = stoi(argv[2]);
  int k = stoi(argv[3]);
  reduction re;
  re.build_recolored_matrix (k, fastqfile , num_reads);
  fastqfile.close();
  cerr<<"Number of unique kmers:  "<<re.pairs.size()<<endl;
  cerr<<"Number of colors before reduction:  "<<num_reads <<endl;
  cerr<<"Number of colors after reduction:  "<<re.num_color <<endl;
  cerr<<"The reduction percentage:  "<<(num_reads - (double)re.num_color) * 100 /num_reads<<"%\n";
  cerr<<"DONE!"<<endl;

}
