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
#include <boost/algorithm/string.hpp>
#include "bubbles.hpp"
#include <sys/timeb.h>

#include <utility>
#include <ctime>

// TCLAP
#include "tclap/CmdLine.h"

 #include <cstdio>

 #include <cstdlib>
#include "kmer.hpp"

int main(){

  ifstream input("filtered_kmc2_list.packed", ios::in|ios::binary|ios::ate);
  debruijn_graph<> dbg = debruijn_graph<>::load_from_packed_edges(input, "$ACGT"/*, &minus_positions*/);
  input.close();
  cerr << "num_nodes()   : " << dbg.num_nodes() << endl;
  cerr << "num_edges()   : " << dbg.num_edges() << endl;

  resistome re;
  vector<size_t> file_size = re.file_size_function();
  size_t number_of_subreads = file_size[0];
  std::map<string,size_t> pairs = re.index_maker();

  sdsl::sd_vector<> b;
  sdsl::load_from_file(b,"output_matrix");

  cerr<<"number of subreads: "<<number_of_subreads<<endl;
  re.find_incomming_ougoing_degree_noposition(dbg);
  re.SNP(dbg,b,number_of_subreads,pairs);
}
