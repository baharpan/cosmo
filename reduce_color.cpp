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
#include "reduce_color.hpp"

int main(){
reduction re;
re.build_recolored_matrix ();
}
