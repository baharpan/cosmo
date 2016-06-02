#include <iostream>
#include <fstream>
#include <utility>
#include <ctime>
// TCLAP
#include "tclap/CmdLine.h"

#include <sdsl/bit_vectors.hpp>

// C STDLIB Headers
#include <cstdio>
#include <cstdlib>

#include <libgen.h>

// Custom Headers
#include "uint128_t.hpp"
#include "debug.h"
#include "kmer.hpp"

using namespace std;
using namespace sdsl;

#include <cstdlib>
#include <sys/timeb.h>

int getMilliCount();
int getMilliCount(){
  timeb tb;
  ftime(&tb);
  int nCount = tb.millitm + (tb.time & 0xfffff) * 1000;
  return nCount;
}

int getMilliSpan(int nTimeStart);
int getMilliSpan(int nTimeStart){
  int nSpan = getMilliCount() - nTimeStart;
  if(nSpan < 0)
    nSpan += 0x100000 * 1000;
  return nSpan;
}

string extension = ".rrr";
typedef struct p
{
  std::string input_filename = "";
  int num_colors;
  std::string output_prefix = "";
} parameters_t;

void parse_arguments(int argc, char **argv, parameters_t & params);
void parse_arguments(int argc, char **argv, parameters_t & params)
{
  TCLAP::CmdLine cmd("Cosmo Copyright (c) Alex Bowe (alexbowe.com) 2014", ' ', VERSION);
  TCLAP::UnlabeledValueArg<std::string> input_filename_arg("input",
            "Input file. Currently only supports DSK's binary format (for k<=64).", true, "", "input_file", cmd);
  TCLAP::UnlabeledValueArg<std::string> num_colors_arg("num_colors",
            "Number of colors", true, "", "num colors", cmd);
  string output_short_form = "output_prefix";
  TCLAP::ValueArg<std::string> output_prefix_arg("o", "output_prefix",
            "Output prefix. Graph will be written to [" + output_short_form + "]" + extension + ". " +
            "Default prefix: basename(input_file).", false, "", output_short_form, cmd);
  cmd.parse( argc, argv );
  params.input_filename  = input_filename_arg.getValue();
  params.num_colors  = atoi(num_colors_arg.getValue().c_str());
  params.output_prefix   = output_prefix_arg.getValue();

}

void deserialize_color_bv(ifstream &colorfile, color_bv &value)
{
      colorfile.read((char *)&value, sizeof(color_bv));
}

int main(int argc, char * argv[]) {
  cout <<"Starting\n";
  parameters_t params;
  parse_arguments(argc, argv, params);

  const char * file_name = params.input_filename.c_str();
  cout << "file name: " << file_name << "\n";

  // Open File
  ifstream colorfile(file_name, ios::in|ios::binary);

  colorfile.seekg(0, colorfile.end);
  size_t end = colorfile.tellg();
  colorfile.seekg(0, colorfile.beg);
  cout << "file size: " << end << std::endl;
  cout << "sizeof(color_bv): " << sizeof(color_bv) << std::endl;
  size_t num_color = params.num_colors;
  size_t num_edges = end / sizeof(color_bv);

  bit_vector b = bit_vector(num_edges*num_color, 0);
  size_t cnt0 = 0;
  size_t cnt1 = 0;
  for (size_t i=0; i < num_edges; i++) {
      color_bv value;
      deserialize_color_bv(colorfile, value);
      for (size_t j=0; j < num_color; j++) {
          b[i*num_color + j] = value[j];
          if (b[i*num_color + j] == 0)
              cnt0++;
          else
              cnt1++;
      }
  }
  cout << "edges: " << num_edges << " colors: " << num_color << " Total: " << num_edges * num_color << endl;
  cout << cnt0  << ":" << cnt1 << endl;

  int sysTime = getMilliCount();
  /*
  bit_vector bv(b);
  sysTime = getMilliCount();
  cout << "BV Creation Time: " << getMilliSpan(sysTime) << endl;
  for (size_t i=0; i < num_edges*num_color; i++) {
    bv[i];
  }
  cout << "BV Access Time: " << getMilliSpan(sysTime) << endl;
  cout << "BV Size (MB): " << size_in_mega_bytes(b) << endl;
  */
  sysTime = getMilliCount();
  rrr_vector<63> rrrb(b);
  cout << "RRR Creation Time: " << getMilliSpan(sysTime) << endl;
  sysTime = getMilliCount();
  for (size_t i=0; i < num_edges*num_color; i++) {
    rrrb[i];
  }
  cout << "RRR AccessTime: " << getMilliSpan(sysTime) << endl;
  cout << "RRR Size (MB): " << size_in_mega_bytes(rrrb) << endl;
  char * base_name = basename(const_cast<char*>(params.input_filename.c_str()));
  string outfilename = ((params.output_prefix == "")? base_name : params.output_prefix) + extension;
  store_to_file(rrrb, outfilename);

  /*
  sysTime = getMilliCount();
  sd_vector<> sdb(b);
  cout << "SD Creation Time: " << getMilliSpan(sysTime) << endl;
  sysTime = getMilliCount();
  for (size_t i=0; i < num_edges*num_color; i++) {
    sdb[i];
  }
  cout << "SD Access Time: " << getMilliSpan(sysTime) << endl;
  cout << "SD Size (MB): " << size_in_mega_bytes(sdb) << endl;

  sysTime = getMilliCount();
  hyb_vector<> hyb(b);
  cout << "Hyb Creation Time: " << getMilliSpan(sysTime) << endl;
  sysTime = getMilliCount();
  for (size_t i=0; i < num_edges*num_color; i++) {
    hyb[i];
  }
  cout << "Hyb Access Time: " << getMilliSpan(sysTime) << endl;
  cout << "Hyb Size (MB): " << size_in_mega_bytes(hyb) << endl;
  */
}
