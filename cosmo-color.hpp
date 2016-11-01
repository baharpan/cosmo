#ifndef COSMO_COLOR_H
#define COSMO_COLOR_H

struct parameters_t {
  std::string input_filename = "";
  std::string color_filename = "";
  std::string output_prefix = "";
  std::string color_mask1 = "";
  std::string color_mask2 = "";
};

int getMilliCount();
int getMilliSpan(int nTimeStart);
void parse_arguments(int argc, char **argv, parameters_t & params);
void test_symmetry(debruijn_graph<> dbg);
void dump_nodes(debruijn_graph<> dbg, uint64_t * colors);
void dump_edges(debruijn_graph<> dbg, uint64_t * colors);
void find_contigs(debruijn_graph<> dbg);
void find_cycles(debruijn_graph<> dbg);











#endif
