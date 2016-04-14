// TODO: fix up pointer and malloc to use smart arrays
// STL Headers
//#include <fstream>
#include <iostream>
//#include <algorithm>
#include <utility>

// TCLAP
#include "tclap/CmdLine.h"

// BOOST
#include <boost/range/adaptor/transformed.hpp>     // Map function to inputs
//#include <boost/range/algorithm/copy.hpp>
//#include <boost/tuple/tuple.hpp>
//#include <boost/tuple/tuple_comparison.hpp> // for uniqued over zipped iterators
//#include <boost/iterator/zip_iterator.hpp>
#define BOOST_RESULT_OF_USE_DECLTYPE // needed to support lambdas in transformed

// C STDLIB Headers
#include <cstdio>
#include <cstdlib>

#include <libgen.h>

// KMC2 Headers
#include "kmc_api/kmc_file.h"

// Custom Headers
#include "uint128_t.hpp"
#include "kmer.hpp"
#include "io.hpp"
#include "sort.hpp"
#include "dummies.hpp"
#include "debug.h"


using namespace std;
using namespace boost::adaptors;

const static string extension = ".packed";

template <typename kmer_t, class Visitor>
void convert(kmer_t * kmers, size_t num_kmers, const uint32_t k, Visitor visit, bool swap, std::vector<color_bv> &colors) {
  // Convert the nucleotide representation to allow tricks
  convert_representation(kmers, kmers, num_kmers, swap);
  //print_kmers(std::cout, kmers, num_kmers, k);
  // Append reverse complements
  #ifdef ADD_REVCOMPS
  size_t revcomp_factor = 2;
  transform(kmers, kmers + num_kmers, kmers + num_kmers, reverse_complement<kmer_t>(k));
//  print_kmers(std::cout, kmers, num_kmers * revcomp_factor, k);
  //memcpy(colors + num_kmers, colors, sizeof(uint64_t) * num_kmers); // copy all of the color masks for the reverse kmers
  std::copy(colors.begin(), colors.begin() + num_kmers, colors.begin() + num_kmers);
  #else
  size_t revcomp_factor = 1;
  #endif

  // NOTE: There might be a way to do this recursively using counting (and not two tables)
  // After the sorting phase, Table A will in <colex(node), edge> order (as required for output)
  // and Table B will be in colex(row) order. Having both tables is helpful for detecting the missing dummy
  // edges with a simple O(N) merge-join algorithm.
  // NOTE: THESE SHOULD NOT BE FREED (the kmers array is freed by the caller)
  kmer_t * table_a = kmers;
  kmer_t * table_b = kmers + num_kmers * revcomp_factor; // x2 because of reverse complements
  std::vector<color_bv>::iterator colors_a = colors.begin();
  std::vector<color_bv>::iterator colors_b = colors.begin() + num_kmers * revcomp_factor;
  // Sort by last column to do the edge-sorted part of our <colex(node), edge>-sorted table
  colex_partial_radix_sort<DNA_RADIX>(table_a, table_b, num_kmers * revcomp_factor, 0, 1,
                                      &table_a, &table_b, get_nt_functor<kmer_t>(),
                                      0, 0, 0, 0,
                                      colors_a, colors_b, &colors_a, &colors_b);
  // Sort from k to last column (not k to 1 - we need to sort by the edge column a second time to get colex(row) table)
  // Note: The output names are swapped (we want table a to be the primary table and b to be aux), because our desired
  // result is the second last iteration (<colex(node), edge>-sorted) but we still have use for the last iteration (colex(row)-sorted).
  // Hence, table_b is the output sorted from [hi-1 to lo], and table_a is the 2nd last iter sorted from (hi-1 to lo]
  colex_partial_radix_sort<DNA_RADIX>(table_a, table_b, num_kmers * revcomp_factor, 0, k,
                                      &table_b, &table_a, get_nt_functor<kmer_t>(),
                                      0, 0, 0, 0,
                                      colors_a, colors_b, &colors_b, &colors_a);

  // outgoing dummy edges are output in correct order while merging, whereas incoming dummy edges are not in the correct
  // position, but are sorted relatively, hence can be merged if collected in a previous pass
  // count dummies (to allocate space)
  size_t num_incoming_dummies = count_incoming_dummy_edges(table_a, table_b, num_kmers * revcomp_factor, k);
  TRACE("num_incoming_dummies: %zu\n", num_incoming_dummies);
  // allocate space for dummies -> we need to generate all the $-prefixed dummies, so can't just use an iterator for the
  // incoming dummies (the few that we get from the set_difference are the ones we apply $x[0:-1] to, so we need (k-1) more for each
  // (to make $...$x[0])... times two because we need to radix sort these bitches.
  #ifdef ALL_DUMMIES
  size_t all_dummies_factor = (k-1);
  size_t dummy_table_factor = 2;
  #else
  size_t all_dummies_factor = 1;
  size_t dummy_table_factor = 1;
  #endif
  // Don't have to alloc if we aren't preparing all dummies, but this option is only used for testing. Usually we want them
  kmer_t * incoming_dummies = (kmer_t*) malloc(num_incoming_dummies * all_dummies_factor * dummy_table_factor * sizeof(kmer_t));
  if (!incoming_dummies) {
    cerr << "Error allocating space for incoming dummies" << endl;
    exit(1);
  }
  // We store lengths because the prefix before the <length> symbols on the right will all be $ signs
  // this is a cheaper way than storing all symbols in 3 bits instead (although it means we need a varlen radix sort)
  uint8_t * incoming_dummy_lengths = (uint8_t*) malloc(num_incoming_dummies * all_dummies_factor * dummy_table_factor * sizeof(uint8_t));
  if (!incoming_dummy_lengths) {
    cerr << "Error allocating space for incoming dummy lengths" << endl;
    exit(1);
  }
  // extract dummies
  printf("num_kmers=%zu; revcomp_factor=%zu; incoming_dummmies=%zu\n", num_kmers, revcomp_factor, num_incoming_dummies);
  find_incoming_dummy_edges(table_a, table_b, num_kmers*revcomp_factor, k, incoming_dummies);
  // add extra dummies
  #ifdef ALL_DUMMIES
  prepare_incoming_dummy_edges(incoming_dummies, incoming_dummy_lengths, num_incoming_dummies, k-1);
  #else
  // Just set the lengths for merging
  memset(incoming_dummy_lengths, k-1, num_incoming_dummies);
  #endif

  kmer_t * dummies_a = incoming_dummies;
  uint8_t * lengths_a = incoming_dummy_lengths;
  // sort dummies (varlen radix)
  // dont need to sort if not adding the extras, since already sorted
  #ifdef ALL_DUMMIES
  kmer_t * dummies_b = incoming_dummies + num_incoming_dummies * (k-1);
  uint8_t * lengths_b = incoming_dummy_lengths + num_incoming_dummies * (k-1);
  colex_partial_radix_sort<DNA_RADIX>(dummies_a, dummies_b, num_incoming_dummies*(k-1), 0, 1,
                                      &dummies_a, &dummies_b, get_nt_functor<kmer_t>(),
                                      lengths_a, lengths_b, &lengths_a, &lengths_b);
  // Don't need the last iteration (i.e. dont need to go to 0) since we arent doing a set difference like above
  colex_partial_radix_sort<DNA_RADIX>(dummies_a, dummies_b, num_incoming_dummies*(k-1), 1, k-1,
                                      &dummies_a, &dummies_b, get_nt_functor<kmer_t>(),
                                      lengths_a, lengths_b, &lengths_a, &lengths_b);
  #endif
  merge_dummies(table_a, table_b, num_kmers*revcomp_factor, k,
                dummies_a, num_incoming_dummies*all_dummies_factor,
                lengths_a,
                // edge_tag needed to distinguish between dummy out edge or not...
                [=](edge_tag tag, const kmer_t & x, const uint32_t x_k, size_t first_start_node, bool first_end_node) {
                  // TODO: this should be factored into a class that prints full kmers in ascii
                  // then add a --full option
                  #ifdef VERBOSE // print each kmer to stderr for testing
                  if (tag == out_dummy)
                    cerr << kmer_to_string(get_start_node(x), k-1, k-1) << "$";
                  else
                    cerr << kmer_to_string(x, k, x_k);
                  cerr << " " << first_start_node << " " << first_end_node << endl;
                  #endif
                  visit(tag, x, x_k, first_start_node, first_end_node);
                });
  // TODO: impl external-merge (for large input. Read in chunk, sort, write temp, ext merge + add dummies to temp, 3-way-merge)
  // TODO: use SSE instructions or CUDA if available (very far horizon)

  free(incoming_dummies);
  free(incoming_dummy_lengths);
}

typedef struct p
{
    //bool ascii = false;
    std::string input_filename = "";
    std::string output_prefix = "";
    bool cortex = false;
    bool kmc = false;
} parameters_t;

void parse_arguments(int argc, char **argv, parameters_t & params);
void parse_arguments(int argc, char **argv, parameters_t & params)
{
  TCLAP::CmdLine cmd("Cosmo Copyright (c) Alex Bowe (alexbowe.com) 2014", ' ', VERSION);
  /* // Add this option after refactoring the visitors (for now just compile with DEBUG if you want printed edges)
  TCLAP::SwitchArg ascii_arg("a", "ascii",
            "Outputs *full* edges (instead of just last nucleotide) as ASCII.",
            cmd, false);
  */
  TCLAP::UnlabeledValueArg<std::string> input_filename_arg("input",
            "Input file. Currently only supports DSK's binary format (for k<=64).", true, "", "input_file", cmd);
  string output_short_form = "output_prefix";
  TCLAP::ValueArg<std::string> output_prefix_arg("o", "output_prefix",
            "Output prefix. Results will be written to [" + output_short_form + "]" + extension + ". " +
            "Default prefix: basename(input_file).", false, "", output_short_form, cmd);
  TCLAP::SwitchArg cortex_arg("c", "cortex",
            "Input file is cortex binary",
            cmd, false);
  TCLAP::SwitchArg kmc_arg("k", "KMC",
            "Input file is a KMC binary",
            cmd, false);
  cmd.parse( argc, argv );
  //params.ascii         = ascii_arg.getValue();
  params.input_filename  = input_filename_arg.getValue();
  params.output_prefix   = output_prefix_arg.getValue();
  params.kmc          = kmc_arg.getValue();
  params.cortex = cortex_arg.getValue();
}
void serialize_color_bv(std::ofstream &cfs, const color_bv &color, unsigned num_colors)//std::vector<color_bv>::iterator &colors, uint64_t index)
{
    for (int i = num_colors -1; i >= 0; --i) {
        cfs.write((char *)&color[i], sizeof(color_bv[i])); //FIXME: Is this the right endianness that SDSL-lite expects
    }
}

#ifdef ADD_REVCOMPS
static const size_t revcomp_factor = 2;
#else
static const size_t revcomp_factor = 1;
#endif


uint64_t* allocate_blocks(const uint32_t kmer_num_bits, const size_t num_kmers, const char* file_name)
{
    if (kmer_num_bits > MAX_BITS_PER_KMER) {
        fprintf(stderr, "ERROR: Kmers larger than %zu bits are not currently supported."
                " %s uses %d bits per kmer (possibly corrupt?).\n",
                MAX_BITS_PER_KMER, file_name, kmer_num_bits);
        exit(EXIT_FAILURE);
    }

    
    if (num_kmers == 0) {
        fprintf(stderr, "ERROR: File %s has no kmers (possibly corrupt?).\n", file_name);
        exit(EXIT_FAILURE);
    }
    
    assert(kmer_num_bits % 64 == 0);
    uint32_t kmer_num_blocks = (kmer_num_bits  / 8 ) / sizeof(uint64_t) ; 
    
    TRACE("kmer_num_bits, k = %d, %d\n", kmer_num_bits, kmer_size);
    TRACE("kmer_num_blocks = %d\n", kmer_num_blocks);


    // ALLOCATE SPACE FOR KMERS (done in one malloc call)
    // x 4 because we need to add reverse complements, and then we have two copies of the table
    size_t kmer_blocks_size = num_kmers * 2 * revcomp_factor * sizeof(uint64_t) * kmer_num_blocks  /* FIXME: CONSERVATIVELY ASSUMES NO REPEAT KMERS ACROSS COUNTS */;
    uint64_t * kmer_blocks = (uint64_t*)malloc(kmer_blocks_size);

    std::cerr << "Allocating " << kmer_blocks_size << " bytes for kmer_blocks."  << std::endl;
    if (!kmer_blocks) {
        cerr << "Error allocating space for kmers" << endl;
        exit(1);
    }
    return kmer_blocks;
}


void load_kmers(const parameters_t& params, std::vector<color_bv>& kmer_colors, uint64_t*& kmer_blocks, uint32_t& kmer_num_bits, uint32_t& kmer_size, size_t& num_kmers, uint32_t& num_colors)
{
    const char* file_name = params.input_filename.c_str();

    int handle = -1;
    if ( (handle = open(file_name, O_RDONLY)) == -1 ) {
        fprintf(stderr, "ERROR: Can't open file: %s\n", file_name);
        exit(EXIT_FAILURE);
    }

    
    if (params.cortex) {
        std::cerr << "Reading cortex file " << file_name << std::endl;
        if ( !cortex_read_header(handle, &kmer_num_bits, &kmer_size) ) {
            fprintf(stderr, "ERROR: Error reading cortex_file %s\n", file_name);
            exit(EXIT_FAILURE);
        }

        if ( cortex_num_records(handle, kmer_num_bits, num_kmers, num_colors) == -1) {
            fprintf(stderr, "Error seeking cortex file %s\n", file_name);
            exit(EXIT_FAILURE);
        }
        if (num_colors > NUM_COLS) {
            fprintf(stderr, "Cortex file %s contains %d colors which exceeds the compile time limit of %d.  Please recompile with NUM_COLS=%d (or larger).\n", file_name, num_colors, NUM_COLS, num_colors);
            exit(EXIT_FAILURE);
        }
        printf("Got num record %zu \n", num_kmers);
        printf("Got num colors (capacity, not occupancy) %d \n", num_colors);
        // printf("NUM_COLS=%zu\n", NUM_COLS);
        // printf("Each entry in .colors file will occupy %d bytes.\n", sizeof(color_bv));
        kmer_blocks = allocate_blocks(kmer_num_bits, num_kmers, file_name);

            ///
        kmer_colors.reserve(num_kmers * 2 * revcomp_factor );
        printf("Reading kmers\n");
        size_t num_records_read = cortex_read_kmers(handle, kmer_num_bits, num_colors, kmer_size, kmer_blocks, kmer_colors);
        if (num_records_read == 0) {
            fprintf(stderr, "Error reading file %s\n", file_name);
            exit(EXIT_FAILURE);
        }
        assert(num_records_read == num_kmers);
        TRACE("num_records_read = %zu\n", num_records_read);
        
        printf("num_kmers = %zu and num_records_read=%zu\n", num_kmers, num_records_read);


    } else if (params.kmc) {
        std::vector<CKMCFile *>kmer_data_bases; //FIXME: move out of global
        std::cerr << "Reading KMC2 database list file " << file_name << std::endl;
        uint64 peak_kmers;
        if ( !kmc_read_header(file_name, kmer_num_bits, kmer_size, peak_kmers, num_colors, kmer_data_bases) ) {
            fprintf(stderr, "ERROR: Error reading databases listed in KMC2 list file '%s'\n", file_name);
            exit(EXIT_FAILURE);
        }

        std::cerr << "Will read a maximum of " << peak_kmers << " " << kmer_size << "-mers from some color, each with " << kmer_num_bits << " bits." << std::endl;

        if (num_colors > NUM_COLS) {
            fprintf(stderr, "KMC file %s contains %d colors which exceeds the compile time limit of %d.  Please recompile with NUM_COLS=%d (or larger).\n", file_name, num_colors, NUM_COLS, num_colors);
            exit(EXIT_FAILURE);
        }
        
        kmer_colors.reserve(peak_kmers);

        // preallocate enough space for the color0, as we'll certainly need at least that much space
        std::vector<uint64_t> kmer_block_buffer;
        kmer_block_buffer.reserve(peak_kmers * /*kmer_num_blocks*/ (kmer_num_bits  / 8 ) / sizeof(uint64_t) );



        num_kmers = kmc_read_kmers(handle, kmer_num_bits, num_colors, kmer_size, kmer_block_buffer, kmer_colors, kmer_data_bases);
        printf("num_kmers = %zu\n", num_kmers);
        TRACE("num_kmers = %zu\n", num_kmers);

        kmer_blocks = allocate_blocks(kmer_num_bits, num_kmers, file_name);

        // COPY COLORS
        // COPY KMERS
        std::copy(kmer_block_buffer.begin(), kmer_block_buffer.end(), kmer_blocks);
        kmer_colors.resize(num_kmers * 2 * revcomp_factor);
            
    } else {
        std::cerr << "Reading DSK file " << file_name << std::endl;
        if ( !dsk_read_header(handle, &kmer_num_bits, &kmer_size) ) {
            fprintf(stderr, "ERROR: Error reading file %s\n", file_name);
            exit(EXIT_FAILURE);
        }
        TRACE(">> READING DSK FILE\n");

        if ( dsk_num_records(handle, kmer_num_bits, &num_kmers) == -1) {
            fprintf(stderr, "Error seeking file %s\n", file_name);
            exit(EXIT_FAILURE);
        }
        kmer_blocks = allocate_blocks(kmer_num_bits, num_kmers, file_name);
        size_t num_records_read = dsk_read_kmers(handle, kmer_num_bits, kmer_blocks, kmer_size);
        if (num_records_read == 0) {
            fprintf(stderr, "Error reading file %s\n", file_name);
            exit(EXIT_FAILURE);
        }
        TRACE("num_records_read = %zu\n", num_records_read);
        assert ( num_records_read == num_kmers);

    }

    close(handle);
}

           

int main(int argc, char * argv[])
{
    parameters_t params;
    parse_arguments(argc, argv, params);


    std::vector<color_bv> kmer_colors;
    uint64_t * kmer_blocks = 0;

    uint32_t kmer_num_bits = 0;
    uint32_t kmer_size = 0;
    size_t num_kmers = 0;
    uint32_t num_colors = 0;
    
    load_kmers(params,  kmer_colors, kmer_blocks, kmer_num_bits, kmer_size, num_kmers, num_colors);

    char * base_name = basename(const_cast<char*>(params.input_filename.c_str()));

    string outfilename = (params.output_prefix == "")? base_name : params.output_prefix;
    ofstream ofs;
#ifdef VAR_ORDER
    ofstream lcs;
#endif
#if 1
    static const size_t BUFFER_LEN = 2*1024*1024;
    char buffer_a[BUFFER_LEN];
    ofs.rdbuf()->pubsetbuf(buffer_a, BUFFER_LEN);
#ifdef VAR_ORDER
    char buffer_b[BUFFER_LEN];
    lcs.rdbuf()->pubsetbuf(buffer_b, BUFFER_LEN);
#endif
#endif
    // TODO: Should probably do checking here when opening the file...
    ofs.open(outfilename + extension, ios::out | ios::binary);
#ifdef VAR_ORDER
    lcs.open(outfilename + extension + ".lcs", ios::out | ios::binary);
#endif
    PackedEdgeOutputer out(ofs);
    ofstream cfs;
    cfs.open(outfilename + ".colors", ios::out | ios::binary);


    color_bv ones;
    // create an 'all ones' color_bv
    for (unsigned int citer=0; citer < num_colors; ++citer)
        ones[citer] = 1;

    std::vector<color_bv>::iterator colors = kmer_colors.begin() + num_kmers * revcomp_factor;

    size_t index = 0; 

    if (kmer_num_bits == 64) {
        typedef uint64_t kmer_t;
        size_t prev_k = 0; // for input, k is always >= 1



        convert(kmer_blocks, num_kmers, kmer_size,
                [&](edge_tag tag, const kmer_t & x, const uint32_t this_k, size_t lcs_len, bool first_end_node) {
#ifdef VAR_ORDER
                    out.write(tag, x, this_k, (lcs_len != k-1), first_end_node);
                    char l(lcs_len);
                    lcs.write((char*)&l, 1);
#else
                    out.write(tag, x, this_k, lcs_len, first_end_node);
#endif
                    if (tag == standard) {
                        // cerr << kmer_to_string(x, k, this_k) << "c" << colors[index] << "\n";
                        serialize_color_bv(cfs, colors[index++], num_colors);
                        //cfs.write((char *)&colors[index++], sizeof(uint64_t));
                    }
                    else {
                        //uint64_t ones = -1;
                        serialize_color_bv(cfs, ones, num_colors);
                            //assert(!"Not converted to color_bv yet!");
                        //cfs.write((char *)&ones, sizeof(uint64_t));
                    }
                    prev_k = this_k;
                }, !(params.cortex || params.kmc), kmer_colors);
    }
    else if (kmer_num_bits == 128) {
        typedef uint128_t kmer_t;
        size_t prev_k = 0;
        kmer_t * kmer_blocks_128 = (kmer_t*)kmer_blocks;
        convert(kmer_blocks_128, num_kmers, kmer_size,
                [&](edge_tag tag, const kmer_t & x, const uint32_t this_k, size_t lcs_len, bool first_end_node) {
#ifdef VAR_ORDER
                    out.write(tag, x, this_k, (lcs_len != k-1), first_end_node);
                    char l(lcs_len);
                    lcs.write((char*)&l, 1);
#else
                    out.write(tag, x, this_k, lcs_len, first_end_node);
#endif
                    if (tag == standard) {
                        // cerr << kmer_to_string(x, k, this_k) << "c" << colors[index] << "\n";
                        serialize_color_bv(cfs, colors[index++], num_colors);
                        //cfs.write((char *)&colors[index++], sizeof(uint64_t));
                    }
                    else {
                        //uint64_t ones = -1;
                        serialize_color_bv(cfs, ones, num_colors);
                            //assert(!"Not converted to color_bv yet!");
                        //cfs.write((char *)&ones, sizeof(uint64_t));
                    }
                    prev_k = this_k;
                }, !(params.cortex || params.kmc), kmer_colors);
    }

    out.close();
#ifdef VAR_ORDER
    lcs.flush();
    lcs.close();
#endif
    uint64_t t_k(kmer_size); // make uint64_t just to make parsing easier
    // (can read them all at once and take the last 6 values)
    ofs.write((char*)&t_k, sizeof(uint64_t));
    ofs.flush();
    ofs.close();

    cfs.close();

    free(kmer_blocks);

    return 0;
}
