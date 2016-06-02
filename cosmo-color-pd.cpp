#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <libgen.h> // basename

#include "tclap/CmdLine.h"

#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>

#include "io.hpp"
#include "debruijn_graph.hpp"
#include "algorithm.hpp"
#include "cosmo-color-pd.hpp"

using namespace std;
using namespace sdsl;

#include <sys/timeb.h>

bool trace = false;

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

string extension = ".dbg";



void parse_arguments(int argc, char **argv, parameters_t & params)
{
  TCLAP::CmdLine cmd("Cosmo Copyright (c) Alex Bowe (alexbowe.com) 2014", ' ', VERSION);
  TCLAP::UnlabeledValueArg<std::string> input_filename_arg("input",
            ".packed edge file (output from pack-edges).", true, "", "input_file", cmd);
  TCLAP::UnlabeledValueArg<std::string> color_filename_arg("color",
            ".color file (output from pack-edges).", true, "", "color_file", cmd);
  string output_short_form = "output_prefix";
  TCLAP::ValueArg<std::string> output_prefix_arg("o", "output_prefix",
            "Output prefix. Graph will be written to [" + output_short_form + "]" + extension + ". " +
            "Default prefix: basename(input_file).", false, "", output_short_form, cmd);
  string ref_color = "ref_color";
  TCLAP::ValueArg<std::string> ref_color_arg("a", "ref_color",
	    "Ref color, ref_color [" + ref_color + "]", false, "", ref_color, cmd);
  string sample_mask = "sample_mask";
  TCLAP::ValueArg<std::string> sample_mask_arg("b", "sample_mask",
	    "Sample mask, sample_mask [" + sample_mask + "]", false, "", sample_mask, cmd);
  string ref_fasta = "ref_fasta";
  TCLAP::ValueArg<std::string> ref_fasta_arg("r", "ref_fasta",
	    "Reference FASTA filename, ref_fasta [" + ref_fasta + "]", false, "", ref_fasta, cmd);

  cmd.parse( argc, argv );

  params.input_filename  = input_filename_arg.getValue();
  params.color_filename  = color_filename_arg.getValue();
  params.output_prefix   = output_prefix_arg.getValue();
  params.ref_color     = ref_color_arg.getValue();
  params.sample_mask     = sample_mask_arg.getValue();
  params.ref_fasta = ref_fasta_arg.getValue();
}

static char base[] = {'?','A','C','G','T'};


void test_symmetry(debruijn_graph<> dbg) {
  for (unsigned long x = 0; x<dbg.sigma+1;x++) {
    ssize_t in = dbg.incoming(43, x);
    if (in == -1)
      continue;
    for (unsigned long y = 0; y<dbg.sigma+1;y++) {
      ssize_t out = dbg.outgoing(in, y);
      if (out == -1)
	continue;
      cout << "Incoming " << in <<  ":" << out <<"\n";
    }
  }
}



void dump_nodes(debruijn_graph<> dbg, uint64_t * colors) {
  for (size_t i = 0; i < dbg.num_nodes(); i++) {
    cout << i << ":" << dbg.node_label(i) << colors[dbg._node_to_edge(i)] << "\n";
  }
}


void dump_edges(debruijn_graph<> dbg, uint64_t * colors) {
  for (size_t i = 0; i < dbg.num_edges(); i++) {
    cout << i << "e:" << dbg.edge_label(i) << colors[i] << "\n";
  }
}


ssize_t get_first_node(debruijn_graph<> dbg, rrr_vector<63> &colors, uint64_t ref_color, std::string& ref_fasta_content)
{

    
    // find the first edge which is colored 'ref_color' and whose label we can find in the ref genome
    // (since we generate revcomps and color them the same as given, the first edge's label may not
    //  be found in the ref genome; we'll just keep trying till we find one)
    int num_colors = colors.size() / dbg.num_edges();
    ssize_t zeroth_rank_edge = (unsigned long long)-1;
    size_t pos =  std::string::npos;
    ssize_t node_num = 0;
    std::string node_label;
    std::string query = ref_fasta_content.substr(0, dbg.node_label(0).size()/*hack to get the edge kmer size*/);
    for (; node_num < dbg.num_edges(); ++ node_num) {
        if (colors[node_num * num_colors + ref_color]) {
            zeroth_rank_edge = node_num;
            node_label = dbg.node_label(zeroth_rank_edge);
            if (node_label == query) {
                //if ((pos = ref_fasta_content.find(edge_label)) == 0 /*!= std::string::npos*/){
                std::cerr << "Found fasta start at edge num " << node_num << std::endl;
                return node_num;
                break;
            }
        }
    }

    // if (pos != std::string::npos) {
    //     std::cout << "Found edge " << edge_label << " number " << edge_num << " in reference genome at position " << pos << "." << std::endl;
    //     std::cout << "ref_genome[" << pos - 1 << ":" << pos + edge_label.size() << " = "
    //               << ref_fasta_content.substr(pos - 1, 1 + edge_label.size()) << std::endl;
    // }else{
    //     std::cerr << "ERROR: Can't find any ref color edges in the ref FASTA" << std::endl;
    //     exit(EXIT_FAILURE);
    // }

    // while (pos != 0) {
        
        

    // now traverse the graph backward till we get back to the start of the ref genome
    //while (pos != 0) {
        

    

}

unsigned dna_ord(char c)
{
    switch(c) {
    case 'A' : return 1;
    case 'C': return 2;
    case 'G': return 3;
    case 'T': return 4;
    default: std::cerr << "ERROR: Unknown base '" << c << "'." << std::endl;
        exit(EXIT_FAILURE);
        
    };
}
// this is meant to roughly mimic strcmp(), will return 0 if the paths starting at s_pos and node_k_pos match for the first L nodes
void advance(debruijn_graph<> dbg, const std::string& ref_fasta_content, const unsigned amount, ssize_t& node_k, ssize_t& node_k_pos)
{
    unsigned node_label_size = dbg.k - 1;  
    for (unsigned i = 0; i < amount; ++i) {
        
        ssize_t edge = dbg.outgoing_edge(node_k, dna_ord(ref_fasta_content[node_k_pos + node_label_size]));
        if (edge == -1) {
            std::cerr << "Reference fasta guides graph traversal through invalid path at position " << node_k_pos + node_label_size << "." << std::endl;
            exit(EXIT_FAILURE);
        }
        node_k = dbg._edge_to_node(edge);
        node_k_pos++;
    }
    
    
}


int colored_outdegree(debruijn_graph<> dbg, ssize_t v, const uint64_t sample_mask, unsigned num_colors, rrr_vector<63> &colors)
{
    unsigned out_count = 0;

    // for each symbol of the alphabet
    for (unsigned long x2 = 1; x2 < dbg.sigma + 1; x2++) {

        // if there exists an outgoing edge for that symbol
        ssize_t next_edge = dbg.outgoing_edge(v, x2);
        if (next_edge != -1) {

            // compute the colors of that edge
            uint64_t color_mask = 0;
            for (int c = 0; c < num_colors; c++)
                color_mask |= colors[next_edge * num_colors + c] << c;

            // and if any colors of that edge match the sample set of colors, increment the out degree counter
            if (color_mask & sample_mask) {
                out_count += 1;

            }
        }
    }
    return out_count;
    
    
}

int colored_indegree(debruijn_graph<> dbg, ssize_t v, const uint64_t sample_mask, unsigned num_colors, rrr_vector<63> &colors)
{
    unsigned in_count = 0;
    
    // for each symbol of the alphabet
    for (unsigned long x2 = 1; x2 < dbg.sigma + 1; x2++) {

        // if there exists an incoming edge for that symbol        
        ssize_t next_edge = dbg.incoming(v, x2);
        if (next_edge != -1) {

            // compute the colors of that edge            
            uint64_t color_mask = 0;
            for (int c = 0; c < num_colors; c++)
                color_mask |= colors[next_edge * num_colors + c] << c;

            // and if any colors of that edge match the sample set of colors, increment the out degree counter
            if (color_mask & sample_mask) {
                in_count += 1;

            }
        }
    }
    return in_count;
    
    
}

//FIXME: add asserts to check the color for the reference genome; it's somewhat annoying user has to specify both the color and the reference genome; we should be able to derive one from the other for the overall flow and having the data replicated could lead to inconsistency errors.
/*FIXME: sample_mask is limited to 64 colors*/
// FIXME: how to handle multiple colors?  Treat them all the same? Or do we have to bookkeep individual color results during one traversal.  What to do about divergence in the latter case?
void get_supernode(debruijn_graph<> dbg, const ssize_t& node_i, const uint64_t sample_mask , std::vector<ssize_t>& s, unsigned num_colors,  rrr_vector<63> &colors, int& node_i_pos_in_supernode)
{
    node_i_pos_in_supernode = 0;

    //FIXME: what happens if ref enters supernode at the side? do we need to backup to the beginning?  supplement seems to imply no
    if (colored_outdegree(dbg, node_i, sample_mask, num_colors, colors) == 1) {
        
        for (unsigned long x = 1; x < dbg.sigma + 1; x++) { // iterate through the alphabet of outgoing edges from node i
            // follow each strand or supernode
            ssize_t edge = dbg.outgoing_edge(node_i, x);
            if (edge == -1)
                continue;
            //branch[branch_num] += base[x];
            // build color mask
            uint64_t color_mask = 0;
            for (int c = 0; c < num_colors; c++)
                color_mask |= colors[edge * num_colors + c] << c;
            if (color_mask & sample_mask) {
                s.push_back(node_i);
                

                // walk along edges until we encounter 
                ssize_t node_pos = dbg._edge_to_node(edge);
                while (colored_indegree(dbg, node_pos, sample_mask, num_colors, colors) <= 1 /*FIXME: should be == 1*/ && colored_outdegree(dbg, node_pos, sample_mask, num_colors, colors) == 1) {

                    ssize_t next_edge = 0;
                    for (unsigned long x2 = 1; x2 < dbg.sigma + 1; x2++) { // iterate through the alphabet
                        next_edge = dbg.outgoing_edge(node_pos, x2);
                        if (next_edge != -1) {
                            uint64_t color_mask = 0;
                            for (int c = 0; c < num_colors; c++)
                                color_mask |= colors[next_edge * num_colors + c] << c;
                            if (color_mask & sample_mask) {
                                s.push_back(node_pos);
                                break;
                                
                            } 
                                                            

                            //branch[branch_num] += base[x2];

                        }
                    }
                    node_pos = dbg._edge_to_node(next_edge);
                    //cout << node_pos << ":" << dbg.node_label(node_pos) << "\n";
                }
                if (trace) {
                    std::cout << "    terminal supernode node label = " << dbg.node_label(node_pos) << std::endl;
                    std::cout << "   terminal supernode node outdegree = " << colored_outdegree(dbg, node_pos, sample_mask, num_colors, colors) << std::endl;
                    std::cout << "   terminal supernode node indegree = " << colored_indegree(dbg, node_pos, sample_mask, num_colors, colors) << std::endl;
                }
            }
        }
    }
}

// int get_divergent(debruijn_graph<> dbg, const std::string& ref_fasta_content, const std::vector<ssize_t>& s, ssize_t node_k, ssize_t node_k_pos)
// {
//     for (unsigned i = 0;; ++i, advance(dbg, ref_fasta_content, 1, node_k, node_k_pos)) {
//         if (node_k != s[i]) return i;
//     }
    
        

// }



// int path_cmp(debruijn_graph<> dbg, const std::string& ref_fasta_content, const std::vector<ssize_t>& s, ssize_t s_pos, ssize_t node_k, ssize_t node_k_pos, const unsigned L)
// {
//     for (unsigned int i = 0; i < L; ++i, advance(dbg, ref_fasta_content, 1, node_k, node_k_pos), ++s_pos) {
//         if (s[s_pos] != node_k) return 1;
//     }

//     return 0;
// }

unsigned int match_length(debruijn_graph<> dbg, const std::string& ref_fasta_content, const std::vector<ssize_t>& s, ssize_t s_pos, ssize_t node_k, ssize_t node_k_pos, unsigned bound = -1)
{
    unsigned int i = 0;
    
    for (; s_pos < s.size() && node_k_pos < ref_fasta_content.size(); ++i, advance(dbg, ref_fasta_content, /* amount */ 1, node_k, node_k_pos), ++s_pos) {
        if (s[s_pos] != node_k || i == bound) return i;
    }
    return i;
    
    
}

void dump_supernode(debruijn_graph<> dbg, const std::vector<ssize_t>& s, ssize_t lflanks, ssize_t lflanke, ssize_t rflanks)
{
    assert(s.size());
    
    //std::cout << "   Divergent supernode matches ref at (k-1)-mers [" <<   lflanks << ", " << lflanke << ") and [" << rflanks << ", " << s.size() << "). Label: " <<  dbg.node_label(s[0]);

    for (unsigned i = lflanke; i < rflanks; ++i) {
        std::string lab = dbg.node_label(s[i]);
        
        std::cout << lab[lab.size()-1];

    }
    std::cout << std::endl;
    

}

// S is a path representing the supernode in the sample graph
// b is an index into S where sample and reference diverage
// L is the desired flank size
// i is the common node
// M is the maximum variant size
// n is the reference sequence
const unsigned L = 1; // number of (k-1)-mers to match in each flank
//const unsigned M = 17000;
const unsigned M = 200000; // "a maximum size M of variant to be searched for" -CORTEX supplement ASSUMED_EQUAL_TO "'--max_var_len INT'  Maximum variant size searched for."  -CORTEX manual

//node_i_pos = 253395. node_i = 411788.
//node_i_pos = 253396. node_i = 2383686.

// for KMC2 k=32, Found fasta start at edge num 7471236
//  and node_i_pos = 0. node_i = 7471236.

void find_divergent_paths(debruijn_graph<> dbg, rrr_vector<63> &colors, uint64_t ref_color, uint64_t sample_mask, std::string& ref_fasta_content)
{
    int variant_num = 0;
    int num_colors = colors.size() / dbg.num_edges();
    ssize_t first_node =  get_first_node(dbg, colors, ref_color, ref_fasta_content);

    // all these are for KMC's K=31
//    ssize_t first_node = 1875943; // get_first_node(dbg, colors, ref_color, ref_fasta_content);
//    ssize_t first_node = 2383686; // get_first_node(dbg, colors, ref_color, ref_fasta_content);
//    ssize_t first_node = 411788; // get_first_node(dbg, colors, ref_color, ref_fasta_content);    
    unsigned node_label_size = dbg.k - 1;
    
    ssize_t node_i = first_node; // cdbg node labeled with a k-mer existing in the reference sequence
    ssize_t node_i_pos = 0;  // starting position in the reference sequence for the above k-mer
//    ssize_t node_i_pos = 253396;  // starting position in the reference sequence for the above k-mer
//    ssize_t node_i_pos = 253395;  // starting position in the reference sequence for the above k-mer    

    
    // for each node in ref_fasta, if it starts a sample supernode,
    // scan through ref fasta looking for the other end of the supernode up to M nodes away
    // /* - 2 in the following because we need min of 3 (k-1)-mers to have a bubble */
    while(node_i_pos < ref_fasta_content.size() - node_label_size - 2 ) {
        std::cerr << "node_i_pos = " << node_i_pos  << ". node_i = " << node_i <<"." << std::endl;
        std::vector<ssize_t> s; // supernode
        int node_i_pos_in_supernode = -1;
        trace = false; //(node_i_pos == 253396); // turn on tracing for this node
        get_supernode(dbg, node_i, sample_mask, s, num_colors, colors, node_i_pos_in_supernode);

        
        if (s.size()) {


            ssize_t overlap_len = match_length(dbg, ref_fasta_content, s, node_i_pos_in_supernode,
                                               node_i, node_i_pos);
            ssize_t b = node_i_pos_in_supernode +  overlap_len;
            std::cerr << "    Got supernode of size " << s.size() << " which matches the ref genome for " << overlap_len << " nodes" << std::endl;
            if (trace) {
                std::cerr << "   node_i label = " << dbg.node_label(node_i) << std::endl;
                std::cerr << "   node_i = " << node_i << std::endl;
                dump_supernode(dbg, s, node_i_pos_in_supernode, b, s.size() - L);
            }
            std::cerr << "    node_i_pos_in_supernode = " << node_i_pos_in_supernode
                      << " overlap_len = " <<  overlap_len << std::endl;
            // if the ref matches the supernode to the end, we can skip ahead
            if (b == s.size()) {
                std::cerr << "    Skipping remaining " << overlap_len
                          << " nodes that match supernode." << std::endl;
                advance(dbg, ref_fasta_content, overlap_len, node_i, node_i_pos);
                continue;
            }
            
        
            if (overlap_len >= L ) {

                std::cerr << "    searching for right flank match..." << std::endl;
                ssize_t node_j = node_i;
                ssize_t node_j_pos = node_i_pos;
                advance(dbg, ref_fasta_content, overlap_len + 1, node_j, node_j_pos);
                for   (; node_j_pos <= node_i_pos + M;
                       advance(dbg, ref_fasta_content, /* amount */ 1, node_j, node_j_pos)) {
                    int right_overlap = match_length(dbg, ref_fasta_content, s, s.size() - L, node_j,
                                                     node_j_pos, L);
                    if (trace) std::cerr << "    overlap search depth: " << node_j_pos - node_i_pos << " found overlap: " << right_overlap << " node_j label: " << dbg.node_label(node_j)
                                         << " s[s.size() - L] label: " << dbg.node_label(s[s.size() - L]) << std::endl;
                    if (L == right_overlap) {
                        variant_num += 1;
                        std::cerr << "Found bubble " << variant_num << " with " << right_overlap << " node rightmost flank overlap starting at reference position " << node_j_pos << std::endl;

                        // 5' flank
                        std::cout << ">var_" << variant_num << "_5p_flank length:" << overlap_len + dbg.k - 1 << " INFO:KMER=" << dbg.k << std::endl;
                        std::cout << ref_fasta_content.substr(node_i_pos, overlap_len + dbg.k - 1) << std::endl;

                        // branch_1 : reference branch
                        unsigned long long branch_1_start = node_i_pos + overlap_len + dbg.k - 1;
                        unsigned long long branch_1_length = node_j_pos - branch_1_start;
                        std::cout << ">var_" << variant_num << "_branch_1 length:" << branch_1_length << " kmer:" << dbg.k << std::endl;
                        std::cout << ref_fasta_content.substr(branch_1_start, branch_1_length) << std::endl;

                        // branch_2
                        // in the first bubble in the ecoli6 experiment, the supernode branch was 2nd in the output file, so we'll do it that way
                        std::cout << ">var_" << variant_num << "_branch_2 length:" << s.size() - L - b << " kmer:" << dbg.k << std::endl;
                        dump_supernode(dbg, s, node_i_pos_in_supernode, b, s.size() - L);

                        // 3' flank
                        std::cout << ">var_" << variant_num << "_3p_flank length:" << right_overlap + dbg.k - 1 << " kmer:" << dbg.k << std::endl;
                        std::cout << ref_fasta_content.substr(node_j_pos, right_overlap + dbg.k - 1) << std::endl;


                        break;
                    }
                }
            }
            advance(dbg, ref_fasta_content, overlap_len, node_i, node_i_pos);
            continue;
        }
        advance(dbg, ref_fasta_content, 1, node_i, node_i_pos);
        
    }
}

const char *const starts[] = {"GCCATACTGCGTCATGTCGCCCTGACGCGC","GCAGGTTCGAATCCTGCACGACCCACCAAT","GCTTAACCTCACAACCCGAAGATGTTTCTT","AAAACCCGCCGAAGCGGGTTTTTACGTAAA","AATCCTGCACGACCCACCAGTTTTAACATC","AGAGTTCCCCGCGCCAGCGGGGATAAACCG","GAATACGTGCGCAACAACCGTCTTCCGGAG"};
    
void find_bubbles(debruijn_graph<> dbg, rrr_vector<63> &colors, uint64_t ref_color, uint64_t sample_mask)
{
    int t = getMilliCount();
    int num_colors = colors.size() / dbg.num_edges();
    //uint64_t combined_mask = ref_color | sample_mask;
    bit_vector visited = bit_vector(dbg.num_nodes(), 0);
    cout << "Starting to look for bubbles\n";
    std::vector<std::string> branch(2);
    bool found_miss = false;
    for (size_t node_i = 0; node_i < dbg.num_nodes(); node_i++) {
        ssize_t start_node = node_i; // place to store start of branch kmer
        std::string start_label(dbg.node_label(start_node));
        found_miss = false;
        // for (int si = 0; si < 7; ++si) {
        //     if (!start_label.compare(starts[si])) {
        //         std::cerr << "Found missing start node " << starts[si] << " outdegree: " << dbg.outdegree(i) << std::endl;
        //         found_miss = true;
        //     }
        // }

        // cout << "Node " << i << ":" << dbg.node_label(i) << " color: " << color_mask << "\n";
        if (!visited[node_i] && dbg.outdegree(node_i) == 2) { //FIXME: why do we only care about outdegree == 2?
            // initialize bubble tracking variables
            int branch_num = 0;
            ssize_t end[2]; // place to store end of branch kmer

            branch[0].clear();
            branch[1].clear();

            int branch_offset = 0;
            uint64_t branch_color[2];


            // start of a bubble handling
            for (unsigned long x = 1; x < dbg.sigma + 1; x++) { // iterate through the alphabet of outgoing edges from node i
                // follow each strand or supernode
                ssize_t edge = dbg.outgoing_edge(node_i, x);
                if (edge == -1)
                    continue;
                branch[branch_num] += base[x];
                // build color mask
                uint64_t color_mask = 0;
                for (int c = 0; c < num_colors; c++)
                    color_mask |= colors[edge * num_colors + c] << c;
                branch_color[branch_num] = color_mask;

                // walk along edges until we encounter 
                ssize_t node_pos = dbg._edge_to_node(edge);
                while (dbg.indegree(node_pos) == 1 && dbg.outdegree(node_pos) == 1) {
                    visited[node_pos] = 1;
                    ssize_t next_edge = 0;
                    for (unsigned long x2 = 1; x2 < dbg.sigma + 1; x2++) { // iterate through the alphabet
                        next_edge = dbg.outgoing_edge(node_pos, x2);
                        if (next_edge != -1) {
                            branch[branch_num] += base[x2];
                            break;
                        }
                    }
                    node_pos = dbg._edge_to_node(next_edge);
                    //cout << node_pos << ":" << dbg.node_label(node_pos) << "\n";
                }
                if (found_miss) {
                    std::cerr << "dbg.indegree(node_pos) = " << dbg.indegree(node_pos)  << " dbg.outdegree(node_pos) = " << dbg.outdegree(node_pos)  << std::endl;
                    ssize_t next_edge = 0;
                    std::cerr << "outgoing bases: ";
                    for (unsigned long x2 = 1; x2 < dbg.sigma + 1; x2++) { // iterate through the alphabet
                        next_edge = dbg.outgoing_edge(node_pos, x2);
                        uint64_t color_mask = 0;

                        if (next_edge != -1) {
                            for (int c = 0; c < num_colors; c++)
                                color_mask |= colors[next_edge * num_colors + c] << c;


                            std::cerr << base[x2] << " (c " << color_mask << ") (p " << dbg._edge_to_node(next_edge) << ")" << std::endl;
;
                        }
                    }
                    std::cerr << std::endl;
                    
                }
                // cout << "Stopped due to : " << dbg.indegree(node_pos) << ":" << dbg.outdegree(node_pos) << ":" << branch_offset << "\n";

                end[branch_num++] =  (dbg.indegree(node_pos) > 1) ? node_pos : 0;
                branch_offset = 0;
            }
            // check if both branches ended on the same kmer and they pass the requested color masks
            //cout << "Trying " << branch_color[0] << ":" << branch_color[1] << " " << end[0] << ":" << end[1] <<"\n";
            //cout << ref_color << ":" << sample_mask << "\n";
            //cout << "PutativeStart flank: " << dbg.node_label(start) << " c: " << branch_color[0] << ":" << branch_color[1] << "\n";
            if (found_miss) {
                std::cerr << "arm sizes: " << branch[0].size() << "  " << branch[1].size() << std::endl;
                std::cerr << branch[0] << std::endl << branch[1] << std::endl;
            }
                                    
            // check same end node
            if ((end[0] && end[0] == end[1]) ) {
                if (found_miss)  std::cerr << "Missing bubble passed end check" << std::endl;
                // check color:
                if (true || ((ref_color & branch_color[0] && !(~ref_color & branch_color[0]) &&
                  sample_mask & branch_color[1] && !(~sample_mask & branch_color[1])) || 
                 (ref_color & branch_color[1] && !(~ref_color & branch_color[1]) &&
                  sample_mask & branch_color[0] && !(~sample_mask & branch_color[0])))) {
                    cout << "\nStart flank: " << dbg.node_label(start_node) << " c: " << branch_color[0] << ":" << branch_color[1] << "\n";
                    cout << "Branch: " << branch[0] << "\n";
                    cout << "Branch: " << branch[1] << "\n";
                    cout << "End flank: " << dbg.node_label(end[0]) << "\n";
                    if (found_miss) std::cerr << "Reported 'missing' bubble" << std::endl;
                }
            }
        }
    }
    cerr << "Find bubbles time: " << getMilliSpan(t) << std::endl;
}

//FIXME : deal with uppercase/lowercase/unspecified NT in FASTA
int parse_fasta(const std::string& ref_fasta_fname, std::string& ref_fasta_content)
{
    std::ifstream f(ref_fasta_fname);
    std::string s;
    if (!std::getline(f, s) || s[0] != '>') {
        std::cerr << "ERROR: File " << ref_fasta_fname << " doesn't seem like a FASTA file." << std::endl;
        exit(EXIT_FAILURE);
    }
    while(std::getline(f, s)) {
        ref_fasta_content += s;
    }
        

}

int main(int argc, char* argv[]) {
  parameters_t p;
  parse_arguments(argc, argv, p);

  ifstream input(p.input_filename, ios::in|ios::binary|ios::ate);
  // Can add this to save a couple seconds off traversal - not really worth it.
  //vector<size_t> minus_positions;
  debruijn_graph<> dbg = debruijn_graph<>::load_from_packed_edges(input, "$ACGT"/*, &minus_positions*/);
  input.close();

  rrr_vector<63> colors;
  load_from_file(colors, p.color_filename);

  cerr << "k             : " << dbg.k << endl;
  cerr << "num_nodes()   : " << dbg.num_nodes() << endl;
  cerr << "num_edges()   : " << dbg.num_edges() << endl;
  cerr << "colors        : " << colors.size() / dbg.num_edges() << endl; 
  cerr << "Total size    : " << size_in_mega_bytes(dbg) << " MB" << endl;
  cerr << "Bits per edge : " << bits_per_element(dbg) << " Bits" << endl;
  cerr << "Color size    : " << size_in_mega_bytes(colors) << " MB" << endl;

  //dump_nodes(dbg, colors);
  //dump_edges(dbg, colors);
  uint64_t mask1 = (p.ref_color.length() > 0) ? atoi(p.ref_color.c_str()) : -1;
  uint64_t mask2 = (p.sample_mask.length() > 0) ? atoi(p.sample_mask.c_str()) : -1;
  std::string ref_fasta_content;
  std::cerr << "Loading reference FASTA file " << p.ref_fasta  << "...";
  parse_fasta(p.ref_fasta, ref_fasta_content);
  std::cerr << " got " << ref_fasta_content.size() << " nucleotides." << std::endl;
  find_divergent_paths(dbg, colors, mask1, mask2, ref_fasta_content);
}
