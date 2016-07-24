#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <libgen.h> // basename
#include <stdio.h>
#include <ctype.h>


#include "tclap/CmdLine.h"

#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>

#include "io.hpp"
#include "debruijn_graph.hpp"
#include "algorithm.hpp"
#include "cosmo-color.hpp"
#include <boost/algorithm/string.hpp>    

using namespace std;
using namespace sdsl;
std::vector<string> permutation;
#include <sys/timeb.h>
void readFile(  )
{
    ifstream file;
    file.open ("/s/fir/c/nobackup/baharpan/git/cosmo/Ecoli-K12/permutation.payload");
    std::string word;
    char x ;
    word.clear();

    while (file>>word )
    {
        x = file.get();

        while ( x != ' ' )
        {
            word = word + x;
            x = file.get();
        }
	
	permutation.push_back(word);
	//cout<<word<<endl;
    }
}


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
  // TCLAP::UnlabeledValueArg<std::string> color_filename_arg("color",
  //  ".color file (output from pack-edges).", true, "", "color_file", cmd);
  string output_short_form = "output_prefix";
  TCLAP::ValueArg<std::string> output_prefix_arg("o", "output_prefix",
            "Output prefix. Graph will be written to [" + output_short_form + "]" + extension + ". " +
            "Default prefix: basename(input_file).", false, "", output_short_form, cmd);
  //  string color_mask1 = "color_mask1";
  // TCLAP::ValueArg<std::string> color_mask1_arg("a", "color_mask1",
  //	    "Color mask 1, color1 [" + color_mask1 + "]", false, "", color_mask1, cmd);
  //  string color_mask2 = "color_mask2";
  //  TCLAP::ValueArg<std::string> color_mask2_arg("b", "color_mask2",
  //	    "Color mask 2, color2 [" + color_mask2 + "]", false, "", color_mask2, cmd);
  cmd.parse( argc, argv );

  params.input_filename  = input_filename_arg.getValue();
  // params.color_filename  = color_filename_arg.getValue();
  params.output_prefix   = output_prefix_arg.getValue();
  // params.color_mask1     = color_mask1_arg.getValue();
  // params.color_mask2     = color_mask2_arg.getValue();
}

static char base[] = {'?','A','C','G','T'};
std::string ReadNthLine(int N)
{
   std::ifstream in("/s/fir/c/nobackup/baharpan/Ecoli/intermediate/cosmo-input-positions-nokmer");

   std::string s;
   
   for(int i = 1; i <= N; ++i)
       std::getline(in, s);

   std::getline(in,s);
   return s; 
   }

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



/*void dump_edges(debruijn_graph<> dbg) {
  for (size_t i = 0; i < dbg.num_nodes(); i++) {
    cout << i << ":" << dbg.node_label(i)<< "\n";
  }
  }*/


/*void dump_edges(debruijn_graph<> dbg, uint64_t * colors) {
  for (size_t i = 0; i < dbg.num_edges(); i++) {
    cout << i << "e:" << dbg.edge_label(i) << colors[i] << "\n";
  }
}*/

/*const char *const starts[] = {"GCCATACTGCGTCATGTCGCCCTGACGCGC","GCAGGTTCGAATCCTGCACGACCCACCAAT","GCTTAACCTCACAACCCGAAGATGTTTCTT","AAAACCCGCCGAAGCGGGTTTTTACGTAAA","AATCCTGCACGACCCACCAGTTTTAACATC","AGAGTTCCCCGCGCCAGCGGGGATAAACCG","GAATACGTGCGCAACAACCGTCTTCCGGAG"};*/
                std::vector<string> tokens;
std::vector<string> tokens_new;
		vector<vector<string> >::iterator tokenss_iterator;
                vector<vector<string> >::iterator tokenss_iterator_next;
		vector<string>::iterator tokens_iterator;
 
void find_contigs(debruijn_graph<> dbg)
{
  
  readFile(  );
 
  std::vector <vector<string>>tokenss(dbg.num_edges());
  /* for (size_t pay=50; pay<=100permutation.size(); ++pay){
    string poss=ReadNthLine(pay);
    istringstream iss(poss);
    copy(istream_iterator<string>(iss),
    istream_iterator<string>(),
    back_inserter(tokens));
    tokenss.insert(tokenss.begin()+pay,tokens);
    
    }*/
  size_t payload;
  std::vector<int> visited (dbg.num_nodes());
  cout << "Starting to look for contigs\n";
   ofstream myfile;
  
  myfile.open ("/s/fir/c/nobackup/baharpan/git/cosmo/Ecoli-K12/contigs.txt");
 for (size_t i =0; i <dbg.num_nodes() ; i++) {
   ssize_t start = i; // place to store start of branch kmer
        std::string start_label(dbg.node_label(start));
	payload=stoi(permutation[dbg._node_to_edge(i)]);	
	//	cout << i << ":" << dbg.node_label(i)<<"and payload"<<payload<<endl;

	 string poss=ReadNthLine(payload);
    istringstream iss(poss);
    /*  copy(istream_iterator<string>(iss),
    istream_iterator<string>(),
    ostream_iterator<string>(cout, "\t"));*/
    
    
   
    
    //back_inserter(tokens));
    // tokenss.insert(tokenss.begin()+payload,tokens);
    // tokenss_iterator = tokenss.begin()+payload;	
    //for(tokens_iterator = (*tokenss_iterator).begin();tokens_iterator!=(*tokenss_iterator).end();++tokens_iterator) {
    // cout<<"Positions are"<<*tokens_iterator<<" ";
    //	cerr<<endl;
	// }
    //	cout<<"outdegree is"<<dbg.outdegree(i)<<endl;
		if (visited[i]!=1){

	string contig[dbg.sigma];
        
	for (unsigned long x = 1; x < dbg.sigma + 1; x++) {
                ssize_t edge = dbg.outgoing_edge(i, x);
                if (edge == -1)
                    continue;
		contig[x] += base[x];
		
		//	cerr<< "first base is"<<base[x]<<endl;
		ssize_t pos = dbg._edge_to_node(edge);
                while (dbg.indegree(pos) == 1 && dbg.outdegree(pos) == 1) {
		   visited[pos] = 1;
                    ssize_t next_edge = 0;
		     int payyload = stoi(permutation[dbg._node_to_edge(pos)]);
		     
		     //	cout << pos << ":" << dbg.node_label(pos)<<"and payload"<<payyload<<endl;
			 string poss=ReadNthLine(payyload);
    istringstream iss(poss);
    /* copy(istream_iterator<string>(iss),
    istream_iterator<string>(),
    ostream_iterator<string>(cout, "\n"));*/
    //back_inserter(tokens_new));
    // tokenss.insert(tokenss.begin()+payyload,tokens_new);
	
		//std::cerr <<"base is"<< base[x] << std::
                    for (unsigned long x2 = 1; x2 < dbg.sigma + 1; x2++) { // iterate through the alphabet
		       next_edge = dbg.outgoing_edge(pos, x2);
		      //	int next_payload= stoi (permutation [dbg._node_to_edge (dbg._edge_to_node(next_edge))]);
		      //	cout<<"next payload is"<<next_payload;
			
                        if (next_edge != -1){ //&& std::find((*tokenss_iterator).begin(), (*tokenss_iterator).end(), position+1) != (*tokenss_iterator).end()){
			  // tokenss_iterator = tokenss.begin()+payyload;
		        
			  //	for(tokens_iterator = (*tokenss_iterator).begin();tokens_iterator!=(*tokenss_iterator).end();++tokens_iterator) {
			  // cout<<"Positions are"<<*tokens_iterator<<" ";
			  //  tokenss_iterator_next = tokenss.begin()+next_payload;
			  //  float position= stof( *tokens_iterator);
			   //  stringstream ss;
			   // ss << position++;

			  //  string position_str = ss.str();
			  //  if(std::find( (*tokenss_iterator_next).begin(),  (*tokenss_iterator_next).end(),position_str ) != (*tokenss_iterator_next).end()){
			       
			  
                            contig[x]+= base[x2];
			    //cerr<<base[x2]<<" ";
                            break;
			    
			}
			 
			}
                    pos = dbg._edge_to_node(next_edge);}
		if (payload!= -1){
	myfile<<">"<<poss<<endl;
	myfile<<contig [x];
	myfile<<endl;
	  
	//	cerr <<"The"<<x-1<<"th contig:"<<contig[x]<<endl;
		}		  
	}
		}
 }
 	myfile.close();
 }

int main(int argc, char* argv[]) {
  parameters_t p;
  parse_arguments(argc, argv, p);
  //std::map<string,int> pairs;
  //int payload;
  ifstream input(p.input_filename, ios::in|ios::binary|ios::ate);
  // Can add this to save a couple seconds off traversal - not really worth it.
  //vector<size_t> minus_positions;
  debruijn_graph<> dbg = debruijn_graph<>::load_from_packed_edges(input, "$ACGT"/*, &minus_positions*/);
  input.close();		
 
	 cerr << "k             : " << dbg.k << endl;
  cerr << "num_nodes()   : " << dbg.num_nodes() << endl;
  cerr << "num_edges()   : " << dbg.num_edges() << endl;
  cerr << "colors        : " << 74319898 << endl; 
  cerr << "Total size    : " << size_in_mega_bytes(dbg) << " MB" << endl;
  cerr << "Bits per edge : " << bits_per_element(dbg) << " Bits" << endl;
  // cerr << "Color size    : " << size_in_mega_bytes(colors) << " MB" << endl;
  cerr<<"Permutation size : "<<permutation.size()<<endl;
  //dump_nodes(dbg, colors);
 // dump_edges(dbg);
  // uint64_t mask1 = (p.color_mask1.length() > 0) ? atoi(p.color_mask1.c_str()) : -1;
   // uint64_t mask2 = (p.color_mask2.length() > 0) ? atoi(p.color_mask2.c_str()) : -1;
	 find_contigs(dbg);
}

