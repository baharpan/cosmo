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

string extension = ".contigs";


void parse_arguments(int argc, char **argv, parameters_t & params)
{
  TCLAP::CmdLine cmd("Cosmo Copyright (c) Alex Bowe (alexbowe.com) 2014", ' ', VERSION);
  TCLAP::UnlabeledValueArg<std::string> input_filename_arg("input",
            ".packed edge file (output from pack-edges).", true, "", "input_file", cmd);
   TCLAP::UnlabeledValueArg<std::string> color_filename_arg("color",
    ".color file (output from cosmo-pack.cpp).", true, "", "color_file", cmd);
  //  TCLAP::UnlabeledValueArg<std::string> positions_filename_arg("positions",
   // ".positions file (input by user)", true, "", "positions_file", cmd);
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
  params.color_filename  = color_filename_arg.getValue();
  params.output_prefix   = output_prefix_arg.getValue();
  // params.color_mask1     = color_mask1_arg.getValue();
  // params.color_mask2     = color_mask2_arg.getValue();
}


static char base[] = {'?','A','C','G','T'};

void readFile(ifstream &color)
{
    std::string word;
    char x ;
    word.clear();
    while (color>>word ){
        x = color.get();
        while ( x != ' ' ){
            word = word + x;
            x = color.get();
        }
    permutation.push_back(word);

    //cout<<word<<endl;
        }
color.close();
cerr << "num_payloads(): " << permutation.size()<<" must be equal to the number of edges"<< endl;
    }

std::string ReadNthLine(int N){
   std::ifstream positions("/s/fir/c/nobackup/baharpan/Ecoli/intermediate/cosmo-input-positions-nokmer");
   std::string s;
   for(int i = 1; i <= N; ++i)
       std::getline(positions,s);
   std::getline(positions,s);
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
vector<string> tokens;
vector<string> tokens_new;
vector<vector<string> >::iterator tokenss_iterator;
vector<vector<string> >::iterator tokenss_iterator_next;
vector<string>::iterator tokens_iterator;
 
void find_contigs(debruijn_graph<> dbg, ofstream &output)
{
//std::vector <vector<string>>tokenss(dbg.num_edges());
/*cerr<<"Constructing The Vector of position vectors"<<endl;
    for  (size_t pay=0; pay<permutation.size(); ++pay){
    string poss=ReadNthLine(pay);
    istringstream iss(poss);
    copy(istream_iterator<string>(iss),
    istream_iterator<string>(),
    back_inserter(tokens));
    tokenss.insert(tokenss.begin()+pay,tokens);}
cerr<<"End of Constructing"<<endl;*/
size_t  payload;
size_t sum=0;
size_t average=0;
  std::vector<int> visited (dbg.num_nodes());
  cout << "Starting to look for contigs\n";
  
 for (size_t i =1; i <=dbg.num_nodes() ; i++) {
   ssize_t start = i; // place to store start of branch kmer
   std::string start_label(dbg.node_label(start));
   payload=stoi(permutation[dbg._node_to_edge(i)]);	
   //cout << i << ":" << dbg.node_label(i)<<"and payload"<<payload<<endl;

   string poss=ReadNthLine(payload);
   istringstream iss(poss);
    /* copy(istream_iterator<string>(iss),
    istream_iterator<string>(),
    ostream_iterator<string>(cout, "\t"));*/
    
    	
     //cerr<<"outdegree is"<<dbg.outdegree(i)<<endl;
	  string contig;
    //tokenss_iterator = tokenss.begin()+payload;	
    //for(tokens_iterator = (*tokenss_iterator).begin();tokens_iterator!=(*tokenss_iterator).end();++tokens_iterator) {
     //cout<<"Positions are"<<*tokens_iterator<<" ";
    //	cerr<<endl;
	// }
     	if (visited[i]!=1){


         
	for (unsigned long x = 1; x < dbg.sigma + 1; x++) {
	  
                ssize_t edge = dbg.outgoing_edge(i, x);
        
                if (edge == -1)
                    continue;
		contig += base[x];
		
		//	cerr<< "first base is"<<base[x]<<endl;
		ssize_t pos = dbg._edge_to_node(edge);
		//	cerr<<"outdegree is"<<dbg.outdegree(pos)<<endl;
                while ( dbg.indegree(pos)&& dbg.outdegree(pos) == 1) {
		  
		   visited[pos] = 1;
                    ssize_t next_edge = 0;
		     int payyload = stoi(permutation[dbg._node_to_edge(pos)]);
		     
	//	     cout << pos << ":" << dbg.node_label(pos)<<"and payload"<<payyload<<endl;
			 string poss=ReadNthLine(payyload);
    istringstream iss(poss);
    /* copy(istream_iterator<string>(iss),
    istream_iterator<string>(),
    ostream_iterator<string>(cout, "\n"));*/
  /* tokenss_iterator_next = tokenss.begin()+payyload;
    for(tokens_iterator = (*tokenss_iterator_next).begin();tokens_iterator!=(*tokenss_iterator_next).end();++tokenss_iterator_next) {
           cout<<"Positions are"<<*tokens_iterator<<" ";}
	
   */
            for (unsigned long x2 = 1; x2 < dbg.sigma + 1; x2++) { // iterate through the alphabet
		       next_edge = dbg.outgoing_edge(pos, x2);
			
                        if (next_edge != -1 ){//&& std::find((*tokenss_iterator).begin(), (*tokenss_iterator).end(), position+1) != (*tokenss_iterator).end()
			  
			  //  float position= stof( *tokens_iterator);
			   //  stringstream ss;
			   // ss << position++;

			  //  string position_str = ss.str();
			  //  if(std::find( (*tokenss_iterator_next).begin(),  (*tokenss_iterator_next).end(),position_str ) != (*tokenss_iterator_next).end()){
			       
			  
                            contig+= base[x2];
			    //   std::cerr <<"base is"<< base[x2]<<endl;
                            break;
			    
			}
			 
			}
                    pos = dbg._edge_to_node(next_edge);}
		    //cerr <<"The"<<x<<"th contig:"<<contig[x]<<endl;}
	}
	//	cerr <<"contig:"<<dbg.node_label(start)<<contig<<"with length"<<contig.size()+dbg.node_label(start).size()<<endl;
	sum=sum+contig.size()+dbg.node_label(start).size();
    average= sum/i;
    
    output<<">"<<poss<<" "<<"length:"<<contig.size()+dbg.node_label(start).size()<<" "<<"average:"<<average<<endl;
	output<<dbg.node_label(start)<<contig<<endl;        }
 }
}

int main(int argc, char* argv[]) {
  parameters_t p;
  parse_arguments(argc, argv, p);
  const char * file_name = p.color_filename.c_str();
  string outfilename = (p.output_prefix == "")? p.input_filename : p.output_prefix;  
  ifstream color(file_name, ios::in);
  ifstream input(p.input_filename, ios::in|ios::binary|ios::ate);
//  ifstream positions (p.positions_filename, ios::in);
  ofstream output;
  output.open(outfilename + extension, ios::out);
  debruijn_graph<> dbg = debruijn_graph<>::load_from_packed_edges(input, "$ACGT"/*, &minus_positions*/);
  input.close();		
  cerr << "k             : " << dbg.k << endl;
  cerr << "num_nodes()   : " << dbg.num_nodes() << endl;
  cerr << "num_edges()   : " << dbg.num_edges() << endl;
  cerr << "Total size    : " << size_in_mega_bytes(dbg) << " MB" << endl;
  cerr << "Bits per edge : " << bits_per_element(dbg) << " Bits" << endl;
  cerr << "colors        : " << 74319898 << endl; 
  readFile(color);
  find_contigs(dbg,output);
}

