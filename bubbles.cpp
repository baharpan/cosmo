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
#include "kmer.hpp"
#include <sys/timeb.h>
#include <future>
#include <boost/algorithm/string.hpp>
#include <sys/timeb.h>
#include <utility>
#include <ctime>
#include "tclap/CmdLine.h"
#include <cstdio>
#include <cstdlib>
using namespace std;
using namespace sdsl;


static char base[] = {'?','A','C','G','T'};


struct parameters {
  size_t k = 0;
  //size_t n = 0;
  size_t b = 0;
  size_t u = 0;
  size_t c = 0;
};


parameters parse_arguments(int argc, char **argv) {
  parameters params;
  TCLAP::CmdLine cmd(" ");
  TCLAP::ValueArg<size_t> k_value("k", "kmer_length", "Length of edges (node is k-1)", true, 32 , "k value", cmd);
  //TCLAP::ValueArg<size_t> num_reads("n", "num_reads", "number of the reads", true, 0,"number of reads", cmd);
  TCLAP::ValueArg<size_t> len_br("b", "length_of_branches", "length of output sequences", true, 1000, "length of output sequences", cmd);
  TCLAP::ValueArg<size_t> num_br("u", "number_of_branches", "number of output sequences", true, 2, "number of output sequences", cmd);
  TCLAP::ValueArg<size_t> cov ("c", "coverage", "minimum coverage of the reads", true, 0, "minimum coverage", cmd);

  cmd.parse( argc, argv );
  params.k = k_value.getValue();
  //params.n = num_reads.getValue();
  params.b = len_br.getValue();
  params.u = num_br.getValue();
  params.c = cov.getValue();
  return params;
}



std::map<string,size_t> index_maker(){
  ifstream f ("pair" , ifstream::in);
  map<string,size_t> pairs;
  string s;
  size_t index = 0;
  while(std::getline(f,s)){
    std::locale loc;
    string ss="";
    for (std::string::size_type i = 0; i < s.length(); ++i)
      ss += std::toupper(s[i],loc);
      pairs.insert(std::pair<string,size_t>(ss,index++));
    }
    return pairs;
  }

  vector<pair<string,size_t>> repeats_maker(){
    ifstream f ("repeats" , ifstream::in);
    vector<pair<string,size_t>> cycles;
    string s;

    while(std::getline(f,s)){
      std::locale loc;
      string ss="";
      for (std::string::size_type i = 0; i < s.length(); ++i)
        ss += s[i];
        size_t len = ss.find(" ");
	      string kmer = ss.substr(0, len);
        string repeat = ss.substr(len+1,len+2);
	      cycles.push_back(make_pair(kmer,stoi(repeat)));

      }
      return (cycles);

    }

  void find_incomming_ougoing_degree_noposition(debruijn_graph<> &dbg){
  size_t  in2 = 0, in3 = 0, in4 = 0, out2 = 0, out3 = 0, out4 = 0;
  std::cerr<<"finding incomming and outgoing degrees nonposition\n";
  for (size_t i=0; i< dbg.num_nodes(); i++){

                  switch(dbg.indegree(i)){
                          case 2 : in2++; break;
                          case 3 : in3++; break;
                          case 4 : in4++; break;
                          }
                  switch(dbg.outdegree(i)){
                          case 2 : out2++; break;
                          case 3 : out3++; break;
                          case 4 : out4++; break;
                          }

  }
  cerr<<"Number of nodes with incoming degree 2: "<<in2<<endl;
  cerr<<"Number of nodes with incoming degree 3: "<<in3<<endl;
  cerr<<"Number of nodes with incoming degree 4: "<<in4<<endl;
  cerr<<"-------------------------------------------"<<endl;
  cerr<<"Number of nodes with outgoing degree 2: "<<out2<<endl;
  cerr<<"Number of nodes with outgoing degree 3: "<<out3<<endl;
  cerr<<"Number of nodes with outgoing degree 4: "<<out4<<endl;
}



  vector<unsigned long long int> Coherency (string kmer ,sd_vector<> &b, size_t number_of_subreads, std::map<string,size_t> &pairs){
    rank_support_sd<1>b_rank(&b);
    select_support_sd<1>b_sel(&b);
    vector<unsigned long long int> temp;
    map<string,size_t>::iterator find = pairs.find(kmer);
    if (find != pairs.end()) {
      size_t key1 = (pairs[kmer]) * number_of_subreads;
      size_t key2 = (pairs[kmer]+1) * number_of_subreads;
      size_t x = b_rank(key1);
      x++;
      size_t y = b_sel(x);
      while(y < key2){
        y = b_sel(x);
        if (y >= key2)  break;
        temp.push_back(y % number_of_subreads);
        x++;
        }
    }
  return temp;
}


bool Intersection (vector<unsigned long long int> &active_color , vector<unsigned long long int> &temp , size_t overlap) {
  if (active_color.empty() || temp.empty()) return false;
  size_t count = 0;
  for (auto it = active_color.begin(); it != active_color.end(); it++) {
    if (find(temp.begin(), temp.end(), *it) != temp.end()) {
      count++;
      if (count >= overlap) return true;
    }
  }

  return false;
}


string back(size_t k, debruijn_graph<> &dbg, size_t i, sd_vector<> &b, size_t number_of_subreads, std::map<string,size_t>& pairs, size_t max_branch_length, size_t overlap){
  bit_vector visited_cycle = bit_vector(dbg.num_nodes(), 0);
  size_t poss = dbg.backward(i);
  char base = dbg.node_label(i)[k-2];
  string kmer = dbg.node_label(poss) += base;
  vector<unsigned long long int> active_color = Coherency(kmer, b, number_of_subreads, pairs);
  string prefix="";
  prefix.insert(0, 1, dbg.node_label(dbg.backward(i))[0]);

  while(dbg.indegree(poss) != 0 && !visited_cycle[poss] && prefix.size() < max_branch_length){
    visited_cycle[poss] = 1;
    if (dbg.node_label(dbg.backward(poss))[0]=='$') break;
    char base2 = dbg.node_label(poss)[k-2];
    string kmer2 = dbg.node_label(dbg.backward(poss)) += base2;
    vector <unsigned long long int> temp = Coherency(kmer2, b, number_of_subreads, pairs);

        if (Intersection(active_color,temp, overlap)){
          prefix.insert(0, 1, dbg.node_label(dbg.backward(poss))[0]);
        }

        else return (prefix);

     poss = dbg.backward(poss);
     active_color = temp;
   }
return (prefix);
}



map <size_t, pair<size_t,vector<string>>> FoundBulge;
size_t bulge_count = 0;


void SNPCall (size_t k, size_t i,  debruijn_graph<>& dbg, sd_vector<>& b, size_t number_of_subreads,
  map<string,size_t>& pairs, ofstream &variations, bit_vector &visited, size_t max_branch_length, size_t max_branch_num, size_t overlap, vector<pair<string,size_t>> cycles){

    bit_vector visited_bulge = bit_vector(dbg.num_nodes(), 0);
    size_t end_node;
    visited[i] = 1;
    visited_bulge[i] = 1;
    size_t branch_num = -1;
    vector<string> branch_labels;
    vector<char> var_base;




for (size_t x = 1; x < dbg.sigma + 1; x++) {

  ssize_t edge = dbg.outgoing_edge(i, x);
  if (edge == -1) continue;
  branch_num++;
  bit_vector visited_cycle = bit_vector(dbg.num_nodes(), 0);
  visited_cycle[i] = 1;
  string branch = "";
  branch_labels.push_back(branch);
  var_base.push_back(base[x]);
  string kmer = dbg.node_label(i)+ base[x];
  branch_labels[branch_num]+= base[x];
  vector<unsigned long long int> active_color = Coherency(kmer,b,number_of_subreads,pairs);
  if (active_color.empty()) continue;
  string left_flank_embed;
  ssize_t pos = dbg._edge_to_node(edge);


  bool bulge = false;
  size_t c = 0;
  bool firstNodeCycle = true;
  vector<size_t> end_cycles;
  bool cycleFine = false;
  size_t check;
  while ( c < max_branch_num && dbg.outdegree(pos) != 0 && branch_labels[branch_num].size() < max_branch_length ) {
      if (visited_bulge[pos]){
      bulge = true;
      end_node = pos;
      }

    if (dbg.outdegree(pos) > 1 ) {
      map<size_t, pair<size_t,vector<string>>>::iterator it = FoundBulge.find(pos);

      if (visited[pos] && it != FoundBulge.end()) {
        pos = it -> second.first ;
        left_flank_embed = branch_labels[branch_num];

        for (vector<string>::iterator v = (it->second.second).begin(); v!= (it->second.second).end(); v++){


                branch_labels[branch_num] += *v;
                if (v == (it->second.second).end()-1) break;
                branch_labels.push_back(branch);
                branch_num++;
                branch_labels[branch_num] = left_flank_embed;



            }
        if (visited_bulge[pos]){
          bulge = true;
          end_node = pos;
        }
        visited_bulge[pos] = 1;

        if (dbg.outdegree(pos) == 0) break;
      }

      if (!visited[pos]) {  SNPCall (k, pos, dbg, b, number_of_subreads, pairs, variations, visited, max_branch_length, max_branch_num, overlap, cycles);}
    }
    vector<unsigned long long int> temp;
    size_t next_edge = 0;
    string kmer2;
    for (size_t x2 = 1; x2 < dbg.sigma + 1; x2++) {
      next_edge = dbg.outgoing_edge(pos, x2);
      if (next_edge != -1) {
          kmer2 = dbg.node_label(pos)+base[x2];
          temp = Coherency(kmer2, b, number_of_subreads, pairs);
          if (temp.empty())  continue;

          if (visited_cycle[pos]){
            if (firstNodeCycle){

              for (size_t i = 0; i < cycles.size(); i++){
                if (cycles[i].first == dbg.node_label(pos)){
                  end_cycles.push_back(cycles[i].second);
                }
              }
            if (end_cycles.empty() ) {c+= 1; break;}
                firstNodeCycle = false; }

            for (size_t i = 0; i < end_cycles.size(); i++){

              if ( find(temp.begin(), temp.end(), end_cycles[i]) != temp.end() ) {cycleFine = true; break;}}
              if (!cycleFine){
                c+= 1;  end_cycles.clear(); firstNodeCycle = true; break;}
            }
            cycleFine = false;
          visited_cycle[pos]=1;
          map<string,size_t>::iterator finddd = pairs.find(kmer2);
          if (Intersection (active_color, temp, overlap)){
              visited_bulge[pos] = 1;
              branch_labels[branch_num] += base[x2];
	      break;
            }
          else { c += 1;  break;}

        }
      }

    pos = dbg._edge_to_node(next_edge);
    active_color = temp;
  }
    if (bulge == true && branch_labels.size() > 1){
      FoundBulge.insert(std::make_pair(i,make_pair(end_node,branch_labels)));
      bulge_count++;


      cout<<"\n>Bulge number "<<bulge_count<<" From node number "<<i<<" variation of ";
      for (size_t i = 0; i < var_base.size(); ++i) cout<<var_base[i];
      cout << "\n\n";
      cout << "Start node: " << dbg.node_label(i) << "\n";
      string prefix = back(k, dbg, i, b, number_of_subreads, pairs,max_branch_length, overlap);
      size_t index = prefix.size() + k;
      char var;

      for (size_t j = 0; j < branch_labels.size(); ++j){
          if (branch_labels[j].size() < k) continue;
          if (max_branch_num == 2 && branch_labels[j][0] == var) continue;
            cout<< "Branch: "<<branch_labels[j]<<endl;
            variations<<">"<<bulge_count<<"."<<j<<"."<<index<<".";
            for (size_t i = 0; i < var_base.size(); ++i) variations<<var_base[i];
              variations<<"\n"<< prefix<<dbg.node_label(i)<<branch_labels[j] << "\n";
              var = branch_labels[j][0];
      }

    }

  }
}


void SNP (size_t k, debruijn_graph<> &dbg, sd_vector<>& b, size_t number_of_subreads, map<string,size_t> &pairs, size_t max_branch_length,
  size_t max_branch_num, size_t overlap){
  vector<pair<string,size_t>> cycles = repeats_maker();
  cout<<"Beginning of SNP calling ..."<<endl;
  cout<<"========================="<<endl;
  bit_vector visited = bit_vector(dbg.num_nodes(), 0);
  ofstream variations;
  variations.open ("variations");
  for (size_t i = 0; i < dbg.num_nodes(); i++) {
     if ( dbg.outdegree(i) > 1 && !visited[i]){
         if (i % 10000 == 0) cout<<i*100/dbg.num_nodes()<<" % of graph is processed"<<endl;
         SNPCall(k, i, dbg, b, number_of_subreads, pairs, variations , visited, max_branch_length, max_branch_num, overlap, cycles);
       }
     }
    cout<<"100% of graph is processed"<<endl;
    cout<<"========================="<<endl;
    cout<<"End of SNP calling.\n";
   }







int main(int argc, char *argv[]){

  auto params = parse_arguments(argc, argv);
  size_t k = params.k;
  //size_t number_of_subreads = params.n;
  size_t max_branch_length = params.b;
  size_t max_branch_num = params.u;
  size_t overlap = params.c;


  ifstream f;
  f.open("subreadInfo");
  string s;
  getline(f,s);
  size_t number_of_subreads = stoi(s);
  ifstream input("filtered_kmc2_list.packed", ios::in|ios::binary|ios::ate);
  debruijn_graph<> dbg = debruijn_graph<>::load_from_packed_edges(input, "$ACGT"/*, &minus_positions*/);
  input.close();
  cerr << "num_nodes()   : " << dbg.num_nodes() << endl;
  cerr << "num_edges()   : " << dbg.num_edges() << endl;
  std::map<string,size_t> pairs = index_maker();
  sdsl::sd_vector<> b;
  sdsl::load_from_file(b,"output_matrix");
  cerr << "Total size of matrix: "<<size_in_mega_bytes(b) <<" Mb"<<endl;
  cerr<<"number of reads: "<<number_of_subreads<<endl;



  //find_incomming_ougoing_degree_noposition(dbg);
  SNP(k, dbg,b,number_of_subreads,pairs, max_branch_length, max_branch_num, overlap);
}
