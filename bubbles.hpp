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
#include <sys/timeb.h>

#include <utility>
#include <ctime>

// TCLAP
#include "tclap/CmdLine.h"

 #include <cstdio>

 #include <cstdlib>
#include "kmer.hpp"
static char base[] = {'?','A','C','G','T'};

class resistome {



public:





  vector<size_t> file_size_function(){
  ifstream f ("subreads_info");
    vector<size_t> f_size; //in case we have more than one sample
    string l;
    while (std::getline(f, l)) {
          f_size.push_back(stol(l)/4);
      }
      f.close();
      return f_size;
  };

  std::map<string,size_t> index_maker(){
    ifstream f ;
     std::map<string,size_t> pairs;
    string s;
    size_t index=0;
    while(std::getline(f,s)){
      std::locale loc;
      string ss="";
      for (std::string::size_type i=0; i<s.length(); ++i)
        ss += std::toupper(s[i],loc);
        pairs.insert(std::pair<string,size_t>(ss,index++));
      }
      std::map<string,size_t>::iterator it = pairs.begin();
      return pairs;
    };

    void find_incomming_ougoing_degree_noposition(debruijn_graph<> dbg){
    size_t  in2=0, in3=0, in4=0, out2=0, out3=0, out4=0;
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
  };


    std::vector<unsigned long long int> Coherency (string &kmer ,sd_vector<> &b, size_t number_of_subreads, std::map<string,size_t> &pairs){

      rank_support_sd<1>  b_rank(&b);
      select_support_sd<1> b_sel (&b);

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
  };


  bool Intersection (vector<unsigned long long int> &active_color , vector<unsigned long long int>&temp) {
    bool status = false;
    size_t in = 0;
    if (!active_color.empty() && !temp.empty()){
      while (in != temp.size()){
         if (std::find(active_color.begin(), active_color.end(), temp[in])!= active_color.end() ){
                status = true;
                 break;
               }
                in++;
              }
      }
    else status = true;
  return status;
};


  string back(debruijn_graph<> dbg, size_t i, sd_vector<> &b, size_t &number_of_subreads, std::map<string,size_t>& pairs){

    size_t poss = dbg.backward(i);
    char base = dbg.node_label(i)[30];
    string kmer = dbg.node_label(poss) += base;
    vector<unsigned long long int> active_color = this->Coherency(kmer, b, number_of_subreads, pairs);
    string prefix="";
    prefix.insert(0, 1, dbg.node_label(dbg.backward(i))[0]);

    while(dbg.indegree(poss) != 0 ){
      if (dbg.node_label(dbg.backward(poss))[0]=='$') break;
      char base2 = dbg.node_label(poss)[30];
      string kmer2 = dbg.node_label(dbg.backward(poss)) += base2;
      bool coherent = false;
      vector <unsigned long long int> temp = this->Coherency(kmer2, b, number_of_subreads, pairs);

          if (this->Intersection(active_color,temp)){

            coherent = true;
            prefix.insert(0, 1, dbg.node_label(dbg.backward(poss))[0]);

          }




       if (coherent == false) return (prefix);

       poss = dbg.backward(poss);
       active_color = temp;
     }
  return (prefix);
};








  map <size_t, pair<size_t,vector<string>>> FoundBulge;
  size_t bulge_count = 0;


  void SNPCall (size_t i,  debruijn_graph<>& dbg, sd_vector<>& b, size_t &number_of_subreads, std::map<string,size_t>& pairs, ofstream &variations, bit_vector &visited){
  bit_vector visited_cycle = bit_vector(dbg.num_nodes(), 0);
  bit_vector visited_bulge = bit_vector(dbg.num_nodes(), 0);
  size_t end_node;
  visited[i] = 1;
  visited_bulge[i] = 1;
  size_t branch_num = -1;
  std::vector<std::string> branch_labels;


  for (unsigned long x = 1; x < dbg.sigma + 1; x++) {

    ssize_t edge = dbg.outgoing_edge(i, x);
    if (edge == -1) continue;
    branch_num ++;
    string branch = "";
    branch_labels.push_back(branch);
    string kmer = dbg.node_label(i)+ base[x];
    branch_labels[branch_num]+= base[x];
    vector<unsigned long long int> active_color = this->Coherency(kmer,b,number_of_subreads,pairs);

    ssize_t pos = dbg._edge_to_node(edge);


    bool bulge = false;
    while ( dbg.outdegree(pos) != 0 && branch_labels[branch_num].size() < 1000 && !visited_cycle[pos]) {
        if (visited_bulge[pos]){
        bulge = true;
        end_node = pos;
        }


      if (dbg.outdegree(pos) > 1 ) {
        map<size_t, pair<size_t,vector<string>>>::iterator it = FoundBulge.find(pos);

        if (visited[pos] && it != FoundBulge.end()) {
          pos = it -> second.first ;
          string left_flank_embed = branch_labels[branch_num];
          size_t num_copy = 0;
          for (vector<string>::iterator v = (it->second.second).begin(); v!= (it->second.second).end(); v++){
              if (num_copy <100){
                  num_copy ++;
                  branch_labels[branch_num] += *v;
                  if (v == (it->second.second).end()-1) break;
                  branch_labels.push_back(branch);
                  branch_num++;
                  branch_labels[branch_num] = left_flank_embed;
                }
              }
          visited_cycle[pos] = 1;
          if (visited_bulge[pos]){
            bulge = true;
            end_node = pos;
          }
          visited_bulge[pos] =1;
          if (dbg.outdegree(pos) == 0) break;
        }

        if (!visited[pos]) {  SNPCall (pos, dbg, b, number_of_subreads, pairs, variations, visited);}
      }


      vector<unsigned long long int> temp;
      size_t next_edge = 0;
      string kmer2;
      for (unsigned long x2 = 1; x2 < dbg.sigma + 1; x2++) {
        next_edge = dbg.outgoing_edge(pos, x2);
        if (next_edge != -1) {
            kmer2 = dbg.node_label(pos)+base[x2];
            temp = this->Coherency(kmer2, b, number_of_subreads, pairs);
            map<string,size_t>::iterator finddd = pairs.find(kmer2);
            if (this->Intersection (active_color, temp)){
                visited_bulge[pos] = 1;
                branch_labels[branch_num] += base[x2];
                break;
              }
            else  break;
          }
        }  

      pos = dbg._edge_to_node(next_edge);
      active_color = temp;
    }
      if (bulge == true){
        FoundBulge.insert(std::make_pair(i,make_pair(end_node,branch_labels)));
        bulge_count++;
        cout<<"\n>Bulge number "<<bulge_count<<" From node number "<<i<<endl;
        cout << "Start node: " << dbg.node_label(i) << "\n";
        string prefix = back(dbg, i, b, number_of_subreads, pairs);
        size_t index = prefix.size() + 32;
        for (size_t j = 0; j < branch_labels.size(); ++j){
          cout<< "Branch: "<<branch_labels[j]<<endl;
          variations<<">"<<bulge_count<<"."<<j<<"."<<index<<".";
          variations<<"\n"<< prefix<<dbg.node_label(i)<<branch_labels[j] << "\n";
        }
      }
    }
  }


  void SNP (debruijn_graph<> &dbg, sd_vector<>& b, size_t &number_of_subreads, std::map<string,size_t> &pairs){
    bit_vector visited = bit_vector(dbg.num_nodes(), 0);
    ofstream variations;
    variations.open ("variations");
    for (size_t i = 0; i < dbg.num_nodes(); i++) {
  	   if ( dbg.outdegree(i) > 1 && !visited[i] ){
  		     SNPCall(i, dbg, b, number_of_subreads, pairs, variations , visited);
         }
       }
     }


};
