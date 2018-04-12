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


string push(vector<string> tr_branches , size_t i){

string a = tr_branches[tr_branches.size()-1];
a.erase (a.begin());
std::size_t found = tr_branches[i].find(a);
if (found!=std::string::npos){
	a=tr_branches[i].substr(found+a.length());}
else a="";

return (a);
}


string back(debruijn_graph<> dbg,size_t i , sd_vector<> &b, size_t &number_of_subreads, std::map<string,size_t>& pairs){
rank_support_sd<1>  b_rank(&b);
select_support_sd<1> b_sel (&b);

bit_vector visited_cycle = bit_vector(dbg.num_nodes(), 0);
visited_cycle[i] = 1;

vector<unsigned long long int> active_color;
bit_vector A = bit_vector(number_of_subreads , 0);
size_t poss = dbg.backward(i);
char base = dbg.node_label(i)[30];
string kmer = dbg.node_label(poss) += base;
size_t x = b_rank((pairs[kmer])*number_of_subreads);
size_t y = 0;
x++;
y = b_sel(x);
while( y < ((pairs[kmer] + 1) * number_of_subreads)){
        y = b_sel(x);
        if (y >= ((pairs[kmer] + 1)*number_of_subreads)) break;
        active_color.push_back(y % number_of_subreads);
        x++;
      }

string prefix="";
prefix.insert(0,1,dbg.node_label(dbg.backward(i))[0]);
bool is = false;


while(dbg.indegree(poss) != 0 && is == false && !visited_cycle[poss]){
  visited_cycle[poss] = 1;
  bit_vector B = bit_vector(number_of_subreads , 0);
  vector <unsigned long long int> temp;
  if (dbg.node_label(dbg.backward(poss))[0]=='$') break;
  char base2 = dbg.node_label(poss)[30];
  string kmer2 = dbg.node_label(dbg.backward(poss))+=base2;
  bool here = false;

  size_t xx = b_rank((pairs[kmer2])*number_of_subreads);
  size_t yy = 0;
  xx++;
  yy = b_sel(xx);
  while(yy < ((pairs[kmer2]+1)*number_of_subreads)){
      yy = b_sel(xx);
      if (yy >= ((pairs[kmer2]+1)*number_of_subreads)) break;
        temp.push_back(yy%number_of_subreads);
        xx++;}


  size_t in=0;
  while (in != temp.size()){
        if (std::find(active_color.begin(), active_color.end(), temp[in])!= active_color.end()){
          here = true;
          prefix.insert(0,1,dbg.node_label(dbg.backward(poss))[0]);
          break;
        }
	      else in++;
}


if(here == false) is = true;
poss = dbg.backward(poss);
active_color = temp;
}


return (prefix);

}


void cycle_remove (std::map<long,vector<std::pair<long,long>>> & cymap , ifstream file){

string s;
std::string delimiter = " ";
string delimiter_2 = ","  ;
std::string index;
string read;
string mult;
vector<std::pair<long,long>> repeats;
while (std::getline(file, s)){
        unsigned long long pos = 0;
	      while ((pos = s.find(delimiter)) != std::string::npos) {

          index = s.substr(0, pos);
          s.erase(0, pos + delimiter.length());
          unsigned long long pos2 = 0;
          while ((pos2 = s.find(delimiter_2)) != std::string::npos) {

              read = s.substr(0,pos2);
              s.erase(0, pos2 + delimiter.length());
              mult = s;

              break;
            }
          }
map<long,vector<std::pair<long,long>>>::iterator it = cymap.find(stol(index));
if(it == cymap.end()){
  repeats.clear();
  repeats.push_back(std::make_pair(stol(read),stol(mult)));
  cymap.insert(std::make_pair(stol(index),repeats));
}
else {
repeats.push_back(std::make_pair(stol(read),stol(mult)));
it->second = repeats;
}

}
}



int substring(string& branch, string& repeat){
int count = 0;
    size_t nPos = branch.find(repeat, 0); // fist occurrence
    while(nPos != string::npos)
    {
        count++;
        nPos = branch.find(repeat, nPos+1);
    }

    return count;
}



bool Intersection (vector<long unsigned int> &active_color , vector<long unsigned int>&temp) {
  size_t in = 0;
  bool status = false;
  if (!active_color.empty() && !temp.empty()){ //dummies are allowed
    while (in != temp.size()){
       if (std::find(active_color.begin(), active_color.end(), temp[in])!=active_color.end() ){
              status = true;
               break;
             }

              else in++;
          }
    }
  else status = true; //dummies are allowed


return status;
}



std::vector<long unsigned int> Coherency (string &kmer2 , sd_vector<>& b ,size_t &number_of_subreads, std::map<string,size_t> &pairs){

        rank_support_sd<1>  b_rank(&b);
        select_support_sd<1> b_sel (&b);
        vector<long unsigned int> temp;
        map<string,size_t>::iterator finddd = pairs.find(kmer2);
        if (finddd != pairs.end()) {
                bit_vector B = bit_vector(number_of_subreads, 0);
                size_t key11 = (pairs[kmer2])*number_of_subreads;
	              size_t key12 = (pairs[kmer2]+1)*number_of_subreads;
                size_t xx = b_rank(key11);
                size_t yy = 0;
                xx++;
                yy=b_sel(xx);
                while(yy < key12){
                    yy = b_sel(xx);
                    if (yy >= key12)  break;
                    temp.push_back(yy%number_of_subreads);
                    xx++;

                  }

            }
return temp;}



bit_vector visited = bit_vector(400000000, 0);
map <size_t, pair<size_t,vector<string>>> FoundBulge;
size_t bulge_count = 0;
std::ofstream variations ("variations");

void SNPCall (size_t i, map<long,vector<std::pair<long,long>>> &cymap,  debruijn_graph<>& dbg, sd_vector<>& b, size_t &number_of_subreads, std::map<string,size_t>& pairs){
bit_vector visited_bulge = bit_vector(dbg.num_nodes(), 0);
size_t end_node;
visited[i] = 1;
visited_bulge[i] = 1;
size_t branch_num = -1;
std::vector<std::string> branch_labels;

for (unsigned long x = 1; x < dbg.sigma + 1; x++) {
	ssize_t edge = dbg.outgoing_edge(i, x);
	vector<string> pattern;
  if (edge == -1) continue;
  bit_vector visited_cycle = bit_vector(dbg.num_nodes(), 0);
  visited_cycle [i] = 1;
	branch_num ++;
	string branch = "";
	branch_labels.push_back(branch);
	string kmer = dbg.node_label(i)+ base[x];
	branch_labels[branch_num]+= base[x];
	vector<long unsigned int> active_color = Coherency(kmer,b,number_of_subreads,pairs);
	ssize_t pos = dbg._edge_to_node(edge);
	bool stop = false;
	bool bulge = false;
	while ( dbg.outdegree(pos) !=0 && stop == false && branch_labels[branch_num].size() < 4000 && !visited_cycle[pos]) {
    if (visited_cycle[pos]) {
                int max = 0;
                string check;
                string bra = dbg.node_label(i) + branch_labels[branch_num];
                for (unsigned i = bra.length()-32; i < bra.length(); ++i) {
                        check += bra.at(i);
                 }

                if  (pattern.empty() || cymap.find(pairs[pattern[0]]) != cymap.end() ) {
                        if (pattern.empty()) pattern.push_back(check);

                        for (size_t m = 0; m < cymap[pairs[pattern[0]]].size(); m++){
                                        if (cymap[pairs[pattern[0]]][m].second > max) max = cymap[pairs[pattern[0]]][m].second ;
                                }
                        if (max==0) { stop = true; break;}
                        if (substring(bra,pattern[0]) > max-1 || branch_labels[branch_num].length()>4000 ){
                                stop=true; break;  }
                       }
                else   stop=true; break;
              }

    visited_cycle[pos] = 1;




    if (visited_bulge[pos]){
			bulge = true;
      end_node = pos;
      //FoundBulge.insert(std::make_pair(i,pos));

					 }


		if (dbg.outdegree(pos) > 1 ) {
      map<size_t, pair<size_t,vector<string>>>::iterator it = FoundBulge.find(pos);
      if (visited[pos] && it != FoundBulge.end()) {
        pos = it -> second.first ;
        string left_flank_embed = branch_labels[branch_num];
        size_t num_copy = 0; //Max number of branches
        for (vector<string>::iterator v = (it->second.second).begin(); v!= (it->second.second).end(); v++){
            if (num_copy <100){
                num_copy ++;
                branch_labels[branch_num] += *v;

                if (v== (it->second.second).end()-1) break;
                branch_labels.push_back(branch);
                branch_num++;
                branch_labels[branch_num] = left_flank_embed;
              }
            }

        if (visited_bulge[pos]){
    			bulge = true;
          end_node = pos;

    			}
          visited_bulge[pos] =1;


        if (dbg.outdegree(pos) == 0) break;



        }



      if (!visited[pos]) SNPCall (pos,cymap, dbg,b,number_of_subreads,pairs);
    }

		vector<long unsigned int> temp;
		size_t next_edge = 0;
		string kmer2;
		for (unsigned long x2 = 1; x2 < dbg.sigma + 1; x2++) {
        next_edge = dbg.outgoing_edge(pos, x2);
        if (next_edge != -1) {
		       kmer2 = dbg.node_label(pos)+base[x2];
	         temp = Coherency(kmer2,b,number_of_subreads,pairs);
		       //map<string,size_t>::iterator finddd = pairs.find(kmer2);
          //if (finddd!=pairs.end())
		       if (Intersection (active_color, temp)){
		           visited_bulge[pos]=1;
               branch_labels[branch_num] += base[x2];  break;
             }
          else{  stop = true; break; }
        }
      }
		pos = dbg._edge_to_node(next_edge);
		active_color = temp;}
		if (bulge == true){
      FoundBulge.insert(std::make_pair(i,make_pair(end_node,branch_labels)));
      bulge_count++;
		  cout<<"\n>Bulge number "<<bulge_count<<" From node number "<<i<<endl;
      cout << "Start node: " << dbg.node_label(i) << "\n";
      string prefix = back(dbg,i,b,number_of_subreads, pairs);
      size_t index = prefix.size()+32;
      for (size_t j = 0; j < branch_labels.size(); ++j){
        cout<< "Branch: "<<branch_labels[j]<<endl;
        variations<<">"<<bulge_count<<"."<<j<<"."<<index<<".";
        variations<<"\n"<< prefix<<dbg.node_label(i)<<branch_labels[j] << "\n";

        }

    }
}
}




void SNP (debruijn_graph<> &dbg,sd_vector<>& b,size_t &number_of_subreads, std::map<string,size_t> &pairs){
std::map<long,vector<std::pair<long,long>>> cymap;
cycle_remove(cymap, ifstream ("repeats.txt"));
for (size_t i = 0; i < dbg.num_nodes(); i++) {
	if ( dbg.outdegree(i)>1 && !visited[i] ){
		SNPCall(i,cymap,dbg,b,number_of_subreads,pairs);
     }
   }
}



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
}




void index_maker(ifstream f , std::map<string,size_t>& pairs){
string s;
size_t index = 0;
while(std::getline(f,s)){
  std::locale loc;
  string ss="";
  for (std::string::size_type i = 0; i <s.length(); ++i)
    ss += std::toupper(s[i],loc);
    pairs.insert(std::pair<string,size_t>(ss,index++));
  }
std::map<string,size_t>::iterator it = pairs.begin();

}




int main(){

vector<unsigned long long> file_size;
ifstream infile("subreads_info");
string line2;
while (std::getline(infile, line2)) {
      file_size.push_back(stol(line2)/4);
  }


ifstream input("filtered_kmc2_list.packed", ios::in|ios::binary|ios::ate);
debruijn_graph<> dbg = debruijn_graph<>::load_from_packed_edges(input, "$ACGT");
input.close();

cerr << "k             : " << dbg.k << endl;
cerr << "num_nodes()   : " << dbg.num_nodes() << endl;
cerr << "num_edges()   : " << dbg.num_edges() << endl;
cerr << "Total size    : " << size_in_mega_bytes(dbg) << " MB" << endl;
cerr << "Bits per edge : " << bits_per_element(dbg) << " Bits" << endl;

std::map<string,size_t> pairs;
cerr<<"Making the pair map"<<endl;
index_maker(ifstream ("pair") ,  pairs);
cerr<<"pair size: "<<pairs.size()<<endl;

sdsl::sd_vector<> b;
sdsl::load_from_file(b,"output_matrix");
cerr << "Total size: 			"<<size_in_mega_bytes(b) <<" Mb"<<endl;

ifstream sub("subreads_info");
string line;
size_t number_of_subreads;
while (std::getline(sub, line)) {
    number_of_subreads = stoi(line)/4;
  }
  cerr<<"number of subreads: "<<number_of_subreads<<endl;
//find_incomming_ougoing_degree_noposition(dbg);
SNP(dbg,b,number_of_subreads,pairs);
}
