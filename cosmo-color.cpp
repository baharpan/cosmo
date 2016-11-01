#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <libgen.h> // basename
#include <stdio.h>
#include <ctype.h>
#include <stdint.h>
#include <sdsl/int_vector.hpp>
#include <sdsl/util.hpp>
#include <algorithm> // for sort


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
std::vector<int> permutation;
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
        permutation.push_back(stoi(word));
        
        //cout<<word<<endl;
    }
    color.close();
    cerr << "num_payloads(): " << permutation.size()<<" must be equal to the number of edges"<< endl;
}

/*std::string ReadNthLine(int N){
    std::ifstream positionss("/s/fir/c/nobackup/baharpan/git/cosmo/cosmo-positions-test");
    std::string s;
    for(int i = 1; i <=N; ++i)
        std::getline(positionss,s);
    std::getline(positionss,s);
    return s;
}*/



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

void find_cycles(debruijn_graph<> dbg){
    
    bit_vector visited = bit_vector(dbg.num_nodes(), 0);
    bit_vector visited_cycle = bit_vector(dbg.num_nodes(), 0);
    
    cout << "Starting to look for cycles\n";
    size_t cycle=0;
    int count=0;
    for (size_t i =1; i <=dbg.num_nodes() ; i++) {
       // cout << i << ":" << dbg.node_label(i)<<"orginal"<<endl;
        ssize_t start = i;
        string contig;
        bit_vector visited_in_contig = bit_vector(dbg.num_nodes(), 0);
        
                    if (visited[i]!=1){
                        count++;
                       
                for (unsigned long x = 1; x < dbg.sigma + 1; x++) {
                    
                    ssize_t edge = dbg.outgoing_edge(i, x);
                    
                    if (edge == -1)
                        continue;
                    
                contig += base[x];
                    ssize_t pos = dbg._edge_to_node(edge);
                   

                    bool nocycle=false;
                
                    while (dbg.outdegree(pos)!=0 && nocycle==false){
                        visited[pos] = 1;
                        //count++;
                        
                        //cout << pos << ":" << dbg.node_label(pos)<<"pos"<<endl;
                        //visited_in_contig[pos]=1;
                         if(visited_in_contig[pos]==1 || pos==i){
                             if (visited_cycle[pos]!=1){
                                 cycle++;
                                 cout<<"cycle in "<<pos<<endl;
                                 //cerr<<cycle<<endl;
                                 visited_cycle[pos]=1;}
                             nocycle=true;}
                        visited_in_contig[pos]=1;
                         ssize_t next_edge = 0;
                         for (unsigned long x2 = 1; x2 < dbg.sigma + 1; x2++) {
                             next_edge = dbg.outgoing_edge(pos, x2);
                             if (next_edge!=-1){
                                
                               contig+= base[x2];
                            
                                 break;}
                        }
                    
                        pos = dbg._edge_to_node(next_edge);}
                    
                    
        
                }
        //                cerr <<"contig:"<<dbg.node_label(start)<<contig<<"with length"<<contig.size()<<endl;
                        cerr<<"cycle is"<<cycle<<endl;
       //                 cerr<<"cont is"<<count<<endl;
    //                    output<<">"<<"length"<<contig.size()+dbg.node_label(start).size()<<"_"<<endl;//a<<endl;//" "<<"average:"<<average<<endl;
      //                              output<<dbg.node_label(start)<<contig<<endl;
       //
                    }
        
    }
    cerr<<"Nember of cycles in original graph is "<<cycle<<endl;
}



/*void dump_edges(debruijn_graph<> dbg) {
 for (size_t i = 0; i < dbg.num_nodes(); i++) {
 cout << i << ":" << dbg.node_label(i)<< "\n";
 }
 }*/


/*void dump_edges(debruijn_graph<> dbg, uint64_t * colors) {
 for (size_t i = 0; i < dbg.num_edges(); i++) {
 cout << i << "e:" << dbg.edge_label(i) << colors[i] << "\n";
 }/{
 }*/

/*const char *const starts[] = {"GCCATACTGCGTCATGTCGCCCTGACGCGC","GCAGGTTCGAATCCTGCACGACCCACCAAT","GCTTAACCTCACAACCCGAAGATGTTTCTT","AAAACCCGCCGAAGCGGGTTTTTACGTAAA","AATCCTGCACGACCCACCAGTTTTAACATC","AGAGTTCCCCGCGCCAGCGGGGATAAACCG","GAATACGTGCGCAACAACCGTCTTCCGGAG"};*/

//size_t branch=0;
void find_contigs(debruijn_graph<> dbg, ofstream &output)
{
    std::ifstream positionss("/s/fir/c/nobackup/baharpan/yeast/lormapositionis");
    string linee;
    vector<int_vector<>> generall;
   // cerr<<"general is generated"<<endl;
    //vector<int_vector<>>general;
    //cout<<sizeof(general)<<"general"<<endl;
    
    //cout<<get_size_in_mega_bytes(vec)<<"vec"<<endl;
    //int j=0;
    while(std::getline(positionss, linee)){
        
        istringstream iss(linee);
        std::string word;
        
        int number_of_positions = 0;
        while (iss >> word){
            number_of_positions++;}
        int n;
        istringstream isss(linee);
  //      vector<int> param;
        int_vector<> vec(number_of_positions);
        //cout<<sizeof(param)<<"param"<<endl;
        //cout<<sizeof(param)<<"paramsize"<<endl;
        int_vector<>::size_type i =0;
    //    cerr<<"vec is generated"<<endl;
        //int_vector<>::size_type i=0
        //for(int_vector<>::size_type i=0;i<vec.size();i++) {
        while( isss >> n ) {
            
            //param.push_back(n);
            //cout<<sizeof(param)<<"befor"<<endl;
            vec[i]=n;
            i++;
           
            
            
        
        
            
            //cout<<sizeof(vec)<<"size of vec"<<endl;
            
        }
    
        
            //cout<<size_in_bytes(vec)<<"before"<<endl;
            util::bit_compress(vec);
        //vlc_vector<> vv(vec);
            //j+=size_in_bytes(vec);
            //cout<<size_in_bytes(vec)<<"after"<<endl;
            //cout<<vec<<endl;
        //cout<<vec[1]<<"element"<<endl;
    
        generall.push_back(vec);}
        
        
        //j++;
        //cerr<<param.size()<<"size of param"<<endl;
        //general.push_back(vec);
        //cerr<<general[j-1].size()<<"size of first element"<<endl;
    //}
    //j+=sizeof(vector<int>);
    cout<<generall.size()<<"size of general"<<endl;//add
    
    
    //cout<<j<<"total_size: Mb"<<endl;
    
    //int size=sizeof(std::vector<int>)*(general.size()+1);
    //int sizee;
    //for (int i=0; i!=general.size(); i++){
        
      //  sizee= size+=(sizeof(int) * general[i].size());
    //}
    
                 
      //  cout<<sizee<<"before compression MB"<<endl;

    int t = getMilliCount();
    int count=0;
   
    
    size_t  payload;
    
    //size_t contignum=0;
    //size_t average=0;
    bit_vector visited = bit_vector(dbg.num_nodes(), 0);
    bit_vector visited_cycle = bit_vector(dbg.num_nodes(), 0);
    bit_vector visited_cycle2 = bit_vector(dbg.num_nodes(), 0);
    
    
   // std::ifstream positions("/s/fir/c/nobackup/baharpan/git/cosmo/Ecoli-narrow/cosmo-positions-second");
    //std::vector <std::string> poslines;
    //std::string line;
    //int andaz=0;
   // while (std::getline(positions, line)){
     //   andaz+=line.size();
       // poslines.push_back(line);
        
    //}
   /* andaz+=sizeof(vector<string>);
    cout<<andaz<<"size of andaz Byte"<<endl;
    cout<<sizeof(vector<string>)<<"string"<<endl;
    cout<<sizeof(vector<int>)<<"int"<<endl;
    cout<<sizeof(int_vector<>)<<"int"<<endl;*/
    
    
    //cout << myLines.size() <<endl;
    cerr << "Finding contigs:\n";
   size_t cycle=0;
    size_t allcycle=0;
    for (size_t i =1; i<=dbg.num_nodes(); i++) {
        //std::vector<string> position;i
        int_vector<> posl;
        int_vector<>::iterator posit;
        ssize_t start = i; // place to store start of branch kmer
        std::string start_label(dbg.node_label(start));
        //int_vector<> posl;
        payload=permutation[dbg._node_to_edge(i)];
        
     //   cout << i << ":" << dbg.node_label(i)<<"and payload orginal"<<payload<<endl;//add
        //cout<<"outdegree is"<<dbg.outdegree(i)<<endl;
        //cout<<"indegree is"<<dbg.indegree(i)<<endl;
        //if (dbg.outdegree(i)>1){
          //  cout<<"wow"<<endl;
        //}
        if(payload==-1){
            continue;}
        //    cout<<"save from $"<<endl;}
        if(payload!=-1){
        //string pos=poslines[payload];
          //  cout<<pos<<endl;
             posl=generall[payload];}
         //   cout<<posl<<"THis is new"<<endl;}
        //istringstream iss(pos);
        //copy(istream_iterator<string>(iss),
          //   istream_iterator<string>(),
            // back_inserter(position));}
        bool empty=false;
       // cout<<"positions in i"<<posl.size()<<endl;
        size_t qnum=0;
        for(size_t q=1;q!=posl.size();q++){
            
            if(posl[q]==0){
                qnum++;}}
            if(qnum==posl.size()){
                empty=true;}
       // if(visited[i]==0){
         //   break;}
        if(empty==false && !visited[i]){
            
            count++;
         //   cout<<i<<"is seen count"<<endl;
       // string contig;    
        posit=posl.begin();
       while(posit!=posl.end()&&*posit==0){posit++;} 
            
        while(posit!=posl.end()&& *posit!=0 ){
           // if(*posit==0){;}//&& *posit!=0){
        //    bool howw=false;
        bit_vector visited_in_contig = bit_vector(dbg.num_nodes(), 0);
    
        bit_vector visited_in_contig2 = bit_vector(dbg.num_nodes(), 0);
           // bit_vector visited_in_position=bit_vector(general[payload].size() , 0);
       string contig;
        
        //if (!visited[i]){
            
            
            //visited[i] = 1;
            //std::vector<string> next_position;
            int_vector<> poss;
            int_vector<>::iterator nextpos;
            
            for (unsigned long x = 1; x < dbg.sigma + 1; x++) {
                
                ssize_t edge = dbg.outgoing_edge(i, x);
                
                if (edge == -1)
                    continue;
                
                //contig += base[x];
                
                //cerr<< "first base is"<<base[x]<<endl;
                
                ssize_t pos = dbg._edge_to_node(edge);
                if (visited[pos]==1){break;}
                //cout<<"outdegree is"<<dbg.outdegree(pos)<<endl;
                int payyload = permutation[dbg._node_to_edge(pos)];
                
                
             //  cout << pos << ":" << dbg.node_label(pos)<<"and payload"<<payyload<<endl;//add
                
                visited_in_contig[pos]=1;
                visited_in_contig2[pos]=1;
                
                
                if (payyload == -1){
                    continue;
                }
                
                
                //string poss=ReadNthLine(payyload);
                poss=generall[payyload];
                int_vector<> posss;
                //cout<<"1:"<<posst<<"2:"<<poss;
             //   cout<<poss<<endl;
                
                
              bool howw=false;  
                
                //istringstream iss(poss);
                //copy(istream_iterator<string>(iss),
                  //   istream_iterator<string>(),
                    // back_inserter(next_position));
               // if (!!!next_position.size()==0){//|| visited[pos]==1){
                 //   break;
                //}
                //if (visited[pos]==1){
                  //  break;
                //}
                //nextpos=next_position.begin();
                //nextpos=next_position/visited2.begin();
                
                //while (nextpos2!=next_position2.end()){
                //if(stoi(*nextpos2)-stoi(*nextpos)>=0){
                //int_vector<> test(next_position.size());
             /*   if (visited[pos]==1){
                    break;
                 
                }*/
                size_t m=0;
                for (nextpos=poss.begin();nextpos!=poss.end();nextpos++){
                    //cout<<"difference"<<(*nextpos -*posit)<<endl;
                    //&&(*nextpos-*posit)>-500
                    if ((*nextpos-*posit)<=500&&  *nextpos!=0 ){
                     //   if(*nextpos==0){continue;}
                        
                        //&& !visited_in_position[*nextpos]){
                        
                        //&& stoi(*nextpos2)-stoi(*nextpos)>=0){
                        //visited_in_position[*nextpos]=1;                        //found+=1;
                        howw=true;
               //         cout<<*posit<<":"<<*nextpos<<endl;
                 //       cout<<"howw is pos"<<howw<<endl;

                        break;}
                    m++;
                    
                }
                if(howw==true){
                contig+= base[x];}
                //break};
                                                else{
                                                                    continue;}
              // cout<<m<<"m is"<<endl;
                if(m<poss.size()&& howw==true){ 
               
                //cout<<poss<<"before"<<endl;
                (generall[payyload])[m]=0;}
              //  cout<<generall[payyload]<<"after"<<endl;}
                size_t lnum=0;
                for(size_t l=0; l!=generall[payyload].size(); l++){    
                    if (generall[payyload][l]==0){
                        lnum++;}
                        }
                   if(lnum==generall[payyload].size()){
                     visited[pos]=1;
                    
                    count++;
                //    cout<<"pos is seen"<<endl;
                    //sdsl::util::clear(general[payyload]);;
                }
                
               
                
                    
                

                //cout<<"size of next_position:"<<next_position.size()<<endl;
                
                //std::size_t indexx = general[payyload].find(*nextpos);
                //cout<<"indexx is"<<indexx<<endl;
               // if (indexx != general[payyload].end()){
                    //poslines[payyload]=poslines[payyload].erase(indexx,(*nextpos).length()+1);
                    //general[payyload]=general[payyload].erase(general[payyload].begin()+indexx);
                    //cout<<"for pos:"<<general[payyload];
                    //if (general[payyload].empty()){
                      //     count++;
                        //    cout<<"count for pos"<<endl;}
                //}
                //else{
                  //  break;
                //}
                //cout<<poslines[payyload]<<endl;
                
                
                //cout<<*posit<<":"<<*nextpos<<endl;                //poslines[payyload]=poslines[payyload].substr(poslines[payyload].find_first_of(" \t")+1);
                int payyload2=0;
                // while ( dbg.indegree(pos)&& dbg.outdegree(pos) == 1/*&& !visited[pos]*/) {
                bool breakk=true;
                //std::vector<int> visited_pos (dbg.num_nodes(), 0);
                
                 bool nocycle=false;
                
                size_t index=0;
                size_t k=0;
                nextpos=poss.begin()+m;
                while (dbg.outdegree(pos)!=0 && breakk==true){// && nocycle==false){ //&& breakk==true){/* && !visited[pos] e& dbg.outdegree(pos) == 1 && breakk==true){//breakk==false){*/
                    //visited[pos] = 1;
                    ssize_t next_edge = 0;
                    //ssize_t next = 0;
                    
                    //std::vector<string>next_position2;
                    int_vector<>::iterator nextpos2;
                    for (unsigned long x2 = 1; x2 < dbg.sigma + 1; x2++) { // iterate through the alphabet
                        next_edge = dbg.outgoing_edge(pos, x2);
                        if(next_edge != -1){
                            ssize_t next = dbg._edge_to_node(next_edge);
                            //if(dbg.outdegree(next)>1){
                              //  branch++;
                                //cout<<"more thaun one"<<next << ":" << dbg.node_label(next)<<endl;
                            if(visited[next]==1){
                            breakk=false; }
                          //  cout<<"Can not seen again"<<endl;}
                            
                            
                            payyload2 = permutation[dbg._node_to_edge(next)];
                           // cout << next << ":" << dbg.node_label(next)<<"and payload hello"<<payyload2<<endl;//add
                         // if (visited_cycle[next]!=1){
                           //     visited_cycle[next]=1;
                            //visited_cycle[next]=1;
 // visited_pos[next]++;
                            // cout<<visited_p/os[next]<<endl;
                            //if (visited_pos[next]>1){
                            //   breakk=false;
                            //}
                            if (payyload2==-1)
                                break;
                       
//add if(visited_in_contig2[next]==1 || next==i ){
           //add                    if (visited_cycle2[next]!=1){
                //add                allcycle++;
                   //             cout<<"cycle in "<<
                     //add                               visited_cycle2[next]=1;}}
                //add   visited_in_contig2[next]=1; 
                            //string posss=ReadNthLine(payyload2);
                            posss=generall[payyload2];
                            //cout<<"1:"<<possst<<"2:"<<posss;
                      //      cout<<posss<<endl;
                            
                            //istringstream isss(posss);
                            //copy(istream_iterator<string>(isss),
                              //   istream_iterator<string>(),
                                // back_inserter(next_position2));
                            //if (next_position2.size()==0){
                              //  break;
                            //}
                         /*   if(visited[next]==1){
                                breakk=false;}*/

                            
                            //cout<<"size of nextpos2 is"<<next_position2.size()<<endl;
                            
                            bool how=false;
                            
                            
                            
                           //for (size_t i=0; i<posss.size();i++){
                             //       for (size_t j=0; j<poss.size();j++){ //
                            
                            
                            //nextpos2=next_position2.begin();
                            
                            //while (nextpos2!=next_position2.end()){
                                //if(stoi(*nextpos2)-stoi(*nextpos)>=0){
                            
                            //std::vector<int> difference;
                            //std::vector<int>::iterator iter;
                            //auto iter=difference.begin();
                            size_t n =0;
                            for (nextpos2=posss.begin();nextpos2!=posss.end();nextpos2++){
                                //&&*nextpos2-*nextpos>-500
                               
                                if (*nextpos2-*nextpos<500  && *nextpos2!=0){
   //                                 if(*nextpos2==0){continue;}
                                    how=true;
                                 //   cout<<*nextpos<<":"<<*nextpos<<endl;
                                    break;
                                    }
                                n++;
                                }
                            k=n;
                           // else{
                                
                          //  breakk=false;
                          //  }
                            if (k+1>=posss.size()&& how==true){
                                visited[next]=1;
                                count++;}
                             //   cout<<"next is seen"<<endl;
                          //  }
                        //    cout<<k<<"n is"<<endl;
                            if(k<posss.size() && how==true){
                            //cout<<k<<"n is"<<endl;
                         //   cout<<"how is"<<how<<endl;
                           // cout<<posss<<"before"<<endl;
                            (generall[payyload2])[k]=0;}
                          //  cout<<generall[payyload2]<<"after"<<endl;}
                             size_t gnum=0;
                             for (size_t g=1; g!=generall[payyload2].size();g++){
                                if (generall[payyload2][g]==0){
                                    gnum++;}}
                            if (gnum==generall[payyload2].size()){
                                visited[next]=1;
                                count++;
                              //  cout<<"next is seen"<<endl;
                            }
                            if(how==false)
                                breakk=false;

                           
                            if(visited_in_contig[next]==1 || next==i){
                               if (visited_cycle[next]!=1){
                                cycle++;
                   //             cout<<"cycle in "<<next<<endl;
                                                                   visited_cycle[next]=1;}}
                   //                                                                             nocycle=true;}
                   //                                                                                                       //  breakk=false;}
                                                                                                                                                    visited_in_contig[next]=1;
//visited_in_contig2[next]=1;
                             
                            contig+= base[x2];
                            break;
                          
                            
                        }
                    
                    }//for x2
                    pos = dbg._edge_to_node(next_edge);
                   poss=posss;
                    nextpos=poss.begin()+k;
                }//main while
                
            }//for x
            
          //  cerr <<"contig:"<<dbg.node_label(start)<<contig<<"with length"<<contig.size()+dbg.node_label(start).size()<<endl;
            //cerr <<"Number of branch is"<<branch<<endl;
        //add    cerr<<"cycle is"<<cycle<<endl;
      //add      cerr<<"allcycleis"<<allcycle<<endl;
            //cerr<<"count is"<<count<<endl;
            //cerr<<"empty is"<<empty<<endl;
            //cerr<<"total traveresed"<<empty+count<<endl;
            //for(int i=0; i<equal+1; ++i){
            //ostringstream convert;
            //convert<<i;
            //string a= convert.str();
            output<<">"<<"length"<<contig.size()+dbg.node_label(start).size()<<"_"<<endl;//a<<endl;//" "<<"average:"<<average<<endl;
            output<<dbg.node_label(start)<<contig<<endl;
            // }
            posit++;
            
        }//while posit
            
        }//if i
        visited[i] = 1;
        
        //general[payload]={0};
  //   cerr<<"count is"<<count<<endl; 
       // cerr <<"contig:"<<dbg.node_label(start)<<contig<<"with length"<<contig.size()+dbg.node_label(start).size()<<endl;
    }//for(i
   // cerr << "Find contigs time: " << getMilliSpan(t) << std::endl;
   // sdsl::util::clear (general);
   cerr<<"cycle is"<<cycle<<endl;
}//void


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
    //cerr << "colors        : " << 74319898 << endl;
    readFile(color);
    find_contigs(dbg,output);
    //find_cycles(dbg);
    }
