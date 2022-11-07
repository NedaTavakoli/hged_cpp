#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector> 
#include <algorithm>
#include <ctime>
#include <cstring>
#include <chrono>
#include <numeric>
#include <cassert>
#include <unordered_map>
#include "ext/prettyprint.hpp"
#include "common.hpp"
#include <boost/graph/graph_traits.hpp>
#include "ext/subprocess.hpp"
#include <boost/tokenizer.hpp>
#include "boost/tuple/tuple.hpp"  // to use tuple, https://stackoverflow.com/questions/6482338/how-to-create-a-list-of-tuples-c
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include <set>

using boost::tuple;

//typedef std::unordered_map<int, std::vector<boost::tuple<std::string, std::string, std::string, std::string>>> tuple_list; 
using value_elemet = boost::tuple<std::string, std::string, std::string, std::string>;
using Value_type = std::vector<value_elemet>;
typedef std::unordered_map<int, Value_type> tuple_list; 

//typedef std::vector<boost::tuple<int, std::string, std::string, std::string, std::string>> tuple_list; // this vector is used for POS, REF, ALT, SAMPLE=GT

/**
 * @brief   Extract POS, REF, ALT, SAMPLe=GT, from the vcf and fasta file 
 *          Input:  vcf file, reference fasta file, alpha, chr
 *          Output: tuple_list contain list of tuples of POS, REF, ALT, SAMPLE=GT  
 *
 */
void  extract_vcf_pos_ref_alt_sample_gt (const std::string &vcf_file, const std::string &fasta_file, const int &chr,\
      int &start_pos, int &last_pos, std::vector<int> &variant_positions, tuple_list &tl)
{
   // seed random generator by time in seconds (this may create issue if two instances are launched at the same time)
    srand(time(0)); int random = rand() % 100000;  
    std::string tmp_file = ".pos_ref_alt_sample_GT_chr_" + std::to_string(chr) + "_" + std::to_string(random) + ".txt";
    std::string cmd = std::string(TOSTRING(BCFTOOLSPATH)) + " query -i \'(TYPE=\"snp\" || TYPE=\"indel\") && GT=\"alt\"\'" + " -f \'%POS\t%REF\t%ALT[\t%SAMPLE\t%GT]\n\'  "  + vcf_file + " >  " + tmp_file;

    std::cout << "INFO, hged::main, extracting pos, ref, alt, non-zero GT from vcf file using command: " << cmd << std::endl;
    std::system(cmd.c_str());

    std::ifstream file (tmp_file);
    std::string line;
    std::string REF, ALT;
    std::string sample;
    std::string gt;
    //tuple_list tl;

    while (std::getline(file, line))
    {
      std::istringstream iss(line);
      int pos;
      
      iss >> pos >> REF >> ALT;
      
      if(line.size() > 0) 
      {
        variant_positions.push_back(pos);
        while (iss >> sample >> gt)
        {
          //tl.push_back({pos, REF, ALT, sample, gt});
          //tl[pos].push_back({ REF, ALT, sample, gt});
          tl[pos].push_back(boost::make_tuple(REF, ALT, sample, gt)); 
        }
      }
    }

     
    start_pos = variant_positions.front();
    std::cout <<  "INFO, hged::main, start position: " << start_pos << std::endl;
    last_pos  = variant_positions.back();
    std::cout <<  "INFO, hged::main, last position: " << last_pos << std::endl;

    // for (const auto& i: tl) {
    //   std::cout << "POS: " << i.first << std::endl;
    //   for (std::vector<value_elemet>::size_type j =0 ; j < i.second.size() ; ++j){
    //     std::cout << "POS: " << i.first << std::endl;
    //     std::cout << "REF: " << boost::get<0> (i.second.at(j))  << std::endl;
    //     std::cout << "ALT: " <<  boost::get<1> (i.second.at(j))  << std::endl;
    //     std::cout << "SAMPLE: " << boost::get<2> (i.second.at(j))  << std::endl;
    //     std::cout << "GT: " << boost::get<3> (i.second.at(j))  << std::endl;
    // }
    // }

}


void get_linear_backbone(const std::string &fasta_file, const int &chr, int &start_pos, int &last_pos, std::string &backbone_seq){

    // Extract linear backbone using samtools
    // seed random generator by time in seconds (this may create issue if two instances are launched at the same time)
    srand(time(0)); int random = rand() % 100000;  
    std::string tmp_file = ".linear_bc_chr" + std::to_string(chr) + "_" + std::to_string(random) + ".txt";
    std::string cmd = std::string(TOSTRING(SAMTOOLSPATH)) + " faidx " +  fasta_file + " " + std::to_string(chr) + ":" + std::to_string(start_pos) + "-" + std::to_string(last_pos) + " >  " + tmp_file;
    std::cout << "INFO, hged::main, get linear backbone genome from fasta file using command: " << cmd << std::endl;
    std::system(cmd.c_str());

    std::ifstream file (tmp_file);
    std::string line;

    // ignore fist line and remove whitespaaces
    std::getline(file, line);  //ignore the first line
   
    // store the rest of file into a string named content
    // https://stackoverflow.com/questions/2912520/read-file-contents-into-a-string-in-c
    std::string content( (std::istreambuf_iterator<char>(file) ),(std::istreambuf_iterator<char>()) );

    // using the erase, remove_if, and ::isspace functions --> remove all the whitespaces from the string
    content.erase(std::remove_if(content.begin(), content.end(), ::isspace),
            content.end());                   
    backbone_seq= content;  
 }

  void obtain_substrings(std::string &backbone_seq, int &start_pos, const int &alpha, tuple_list &tl ){


    // seed random generator by time in seconds (this may create issue if two instances are launched at the same time)
    srand(time(0)); int random = rand() % 100000;  
    std::string pos_sub = ".pos_sub_alpha_" +std::to_string(alpha) + "_" + std::to_string(random) + ".txt";
    std::cout << "INFO, hged::obtain_substring, extracting substrings per position and save to the file named: "  << pos_sub << std::endl;

    std::ofstream pos_sub_file (pos_sub);
    
    std::unordered_map<int, std::set<std::string>> pos_substrings;
    
    int final_pos = start_pos + backbone_seq.size();
    std::cout << "final pos: " << final_pos << std::endl;
    std::string sample;
    std::string s;

    for (const auto& i: tl) {
      std::cout << "POS: " << i.first << std::endl;
      for (std::vector<value_elemet>::size_type j =0 ; j < i.second.size() ; ++j){
        // std::cout << "POS: " << i.first << std::endl;
        // std::cout << "REF: " << boost::get<0> (i.second.at(j))  << std::endl;
        // std::cout << "ALT: " <<  boost::get<1> (i.second.at(j))  << std::endl;
        // std::cout << "SAMPLE: " << boost::get<2> (i.second.at(j))  << std::endl;
        // std::cout << "GT: " << boost::get<3> (i.second.at(j))  << std::endl;

        sample = boost::get<2> (i.second.at(j));

        std::string alt = boost::get<3> (i.second.at(j));
        if( (alt.compare(0, 1, "0") != 0) || (alt.compare(2, 1, "0") != 0) ){

        
          int current_pos = i.first ;
          std::string s = " ";
          bool sample_found_at_pos;
          std::string alt_choice;
        

          while (current_pos < final_pos && s.size() < alpha){
            if (tl.find(current_pos) != tl.end())  { // the position exists in the unordered map
                sample_found_at_pos = false;
                std::string sam = boost::get<2> (i.second.at(j));
                if (sam.compare(sample) == 0) {
                  std::string gt = boost::get<3> (i.second.at(j));
                  if (gt.compare(0,1,"0") != 0) 
                    alt_choice = gt.substr(0,1);
                  if (gt.compare(2,1,"0") != 0) 
                    alt_choice = gt.substr(2,1); 
                   if (alt_choice.compare("0") !=0){
                    std::vector<std::string> alt_splited;
                    std::vector<std::string>::iterator it;
                    boost::split(alt_splited, boost::get<1> (i.second.at(j)), boost::is_any_of(",")); //split ALT over comma, save to a vector
                    s += alt_splited[atoi(alt_choice.c_str())-1];
                    current_pos +=  (boost::get<0> (i.second.at(j))).size(); // add reference length
                    sample_found_at_pos = true;
                   }  
                 }

                if (sample_found_at_pos == false){
                  s += backbone_seq.substr(current_pos-start_pos, 1);
                  current_pos +=1;
               }
            }
            else{
              s += backbone_seq.substr(current_pos-start_pos, 1);
              current_pos +=1;
           }

           // pos_substrings[i.first].insert(s);

          }
            pos_substrings[i.first].insert(s);
            
        }
       // pos_substrings[i.first].insert(s);
      
    }
        pos_sub_file << std::to_string(i.first) + " " << pos_substrings[i.first] << std::endl;
        std::cout << "POS_substring for POS " << i.first << " set of substring" << pos_substrings[i.first] <<std::endl;
    }
   
  }

  void construct_graph (std::string &backbone_seq, int &start_pos, int &last_pos, tuple_list &tl){
    
    // seed random generator by time in seconds (this may create issue if two instances are launched at the same time)
    srand(time(0)); int random = rand() % 100000;  
    std::string graph_file = ".construct_graph_" + std::to_string(random) + ".txt";
    std::cout << "INFO, hged::construct graph, constructing graph file: "  << graph_file << std::endl;

    std::ofstream graph_file_name (graph_file);
    

    for (int i = 0 ; i < backbone_seq.size() ; ++i){

      // add backbone edges
      graph_file_name << std::to_string(start_pos + i) + " " + std::to_string(start_pos + i+ 1) + " " +backbone_seq[i] \
      + " " + "-" << std::endl;
    }

    // add alt paths
     int start, end, new_vertex;
     new_vertex = last_pos + 2;
     int index = 0 ;

     for (const auto& i: tl) {
      std::cout << "POS: " << i.first << std::endl;
      for (std::vector<value_elemet>::size_type j =0 ; j < i.second.size() ; ++j){
        // std::cout << "POS: " << i.first << std::endl;
        // std::cout << "REF: " << boost::get<0> (i.second.at(j))  << std::endl;
        // std::cout << "ALT: " <<  boost::get<1> (i.second.at(j))  << std::endl;
        // std::cout << "SAMPLE: " << boost::get<2> (i.second.at(j))  << std::endl;
        // std::cout << "GT: " << boost::get<3> (i.second.at(j))  << std::endl;

        std::string alt = boost::get<1> (i.second.at(j));
        std::cout << "POS " << i.first << " ALT " << alt << std::endl;
        std::string ref = boost::get<0> (i.second.at(j));
        std::vector<std::string> alt_splited;
        std::vector<std::string>::iterator it;
        boost::split(alt_splited, alt, boost::is_any_of(",")); //split ALT over comma, save to a vector
        for (size_t k = 0; k < alt_splited.size() ; ++k ){ // for each element in alt
          if (k == 0){
            start = i.first;
          }else{
            start = new_vertex;
          }   
          
          if (k == (alt_splited[k].size()-1) ){
             end = i.first + ref.size();

          }else{
            new_vertex += 1;
            end = new_vertex;

          }

          // TODO : remove duplicates
          graph_file_name << std::to_string(start) + " "  + std::to_string(end) + " " + alt_splited[k] + " " + std::to_string(i.first) << std::endl;
          //graph_file_name << std::to_string(start) + " "  + std::to_string(end) + " " + alt_splited[k] + " " + std::to_string(index) << std::endl;

        }
        
    }
      index +=1;
      

    }

    // remove duplicate lines


  }



int main(int argc, char **argv) {

  //parse command line arguments
  Parameters parameters;
  parseandSave_ILP(argc, argv, parameters);

  tuple_list tl;
  std::ofstream pos_sub_file, graph_file_name;
  int start_pos, last_pos, pos, index;
  std::vector<std::string>  samples;
  std::unordered_map<int, std::string> Sample_GT;
  std::vector<int> variant_positions;
  std::string backbone_seq;


  extract_vcf_pos_ref_alt_sample_gt (parameters.vcffile, parameters.fasta_ref_file, parameters.chr, start_pos,\
   last_pos, variant_positions, tl);
  get_linear_backbone (parameters.fasta_ref_file, parameters.chr, start_pos, last_pos, backbone_seq);
  obtain_substrings (backbone_seq, start_pos, parameters.alpha, tl);
  construct_graph (backbone_seq, start_pos, last_pos, tl);


  return 0;
}




