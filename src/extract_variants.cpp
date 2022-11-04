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

using boost::tuple;

typedef std::unordered_map<int, std::vector<boost::tuple<std::string, std::string, std::string, std::string>>> tuple_list; 
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
          tl[pos].push_back({ REF, ALT, sample, gt}); 
        }
      }
    }

     
    start_pos = variant_positions.front();
    std::cout <<  "INFO, hged::main, start position: " << start_pos << std::endl;
    last_pos  = variant_positions.back();
    std::cout <<  "INFO, hged::main, last position: " << last_pos << std::endl;

    for (tuple_list::const_iterator i = tl.begin(); i != tl.end(); ++i) {
        std::cout << "Pos: " << i->get<0>() << std::endl;
        std::cout << "REF: " << i->get<1>() << std::endl;
        std::cout << "ALT: " << i->get<2>() << std::endl;
        std::cout << "SAMPLE: " << i->get<3>() << std::endl;
        std::cout << "GT: " << i->get<4>() << std::endl;
    }

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

  void obtain_substrings(std::string &backbone_seq, int &start_pos, const int $alpha, tuple_list &tl){

  }



int main(int argc, char **argv) {

  //parse command line arguments
  Parameters parameters;
  parseandSave_ILP(argc, argv, parameters);

  tuple_list tl;
  int start_pos, last_pos, pos, index;
  std::vector<std::string>  samples;
  std::unordered_map<int, std::string> Sample_GT;
  std::vector<int> variant_positions;
  std::string backbone_seq;


  extract_vcf_pos_ref_alt_sample_gt (parameters.vcffile, parameters.fasta_ref_file, parameters.chr, start_pos,\
   last_pos, variant_positions, tl);
  get_linear_backbone (parameters.fasta_ref_file, parameters.chr, start_pos, last_pos, backbone_seq);
  obtain_substrings (backbone_seq, start_pos, parameters.alpha, tl);


  return 0;
}




