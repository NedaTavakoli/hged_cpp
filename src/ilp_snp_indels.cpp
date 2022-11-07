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
#include "gurobi_c++.h"
#include "ext/subprocess.hpp"
#include <boost/tokenizer.hpp>
#include "boost/tuple/tuple.hpp"

using boost::tuple;

/**
 * @brief   Extract POS, REF, ALT, SAMPLe=GT, from the vcf and fasta file 
 *          Input:  vcf file, reference fasta file, alpha, chr
 *          Output:  start_pos, last_pos, varriant_positions, REFs, ALTs, Sample_GT
 *
 */
void  extract_vcf_pos_ref_alt_sample_gt (const std::string &vcf_file, const std::string &fasta_file, const int &chr,\
      int &start_pos, int &last_pos, std::vector<int> &variant_positions, std::vector<std::string> &REFs,\
      std::vector<std::string> &ALTs, std::unordered_map<int, std::string> &Sample_GT)
{
   // seed random generator by time in seconds (this may create issue if two instances are launched at the same time)
    srand(time(0)); int random = rand() % 100000;  
    std::string tmp_file = ".pos_ref_alt_sample_GT_chr_" + std::to_string(chr) + "_" + std::to_string(random) + ".txt";
    std::string cmd = std::string(TOSTRING(BCFTOOLSPATH)) + " query -i \'(TYPE=\"snp\" || TYPE=\"indel\") && GT=\"alt\"\'" + " -f \'%POS\t%REF\t%ALT[\t%SAMPLE=%GT]\n\'  "  + vcf_file + " >  " + tmp_file;

    std::cout << "INFO, hged::main, extracting pos, ref, alt, non-zero GT from vcf file using command: " << cmd << std::endl;
    std::system(cmd.c_str());

    std::ifstream file (tmp_file);
    std::string line;
    std::string REF, ALT;
    std::string sample_gt;
    int count = 0; 

    while (std::getline(file, line))
    {
      std::istringstream iss(line);
      int pos;
      
      iss >> pos >> REF >> ALT;
      if(line.size() > 0) 
      {
        variant_positions.push_back(pos);
        REFs.push_back(REF);
        ALTs.push_back(ALT);
      
      while (iss >> sample_gt)
      {
        Sample_GT[count]  += sample_gt + " "; 
      }
      }
      count += 1;
    }

     
    start_pos = variant_positions.front();
    std::cout <<  "INFO, hged::main, start position: " << start_pos << std::endl;
    last_pos  = variant_positions.back();
    std::cout <<  "INFO, hged::main, last position: " << last_pos << std::endl;
 }

void get_index_for_variant_position (std::vector<int> &variant_positions, const int &pos, int &index){
    auto it = find(variant_positions.begin(), variant_positions.end(), pos);
    index = it - variant_positions.begin();
 }

 void get_linear_backbone(const std::string &fasta_file, const int &chr, int &start_pos, int &last_pos, std::string &linear_backbone){

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
    linear_backbone = content;  
 }

void get_samples(const std::string &vcf_file, std::vector<std::string> &samples){

  // Extract samples from vcf file
    srand(time(0)); int random = rand() % 100000;  
    std::string tmp_file = "sample." + std::to_string(random) + ".txt";
    std::string cmd = std::string(TOSTRING(BCFTOOLSPATH)) + " query -l " + vcf_file + " >  " + tmp_file;
    std::cout << "INFO, hged::main, extracting samples from vcf file using command: " << cmd << std::endl;
    std::system(cmd.c_str());

    std::ifstream file (tmp_file);
    std::string line;
    while (std::getline(file, line))
    {
      std::string col;
      std::istringstream iss(line);
      iss >> col;
      samples.push_back(col);
    }
}

 void get_variants_in_range_alpha(std::vector<int> &variant_positions, const int &alpha,\
     std::unordered_map<int, std::vector<int>> &variants_in_range ){
    
    std::cout << "INFO, hged::main, get variants in range alpha " << std::endl;
    for (std::size_t i = 0; i < variant_positions.size(); i++){
      for (std::size_t j = i+1; j < variant_positions.size(); j++)
      {
        if (variant_positions[j]- variant_positions[i] <= alpha) variants_in_range[i].push_back(variant_positions[j]);
      }
    }
  }
  

 void extract_pos_substring(std::string &linear_backbone, const int &alpha,\
     int &start_pos, std::vector<std::string> &samples, std::vector<int> &variant_positions, std::unordered_map<int,std::vector<int>> &variants_in_range, std::vector<std::string> &ALTs, std::unordered_map<int, std::string> &Sample_GT)  
 {

  srand(time(0)); int random = rand() % 100000;  
 // std::string output_file = ".pos_substrings_chr22_" +  std::to_string(random) + ".txt";
   std::ofstream output_file(".pos_substrings_chr22_" +  std::to_string(random) + ".txt");

    std::cout << "INFO, hged::main, extracting pos substring " << std::endl;
    
   //for (int w= 0; w <100; w++)
   for (int w= 0; w <variant_positions.size(); w++)
   {

    std::vector<std::string> substrings;
    std::vector<std::string> substrings_u;
    //int v = 16051453; // commet this if the above is uncomment
    int v = variant_positions[w];

    int index;
    get_index_for_variant_position (variant_positions, v, index); // i is the variant position index, should be changed with the above one
    
    std::string ref_haplotype = linear_backbone.substr(variant_positions[index] - start_pos, alpha);
    std::string h1 = ref_haplotype;
    std::string h2 = ref_haplotype;
    
    std::vector<int> list_variants_in_range;
    list_variants_in_range = variants_in_range[index];
    std::vector<int>::iterator it;
    // add variant i to the list_vatiants_in_range: List of variants needs to be addressed
    list_variants_in_range.insert (list_variants_in_range.begin(), variant_positions[index]);
    
    std::string  gt_h1, gt_h2;
    int idx;


    for (std::size_t k = 0 ; k < list_variants_in_range.size(); k++)
    {

          get_index_for_variant_position (variant_positions, list_variants_in_range[k] , idx); // i is the variant position index, should be changed with the above one
    
          std::size_t found;
          for (std::size_t j = 0 ; j < samples.size(); j++)
          {
            found = Sample_GT[idx].find(samples[j]);
            if(found != std::string::npos) // if the sample found in the variant position
            {
              gt_h1 = Sample_GT[idx].substr(found+8, 1);  // to extract from SAMPLE=GT
              gt_h2 = Sample_GT[idx].substr(found+10, 1);

              // split ALT for that variant position over commas 
              std::vector<std::string> result;
              std::stringstream s_stream(ALTs[idx]); //create string stream from the string
              while(s_stream.good()) {
                  std::string substr;
                  getline(s_stream, substr, ','); //get first string delimited by comma
                  result.push_back(substr); // ALT1/ ALT2, ALT3
              }


              if (std::stoi(gt_h1) !=0)
              {
                if (std::stoi(gt_h1) ==1) h1.replace(list_variants_in_range[k] - list_variants_in_range[0],  result[0].length(), result[0]);  //ALT1
                if (std::stoi(gt_h1) ==2) h1.replace(list_variants_in_range[k] - list_variants_in_range[0],  result[1].length(), result[1]);  //ALT2
                if (std::stoi(gt_h1) ==3) h1.replace(list_variants_in_range[k] - list_variants_in_range[0],  result[2].length(), result[2]);  //ALT3
                substrings.push_back(h1);
              }

              if (std::stoi(gt_h2) !=0){
                if (std::stoi(gt_h2) ==1) h2.replace(list_variants_in_range[k] - list_variants_in_range[0],   result[0].length(), result[0]);  //ALT1
                if (std::stoi(gt_h2) ==2) h2.replace(list_variants_in_range[k] - list_variants_in_range[0],   result[1].length(), result[1]);  //ALT2
                if (std::stoi(gt_h2) ==3) h2.replace(list_variants_in_range[k] - list_variants_in_range[0],  result[2].length(), result[2]);  //ALT3
                substrings.push_back(h2);
                
              }
            }
          }
      }

    
      //substrings (unique values)
      substrings_u.insert (substrings_u.end(), substrings.begin(), substrings.end());
      std::sort (substrings_u.begin(), substrings_u.end());
      substrings_u.erase(std::unique(substrings_u.begin(), substrings_u.end()), substrings_u.end() );
      //std::cout << "substrings for varinat position " << variant_positions[w] << " is " << substrings_u << std::endl; 
      output_file << substrings_u << "\n"; 

  }  
 


 }   



int main(int argc, char **argv) {

  //parse command line arguments
  Parameters parameters;
  parseandSave_ILP(argc, argv, parameters);

  int start_pos, last_pos, pos, index;
  std::vector<std::string>  samples;
  std::unordered_map<int, std::vector<int>> variants_in_range;
  std::unordered_map<int, std::string> Sample_GT;
  //std::unordered_map<int, std::unordered_map<std::string, std::string>> Sample_GT;
  std::vector<int> variant_positions;
  std::vector<std::string> REFs;
  std::vector<std::string> ALTs;
  //std::unordered_map<std::string, std::string> Sample_GT;  // dictionary sample:GT
  std::string linear_backbone;

  extract_vcf_pos_ref_alt_sample_gt (parameters.vcffile, parameters.fasta_ref_file, parameters.chr, start_pos,\
   last_pos, variant_positions, REFs, ALTs, Sample_GT);
  get_index_for_variant_position(variant_positions, pos, index);
  get_linear_backbone (parameters.fasta_ref_file, parameters.chr, start_pos, last_pos, linear_backbone);
  get_samples (parameters.vcffile, samples);
  get_variants_in_range_alpha (variant_positions,parameters.alpha, variants_in_range);
  extract_pos_substring (linear_backbone, parameters.alpha, start_pos, samples, variant_positions, variants_in_range, ALTs,\
   Sample_GT);





  return 0;
}

  
     
  
  
  
  