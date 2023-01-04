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
#include "boost/tuple/tuple.hpp"
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include <set>
#include <chrono>

using boost::tuple;

using value_elemet = boost::tuple<std::string, std::string, std::string, std::string>;
using Value_type = std::vector<value_elemet>;
typedef std::map<int, Value_type> tuple_list;

/**
 * @brief   Extract POS, REF, ALT, SAMPLe, GT, from the vcf and fasta file
 *          Input:  vcf file, reference fasta file, chr, first variant position, last variant position, list of variant positions
 *          Output: tuple_list contains is a dictionary with key as variant position and value is a list of tuples containing (REF, ALT, SAMPLE, GT)
 *          the form of the output is {POS_1:[(REF, ALT, SAMPLE, GT), (REF, ALT, SAMPLE, GT),.... ], POS_2, [(REF, ALT, SAMPLE, GT), (REF, ALT, SAMPLE, GT),.... ],....}
 *
 */
void extract_vcf_pos_ref_alt_sample_gt(const std::string &vcf_file, 
                                const std::string &fasta_file, 
                                const int &chr,
                                int &start_pos,
                                int &last_pos, 
                                std::vector<int> &variant_positions, 
                                tuple_list &tl)
{
    // seed random generator by time in seconds (this may create issue if two instances are launched at the same time)
    srand(time(0));
    int random = rand() % 100000;
    std::string tmp_file = "pos_ref_alt_sample_GT_chr_" + std::to_string(chr) + "_" + std::to_string(random) + ".txt";
    std::string cmd = std::string(TOSTRING(BCFTOOLSPATH)) + " query -i \'(TYPE=\"snp\" || TYPE=\"indel\") && GT=\"alt\"\'" + " -f \'%POS\t%REF\t%ALT[\t%SAMPLE\t%GT]\n\'  " + vcf_file + " >  " + tmp_file;

    std::cout << "INFO, hged::main, extracting pos, ref, alt, non-zero GT from vcf file using command: " << cmd << std::endl;
    std::system(cmd.c_str());

    std::ifstream file(tmp_file);
    std::string line;
    std::string REF, ALT;
    std::string sample;
    std::string gt;

    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        int pos;
        iss >> pos >> REF >> ALT;

        if (line.size() > 0)
        {
            variant_positions.push_back(pos);
            while (iss >> sample >> gt)
            {
                tl[pos].push_back(boost::make_tuple(REF, ALT, sample, gt));
            }
        }
    }

    start_pos = variant_positions.front();
    std::cout << "INFO, hged::main, start position: " << start_pos << std::endl;
    last_pos = variant_positions.back();
    std::cout << "INFO, hged::main, last position: " << last_pos << std::endl;
    last_pos = variant_positions.back();
    std::cout << "INFO, hged::main, total variant position: " << variant_positions.size() << std::endl;
}

void get_linear_backbone(const std::string &fasta_file, const int &chr, int &start_pos, int &last_pos, std::string &backbone_seq)
{

    // Extract linear backbone using samtools
    // seed random generator by time in seconds (this may create issue if two instances are launched at the same time)
    srand(time(0));
    int random = rand() % 100000;
    std::string tmp_file = "linear_bc_chr" + std::to_string(chr) + "_" + std::to_string(random) + ".txt";
    std::string cmd = std::string(TOSTRING(SAMTOOLSPATH)) + " faidx " + fasta_file + " " + std::to_string(chr) + ":" + std::to_string(start_pos) + "-" + std::to_string(last_pos) + " >  " + tmp_file;
    std::cout << "INFO, hged::main, get linear backbone genome from fasta file using command: " << cmd << std::endl;
    std::system(cmd.c_str());

    std::ifstream file(tmp_file);
    std::string line;

    // ignore fist line and remove whitespaaces
    std::getline(file, line); // ignore the first line

    // store the rest of file into a string named content
    std::string content((std::istreambuf_iterator<char>(file)), (std::istreambuf_iterator<char>()));

    // using the erase, remove_if, and ::isspace functions --> remove all the whitespaces from the string
    content.erase(std::remove_if(content.begin(), content.end(), ::isspace),
                  content.end());
    backbone_seq = content;
}

void obtain_substrings(std::string &backbone_seq, int &start_pos, const int &alpha, tuple_list &tl, const int &chr)
{

    // seed random generator by time in seconds (this may create issue if two instances are launched at the same time)
    srand(time(0));
    int random = rand() % 100000;
    std::string pos_sub = "pos_sub_alpha_" + std::to_string(alpha) + "_chr" + std::to_string(chr) + "_" + std::to_string(random) + ".txt";
    std::cout << "INFO, hged::obtain_substring, extracting substrings per position and save to the file named: " << pos_sub << std::endl;

    std::ofstream pos_sub_file(pos_sub);

    std::map<int, std::set<std::string>> pos_substrings;

    int final_pos = start_pos + backbone_seq.size();
    std::string sample;
    std::string s;

    for (const auto &i : tl)
    {
        for (std::vector<value_elemet>::size_type j = 0; j < i.second.size(); ++j)
        {
            sample = boost::get<2>(i.second.at(j));
            std::string GT = boost::get<3>(i.second.at(j));

            for (int hap = 0; hap < 2; hap++)
            {
                if ((hap == 0 && (GT.compare(0, 1, "0") != 0)) || (hap == 1 && (GT.compare(2, 1, "0")) != 0))
                {
                    int current_pos = i.first;
                    std::string s = " ";
                    bool sample_found_at_pos;
                    std::string alt_choice;

                    while (current_pos < final_pos + 1 && s.size() < alpha + 1)
                    {   
                        std::map<int, Value_type>::iterator i2 = tl.find(current_pos); // find the index of the current pos in map
                        if (i2 != tl.end()) // the position exists in the  map
                        {
                            sample_found_at_pos = false;
                            for (std::vector<value_elemet>::size_type t = 0; t < i2->second.size(); ++t)
                            {
                                std::string sam = boost::get<2>(i2->second.at(t));
                                if (sam.compare(sample) == 0)
                                {
                                    std::string gt = boost::get<3>(i2->second.at(t));
                                    if (hap == 0) alt_choice = gt.substr(0, 1);
                                    if (hap == 1) alt_choice = gt.substr(2, 1);
                                    if (alt_choice.compare("0") != 0)
                                    {
                                        std::vector<std::string> alt_splited;
                                        boost::split(alt_splited, boost::get<1>(i2->second.at(t)), boost::is_any_of(",")); // split ALT over comma, save to a vector
                                        s += alt_splited[atoi(alt_choice.c_str()) - 1];
                                        current_pos += (boost::get<0>(i2->second.at(t))).size(); // add reference length
                                        sample_found_at_pos = true;
                                    }
                                }
                            }

                            if (sample_found_at_pos == false)
                            { // current pos is a variant position but sample not found
                                s += backbone_seq.substr(current_pos - start_pos, 1);
                                current_pos += 1;
                            }
                        }
                        else // current pos is not a variant position
                        {
                            s += backbone_seq.substr(current_pos - start_pos, 1);
                            current_pos += 1;
                        }
                    }
                    // add to this positions set
                    if (s.size() < alpha + 1)
                    {
                        pos_substrings[i.first].insert(s);
                    }
                    else
                    {
                        pos_substrings[i.first].insert(s.substr(0, alpha + 1));
                    }
                }
            }
        }

        // write to file
        pos_sub_file << std::to_string(i.first - start_pos);
        for (const auto &element : pos_substrings[i.first])
            pos_sub_file << +" " << element;
        pos_sub_file << std::endl;
    }
}

void construct_graph(std::string &backbone_seq, int &start_pos, int &last_pos, tuple_list &tl, std::vector<int> &variant_positions, const int &chr)
{

    // seed random generator by time in seconds (this may create issue if two instances are launched at the same time)
    srand(time(0));
    int random = rand() % 100000;
    std::string graph_file = "graph_chr" + std::to_string(chr) + "_" + std::to_string(random) + ".txt";
    std::cout << "INFO, hged::construct graph, constructing graph file: " << graph_file << std::endl;

    std::ofstream graph_file_name(graph_file);

    for (int i = 0; i < backbone_seq.size(); ++i)
    {
        // add backbone edges
        graph_file_name << std::to_string(start_pos + i - start_pos) + " " + std::to_string(start_pos + i + 1 - start_pos) + " " + backbone_seq[i] + " " + "-" << std::endl;
    }

    // add alt paths
    int start, end, new_vertex;
    new_vertex = last_pos + 2;
    int index = 0;
    std::string alt;
    std::string ref;
    std::map<int, Value_type>::iterator it = tl.begin();

    for (auto i = tl.begin(); i != tl.end(); i++)
    {

        // get the index of the variant position POS i.first
        auto it = find(variant_positions.begin(), variant_positions.end(), i->first);
        int index = it - variant_positions.begin();

        alt = boost::get<1>(i->second.at(0));
        ref = boost::get<0>(i->second.at(0));

        std::vector<std::string> alt_split_by_comma;
        boost::split(alt_split_by_comma, alt, boost::is_any_of(","));

        for (int elm = 0; elm < alt_split_by_comma.size(); elm++)
        {
            const char *alt_splited = alt_split_by_comma[elm].c_str(); // get characters of each alt

            for (size_t k = 0; k < strlen(alt_splited); ++k) // for each element in alt
            {
                if (k == 0){
                    start = i->first;
                }
                else{
                    start = new_vertex;
                }

                if (k == (strlen(alt_splited) - 1)) {
                    end = i->first + ref.size();
                }
                else{
                    new_vertex += 1;
                    end = new_vertex;
                }
                graph_file_name << std::to_string(start - start_pos) + " " + std::to_string(end - start_pos) + " " + alt_splited[k] + " " + std::to_string(index) << std::endl;
            }
        }
    }
}

int main(int argc, char **argv)
{

    // parse command line arguments
    Parameters parameters;
    parseandSave_ILP(argc, argv, parameters);

    tuple_list tl;
    std::ofstream pos_sub_file, graph_file_name;
    int start_pos, last_pos, pos, index;
    std::vector<std::string> samples;
    std::map<int, std::string> Sample_GT;
    std::vector<int> variant_positions;
    std::string backbone_seq;

    extract_vcf_pos_ref_alt_sample_gt(parameters.vcffile, parameters.fasta_ref_file, parameters.chr, start_pos,
                                      last_pos, variant_positions, tl);

    last_pos = last_pos + (boost::get<0>(tl[last_pos].at(0))).size();
    get_linear_backbone(parameters.fasta_ref_file, parameters.chr, start_pos, last_pos, backbone_seq);
    obtain_substrings(backbone_seq, start_pos, parameters.alpha, tl, parameters.chr);
    construct_graph(backbone_seq, start_pos, last_pos, tl, variant_positions, parameters.chr);

    return 0;
}
