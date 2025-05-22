#include "constraint_large_dollo_parsimony.h"
#include<cstdio>
#include<iostream>
#include<fstream>
#include<vector>
#include <cstdint>
#include<cstdlib>
#include<unistd.h>
#include "binary_character_matrix.hpp"
using namespace std;

const std::string help = 
"===================================== Dollo-CDP =====================================\n"
"Dollo-CDP is a program that solves the large Dollo parsimony problem for binary\n"
"characters (missing values allowed) within a clade-constrained version of tree space.\n\n" 
"USAGE for Large Dollo problem:\n"
"./dollo-cdp -i <input characters file> -l <character BED file> -d <character linking threshold> -g <outgroup name> -o <output file>\n\n"
"USAGE for small Dollo parsimony problem:\n"
"./dollo-cdp -i <input characters file> -l <character BED file> -d <character linking threshold> -q <input species tree>\n\n"
"OPTIONS:\n"
"[-h|--help]\n"
"        Prints this help message.\n"
"(-i|--input) <input characters file>\n"
"        Name of file containing input characters in nexus format (required)\n"
"(-l|--location) <character BED file>\n"
"	 Name of the file containing BED entries for each character\n"
"(-d|--distance) <character linking threshold>\n"
"	 Maximum distance between characters for loss linking\n"
"(-x|-g|--outgroup) <outgroup name>\n"
"        Comma separated list of outgroup taxa used to root solution space\n"
"[(-t|--trees) <input trees file>]\n"
"        Name of file containing trees in newick format for constructing solution\n" 
"        space with ASTRAL. If no file is specified, characters will be treated as\n"
"        rooted trees with at most one internal branch.\n"
"[(-q) <input species file>]\n"
"        Name of file containing species trees in newick format\n"
"[(-k)]\n"
"        Write characters as rooted trees with one internal branch and exit\n"
"[(-o|--output) <output file>]\n"
"        Name of file for writing output species tree (default: stdout)\n\n"
"Contact: Post issue to Github (https://github.com/molloy-lab/Dollo-CDP)\n"
"         or email Junyan Dai (jdai1234@umd.edu) & Erin Molloy (ekmolloy@umd.edu)\n\n"
"If you use Dollo-CDP in your work, please cite:\n"
"  Dai, Rubel, Han, Molloy, 2024, \"Dollo-CDP: a polynomial-time algorithm for\n"
"  the clade-constrained large Dollo parsimony problem\", Algorithms for Molecular\n"
"  Biology, https://doi.org/10.1186/s13015-023-00249-9.\n"
"====================================================================================\n\n";





std::vector<std::string> get_outgroup(std::string oname) {
  std::vector<std::string> lat;
  boost::split(lat, oname, boost::is_any_of(","));
  return lat;
}


// For test reading clades purpose
void print_clades_set(clades_set cs, vector<string> &labels) {
  for (const auto& e: cs) {
    cout << e.to_labels(labels) << endl;
  }
}

std::tuple<int,int,int,std::string> split(std::string &line) {
  std::stringstream stream(line);
  std::string chr;
  std::string begin;
  std::string end;
  std::string name; 


  std::getline(stream,chr,'\t');
  std::getline(stream,begin,'\t');
  std::getline(stream,end,'\t');
  std::getline(stream,name,'\t');

  int num_chr = 0;
  int num_begin = stoi(begin);
  int num_end = stoi(end);
  if(chr == "X"){
    num_chr = 23;
  }
  else if(chr == "Y"){
    num_chr = 24;
  }
  else{
    num_chr = stoi(chr);
  }

  return make_tuple(num_chr,num_begin,num_end,name);
}


std::vector<std::tuple<int,int,int,std::string>> get_loc_vector(std::string &filename, int num_chars) {
	
  std::vector<std::tuple<int,int,int,std::string>> loc(num_chars);
  
  std::ifstream file(filename);
  
  if(!file) {
    std::cout << "Big problem" << std::endl;
    return {}; 
  }

  std::string line;
  int line_num = 0;

  while(std::getline(file,line)){
	  
    std::tuple<int,int,int,std::string> parts = split(line);
    loc[line_num] = (make_tuple(std::get<0>(parts),std::get<1>(parts),std::get<2>(parts),std::get<3>(parts)));
  
    line_num++;
  }
  
  return loc;
}


int main(int argc, char** argv) {
  std::cout << "Dollo-CDP version 1.0.1\nCOMMAND: ";
  for (int j = 0; j < argc; j ++) 
    std::cout << argv[j] << " ";
  std::cout << std::endl;


  if (argc == 1) {std::cout << help; return 0;}
  auto start = std::chrono::high_resolution_clock::now();
  string filename1 = "";
  string filename2 = "";
  string filename3 = "";
  string filename_bed = "";
  string input_tree = "";
  string outname;
  bool large_Dollo = true;
  bool write_bptrees_and_exit = false;
  bool user_defined_search_space = false;
  unsigned int d = 0;

  for (int i = 0; i < argc; i++) {
    string opt(argv[i]);
    if (opt == "-h" || opt == "--help") {std::cout << help; return 0;}
    if (opt == "-i" || opt == "--input" && i < argc - 1) filename1 = argv[++ i];
    if (opt == "-k") {write_bptrees_and_exit = true;} 
    if (opt == "-t" || opt == "--trees" && i < argc - 1) {user_defined_search_space = true; filename2 = argv[++ i];}
    if (opt == "-l" || opt == "--locations" && i < argc - 1) filename_bed = argv[++ i];
    if (opt == "-d" || opt == "--distance" && i < argc - 1) d = std::stoi(argv[++ i]);
    if (opt == "-q" && i < argc - 1) {large_Dollo = false; input_tree = argv[++ i];}
    if (opt == "-x" || opt == "--outgroup" || opt == "-g" && i < argc - 1) outname = argv[++ i];
    if (opt == "-o" || opt == "--output"&& i < argc - 1) filename3 = argv[++ i];
  }

  if (argc < 3) {
    std::cout << "Not enough arguments given!\n\n";
    std::cout << help;
    return 1;
  }

  if (filename1 == "") {
    std::cout << "Need to provide file name for input characters!\n\n";
    std::cout << help;
    return 1;
  }
  if (filename1[0] == '-') {
     std::cout << "Warning: May not have correctly specified file name for -i option!\n\n";
  }
  if (filename_bed == "" && large_Dollo) {
     std::cout << "Need to provide character BED file!\n\n";
     std::cout << help;
     return 1;
  }
  if (d == 0 && large_Dollo) {
     std::cout << "Need to provide a valid character linking distance!\n\n";
     std::cout << help;
     return 1;
  }
  if (outname == "" && large_Dollo) {
    std::cout << "Need to provide outgroup name!\n\n";
    std::cout << help;
    return 1;
  }
  if (input_tree == "" && !large_Dollo) {
    std::cout << "Need to provide file name for input species tree when using -q option!\n\n";
    std::cout << help;
    return 1;
  }
  if (input_tree[0] == '-') {
     std::cout << "Warning: May not have correctly specified file name for -q option!\n\n";
  }

  unsigned int k;

  boost::unordered_map<string, unsigned int> label2index;

  vector<string> labels;

  boost::unordered_map<std::string, std::string> record;

  uint8_t** C = read_characters(filename1, k, label2index, labels, record);
  
  auto read_characters_end = std::chrono::high_resolution_clock::now();
 
  std::vector<std::tuple<int,int,int,std::string>> loc = get_loc_vector(filename_bed, k);

  if (!large_Dollo) {
    
    //add new parameters here
    //loc table
    //d threshold
    unsigned int score = Dollo_parsimony_score(C, k, input_tree, label2index, loc, d);
    std::ofstream fout(filename3);
    fout << score << endl;
    std::cout << "The Dollo parsimony score: " << score << std::endl;
    auto small_dollo_end = std::chrono::high_resolution_clock::now();
    auto duration_of_small_dollo =  std::chrono::duration_cast<std::chrono::milliseconds>(small_dollo_end - start);
    std::cout << "execution time of small Dollo Parsimony: " << duration_of_small_dollo.count()  << "ms" << std::endl;
    return 0;
  }
   
  vector<string> outgroup = get_outgroup(outname);



  if (!user_defined_search_space && large_Dollo) {
    phylotools::BinaryCharacterMatrix S = phylotools::BinaryCharacterMatrix();
    
    std::ifstream matrix_file_stream(filename1);

    if (!matrix_file_stream.is_open()) {
      std::cerr << "Failed to open character matrix file ";
      return 1;
    }
    
    filename2 = filename1 + "-auto-generated.bptrees";

    std::ofstream generated_clades_file(filename2);
    
    if (!generated_clades_file.is_open()) {
      std::cerr << "Error: could not create generated_clades_file " << std::endl;
      return 1;
    }
    std::string format;
    if (filename1.length() > 4) format = filename1.substr(filename1.length() - 4);

    if (format == ".nex") {
      write_newick_from_C(generated_clades_file, k, labels, C);
    } else {
      std::cout << "Input character matrix must be in NEXUS format with extension .nex!\n\n";
      std::cout << help;
      return 1;
      //S.readMatrix(matrix_file_stream);
      //S.writeNewick(generated_clades_file);
    }
    matrix_file_stream.close();
    generated_clades_file.close();
  }

  if (write_bptrees_and_exit) return 0;

  clades_set X = read_search_space(filename2, label2index, labels, outgroup);

  if (!user_defined_search_space) {
    std::remove(filename2.c_str());
  }
  auto search_space_end = std::chrono::high_resolution_clock::now();
  
  opt f;
  states_map St;
  state_record st_r, st_l;
  child_record g;
  
  auto duration_of_characters =  std::chrono::duration_cast<std::chrono::milliseconds>(read_characters_end - start);
  std::cout << "execution time of reading characters matrix: " << duration_of_characters.count()  << "ms" << std::endl;
  auto duration_of_search_space =  std::chrono::duration_cast<std::chrono::milliseconds>(search_space_end - read_characters_end);
  std::cout << "execution time of computing search space: " << duration_of_search_space.count() << "ms" << std::endl;
  
 
  std::ofstream fout(filename3);
  Tree reso = constraint_large_dollo_parsimony(X, record, C, k, labels, label2index, loc, d,filename3);
  cout << reso.newick() << endl;
  fout << reso.newick() << endl;

  auto end = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  std::cout << "execution time: " << duration.count() << "ms" << std::endl;
  return 0;

  
  
}
