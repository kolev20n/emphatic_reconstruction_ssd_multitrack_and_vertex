#ifndef histogram_manager_h
#define histogram_manager_h

#include <string>

#include "common_data_structures_and_functions.h"


class TTree;
class TFile;
class TH1D;

class histogram_manager
{
public:
  histogram_manager(std::string an_output_file_name, output_data_structure* an_output_data_structure);
  ~histogram_manager();
  bool book();
  bool fill(input_root_data_structure* the_input_root_data_structure);
  bool save();

private:
  std::string output_message;
  
  std::string output_root_file_name;
  
  TH1D* h101;

  TFile* out_file;
  output_data_structure* the_output_data_structure;
  TTree* output_tree;
  output_root_data_structure* the_output_root_data_structure;
};

#endif
