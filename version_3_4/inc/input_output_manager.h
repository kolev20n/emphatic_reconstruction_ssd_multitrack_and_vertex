#ifndef input_output_manager_h
#define input_output_manager_h

#include <string>

#include "common_data_structures_and_functions.h"

class TFile;
class TTree;
class partial_geometry;

//class histogram_manager;

class input_output_manager
{
public:
  input_output_manager();
  ~input_output_manager();
  bool initialize_files(int argc, char** argv);
  bool initialize_geometry(partial_geometry* the_geometry);
  bool initialize_root_trees(TFile*& an_input_root_file, TTree*& an_input_ssd_energy_tree, TTree*& an_input_particle_truth_tree, TTree*& an_input_steps_truth_tree, input_root_data_structure*& the_input_root_data_structure);
  int get_number_of_events_to_process();
  // int get_number_of_plates();
  reconstruction_options get_reconstruction_options();
  
private:
  std::string output_message;
  
  // bool use_magnetic_field;
  std::string magnetic_field_file_name;
  std::string geometry_file_name;
  std::string input_root_file_name;
  TFile* the_input_root_file;
  TTree* the_input_ssd_energy_tree;
  TTree* the_input_particle_truth_tree;
  TTree* the_input_steps_truth_tree;
  
  std::string output_root_file_name;

  std::string input_root_ssd_energy_tree_name;
  std::string input_root_particle_truth_tree_name;
  std::string input_root_steps_truth_tree_name;

  int number_of_events_to_process;
  
  reconstruction_options* the_reconstruction_options;
};

#endif
