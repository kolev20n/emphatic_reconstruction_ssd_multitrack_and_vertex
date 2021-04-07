#ifndef reconstruction_manager_h
#define reconstruction_manager_h

#include "common_data_structures_and_functions.h"

class TTree;
class TFile;

class input_output_manager;
class cluster_manager;
class transformation_manager;
class track_finding_manager;
class partial_geometry;

class reconstruction_manager
{
public:
  reconstruction_manager(input_output_manager* an_input_output_manager);
  void process_events();

private:
  input_output_manager* the_input_output_manager;
  cluster_manager* the_cluster_manager;
  transformation_manager* the_transformation_manager;
  partial_geometry* the_geometry;
  TFile* the_input_root_file;
  TTree* the_input_ssd_energy_tree;
  TTree* the_input_particle_truth_tree;
  TTree* the_input_steps_truth_tree;
  input_root_data_structure* the_input_root_data_structure;
  run_output_data_structure* the_run_output_data_structure;
  event_mc_output_data_structure* the_event_mc_output_data_structure;
  event_reco_output_data_structure* the_event_reco_output_data_structure;

  track_finding_manager* the_track_finding_manager;
  event_characteristics* the_event_characteristics;
  
  int the_algorithm;
  // 1: no field; 1 triplet group (should be in region 2, as being
  // in region 1 makes no sense); 2 regions
};

#endif
