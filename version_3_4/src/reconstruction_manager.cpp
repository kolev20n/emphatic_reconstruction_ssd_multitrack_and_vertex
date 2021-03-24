#include <iostream>
#include <vector>
//#include <map>
//#include <string>
//#include <sstream>
//#include <stdlib.h>
//#include <cmath>

#include <TFile.h>
#include <TTree.h>

//#include <TMath.h>
//#include <Math/Vector3D.h>
//#include <Math/GenVector/Rotation3D.h>
//#include <Math/GenVector/RotationX.h>
//#include <Math/GenVector/RotationY.h>
//#include <Math/GenVector/RotationZ.h>

#include "reconstruction_manager.h"
#include "input_output_manager.h"
#include "cluster_manager.h"
#include "transformation_manager.h"
#include "track_finding_manager.h"
//#include "histogram_manager.h"
#include "global_data_dispatcher.h"
#include "partial_geometry.h"

using namespace std;
using namespace ROOT::Math;

reconstruction_manager::reconstruction_manager(input_output_manager* an_input_output_manager)
:the_input_output_manager(an_input_output_manager), the_cluster_manager(0), the_transformation_manager(0), the_geometry(0), the_input_root_file(0), the_input_ssd_energy_tree(0), the_input_particle_truth_tree(0), the_input_steps_truth_tree(0), the_input_root_data_structure(0),
    the_track_finding_manager(0), the_event_characteristics(0)
{
  cout << "Created the_reconstruction_manager." << endl;
  
  the_geometry = the_global_data_dispatcher->get_partial_geometry();
  
  the_input_root_data_structure = new input_root_data_structure();
  the_global_data_dispatcher->register_input_root_data_structure(the_input_root_data_structure);
  
  the_event_characteristics = new event_characteristics();
  the_global_data_dispatcher->register_event_characteristics(the_event_characteristics);
  
  the_input_output_manager->initialize_root_trees(the_input_root_file, the_input_ssd_energy_tree, the_input_particle_truth_tree, the_input_steps_truth_tree, the_input_root_data_structure);
  
  the_cluster_manager = new cluster_manager();
  
  the_transformation_manager = new transformation_manager(the_geometry);
  
  if (!the_global_data_dispatcher->get_reconstruction_options()->use_magnetic_field &&
      the_global_data_dispatcher->get_partial_geometry()->get_number_of_tracking_regions() == 2 &&
      the_global_data_dispatcher->get_partial_geometry()->get_number_of_triple_groups() == 1)
  {
    the_algorithm = 1;
  }
  
  the_track_finding_manager = new track_finding_manager(the_geometry, the_algorithm);
}

void reconstruction_manager::process_events()
{
  bool clustering_success;
  bool track_finding_success;
  
  int number_of_events_to_process = the_input_output_manager->get_number_of_events_to_process();
  int actual_number_of_plates = the_geometry->get_number_of_plates();;

  hit_pair dummy_hit_pair;
  vector<hit_pair> hit_pair_list[max_number_of_plates];
  vector<double> clusters[max_number_of_plates];
  
  vector<line_3d> cluster_lines[max_number_of_plates];
  line_3d dummy_line;
  
  vector<track> the_tracks;
  
  cout << "Reconstruction options: " << endl;
  cout << "max_clusters_on_upstream_plates_with_multiple_clusters: " << the_global_data_dispatcher->get_reconstruction_options()->max_clusters_on_upstream_plates_with_multiple_clusters << endl;
  cout << "max_upstream_plates_with_multiple_clusters: " << the_global_data_dispatcher->get_reconstruction_options()->max_upstream_plates_with_multiple_clusters << endl;
  cout << "max_midstream_plus_downstream_plates_with_multiple_clusters: " << the_global_data_dispatcher->get_reconstruction_options()->max_midstream_plus_downstream_plates_with_multiple_clusters << endl;
  cout << "max_midstream_plus_downstream_plates_with_zero_clusters: " << the_global_data_dispatcher->get_reconstruction_options()->max_midstream_plus_downstream_plates_with_zero_clusters << endl;
  cout << "max_midstream_plus_downstream_plates_with_zero_or_multiple_clusters: " << the_global_data_dispatcher->get_reconstruction_options()->max_midstream_plus_downstream_plates_with_zero_or_multiple_clusters << endl;
  cout << "max_midstream_plus_downstream_plates_with_less_than_2_clusters: " << the_global_data_dispatcher->get_reconstruction_options()->max_midstream_plus_downstream_plates_with_less_than_2_clusters << endl;
  
  cout << "Processing events..." << endl;

  for (int i = 0; i < number_of_events_to_process; i++)
    //for (int i = 0; i < 1; i++)
  {
    for (int j = 0; j < max_number_of_plates; j++)
    {
      cluster_lines[j].clear();
    }
    
    for (int j = 0; j < max_number_of_plates; j++)
    {
      hit_pair_list[j].clear();
      clusters[j].clear();
    }
    
    if ((i + 1) % 1000 == 0)
    {
      cout << i + 1 << " events processed." << endl;
    }
    
    the_input_ssd_energy_tree->GetEntry(i);
    the_input_particle_truth_tree->GetEntry(i);
    the_input_steps_truth_tree->GetEntry(i);
    
    // remove this if geant gives cm
    for (int j = 0; j < the_input_root_data_structure->x_vertex->size(); j++)
    {
      the_input_root_data_structure->x_vertex->at(j) *= mm;
      the_input_root_data_structure->y_vertex->at(j) *= mm;
      the_input_root_data_structure->z_vertex->at(j) *= mm;
    }
    
    // ditto
    for (int j = 0; j < the_input_root_data_structure->step_track_id->size(); j++)
    {
      the_input_root_data_structure->step_x_i->at(j) *= mm;
      the_input_root_data_structure->step_y_i->at(j) *= mm;
      the_input_root_data_structure->step_z_i->at(j) *= mm;
      
      the_input_root_data_structure->step_x_f->at(j) *= mm;
      the_input_root_data_structure->step_y_f->at(j) *= mm;
      the_input_root_data_structure->step_z_f->at(j) *= mm;
    }
    
    for (int j = 0; j < the_input_root_data_structure->plate_number->size(); j++)
    {
      dummy_hit_pair.strip_number = the_input_root_data_structure->strip_number->at(j);
      dummy_hit_pair.energy_value = the_input_root_data_structure->total_energy_values->at(j);
      
      //if (dummy_hit_pair.energy_value > 0.05)
      {
        hit_pair_list[the_input_root_data_structure->plate_number->at(j)].push_back(dummy_hit_pair);
      }
    }
    
    clustering_success = the_cluster_manager->cluster(hit_pair_list, clusters, actual_number_of_plates, the_global_data_dispatcher->get_reconstruction_options()->cluster_interpretation_method);
    
    if (!clustering_success)
    {
      cout << "Clustering failed. This is a serious problem and no reasonable reconstruction can be done. Quitting..." << endl;
      exit(EXIT_FAILURE);
    }
    
    the_transformation_manager->transform_to_global_tracking_frame(clusters, cluster_lines);

    the_tracks.clear();

    track_finding_success = the_track_finding_manager->find_tracks(cluster_lines, the_tracks);
    
  } // end loop over events
  
  return;
}
