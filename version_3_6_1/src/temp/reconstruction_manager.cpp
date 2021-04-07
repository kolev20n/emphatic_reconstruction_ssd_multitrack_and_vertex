#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <cmath>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <Math/Vector3D.h>
#include <Math/GenVector/Rotation3D.h>
#include <Math/GenVector/RotationX.h>
#include <Math/GenVector/RotationY.h>
#include <Math/GenVector/RotationZ.h>

#include "reconstruction_manager.h"
#include "input_output_manager.h"
#include "cluster_manager.h"
#include "transformation_manager.h"
#include "track_finding_manager.h"
#include "histogram_manager.h"
#include "global_data_dispatcher.h"
#include "partial_geometry.h"

using namespace std;
using namespace ROOT::Math;

reconstruction_manager::reconstruction_manager(input_output_manager* an_input_output_manager)
:the_input_output_manager(an_input_output_manager), the_cluster_manager(0), the_transformation_manager(0), the_geometry(0), the_input_root_file(0), the_input_ssd_energy_tree(0), the_input_particle_truth_tree(0), the_input_steps_truth_tree(0), the_input_root_data_structure(0),
    the_track_finding_manager(0)
{
  cout << "Created the_reconstruction_manager." << endl;
  
  the_geometry = new partial_geometry();
  
  the_global_data_dispatcher->register_geometry(the_geometry);
  
  the_input_output_manager->initialize_geometry(the_geometry);
  
  cout << "Info on geometry:????????????" << endl;
  
  the_input_root_data_structure = new input_root_data_structure();
  the_global_data_dispatcher->register_input_root_data_structure(the_input_root_data_structure);
  
  the_input_output_manager->initialize_root_trees(the_input_root_file, the_input_ssd_energy_tree, the_input_particle_truth_tree, the_input_steps_truth_tree, the_input_root_data_structure);
  
  the_cluster_manager = new cluster_manager();
  
  the_transformation_manager = new transformation_manager(the_geometry);
  
  the_track_finding_manager = new track_finding_manager(the_geometry);
}

void reconstruction_manager::process_events()
{
  int number_of_events_to_process = the_input_output_manager->get_number_of_events_to_process();
  
  int cluster_multiplicity[max_number_of_plates];
  
  hit_pair dummy_hit_pair;
  vector<hit_pair> hit_pair_list[max_number_of_plates];
  vector<double> clusters[max_number_of_plates];
  
  int dummy_int;
  
  double dummy_double;
  
  // every cluster is determined by the center of the line segment
  // parallel to the strips where the cluster is (a point) and the direction vector
  // these need to be set to lab coordinates by the transformation_manager
  
  vector<silicon_strip> cluster_lines[max_number_of_plates];
  silicon_strip dummy_strip;
  
  bool clustering_success;
  bool track_finding_success;
  
  int actual_number_of_plates;
 
  string dummy_string;
  stringstream dummy_stream;
  
  map<vector<int>, int> cluster_distribution_map;
  vector<int> cluster_distribution_vector;
  vector<pair<string, int>> sorted_cluster_distributions;
  
  map<vector<int>, int> upstream_cluster_distribution_map;
  vector<int> upstream_cluster_distribution_vector;
  vector<pair<string, int>> sorted_upstream_cluster_distributions;
  
  map<vector<int>, int> midstream_cluster_distribution_map;
  vector<int> midstream_cluster_distribution_vector;
  vector<pair<string, int>> sorted_midstream_cluster_distributions;
  
  map<vector<int>, int> downstream_cluster_distribution_map;
  vector<int> downstream_cluster_distribution_vector;
  vector<pair<string, int>> sorted_downstream_cluster_distributions;
  
  
  
  
  
  int min_clusters;
  int plate_index;
  map<vector<int>, int> group_cluster_distribution_map;
  vector<int> group_cluster_multiplicity_vector;
  vector<pair<string, int>> sorted_group_cluster_distributions;
  
  
  
  
  
  
  vector<track> the_tracks;
  the_tracks.clear();
  
  actual_number_of_plates = the_geometry->get_number_of_plates();
  the_output_data_structure->actual_number_of_plates = actual_number_of_plates;
  
  if (!determine_algorithm_type(the_geometry))
  {
    cout << "The geometry does not allow any know reconstruction algorithm. Quitting..." << endl;
    exit(EXIT_FAILURE);
  }
  else
  {
    cout << "Will use reconstruction algorithm: ";
    
    dummy_int = the_global_data_dispatcher->get_reconstruction_options()->algorithm_type;
    
    if (dummy_int == 1)
    {
      cout << "Single Triplet of Plates" << endl;
    }
  }
  
  cout << "Processing events..." << endl;
  
  cluster_distribution_map.clear();
  sorted_cluster_distributions.clear();
  
  upstream_cluster_distribution_map.clear();
  sorted_upstream_cluster_distributions.clear();
  
  midstream_cluster_distribution_map.clear();
  sorted_midstream_cluster_distributions.clear();
  
  downstream_cluster_distribution_map.clear();
  sorted_downstream_cluster_distributions.clear();
  
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
    
    cluster_distribution_vector.clear();
    upstream_cluster_distribution_vector.clear();
    midstream_cluster_distribution_vector.clear();
    downstream_cluster_distribution_vector.clear();
    group_cluster_multiplicity_vector.clear();
    
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
    
    the_output_data_structure->same_vertex_for_high_energy_particles = true;
    
    // cout << "----------------" << endl;
    for (int j = 0; j < the_input_root_data_structure->track_id->size(); j++)
    {
      dummy_double = -1000.;
      
      if (the_input_root_data_structure->parent_track_id->at(j) == 1 && the_input_root_data_structure->ekin_vertex->at(j) > 500.)
      {
        /*
        cout << the_input_root_data_structure->track_id->at(j) << " "
             << the_input_root_data_structure->particle_code->at(j) << " "
             << the_input_root_data_structure->ekin_vertex->at(j) << " "
             << the_input_root_data_structure->x_vertex->at(j) << " "
             << the_input_root_data_structure->y_vertex->at(j) << " "
             << the_input_root_data_structure->z_vertex->at(j) << endl;
        */
        
        if (dummy_double < -100.)
        {
          dummy_double = the_input_root_data_structure->z_vertex->at(j);
          the_output_data_structure->mc_guessed_vertex.SetCoordinates(the_input_root_data_structure->x_vertex->at(j), the_input_root_data_structure->y_vertex->at(j), the_input_root_data_structure->z_vertex->at(j));
        }
        else
        {
          if (the_input_root_data_structure->z_vertex->at(j) != dummy_double)
          {
            // cout << dummy_double << " " << the_input_root_data_structure->z_vertex->at(j) << endl;
            // cout << "aaaaa" << endl;
            the_output_data_structure->same_vertex_for_high_energy_particles = false;
            
          }
        }
      }
    }
    
    cout << "MC guessed vertex: " << the_output_data_structure->mc_guessed_vertex << endl;
    
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
    
    for (int j = 0; j < actual_number_of_plates; j++)
    {
      cluster_distribution_vector.push_back(clusters[j].size());
    }
    
    for (int j = 0; j < the_geometry->get_number_of_upstream_plates(); j++)
    {
      upstream_cluster_distribution_vector.push_back(clusters[j].size());
    }
    
    for (int j = the_geometry->get_number_of_upstream_plates(); j < (the_geometry->get_number_of_upstream_plates() + the_geometry->get_number_of_midstream_plates()); j++)
    {
      midstream_cluster_distribution_vector.push_back(clusters[j].size());
    }
    
    for (int j = (the_geometry->get_number_of_upstream_plates() + the_geometry->get_number_of_midstream_plates()); j < actual_number_of_plates; j++)
    {
      downstream_cluster_distribution_vector.push_back(clusters[j].size());
    }
    
    cluster_distribution_map[cluster_distribution_vector]++;
    
    upstream_cluster_distribution_map[upstream_cluster_distribution_vector]++;
    
    midstream_cluster_distribution_map[midstream_cluster_distribution_vector]++;
    
    downstream_cluster_distribution_map[downstream_cluster_distribution_vector]++;
    
    the_transformation_manager->transform_to_lab_frame(clusters, cluster_lines);
    
    for (int j = 0; j < actual_number_of_plates; j++)
    {
      the_output_data_structure->cluster_lines[j].clear();

      for (int k = 0; k < cluster_lines[j].size(); k++)
      {
        the_output_data_structure->cluster_lines[j].push_back(cluster_lines[j].at(k));
      }
    }

    
    
    
    
    for (int j = 0; j < the_geometry->the_plate_groups.size(); j++)
    {
      min_clusters = 1000;
      
      for (int k = 0; k < the_geometry->the_plate_groups.at(j).plate_index_list.size(); k++)
      {
        plate_index = the_geometry->the_plate_groups.at(j).plate_index_list.at(k);
        
        if (cluster_lines[plate_index].size() < min_clusters)
        {
          min_clusters = cluster_lines[plate_index].size();
        }
      }
      
      if (min_clusters == 1000) min_clusters = 0;
      
      group_cluster_multiplicity_vector.push_back(min_clusters);
      
      // cout << min_clusters << " ";
    }
    
    group_cluster_distribution_map[group_cluster_multiplicity_vector]++;
    
    
    
    /*
    // temporary debugging spurious clusters
    //
    cout << "----------------------------" << endl;
    cout << "Clusters:" << endl;
    for (int j = 0; j < actual_number_of_plates; j++)
    {
      for (int k = 0; k < cluster_lines[j].size(); k++)
      {
        cout << j << " " << cluster_lines[j].at(k).line_center << " " << cluster_lines[j].at(k).line_direction << endl;
      }
    }
    
    cout << "-----" << endl;
    cout << "Contributing tracks:" << endl;
    for (int j = 0; j < the_input_root_data_structure->plate_number->size(); j++)
    {
      cout << the_input_root_data_structure->plate_number->at(j) << " ";
      for (int k = 0; k < the_input_root_data_structure->contributing_tracks->at(j).size(); k++)
      {
        cout << the_input_root_data_structure->contributing_tracks->at(j).at(k) << " ";
      }
      
      cout << endl;
    }
    
    cout << "-----" << endl;
    cout << "Tracks:" << endl;
    for (int j = 0; j < the_input_root_data_structure->track_id->size(); j++)
    {
      cout << the_input_root_data_structure->track_id->at(j) << " "
           << the_input_root_data_structure->particle_code->at(j) << " "
      << the_input_root_data_structure->ekin_vertex->at(j) << " "
      << the_input_root_data_structure->creator_process->at(j) << endl;
    }
    
    // end temporary debugging spurious clusters
    
    */
    
    
    the_tracks.clear();
    
    track_finding_success = the_track_finding_manager->find_tracks(cluster_lines, the_tracks);
    
    the_histogram_manager->fill(the_input_root_data_structure);
    
  } // end loop over events
  
  the_histogram_manager->save();
  
  for (map<vector<int>, int>::iterator it = cluster_distribution_map.begin(); it != cluster_distribution_map.end(); ++it)
  {
    dummy_string = "";
    dummy_stream.str("");
    
    //if (it->second >= 5)
    {
      for (int j = 0; j < actual_number_of_plates; j++)
      {
        dummy_stream << it->first.at(j) << " ";
      }
      dummy_string = dummy_stream.str();
      
      sorted_cluster_distributions.push_back(make_pair(dummy_string, it->second));
    }
  }
  
  sort(sorted_cluster_distributions.begin(), sorted_cluster_distributions.end(), compare_map);
  
  cout << "All plates: " << endl;
  for (int i = 0; i < sorted_cluster_distributions.size(); i++)
  {
    cout << sorted_cluster_distributions.at(i).first << "=> " << sorted_cluster_distributions.at(i).second << endl;
  }
  
  for (map<vector<int>, int>::iterator it = upstream_cluster_distribution_map.begin(); it != upstream_cluster_distribution_map.end(); ++it)
  {
    dummy_string = "";
    dummy_stream.str("");
    
    //if (it->second >= 5)
    {
      for (int j = 0; j < the_geometry->get_number_of_upstream_plates(); j++)
      {
        dummy_stream << it->first.at(j) << " ";
      }
      dummy_string = dummy_stream.str();
      
      sorted_upstream_cluster_distributions.push_back(make_pair(dummy_string, it->second));
    }
  }
  
  sort(sorted_upstream_cluster_distributions.begin(), sorted_upstream_cluster_distributions.end(), compare_map);
  
  cout << "Upstream plates: " << endl;
  
  for (int i = 0; i < sorted_upstream_cluster_distributions.size(); i++)
  {
    cout << sorted_upstream_cluster_distributions.at(i).first << "=> " << sorted_upstream_cluster_distributions.at(i).second << endl;
  }
  
  cout << the_geometry->get_number_of_upstream_plates() << " "
  << the_geometry->get_number_of_midstream_plates() << " "
  << the_geometry->get_number_of_downstream_plates() << endl;
  
  for (map<vector<int>, int>::iterator it = midstream_cluster_distribution_map.begin(); it != midstream_cluster_distribution_map.end(); ++it)
  {
    dummy_string = "";
    dummy_stream.str("");
    
    //if (it->second >= 5)
    {
      for (int j = 0; j < the_geometry->get_number_of_midstream_plates(); j++)
      {
        dummy_stream << it->first.at(j) << " ";
      }
      dummy_string = dummy_stream.str();
      
      sorted_midstream_cluster_distributions.push_back(make_pair(dummy_string, it->second));
    }
  }
  
  sort(sorted_midstream_cluster_distributions.begin(), sorted_midstream_cluster_distributions.end(), compare_map);
  
  cout << "Midstream plates: " << endl;
  
  for (int i = 0; i < sorted_midstream_cluster_distributions.size(); i++)
  {
    cout << sorted_midstream_cluster_distributions.at(i).first << "=> " << sorted_midstream_cluster_distributions.at(i).second << endl;
  }
  
  for (map<vector<int>, int>::iterator it = downstream_cluster_distribution_map.begin(); it != downstream_cluster_distribution_map.end(); ++it)
  {
    dummy_string = "";
    dummy_stream.str("");
    
    //if (it->second >= 5)
    {
      for (int j = 0; j < the_geometry->get_number_of_downstream_plates(); j++)
      {
        dummy_stream << it->first.at(j) << " ";
      }
      dummy_string = dummy_stream.str();
      
      sorted_downstream_cluster_distributions.push_back(make_pair(dummy_string, it->second));
    }
  }
  
  sort(sorted_downstream_cluster_distributions.begin(), sorted_downstream_cluster_distributions.end(), compare_map);
  
  cout << "Downstream plates: " << endl;
  
  for (int i = 0; i < sorted_downstream_cluster_distributions.size(); i++)
  {
    cout << sorted_downstream_cluster_distributions.at(i).first << "=> " << sorted_downstream_cluster_distributions.at(i).second << endl;
  }
  
  
  
  
  
  
  for (map<vector<int>, int>::iterator it = group_cluster_distribution_map.begin(); it != group_cluster_distribution_map.end(); ++it)
  {
    dummy_string = "";
    dummy_stream.str("");
    
    //if (it->second >= 5)
    {
      for (int j = 0; j < the_geometry->the_plate_groups.size(); j++)
      {
        dummy_stream << it->first.at(j) << " ";
      }
      dummy_string = dummy_stream.str();
      
      sorted_group_cluster_distributions.push_back(make_pair(dummy_string, it->second));
    }
  }
  
  sort(sorted_group_cluster_distributions.begin(), sorted_group_cluster_distributions.end(), compare_map);
  
  cout << "Groups: " << endl;
  for (int i = 0; i < sorted_group_cluster_distributions.size(); i++)
  {
    cout << sorted_group_cluster_distributions.at(i).first << "=> " << sorted_group_cluster_distributions.at(i).second << endl;
  }
  
  
  
  
  
  
  
  return;
}

bool reconstruction_manager::determine_algorithm_type(partial_geometry* the_geometry)
{
  the_global_data_dispatcher->get_reconstruction_options()->triplet_groups.clear();
  
  the_global_data_dispatcher->get_reconstruction_options()->algorithm_type = 0;
  
  for (int i = 0; i < the_geometry->the_plate_groups.size(); i++)
  {
    if (the_geometry->the_plate_groups.at(i).group_type == 3 || the_geometry->the_plate_groups.at(i).group_type == 4)
    {
      the_global_data_dispatcher->get_reconstruction_options()->triplet_groups.push_back(i);
    }
  }
  
  if (the_global_data_dispatcher->get_reconstruction_options()->triplet_groups.size() == 1)
  {
    // 0 no algorithm found
    // 1 single triplet group (normally with diagonal plate)
    //
    the_global_data_dispatcher->get_reconstruction_options()->algorithm_type = 1;
  }
  
  if (the_global_data_dispatcher->get_reconstruction_options()->algorithm_type == 0)
  {
    cout << "No viable reconstruction algorithm was found. Cannot proceed with any reasonalbe reconstruciton. Quitting..." << endl;
    exit(EXIT_FAILURE);
  }
  
  return true;
}

