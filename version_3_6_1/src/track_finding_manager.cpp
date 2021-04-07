#include <iostream>
#include <utility>
#include <algorithm>

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnMigrad.h"

#include "track_finding_manager.h"
#include "partial_geometry.h"
#include "global_data_dispatcher.h"
#include "common_data_structures_and_functions.h"
#include "my_straight_line_fcn.h"
#include "my_straight_line_fcn_with_point.h"

using namespace std;
using namespace ROOT::Math;
using namespace ROOT::Minuit2;

track_finding_manager::track_finding_manager(partial_geometry* a_geometry, int an_algorithm)
:the_geometry(a_geometry), the_algorithm(an_algorithm)
{
  cout << "Created the_track_finding_manager." << endl;
}

bool track_finding_manager::find_tracks(vector<line_3d> cluster_lines[], vector<track>& the_tracks)
{
  bool track_finding_success = true;
  bool classify_event_sucess;
  bool fitting_success;
  
  event_characteristics the_event_characteristics;
  
  int event_class_by_cluster_multiplicity;
  
  double one_fit_goodness;
  double fit_goodness[100];
  track upstream_track, downstream_track, total_track;
  
  vector<track> downstream_tracks;
  
  if (the_algorithm == 1)
  {
    classify_event_sucess = classify_event(cluster_lines);
    
    if (!classify_event_sucess)
    {
      cout << "Couldn't classify event. This should not normally happen..." << endl;
      return false;
    }
    else
    {
      event_class_by_cluster_multiplicity = the_global_data_dispatcher->get_event_characteristics()->event_class_by_cluster_multiplicity;
    }
    
    if (event_class_by_cluster_multiplicity == 11 || event_class_by_cluster_multiplicity == 12)
    {
      if (event_class_by_cluster_multiplicity == 12)
      {
        assign_cluster_lines_single_track_upstream(cluster_lines);
        assign_cluster_lines_single_track_downstream(cluster_lines);
      }
      upstream_track.cluster_lines.clear();
      upstream_track.plate_index_of_cluster_line.clear();
      
      upstream_track.assumed_start_group = 0;
      upstream_track.assumed_end_group = the_global_data_dispatcher->get_partial_geometry()->get_number_of_groups_in_tracking_region(0) - 1;
      
      for (int i = 0; i < the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates_in_tracking_region(0); i++)
      {
        if (cluster_lines[i].size() > 0)
        {
          upstream_track.cluster_lines.push_back(cluster_lines[i].at(0));
          upstream_track.plate_index_of_cluster_line.push_back(i);
        }
      }

      fitting_success = fit_track(upstream_track);
      fit_goodness[0] = upstream_track.fit_goodness;
      
      downstream_track.cluster_lines.clear();
      downstream_track.plate_index_of_cluster_line.clear();

      downstream_track.assumed_start_group = the_global_data_dispatcher->get_partial_geometry()->get_number_of_groups_in_tracking_region(0);
      downstream_track.assumed_end_group = the_global_data_dispatcher->get_partial_geometry()->get_number_of_groups() - 1;
      
      for (int i = the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates_in_tracking_region(0); i < the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates(); i++)
      {
        if (cluster_lines[i].size() > 0)
        {
          downstream_track.cluster_lines.push_back(cluster_lines[i].at(0));
          downstream_track.plate_index_of_cluster_line.push_back(i);
        }
      }

      fitting_success = fit_track(downstream_track);
      fit_goodness[1] = downstream_track.fit_goodness;
      
      total_track.cluster_lines.clear();
      total_track.plate_index_of_cluster_line.clear();
      
      total_track.assumed_start_group = 0;
      total_track.assumed_end_group = the_global_data_dispatcher->get_partial_geometry()->get_number_of_groups() - 1;
      
      for (int i = 0; i < the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates(); i++)
      {
        if (cluster_lines[i].size() > 0)
        {
          total_track.cluster_lines.push_back(cluster_lines[i].at(0));
          total_track.plate_index_of_cluster_line.push_back(i);
        }
      }

      fitting_success = fit_track(downstream_track);
      fit_goodness[2] = downstream_track.fit_goodness;
      
      // the criterion here has to be worked out
      if (fit_goodness[2] > fit_goodness[0] + fit_goodness[1])
      {
        the_tracks.push_back(total_track);
        the_event_characteristics.event_class_by_reconstructed_tracks = 11;
      }
      else
      {
        the_tracks.push_back(upstream_track);
        the_tracks.push_back(downstream_track);
        the_event_characteristics.event_class_by_reconstructed_tracks = 10; // temporary
      }

      if (the_event_characteristics.event_class_by_reconstructed_tracks != 11)
      {
        find_vertex(the_tracks, the_global_data_dispatcher->get_event_characteristics()->the_vertex);
      }
      
      if (the_global_data_dispatcher->get_event_characteristics()->the_vertex.z() < the_global_data_dispatcher->get_partial_geometry()->get_target()->position.z() + the_global_data_dispatcher->get_partial_geometry()->get_target()->size.z() / 2. && the_global_data_dispatcher->get_event_characteristics()->the_vertex.z() > the_global_data_dispatcher->get_partial_geometry()->get_target()->position.z() - the_global_data_dispatcher->get_partial_geometry()->get_target()->size.z() / 2.)
      {
        the_event_characteristics.event_class_by_reconstructed_tracks = 12;
      }
      else
      {
        the_event_characteristics.event_class_by_reconstructed_tracks = 13;
      }
      
      the_global_data_dispatcher->get_event_characteristics()->event_class_by_reconstructed_tracks = the_event_characteristics.event_class_by_reconstructed_tracks;
    }
    else if (event_class_by_cluster_multiplicity == 21 || event_class_by_cluster_multiplicity == 22)
    {
      if (event_class_by_cluster_multiplicity == 22)
      {
        assign_cluster_lines_single_track_upstream(cluster_lines);
      }
      
      upstream_track.cluster_lines.clear();
      upstream_track.plate_index_of_cluster_line.clear();
      
      upstream_track.assumed_start_group = 0;
      upstream_track.assumed_end_group = the_global_data_dispatcher->get_partial_geometry()->get_number_of_groups_in_tracking_region(0) - 1;
      
      for (int i = 0; i < the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates_in_tracking_region(0); i++)
      {
        if (cluster_lines[i].size() > 0)
        {
          upstream_track.cluster_lines.push_back(cluster_lines[i].at(0));
          upstream_track.plate_index_of_cluster_line.push_back(i);
        }
      }
      
      fitting_success = fit_track(upstream_track);
      one_fit_goodness = upstream_track.fit_goodness;
      
      the_tracks.push_back(upstream_track);
      
      find_guessed_vertex(upstream_track.track_line);
      
      assign_cluster_lines_multiple_track_downstream(cluster_lines, downstream_tracks);
      
      for (int i = 0; i < downstream_tracks.size(); i++)
      {
        fitting_success = fit_track(downstream_tracks.at(i));
        fit_goodness[i] = downstream_tracks.at(i).fit_goodness;
        
        if (true) // !!!! add condition to accept track
        {
          the_tracks.push_back(downstream_tracks.at(i));
        }
      }

      // !!!! do stuff with remaining unassigned cluster lines - partial tracks etc.
    }
    else
    {
      track_finding_success = false;
    }
  }
    
  return track_finding_success;
}

bool track_finding_manager::classify_event(vector<line_3d> cluster_lines[])
{
  int dummy_int, dummy_int_2;
  int min_clusters;
  
  int zero_cluster_counter_all_regions = 0;
  int zero_cluster_counter_for_region[max_number_of_tracking_regions];
  int single_cluster_counter_all_regions = 0;
  int single_cluster_counter_for_region[max_number_of_tracking_regions];
  int double_cluster_counter_all_regions = 0;
  int double_cluster_counter_for_region[max_number_of_tracking_regions];
  int multiple_cluster_counter_all_regions = 0;
  int multiple_cluster_counter_for_region[max_number_of_tracking_regions];
  
  for (int i = 0; i < max_number_of_tracking_regions; i++)
  {
    zero_cluster_counter_for_region[i] = 0;
    single_cluster_counter_for_region[i] = 0;
    double_cluster_counter_for_region[i] = 0;
    multiple_cluster_counter_for_region[i] = 0;
  }
  
  // means two tracking regions only
  if (the_algorithm == 1)
  {
    for (int i = 0; i < the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates(); i++)
    {
      if (cluster_lines[i].size() == 0)
      {
        zero_cluster_counter_all_regions++;
        
        if (i < the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates_in_tracking_region(0))
        {
          zero_cluster_counter_for_region[0]++;
        }
        else
        {
          zero_cluster_counter_for_region[1]++;
        }
      }
      else if (cluster_lines[i].size() == 1)
      {
        single_cluster_counter_all_regions++;
        
        if (i < the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates_in_tracking_region(0))
        {
          single_cluster_counter_for_region[0]++;
        }
        else
        {
          single_cluster_counter_for_region[1]++;
        }
      }
      else if (cluster_lines[i].size() == 2)
      {
        double_cluster_counter_all_regions++;
        
        if (i < the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates_in_tracking_region(0))
        {
          double_cluster_counter_for_region[0]++;
        }
        else
        {
          double_cluster_counter_for_region[1]++;
        }
      }
      else if (cluster_lines[i].size() > 2)
      {
        multiple_cluster_counter_all_regions++;
        
        if (i < the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates_in_tracking_region(0))
        {
          multiple_cluster_counter_for_region[0]++;
        }
        else
        {
          multiple_cluster_counter_for_region[1]++;
        }
      }
    }
    
    if (single_cluster_counter_all_regions == the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates())
    {
      // pure single track, but can be two pieces (from scattering e.g. in target)
      the_global_data_dispatcher->get_event_characteristics()->event_class_by_cluster_multiplicity = 11;
    }
    else if (single_cluster_counter_for_region[0] < the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates_in_tracking_region(0) - the_global_data_dispatcher->get_reconstruction_options()->max_upstream_plates_with_multiple_clusters)
    {
      the_global_data_dispatcher->get_event_characteristics()->event_class_by_cluster_multiplicity = 71;
    }
    // don't use the max_clusters_on_upstream_plates_with_multiple_clusters for now
    else if (single_cluster_counter_for_region[0] == the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates_in_tracking_region(0) - the_global_data_dispatcher->get_reconstruction_options()->max_upstream_plates_with_multiple_clusters && zero_cluster_counter_for_region[0] == 0 && single_cluster_counter_for_region[1] > the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates_in_tracking_region(1) - the_global_data_dispatcher->get_reconstruction_options()->max_midstream_plus_downstream_plates_with_zero_or_multiple_clusters)
    {
      dummy_int = the_global_data_dispatcher->get_partial_geometry()->get_plate_group(the_global_data_dispatcher->get_partial_geometry()->get_triple_plate_group_index(0))->plate_index_list[0];
      
      the_global_data_dispatcher->get_event_characteristics()->event_class_by_cluster_multiplicity = 12;
      
      for (int i = dummy_int; i < dummy_int + 3; i++)
      {
        if (cluster_lines[i].size() != 1)
        {
          the_global_data_dispatcher->get_event_characteristics()->event_class_by_cluster_multiplicity = 13;
        }
      }
    }
    else if (single_cluster_counter_for_region[0] == the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates_in_tracking_region(0) && double_cluster_counter_for_region[1] + multiple_cluster_counter_for_region[1] >= the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates_in_tracking_region(1) - the_global_data_dispatcher->get_reconstruction_options()->max_midstream_plus_downstream_plates_with_less_than_2_clusters)
    {
      the_global_data_dispatcher->get_event_characteristics()->event_class_by_cluster_multiplicity = 21;
    }
    else if (single_cluster_counter_for_region[0] >= the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates_in_tracking_region(0) - the_global_data_dispatcher->get_reconstruction_options()->max_upstream_plates_with_multiple_clusters && double_cluster_counter_for_region[1] + multiple_cluster_counter_for_region[1] >= the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates_in_tracking_region(1) - the_global_data_dispatcher->get_reconstruction_options()->max_midstream_plus_downstream_plates_with_less_than_2_clusters)
    {
      the_global_data_dispatcher->get_event_characteristics()->event_class_by_cluster_multiplicity = 22;
    }
    else
    {
      the_global_data_dispatcher->get_event_characteristics()->event_class_by_cluster_multiplicity = 99;
    }
    
    the_global_data_dispatcher->get_event_reco_output_data_structure()->event_class_by_cluster_multiplicity = the_global_data_dispatcher->get_event_characteristics()->event_class_by_cluster_multiplicity;
    
    for (int i = 0; i < the_global_data_dispatcher->get_partial_geometry()->get_number_of_groups(); i++)
    {
      min_clusters = 1000;
      
      dummy_int = the_global_data_dispatcher->get_partial_geometry()->get_plate_group(i)->number_of_plates_in_group;
      
      for (int j = 0; j < dummy_int; j++)
      {
        dummy_int_2 = the_global_data_dispatcher->get_partial_geometry()->get_plate_group(i)->plate_index_list[j];
        
        if (cluster_lines[dummy_int_2].size() < min_clusters)
        {
          min_clusters = cluster_lines[dummy_int_2].size();
        }
      }

      the_global_data_dispatcher->get_event_characteristics()->group_cluster_multiplicity[i] = min_clusters;
    }
  }
  
  if (PRINT_DEBUG == 3)
  {
    for (int i = 0; i < the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates(); i++)
    {
      cout << cluster_lines[i].size() << " ";
    }
    
    cout << endl;
    
    for (int i = 0; i < the_global_data_dispatcher->get_partial_geometry()->get_number_of_groups(); i++)
    {
      cout << the_global_data_dispatcher->get_event_characteristics()->group_cluster_multiplicity[i] << " ";
    }
    
    cout << endl;
    
    cout << "Class by cluster multiplicity: " << the_global_data_dispatcher->get_event_characteristics()->event_class_by_cluster_multiplicity << endl;
  }
  
  return true;
}

bool track_finding_manager::assign_cluster_lines_single_track_upstream(vector<line_3d> cluster_lines[])
{
  // this function is relying on single plate with more than 1 clusters;
  // if there are more allowed, this function should be made a particular case
  //
 
  bool assignment_success = true;
  
  XYZVector direction_of_interest;
  XYZVector temp_vector;
  
  const int max_number_of_clusters_on_plate = 100;
  double distances[max_number_of_clusters_on_plate];
  double common_distance = 0.;
  double min_distance = 1000.;
  int min_distance_cluster_index = -1;
  
  int mutlicluster_plate;
  int number_of_upstream_plates = the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates_in_tracking_region(0);
  string multicluster_plate_type;
  
  int number_of_plates_of_same_type = 0;
  int plates_of_same_type_index[max_number_of_plates_in_tracking_region];
  
  line_3d temp_cluster_line;
  
  for (int i = 0; i < number_of_upstream_plates; i++)
  {
    if (cluster_lines[i].size() > 1)
    {
      mutlicluster_plate = i;
    }
  }
  
  multicluster_plate_type = the_global_data_dispatcher->get_partial_geometry()->get_plate(mutlicluster_plate)->plate_type;
  
  if (multicluster_plate_type == the_global_data_dispatcher->get_partial_geometry()->get_primary_direction_1_name())
  {
    direction_of_interest = the_global_data_dispatcher->get_partial_geometry()->get_primary_direction_1();
  }
  else if (multicluster_plate_type == the_global_data_dispatcher->get_partial_geometry()->get_primary_direction_2_name())
  {
    direction_of_interest = the_global_data_dispatcher->get_partial_geometry()->get_primary_direction_2();
  }
  else if (multicluster_plate_type == the_global_data_dispatcher->get_partial_geometry()->get_secondary_direction_1_name())
  {
    direction_of_interest = the_global_data_dispatcher->get_partial_geometry()->get_secondary_direction_1();
  }
  else if (multicluster_plate_type == the_global_data_dispatcher->get_partial_geometry()->get_secondary_direction_2_name())
  {
    direction_of_interest = the_global_data_dispatcher->get_partial_geometry()->get_secondary_direction_2();
  }
  
  for (int i = 0; i < number_of_upstream_plates; i++)
  {
    if (the_global_data_dispatcher->get_partial_geometry()->get_plate(i)->plate_type == multicluster_plate_type)
    {
      plates_of_same_type_index[number_of_plates_of_same_type] = i;
      number_of_plates_of_same_type++;
    }
  }
  
  for (int i = 0; i < cluster_lines[mutlicluster_plate].size(); i++)
  {
    distances[i] = 0.;
  }
  for (int i = 1; i < number_of_plates_of_same_type; i++)
  {
    if (i - 1 != mutlicluster_plate && i != mutlicluster_plate)
    {
      temp_vector = cluster_lines[i].at(0).point - cluster_lines[i - 1].at(0).point;
      common_distance += fabs(temp_vector.Dot(direction_of_interest));
    }
    else if (i - 1 == mutlicluster_plate)
    {
      for (int j = 0; j < cluster_lines[mutlicluster_plate].size(); j++)
      {
        temp_vector = cluster_lines[i].at(0).point - cluster_lines[i - 1].at(j).point;
        distances[j] += fabs(temp_vector.Dot(direction_of_interest));
      }
    }
    else
    {
      for (int j = 0; j < cluster_lines[mutlicluster_plate].size(); j++)
      {
        temp_vector = cluster_lines[i - 1].at(0).point - cluster_lines[i].at(j).point;
        distances[j] += fabs(temp_vector.Dot(direction_of_interest));
      }
    }
  }
  
  for (int i = 0; i < cluster_lines[mutlicluster_plate].size(); i++)
  {
    distances[i] += common_distance;
  }

  for (int i = 0; i < cluster_lines[mutlicluster_plate].size(); i++)
  {
    if (distances[i] < min_distance)
    {
      min_distance = distances[i];
      min_distance_cluster_index = i;
    }
  }

  for (int i = 0; i < cluster_lines[mutlicluster_plate].size(); i++)
  {
    if (i == min_distance_cluster_index)
    {
      temp_cluster_line = cluster_lines[mutlicluster_plate].at(min_distance_cluster_index);
      cluster_lines[mutlicluster_plate].clear();
      cluster_lines[mutlicluster_plate].push_back(temp_cluster_line);
    }
  }

  return assignment_success;
}

bool track_finding_manager::assign_cluster_lines_single_track_downstream(vector<line_3d> cluster_lines[])
{
  bool assignment_sucess = true;
  
  int multi_line_plates[max_number_of_plates];
  int number_of_multi_line_plates = 0;
  
  track the_track;
  
  XYZVector track_plate_intersection;
  
  plane_3d plate_plane;
  
  double distance, shortest_distance;
  int index_shortest_distance; // the cluster_line index in the vector that is shortest to the intersection point
  
  for (int i = the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates_in_tracking_region(0); i < the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates(); i++)
  {
    if (cluster_lines[i].size() == 1)
    {
      the_track.cluster_lines.push_back(cluster_lines[i].at(0));
      the_track.plate_index_of_cluster_line.push_back(i);
    }
    else if (cluster_lines[i].size() > 1)
    {
      multi_line_plates[number_of_multi_line_plates] = i;
      number_of_multi_line_plates++;
    }
  }
  
  if (fit_track(the_track))
  {
    for (int i = 0; i < number_of_multi_line_plates; i++)
    {
      shortest_distance = 1000.;
      
      plate_plane = the_global_data_dispatcher->get_partial_geometry()->get_plate(multi_line_plates[i])->plate_plane_equation;
      
      if (!find_line_plate_intersection_point(the_track.track_line, plate_plane, track_plate_intersection))
      {
        assignment_sucess = false;
      }

      if (assignment_sucess)
      {
        for (int j = 0; j < cluster_lines[multi_line_plates[i]].size(); j++)
        {
          distance = find_point_to_line_3d_distance(track_plate_intersection, cluster_lines[multi_line_plates[i]].at(j));
          
          if (distance < shortest_distance)
          {
            shortest_distance = distance;
            index_shortest_distance = j;
          }
        }

        the_track.cluster_lines.push_back(cluster_lines[multi_line_plates[i]].at(index_shortest_distance));
      }
    }
  }
  else
  {
    assignment_sucess = false;
  }

  return assignment_sucess;
}


bool track_finding_manager::assign_cluster_lines_multiple_track_downstream(vector<line_3d> cluster_lines[], vector<track>& downstream_tracks)
{
  bool assignment_sucess = true;
  
  int triplet_group_index;
  ssd_plate triplet_plates[3];
  int plate_index[3], temp_index[3];
  string triplet_group_type;
  int secondary_index;

  int lowest_clusters = 1000;
  int lowest_clusters_index;
  
  XYZVector midpoint;
  XYZVector temp_direction;
  vector<line_3d> lines_to_fit;
  XYZVector closest_point, closest_point_2;

  vector<double> initial_values_of_parameters;
  vector<double> errors_in_initial_values_of_parameters;
  
  XYZVector intersection_point;
  XYZVector point_on_third_plate;
  
  track dummy_track;
  track temp_tracks[40000];
  int number_of_temp_tracks;
  
  vector<track_option> the_track_options;
  track_option dummy_track_option;
  
  double largest_z;
  
  double a_distance, smallest_distance;
  int smallest_distance_cluster_line_index;

  double cut_max_distance_to_accept_in_track = the_global_data_dispatcher->get_reconstruction_options()->cut_max_distance_to_accept_in_track;
  
  // assuming the algorithm is 1
  if (the_algorithm != 1)
  {
    return false;
  }
  
  triplet_group_index = the_global_data_dispatcher->get_partial_geometry()->get_triple_plate_group_index(0);
  triplet_group_type = the_global_data_dispatcher->get_partial_geometry()->get_plate_group(triplet_group_index)->group_type;
  
  if (the_global_data_dispatcher->get_partial_geometry()->get_plate_group(triplet_group_index)->number_of_plates_in_group != 3)
  {
    cout << "Something is wrong with the geometry..." << endl;
    return false;
  }
  
  temp_index[0] = the_global_data_dispatcher->get_partial_geometry()->get_plate_group(triplet_group_index)->plate_index_list[0];
  temp_index[1] = the_global_data_dispatcher->get_partial_geometry()->get_plate_group(triplet_group_index)->plate_index_list[1];
  temp_index[2] = the_global_data_dispatcher->get_partial_geometry()->get_plate_group(triplet_group_index)->plate_index_list[2];
  
  // check if triplet group is OK: X_Y_D1 etc.
  if (triplet_group_type == "X_Y_D1" || triplet_group_type == "Y_X_D1" || triplet_group_type == "X_Y_D2" || triplet_group_type == "Y_X_D2")
  {
    secondary_index = 2;
  }
  else if (triplet_group_type == "D1_X_Y" || triplet_group_type == "D1_Y_X" || triplet_group_type == "D2_X_Y" || triplet_group_type == "D2_Y_X")
  {
    secondary_index = 0;
  }
  // the next doesn't make much sense
  else if (triplet_group_type == "X_D1_Y" || triplet_group_type == "Y_D1_X" || triplet_group_type == "X_D2_Y" || triplet_group_type == "Y_D2_X")
  {
    secondary_index = 1;
  }
  else
  {
    cout << "Unknown triplet group type: " << triplet_group_type << endl;
    return false;
  }
  
  for (int i = 0; i < 3; i++)
  {
    if (cluster_lines[temp_index[i]].size() < lowest_clusters)
    {
      lowest_clusters = cluster_lines[temp_index[i]].size();
    }
  }
  
  // 0 and 1 for primary; 2 for secondary
  if (secondary_index == 0)
  {
    plate_index[0] = temp_index[1];
    plate_index[1] = temp_index[2];
  }
  else if (secondary_index == 1)
  {
    plate_index[0] = temp_index[0];
    plate_index[1] = temp_index[2];
  }
  else if (secondary_index == 2)
  {
    plate_index[0] = temp_index[0];
    plate_index[1] = temp_index[1];
  }
  
  plate_index[2] = temp_index[secondary_index];
  
  triplet_plates[0] = *(the_global_data_dispatcher->get_partial_geometry()->get_plate(plate_index[0]));
  triplet_plates[1] = *(the_global_data_dispatcher->get_partial_geometry()->get_plate(plate_index[1]));
  triplet_plates[2] = *(the_global_data_dispatcher->get_partial_geometry()->get_plate(plate_index[2]));

  // cout << "Lowest clusters: " << lowest_clusters << endl;
  // cout << "Tiplet plate indices: " << plate_index[0] << " " << plate_index[1] << " " << plate_index[2] << endl;
  
  number_of_temp_tracks = 0;
  the_track_options.clear();
  
  for (int i = 0; i < cluster_lines[plate_index[0]].size(); i++)
  {
    for (int j = 0; j < cluster_lines[plate_index[1]].size(); j++)
    {
      if (!find_closest_point(cluster_lines[plate_index[0]].at(i), cluster_lines[plate_index[1]].at(j), closest_point))
      {
        cout << "Could not find point between the two primary lines in triplet group." << endl;
        return false;
      }
  
      initial_values_of_parameters.clear();
      errors_in_initial_values_of_parameters.clear();
      lines_to_fit.clear();
      
      lines_to_fit.push_back(cluster_lines[plate_index[0]].at(i));
      lines_to_fit.push_back(cluster_lines[plate_index[1]].at(j));

      midpoint = 0.5 * (closest_point + guessed_vertex.point);
      temp_direction = (closest_point - guessed_vertex.point).Unit();
  
      initial_values_of_parameters.push_back(midpoint.X());
      initial_values_of_parameters.push_back(midpoint.Y());
      initial_values_of_parameters.push_back(temp_direction.X() / temp_direction.Z());
      initial_values_of_parameters.push_back(temp_direction.Y() / temp_direction.Z());
      
      // need better estimates for these
      errors_in_initial_values_of_parameters.push_back(0.1);
      errors_in_initial_values_of_parameters.push_back(0.1);
      errors_in_initial_values_of_parameters.push_back(0.01);
      errors_in_initial_values_of_parameters.push_back(0.01);

      my_straight_line_fcn_with_point the_fcn(lines_to_fit, guessed_vertex, the_global_data_dispatcher->get_reconstruction_options()->guessed_vertex_weight, midpoint.Z());
      
      VariableMetricMinimizer the_minimizer;
      
      FunctionMinimum min = the_minimizer.Minimize(the_fcn, initial_values_of_parameters, errors_in_initial_values_of_parameters);
  
      temp_tracks[number_of_temp_tracks].track_line.point.SetCoordinates(min.UserState().Value(0), min.UserState().Value(1), midpoint.Z());
      temp_tracks[number_of_temp_tracks].track_line.direction.SetCoordinates(min.UserState().Value(2), min.UserState().Value(3), 1.);
      temp_tracks[number_of_temp_tracks].track_line.direction = temp_tracks[number_of_temp_tracks].track_line.direction.Unit();
      
      find_line_plate_intersection_point(temp_tracks[number_of_temp_tracks].track_line, triplet_plates[2].plate_plane_equation, point_on_third_plate);
      number_of_temp_tracks++;
      
      for (int k = 0; k < cluster_lines[plate_index[2]].size(); k++)
      {
        dummy_track_option.track_index = number_of_temp_tracks - 1;
        dummy_track_option.cluster_index_on_triplet_plates[0] = i;
        dummy_track_option.cluster_index_on_triplet_plates[1] = j;
        dummy_track_option.cluster_index_on_triplet_plates[2] = k;
        dummy_track_option.distance = find_point_to_line_3d_distance(point_on_third_plate, cluster_lines[plate_index[2]].at(k));
        dummy_track_option.strip_on_third_plate = temp_point_to_strip_diagonal_plate(point_on_third_plate);
        the_track_options.push_back(dummy_track_option);
      }
    }
  }
      
  sort(the_track_options.begin(), the_track_options.end(), compare_track_options);
  
/*
  for (int i = 0; i < the_track_options.size(); i++)
  {
    cout << "-> " << the_track_options.at(i).cluster_index_on_triplet_plates[0] << " " << the_track_options.at(i).cluster_index_on_triplet_plates[1] << " " << the_track_options.at(i).cluster_index_on_triplet_plates[2] << " " << the_track_options.at(i).strip_on_third_plate << " " << temp_point_to_strip_diagonal_plate(cluster_lines[plate_index[2]].at(the_track_options.at(i).cluster_index_on_triplet_plates[2]).point) << " " << the_track_options.at(i).distance <<
        " " << cluster_lines[plate_index[0]].at(the_track_options.at(i).cluster_index_on_triplet_plates[0]).point.Y() << " " << cluster_lines[plate_index[1]].at(the_track_options.at(i).cluster_index_on_triplet_plates[1]).point.X() << endl;
  }
*/
  
  dummy_track.cluster_lines.clear();
  dummy_track.plate_index_of_cluster_line.clear();
  
  for (int i = 0; i < lowest_clusters; i++)
  {
    dummy_track.cluster_lines.push_back(cluster_lines[plate_index[2]].at(the_track_options.at(i).cluster_index_on_triplet_plates[2]));
    dummy_track.cluster_lines.push_back(cluster_lines[plate_index[0]].at(the_track_options.at(i).cluster_index_on_triplet_plates[0]));
    dummy_track.cluster_lines.push_back(cluster_lines[plate_index[1]].at(the_track_options.at(i).cluster_index_on_triplet_plates[1]));

    dummy_track.plate_index_of_cluster_line.push_back(plate_index[2]);
    dummy_track.plate_index_of_cluster_line.push_back(plate_index[0]);
    dummy_track.plate_index_of_cluster_line.push_back(plate_index[1]);
    
    dummy_track.track_line = temp_tracks[the_track_options.at(i).track_index].track_line;
    
    dummy_track.assumed_start_group = guessed_vertex.point.Z();
    
    largest_z = -100.;
    
    for (int j = 0; j < dummy_track.cluster_lines.size(); j++)
    {
      if (dummy_track.cluster_lines.at(j).point.Z() > largest_z)
      {
        largest_z = dummy_track.cluster_lines.at(j).point.Z();
      }
    }
    
    dummy_track.assumed_end_group = largest_z;
    dummy_track.sort_cluster_lines_by_plate_index();
    
    downstream_tracks.push_back(dummy_track);
  }
  
  // !!!! remove used cluster lines and do extra stuff with what's left
  
  for (int i = the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates_in_tracking_region(0); i < the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates(); i++)
  {
    if (i != plate_index[0] && i != plate_index[1] && i != plate_index[2])
    {
      for (int j = 0; j < downstream_tracks.size(); j++)
      {
        find_line_plate_intersection_point(downstream_tracks[j].track_line, the_global_data_dispatcher->get_partial_geometry()->get_plate(i)->plate_plane_equation, intersection_point);
        
        smallest_distance = 1000.;

        for (int k = 0; k < cluster_lines[i].size(); k++)
        {
          a_distance = find_point_to_line_3d_distance(intersection_point, cluster_lines[i].at(k));
          
          if (a_distance < smallest_distance)
          {
            smallest_distance = a_distance;
            smallest_distance_cluster_line_index = k;
          }
        }
        
        if (smallest_distance < cut_max_distance_to_accept_in_track)
        {
          downstream_tracks.at(j).cluster_lines.push_back(cluster_lines[i].at(smallest_distance_cluster_line_index));
          downstream_tracks.at(j).plate_index_of_cluster_line.push_back(i);
          
          if (cluster_lines[i].at(smallest_distance_cluster_line_index).point.Z() > downstream_tracks.at(j).assumed_end_group)
          {
            downstream_tracks.at(j).assumed_end_group = cluster_lines[i].at(smallest_distance_cluster_line_index).point.Z();
          }
          
          if (cluster_lines[i].at(smallest_distance_cluster_line_index).point.Z() < downstream_tracks.at(j).assumed_start_group)
          {
            downstream_tracks.at(j).assumed_start_group = cluster_lines[i].at(smallest_distance_cluster_line_index).point.Z();
          }
        }
      }
    }
  }
  
  for (int i = 0; i < downstream_tracks.size(); i++)
  {
    downstream_tracks.at(i).sort_cluster_lines_by_plate_index();
   
    initial_values_of_parameters.clear();
    errors_in_initial_values_of_parameters.clear();
  
    midpoint.SetZ(0.5 * (downstream_tracks.at(i).assumed_end_group + downstream_tracks.at(i).assumed_start_group));
    midpoint.SetX(downstream_tracks.at(i).track_line.point.X() + downstream_tracks.at(i).track_line.direction.X() *  (downstream_tracks.at(i).track_line.point.Z() - midpoint.Z()) / downstream_tracks.at(i).track_line.direction.Z());
    midpoint.SetY(downstream_tracks.at(i).track_line.point.Y() + downstream_tracks.at(i).track_line.direction.Y() *  (downstream_tracks.at(i).track_line.point.Z() - midpoint.Z()) / downstream_tracks.at(i).track_line.direction.Z());
    
    initial_values_of_parameters.push_back(midpoint.X());
    initial_values_of_parameters.push_back(midpoint.Y());
    initial_values_of_parameters.push_back(downstream_tracks.at(i).track_line.direction.X() / downstream_tracks.at(i).track_line.direction.Z());
    initial_values_of_parameters.push_back(downstream_tracks.at(i).track_line.direction.Y() / downstream_tracks.at(i).track_line.direction.Z());
  
    // need better estimates for these
    errors_in_initial_values_of_parameters.push_back(0.1);
    errors_in_initial_values_of_parameters.push_back(0.1);
    errors_in_initial_values_of_parameters.push_back(0.01);
    errors_in_initial_values_of_parameters.push_back(0.01);
  
    my_straight_line_fcn the_fcn(downstream_tracks.at(i).cluster_lines, midpoint.Z());
  
    VariableMetricMinimizer the_minimizer;
  
    FunctionMinimum min = the_minimizer.Minimize(the_fcn, initial_values_of_parameters, errors_in_initial_values_of_parameters);
  
    downstream_tracks.at(i).track_line.point.SetCoordinates(min.UserState().Value(0), min.UserState().Value(1), midpoint.Z());
    downstream_tracks.at(i).track_line.direction.SetCoordinates(min.UserState().Value(2), min.UserState().Value(3), 1.);
    downstream_tracks.at(i).track_line.direction = downstream_tracks.at(i).track_line.direction.Unit();
  }
  
  return assignment_sucess;
}

bool track_finding_manager::find_vertex(vector<track>& the_tracks, XYZVector the_vertex)
{
  XYZVector temp_vector;
  double distance;
  int counter = 0;
  
  for (int i = 0; i < the_tracks.size() - 1; i++)
  {
    for (int j = i + 1; j < the_tracks.size(); j++)
    {
      if (find_line_to_line_distance(the_tracks.at(i).track_line, the_tracks.at(j).track_line, distance))
      {
        if (distance < the_global_data_dispatcher->get_reconstruction_options()->cut_max_distance_to_count_to_vertex)
        {
          find_closest_point(the_tracks.at(i).track_line, the_tracks.at(j).track_line, temp_vector);
          the_vertex += temp_vector;
          counter++;
        }
      }
    }
  }
  
  if (counter > 0)
  {
    the_vertex /= (double)(counter);
    return true;
  }
  else
  {
    return false;
  }
}

bool track_finding_manager::fit_track(track& a_track)
{
  bool fitting_success = true;
  
  double average_x = 0., average_y = 0., average_z = 0.;
  double dx_dz, dy_dz;
  
  string primary_direciton_names[2];
  int primary_direction_plates[2][max_number_of_plates];
  int number_of_primary_direction_plates[2] = {0, 0};
  
  vector<double> initial_values_of_parameters;
  vector<double> errors_in_initial_values_of_parameters;
  
  for (int i = 0; i < 2; i++)
  {
    for (int j = 0; j < max_number_of_plates; j++)
    {
      primary_direction_plates[i][j] = -1;
    }
  }
  
  primary_direciton_names[0] = the_global_data_dispatcher->get_partial_geometry()->get_primary_direction_1_name();
  primary_direciton_names[1] = the_global_data_dispatcher->get_partial_geometry()->get_primary_direction_2_name();
  
  for (int i = 0; i < a_track.cluster_lines.size(); i++)
  {
    if (the_global_data_dispatcher->get_partial_geometry()->get_plate(a_track.plate_index_of_cluster_line.at(i))->plate_type == primary_direciton_names[0])
    {
      primary_direction_plates[0][number_of_primary_direction_plates[0]] = i;
      number_of_primary_direction_plates[0]++;
    }
    else if (the_global_data_dispatcher->get_partial_geometry()->get_plate(a_track.plate_index_of_cluster_line.at(i))->plate_type == primary_direciton_names[1])
    {
      primary_direction_plates[1][number_of_primary_direction_plates[1]] = i;
      number_of_primary_direction_plates[1]++;
    }
  }
  
  if (number_of_primary_direction_plates[0] < 2 || number_of_primary_direction_plates[1] < 2)
  {
    cout << "Too few plates in a given direction to fit a track..." << endl;
    return false;
  }
  
  // also only implemented for X and Y as primary
  for (int i = 0; i < number_of_primary_direction_plates[0]; i++)
  {
    average_z += a_track.cluster_lines.at(primary_direction_plates[0][i]).point.z();
    
    if (primary_direciton_names[0] == "X")
    {
      average_x += a_track.cluster_lines.at(primary_direction_plates[0][i]).point.x();
    }
    else if (primary_direciton_names[0] == "Y")
    {
      average_y += a_track.cluster_lines.at(primary_direction_plates[0][i]).point.y();
    }
    else
    {
      cout << "These primary directions have not been implemented yet" << endl;
      return false;
    }
  }
  
  if (primary_direciton_names[0] == "X")
  {
    average_x /= (double) number_of_primary_direction_plates[0];
  }
  else if (primary_direciton_names[0] == "Y")
  {
    average_y /= (double) number_of_primary_direction_plates[0];
  }
  
  for (int i = 0; i < number_of_primary_direction_plates[1]; i++)
  {
    average_z += a_track.cluster_lines.at(primary_direction_plates[1][i]).point.z();
    
    if (primary_direciton_names[1] == "X")
    {
      average_x += a_track.cluster_lines.at(primary_direction_plates[1][i]).point.x();
    }
    else if (primary_direciton_names[1] == "Y")
    {
      average_y += a_track.cluster_lines.at(primary_direction_plates[1][i]).point.y();
    }
    else
    {
      cout << "These primary directions have not been implemented yet" << endl;
      return false;
    }
  }
  
  if (primary_direciton_names[1] == "X")
  {
    average_x /= (double) number_of_primary_direction_plates[1];
  }
  else if (primary_direciton_names[1] == "Y")
  {
    average_y /= (double) number_of_primary_direction_plates[1];
  }
  
  average_z /= (double)(number_of_primary_direction_plates[0] + number_of_primary_direction_plates[1]);
  
  // this is a particular case; work out the general case later
  //
  if (primary_direciton_names[0] == "X" && primary_direciton_names[1] == "Y")
  {
    dy_dz = (a_track.cluster_lines.at(primary_direction_plates[1][number_of_primary_direction_plates[1] - 1]).point.y() - a_track.cluster_lines.at(primary_direction_plates[1][0]).point.y()) / (a_track.cluster_lines.at(primary_direction_plates[1][number_of_primary_direction_plates[1] - 1]).point.z() - a_track.cluster_lines.at(primary_direction_plates[1][0]).point.z());
    dx_dz = (a_track.cluster_lines.at(primary_direction_plates[0][number_of_primary_direction_plates[0] - 1]).point.x() - a_track.cluster_lines.at(primary_direction_plates[0][0]).point.x()) / (a_track.cluster_lines.at(primary_direction_plates[0][number_of_primary_direction_plates[0] - 1]).point.z() - a_track.cluster_lines.at(primary_direction_plates[0][0]).point.z());
  }
  else if (primary_direciton_names[0] == "Y" && primary_direciton_names[1] == "X")
  {
    dy_dz = (a_track.cluster_lines.at(primary_direction_plates[0][number_of_primary_direction_plates[0] - 1]).point.y() - a_track.cluster_lines.at(primary_direction_plates[0][0]).point.y()) / (a_track.cluster_lines.at(primary_direction_plates[0][number_of_primary_direction_plates[0] - 1]).point.z() - a_track.cluster_lines.at(primary_direction_plates[0][0]).point.z());
    dx_dz = (a_track.cluster_lines.at(primary_direction_plates[1][number_of_primary_direction_plates[1] - 1]).point.x() - a_track.cluster_lines.at(primary_direction_plates[1][0]).point.x()) / (a_track.cluster_lines.at(primary_direction_plates[1][number_of_primary_direction_plates[1] - 1]).point.z() - a_track.cluster_lines.at(primary_direction_plates[1][0]).point.z());
  }
  else
  {
    cout << "This configuration of primary directions has not been implemented yet..." << endl;
    return false;
  }
  
  initial_values_of_parameters.clear();
  errors_in_initial_values_of_parameters.clear();
  
  initial_values_of_parameters.push_back(average_x);
  initial_values_of_parameters.push_back(average_y);
  initial_values_of_parameters.push_back(dx_dz);
  initial_values_of_parameters.push_back(dy_dz);

  // need better estimates for these
  errors_in_initial_values_of_parameters.push_back(0.1);
  errors_in_initial_values_of_parameters.push_back(0.1);
  errors_in_initial_values_of_parameters.push_back(0.01);
  errors_in_initial_values_of_parameters.push_back(0.01);
  
  my_straight_line_fcn the_fcn(a_track.cluster_lines, average_z);
  
  VariableMetricMinimizer the_minimizer;
  
  FunctionMinimum min = the_minimizer.Minimize(the_fcn, initial_values_of_parameters, errors_in_initial_values_of_parameters);
  
  a_track.track_line.point.SetCoordinates(min.UserState().Value(0), min.UserState().Value(1), average_z);
  a_track.track_line.direction.SetCoordinates(min.UserState().Value(2), min.UserState().Value(3), 1.);
  a_track.track_line.direction = a_track.track_line.direction.Unit();
  
  a_track.fit_goodness = min.Fval();
  
  return fitting_success;
}

bool track_finding_manager::find_guessed_vertex(line_3d upstream_track_line)
{
  // assuming the target is in the xy plane
  double target_z, t;
  
  target_z = the_global_data_dispatcher->get_partial_geometry()->get_target()->position.Z();
  
  if (upstream_track_line.direction.Z() != 0)
  {
    t = (target_z - upstream_track_line.point.Z()) / upstream_track_line.direction.Z();
    guessed_vertex.point.SetCoordinates(upstream_track_line.point.X() + upstream_track_line.direction.X() * t, upstream_track_line.point.Y() + upstream_track_line.direction.Y() * t, target_z);
  }
  else
  {
    return false;
  }
  
  guessed_vertex.direction = upstream_track_line.direction;
  
  return true;
}

