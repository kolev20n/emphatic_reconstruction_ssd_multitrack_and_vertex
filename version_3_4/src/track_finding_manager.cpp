#include <iostream>

#include "track_finding_manager.h"
#include "partial_geometry.h"
#include "global_data_dispatcher.h"
#include "common_data_structures_and_functions.h"

using namespace std;
using namespace ROOT::Math;

track_finding_manager::track_finding_manager(partial_geometry* a_geometry, int an_algorithm)
:the_geometry(a_geometry), the_algorithm(an_algorithm)
{
  cout << "Created the_track_finding_manager." << endl;
}

bool track_finding_manager::find_tracks(vector<line_3d> cluster_lines[], vector<track>& the_tracks)
{
  bool classify_event_sucess = true;
  bool track_finding_success = true;
  
  int event_class_by_cluster_multiplicity;
  
  vector<vector<track>> all_tracks_options;
  vector<track> all_tracks;
  track temp_track;
  
  vector<double> fit_quality_options;
  
  vector<line_3d> assigned_cluster_lines[max_number_of_plates];
  
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
        assign_cluster_lines_upstream(cluster_lines, assigned_cluster_lines);
      }
      else
      {
        for (int i = 0; i < the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates(); i++)
        {
          assigned_cluster_lines[i] = cluster_lines[i];
        }
      }
      
for (int i = 0; i < the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates(); i++)
{
  cout << cluster_lines[i].size() << " " << assigned_cluster_lines[i].size() << endl;
}

      all_tracks_options.clear();
      fit_quality_options.clear();
      
      all_tracks.clear();
      
      fit_quality_options.push_back(fit_single_track(assigned_cluster_lines, 0, the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates(), temp_track));
      
      all_tracks.push_back(temp_track);
      all_tracks_options.push_back(all_tracks);
      
      all_tracks.clear();

      fit_quality_options.push_back(fit_single_track(assigned_cluster_lines, 0, the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates_in_tracking_region(0), temp_track));
      
      all_tracks.push_back(temp_track);
      
      fit_quality_options.push_back(fit_single_track(assigned_cluster_lines, the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates_in_tracking_region(0), the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates(), temp_track));
      
      all_tracks.push_back(temp_track);
      all_tracks_options.push_back(all_tracks);
      
      if (fit_quality_options.at(0) < fit_quality_options.at(1) + fit_quality_options.at(2))
      {
        the_tracks.push_back(all_tracks_options.at(0).at(0));
      }
      else
      {
        the_tracks.push_back(all_tracks_options.at(1).at(0));
        the_tracks.push_back(all_tracks_options.at(1).at(1));
      }
    }
    else if (event_class_by_cluster_multiplicity == 21)
    {
      
    }
    else if (event_class_by_cluster_multiplicity == 71 || event_class_by_cluster_multiplicity == 99)
    {
      track_finding_success = false;
    }
  }
  
  return track_finding_success;
}

bool track_finding_manager::classify_event(vector<line_3d> cluster_lines[])
{
  int dummy_int;
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
      the_global_data_dispatcher->get_event_characteristics()->group_cluster_multiplicity[i] = the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates();
      
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
    else if (single_cluster_counter_for_region[0] >= the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates_in_tracking_region(0) - the_global_data_dispatcher->get_reconstruction_options()->max_upstream_plates_with_multiple_clusters && double_cluster_counter_for_region[1] + multiple_cluster_counter_for_region[1] >= the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates_in_tracking_region(1) - the_global_data_dispatcher->get_reconstruction_options()->max_midstream_plus_downstream_plates_with_less_than_2_clusters)
    {
      the_global_data_dispatcher->get_event_characteristics()->event_class_by_cluster_multiplicity = 21;
    }
    else
    {
      the_global_data_dispatcher->get_event_characteristics()->event_class_by_cluster_multiplicity = 99;
    }
    
    for (int i = 0; i < the_global_data_dispatcher->get_partial_geometry()->get_number_of_groups(); i++)
    {
      min_clusters = 1000;
      
      dummy_int = the_global_data_dispatcher->get_partial_geometry()->get_plate_group(i)->number_of_plates_in_group;
      
      for (int j = 0; j < dummy_int; j++)
      {
        dummy_int = the_global_data_dispatcher->get_partial_geometry()->get_plate_group(i)->plate_index_list[j];
        
        if (cluster_lines[dummy_int].size() < min_clusters)
        {
          min_clusters = cluster_lines[dummy_int].size();
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

bool track_finding_manager::fit_single_track(vector<line_3d> some_cluster_lines[], int first_cluster_line, int last_cluster_line, track& a_track)
{
  return true;
}

bool track_finding_manager::assign_cluster_lines_upstream(vector<line_3d> cluster_lines[], vector<line_3d> assigned_cluster_lines[])
{
  // this function is relying on single plate with more than 1 clusters;
  // if there are more allowed, this function should be made a particular case
  //
  
  XYZVector direction_of_interest;
  XYZVector temp_vector;
  
  const int max_number_of_clusters_on_plate = 100;
  double distances[max_number_of_clusters_on_plate];
  double common_distance = 0;
  double min_distance = 1000.;
  int min_distance_cluster_index = -1;
  
  int mutlicluster_plate;
  int number_of_upstream_plates = the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates_in_tracking_region(0);
  string multicluster_plate_type;
  
  int number_of_plates_of_same_type = 0;
  int plates_of_same_type_index[max_number_of_plates_in_tracking_region];
  
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
  
  for (int i = 0; i < the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates(); i++)
  {
    if (i != mutlicluster_plate)
    {
      assigned_cluster_lines[i] = cluster_lines[i];
    }
    else
    {
      assigned_cluster_lines[i].push_back(cluster_lines[i].at(min_distance_cluster_index));
    }
  }
  
  return true;
}
