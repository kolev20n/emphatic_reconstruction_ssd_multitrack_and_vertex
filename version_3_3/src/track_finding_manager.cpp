#include <iostream>
#include <iomanip>

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnMigrad.h"

#include "track_finding_manager.h"
#include "common_data_structures_and_functions.h"
#include "histogram_manager.h"
#include "global_data_dispatcher.h"
#include "my_straight_line_fcn.h"

using namespace std;
using namespace ROOT::Math;
using namespace ROOT::Minuit2;

track_finding_manager::track_finding_manager(histogram_manager* a_histogram_manager, output_data_structure* an_output_data_structure, partial_geometry* a_geometry)
:the_histogram_manager(a_histogram_manager), the_geometry(a_geometry)
{
  cout << "The track finding manager has been created..." << endl;
  
  the_output_data_structure = an_output_data_structure;
  
  //cout << "-.-.-.-.-.-.-.-" << (the_global_data_dispatcher->get_partial_geometry())->get_number_of_plates() << endl;
  
}

bool track_finding_manager::find_preliminary_event_type(vector<silicon_strip> cluster_lines[])
{
  int min_clusters;
  
  int plate_index;
  
  // count how many groups have given multiplicities
  // count_1 groups that have a single minimum clusters
  // count_upstream_1 upstream groups with a single minimum clusters
  // count_increases number of time cluster increase (1 is a normal vertex, if between upstream and midstream)
  // count_decreases to formally recognize messy events (otherwise it is the normal thing to happen)
  //
  int count_1 = 0;
  int count_upstream_1 = 0;
  int count_increases = 0;
  int count_decreases = 0;
  int previous;
  
  bool has_target_vertex = false;

  group_cluster_multiplicity.clear();
  
  
  // cout << "Group cluster mulltiplicities: ";
  
  for (int i = 0; i < the_geometry->the_plate_groups.size(); i++)
  {
    min_clusters = 1000;
    
    for (int j = 0; j < the_geometry->the_plate_groups.at(i).plate_index_list.size(); j++)
    {
      plate_index = the_geometry->the_plate_groups.at(i).plate_index_list.at(j);
      
      if (cluster_lines[plate_index].size() < min_clusters)
      {
        min_clusters = cluster_lines[plate_index].size();
      }
    }
    
    if (min_clusters == 1000) min_clusters = 0;
    
    group_cluster_multiplicity.push_back(min_clusters);
    
    // cout << min_clusters << " ";
  }
  
  // cout << endl;
  
  // cout << "Plate cluster mulltiplicities: ";
  
  for (int i = 0; i < the_geometry->get_number_of_plates(); i++)
  {
    // cout << cluster_lines[i].size() << " ";
  }
  // cout << endl;
  
  previous = -1;
  
  for (int i = 0; i < the_geometry->the_plate_groups.size(); i++)
  {
    if (group_cluster_multiplicity.at(i) == 1)
    {
      count_1++;
      
      if (the_geometry->the_plate_groups.at(i).group_stream == "upstream")
      {
        count_upstream_1++;
      }
    }
    
    if (i > 0)
    {
      if (group_cluster_multiplicity.at(i) > previous)
      {
        count_increases++;
        
        if (i == the_geometry->get_number_of_upstream_groups())
        {
          has_target_vertex = true;
        }
        
        potential_vertex_before_group = i;
      }
      else if (group_cluster_multiplicity.at(i) < previous)
      {
        count_decreases++;
      }
    }
    
    previous = group_cluster_multiplicity.at(i);
  }
  
  /*
  cout << "Has target vertex: " << has_target_vertex << endl;
  cout << "Groups with 1 minimum cluster: " << count_1 << endl;
  cout << "Upstream groups with 1 minimum cluster: " << count_upstream_1 << endl;
  cout << "Increases: " << count_increases << endl;
  cout << "Decreases: " << count_decreases << endl;
  */
  
  if (count_1 == the_geometry->the_plate_groups.size())
  {
    preliminary_event_type = "single_track";
    
    /*
    cout << "????????" << endl;
    
    for (int i = 0; i < the_global_data_dispatcher->get_input_root_data_structure()->track_id->size(); i++)
    {
      if (the_global_data_dispatcher->get_input_root_data_structure()->ekin_vertex->at(i) > 500.)
      {
        cout << the_global_data_dispatcher->get_input_root_data_structure()->track_id->at(i) << " "
             << the_global_data_dispatcher->get_input_root_data_structure()->particle_code->at(i) << " "
             << the_global_data_dispatcher->get_input_root_data_structure()->parent_track_id->at(i) << " "
             << the_global_data_dispatcher->get_input_root_data_structure()->z_vertex->at(i) << " "
             << the_global_data_dispatcher->get_input_root_data_structure()->ekin_vertex->at(i) << " "
             << the_global_data_dispatcher->get_input_root_data_structure()->creator_process->at(i) << " "
             << the_global_data_dispatcher->get_input_root_data_structure()->px_vertex->at(i) << " "
             << the_global_data_dispatcher->get_input_root_data_structure()->py_vertex->at(i) << " "
             << the_global_data_dispatcher->get_input_root_data_structure()->pz_vertex->at(i) << endl;
      }
    }
    
    cout << "????????" << endl;
    */
  }
  else if (count_upstream_1 == the_geometry->get_number_of_upstream_groups())
  {
    if (count_increases == 1)
    {
      if (has_target_vertex)
      {
        preliminary_event_type = "normal_vertex";
      }
      else
      {
        preliminary_event_type = "abnormal_vertex";
      }
    }
    else if (count_increases > 1)
    {
      preliminary_event_type = "secondary_vertex";
    }
  }
  else
  {
    preliminary_event_type = "unknown";
  }
  
  /*
  // should define messy
  if (count_increases > 1 && count_decreases > 1)
  {
    preliminary_event_type = "messy";
  }
  */
  
  // cout << preliminary_event_type << endl;
  
  
  
  if (preliminary_event_type == "unknown")
  {
    /*
    cout << "Start Detailed Info:-------------->" << endl;
    
    cout << "Silicon Strip Hits:" << endl;
    for (int i = 0; i < the_global_data_dispatcher->get_input_root_data_structure()->plate_number->size(); i++)
    {
      cout << the_global_data_dispatcher->get_input_root_data_structure()->plate_number->at(i) << " " << the_global_data_dispatcher->get_input_root_data_structure()->strip_number->at(i) << " "  << the_global_data_dispatcher->get_input_root_data_structure()->total_energy_values->at(i) << " ";
      for (int j = 0; j < the_global_data_dispatcher->get_input_root_data_structure()->contributing_tracks->at(i).size(); j++)
      {
        cout << the_global_data_dispatcher->get_input_root_data_structure()->contributing_tracks->at(i).at(j) << " ";
      }
      cout << endl;
    }
    */
     
    vector<sorted_track> mc_tracks;
    sorted_track dummy_sorted_track;
    
    for (int i = 0; i < the_global_data_dispatcher->get_input_root_data_structure()->track_id->size(); i++)
    {
      dummy_sorted_track.track_id = the_global_data_dispatcher->get_input_root_data_structure()->track_id->at(i);
      dummy_sorted_track.particle_code = the_global_data_dispatcher->get_input_root_data_structure()->particle_code->at(i);
      dummy_sorted_track.parent_track_id = the_global_data_dispatcher->get_input_root_data_structure()->parent_track_id->at(i);
      dummy_sorted_track.ekin_vertex = the_global_data_dispatcher->get_input_root_data_structure()->ekin_vertex->at(i);
      dummy_sorted_track.creator_process = the_global_data_dispatcher->get_input_root_data_structure()->creator_process->at(i);
      
      mc_tracks.push_back(dummy_sorted_track);
    }
    
    sort(mc_tracks.begin(), mc_tracks.end(), compare_mc_tracks);
    
    /*
    cout << "Tracks:" << endl;
    for (int i = 0; i < mc_tracks.size(); i++)
    {
      cout << mc_tracks.at(i).track_id << " " << mc_tracks.at(i).particle_code << " " << mc_tracks.at(i).parent_track_id << " " << mc_tracks.at(i).ekin_vertex << " " << mc_tracks.at(i).creator_process << endl;
    }
    
    cout << "End Detailed Info:<--------------" << endl;
    */
  }
  
  return true;
}

bool track_finding_manager::find_tracks(vector<silicon_strip> cluster_lines[], vector<track>& the_tracks)
{
  track dummy_track;
  
  int current, previous;
  
  bool good_event = true;
  
  /*
  cout << "Track finder. Partial geometry: " << endl;
  
  for (int i = 0; i < the_geometry->get_number_of_plates(); i++)
  {
    cout << the_geometry->get_plate(i).strip_direction << " "
         << the_geometry->get_plate(i).normal_direction << " "
         << the_geometry->get_plate(i).intended_strip_direction << " "
         << the_geometry->get_plate(i).plate_type << " "
         << the_geometry->get_plate(i).plate_group << endl;
  }
  
  cout << "Intended strip directions: " << endl;
  
  for (int i = 0; i < the_geometry->the_intended_strip_directions.size(); i++)
  {
    cout << the_geometry->the_intended_strip_directions.at(i).intended_direction_type << " ";
    
    for (int j = 0; j < the_geometry->the_intended_strip_directions.at(i).plate_index_list.size(); j++)
    {
      cout << the_geometry->the_intended_strip_directions.at(i).plate_index_list.at(j) << " ";
    }
    
    cout << endl;
  }
  
  for (int i = 0; i < the_geometry->the_plate_groups.size(); i++)
  {
    cout << "Group " << i << ", type " << the_geometry->the_plate_groups.at(i).group_type << ", stream " << the_geometry->the_plate_groups.at(i).group_stream << ": ";
    
    for (int j = 0; j < the_geometry->the_plate_groups.at(i).plate_index_list.size(); j++)
    {
      cout << the_geometry->the_plate_groups.at(i).plate_index_list.at(j) << " ";
    }
    
    cout << endl;
  }
  */
   
  find_preliminary_event_type(cluster_lines);
  
  // according to the preliminary_event_type will look for:
  // single_track: 1 track from group 0 to last group
  // normal_vertex: 1 track from group 0 to last upstream; since there are no more increases,
  //   whatever tracks are in the last group from first mid to last, then go back to first mid group
  //   and find the partial tracks there to look for, including tracks with only vertex and first mid group
  // abnormal_vertex: same, but to where the vertex is
  // secondary_vertex: needs another look, because it is rarely clear what to do
  // unknown: ditto
  
  if (preliminary_event_type == "single_track")
  {
    dummy_track.assumed_start_group = 0;
    dummy_track.assumed_end_group = the_geometry->get_number_of_groups() - 1;
    the_tracks.push_back(dummy_track);
  }
  else if (preliminary_event_type == "normal_vertex" || preliminary_event_type == "abnormal_vertex")
  {
    dummy_track.assumed_start_group = 0;
    // dummy_track.assumed_end_group = the_geometry->get_number_of_upstream_groups() - 1;
    dummy_track.assumed_end_group = potential_vertex_before_group - 1;
    the_tracks.push_back(dummy_track);
    
    previous = group_cluster_multiplicity.at(the_geometry->get_number_of_groups() - 1);

    for (int i = 0; i < previous; i++)
    {
      // dummy_track.assumed_start_group = the_geometry->get_number_of_upstream_groups();
      dummy_track.assumed_start_group = potential_vertex_before_group;
      dummy_track.assumed_end_group = the_geometry->get_number_of_groups() - 1;
      the_tracks.push_back(dummy_track);
    }
    
    
    // for (int i = the_geometry->get_number_of_groups() - 2; i > the_geometry->get_number_of_upstream_groups() - 1; i--)
    for (int i = the_geometry->get_number_of_groups() - 2; i > potential_vertex_before_group - 1; i--)
    {
      current = group_cluster_multiplicity.at(i);
      
      for (int j = 0; j < current - previous; j++)
      {
        // dummy_track.assumed_start_group = the_geometry->get_number_of_upstream_groups();
        dummy_track.assumed_start_group = potential_vertex_before_group;
        dummy_track.assumed_end_group = i;
        // for the moment use tracks with length of one group, but investigate more;
        // these seem to be only useful when its the first midstream group
        the_tracks.push_back(dummy_track);
      }
      
      previous = current;
    }
  }
  else
  {
    good_event = false;
  }

  if (good_event)
  {
    good_event = assign_clusters_to_tracks(cluster_lines, the_tracks);
  }
  
  if (good_event)
  {
    /*
    cout << "Tracks: " << endl;
    for (int i = 0; i < the_tracks.size(); i++)
    {
      cout << "track " << i << ": " << the_tracks.at(i).assumed_start_group << " " << the_tracks.at(i).assumed_end_group << endl;
    }
    */
  }
  
  return good_event;
}

bool track_finding_manager::assign_clusters_to_tracks(vector<silicon_strip> cluster_lines[], vector<track>& the_tracks)
{
  bool good_event = true;
  silicon_strip a_vertex;
  
  // Single Triplet of Plates
  //
  
  int triplet_group_index;
  
  XYZVector midpoint;
  
  vector<silicon_strip> triplet_vertex_lines;
 
  XYZVector track_direction;
  XYZVector dummy_vertex;

  vector<double> initial_values_of_parameters;
  vector<double> errors_in_initial_values_of_parameters;
  
  int x_plate = -1;
  int y_plate = -1;
  int d_plate = -1;
  int downstreamest_triplet_plate = -1;
  string downstreamest_triplet_plate_type;
  
  if (find_guessed_vertex(cluster_lines, a_vertex))
  {
    the_output_data_structure->reco_guessed_vertex = a_vertex;
  }
  else
  {
    good_event = false;
  }
  
  if (the_global_data_dispatcher->get_reconstruction_options()->algorithm_type == 1)
  {
    triplet_group_index = the_global_data_dispatcher->get_reconstruction_options()->triplet_groups.at(0);
  }
  
  // these things have to be done at geometry initialization; here the values only should be read
  
  for (int i = 0; i < the_geometry->the_plate_groups.at(triplet_group_index).plate_index_list.size(); i++)
  {
    // cout << the_geometry->get_plate(the_geometry->the_plate_groups.at(triplet_group_index).plate_index_list.at(i)).plate_type << " ";
    
    if (the_geometry->get_plate(the_geometry->the_plate_groups.at(triplet_group_index).plate_index_list.at(i)).plate_type == "x")
    {
      x_plate = i;
    }
    else if (the_geometry->get_plate(the_geometry->the_plate_groups.at(triplet_group_index).plate_index_list.at(i)).plate_type == "y")
    {
      y_plate = i;
    }
    else if (the_geometry->get_plate(the_geometry->the_plate_groups.at(triplet_group_index).plate_index_list.at(i)).plate_type == "d1" || the_geometry->get_plate(the_geometry->the_plate_groups.at(triplet_group_index).plate_index_list.at(i)).plate_type == "d2")
    {
      d_plate = i;
    }
  }
    
  if (x_plate == -1 || y_plate == -1 || d_plate == -1)
  {
    good_event = false;
  }
    
  downstreamest_triplet_plate = y_plate;
  downstreamest_triplet_plate_type = "y";
    
  if (x_plate > downstreamest_triplet_plate)
  {
    downstreamest_triplet_plate = x_plate;
    downstreamest_triplet_plate_type = "x";
  }
  // cout << endl;
  
  if (good_event)
  {
    for (int i = 0; i < cluster_lines[x_plate].size(); i++)
    {
      for (int j = 0; j < cluster_lines[y_plate].size(); j++)
      {
        triplet_vertex_lines.clear();
        errors_in_initial_values_of_parameters.clear();
        
        triplet_vertex_lines.push_back(a_vertex);
        triplet_vertex_lines.push_back(cluster_lines[x_plate].at(i));
        triplet_vertex_lines.push_back(cluster_lines[y_plate].at(j));

        if (downstreamest_triplet_plate_type == "x")
        {
          midpoint = (the_geometry->get_target_position() + cluster_lines[downstreamest_triplet_plate].at(i).line_center) / 2.;
        }
        else if (downstreamest_triplet_plate_type == "y")
        {
          midpoint = (the_geometry->get_target_position() + cluster_lines[downstreamest_triplet_plate].at(j).line_center) / 2.;
        }
        else
        {
          cout << "Unknown downstreamest plate type." << endl;
        }
        
        my_straight_line_fcn the_fcn(triplet_vertex_lines, midpoint.z());
                    
        dummy_vertex.SetCoordinates((a_vertex.line_direction.x() / a_vertex.line_direction.z()) * a_vertex.line_center.z() + a_vertex.line_center.x(), (a_vertex.line_direction.y() / a_vertex.line_direction.z()) * a_vertex.line_center.z() + a_vertex.line_center.y(), 0.);
                    
        initial_values_of_parameters.clear();
        midpoint.SetCoordinates((dummy_vertex.x() + cluster_lines[y_plate].at(j).line_center.x()) / 2., (dummy_vertex.y() + cluster_lines[x_plate].at(i).line_center.y()) / 2., midpoint.z());
        initial_values_of_parameters.push_back(midpoint.x());
        initial_values_of_parameters.push_back(midpoint.y());
        
        if (downstreamest_triplet_plate_type == "x")
        {
          track_direction.SetCoordinates(dummy_vertex.x() - cluster_lines[y_plate].at(j).line_center.x(), dummy_vertex.x() - cluster_lines[x_plate].at(i).line_center.y(), dummy_vertex.z() - cluster_lines[downstreamest_triplet_plate].at(i).line_center.z());
        }
        else if (downstreamest_triplet_plate_type == "y")
        {
          track_direction.SetCoordinates(dummy_vertex.x() - cluster_lines[y_plate].at(j).line_center.x(), dummy_vertex.x() - cluster_lines[x_plate].at(i).line_center.y(), dummy_vertex.z() - cluster_lines[downstreamest_triplet_plate].at(j).line_center.z());
        }
        else
        {
          cout << "Unknown downstreamest plate type." << endl;
        }
        
        initial_values_of_parameters.push_back(track_direction.x() / track_direction.z());
        initial_values_of_parameters.push_back(track_direction.y() / track_direction.z());
        
        errors_in_initial_values_of_parameters.push_back(0.1);
        errors_in_initial_values_of_parameters.push_back(0.1);
        errors_in_initial_values_of_parameters.push_back(0.01);
        errors_in_initial_values_of_parameters.push_back(0.01);
        
        VariableMetricMinimizer the_minimizer;
        
        FunctionMinimum min = the_minimizer.Minimize(the_fcn, initial_values_of_parameters, errors_in_initial_values_of_parameters);
        
        cout << "Chi Square: " << min.Fval() << endl;
        cout << min.UserState().Value(0) << " " << min.UserState().Error(0) << endl;
        cout << min.UserState().Value(1) << " " << min.UserState().Error(1) << endl;
        cout << min.UserState().Value(2) << " " << min.UserState().Error(2) << endl;
        cout << min.UserState().Value(3) << " " << min.UserState().Error(3) << endl;
        
        // calculate the intersection of this line with the d plate
        
        // calculate the strip fired on the d plate
        
        // from the list of strip hits on the d plate select the closest and mark it as used
      }
    }
  }

  
  
  return good_event;
}


bool track_finding_manager::find_guessed_vertex(vector<silicon_strip> cluster_lines[], silicon_strip& a_vertex)
{
  bool good_event = true;
  bool single_cluster_on_all_upstream_plates = true;
 
  XYZVector dummy_vector;
  
  XYZVector midpoint;
  XYZVector track_direction;
  
  vector<int> plate_strips_along_x;
  vector<int> plate_strips_along_y;
  
  vector<silicon_strip> upstream_cluster_lines;
  
  upstream_cluster_lines.clear();
  
  plate_strips_along_x.clear();
  plate_strips_along_y.clear();
  
  for (int i = 0; i < the_geometry->get_number_of_upstream_plates(); i++)
  {
    if (cluster_lines[i].size() != 1)
    {
      single_cluster_on_all_upstream_plates = false;
    }
  }
  
  if (single_cluster_on_all_upstream_plates)
  {
    for (int i = 0; i < the_geometry->get_number_of_upstream_plates(); i++)
    {
      upstream_cluster_lines.push_back(cluster_lines[i].at(0));
    
    
      if (the_geometry->get_plate(i).plate_type == "x")
      {
        plate_strips_along_x.push_back(i);
      }
      else if (the_geometry->get_plate(i).plate_type == "y")
      {
        plate_strips_along_y.push_back(i);
      }
    }
  }
  
  if (single_cluster_on_all_upstream_plates)
  {
    for (int i = 0; i < upstream_cluster_lines.size(); i++)
    {
      // cout << upstream_cluster_lines.at(i).line_center << " " << upstream_cluster_lines.at(i).line_direction << endl;
    }
    
  
    // if (single_cluster_on_all_upstream_plates) cout << "!!!!!!!SINGLE" << endl;
  
    midpoint = (upstream_cluster_lines.at(0).line_center + upstream_cluster_lines.at(upstream_cluster_lines.size() - 1).line_center) / 2.;
    
    my_straight_line_fcn the_fcn(upstream_cluster_lines, midpoint.z());

    vector<double> initial_values_of_parameters;
    initial_values_of_parameters.push_back(midpoint.x());
    initial_values_of_parameters.push_back(midpoint.y());
    
    track_direction = upstream_cluster_lines.at(plate_strips_along_y.at(0)).line_center - upstream_cluster_lines.at(plate_strips_along_y.at(plate_strips_along_y.size() - 1)).line_center;
    initial_values_of_parameters.push_back(track_direction.x() / track_direction.z());
    
    track_direction = upstream_cluster_lines.at(plate_strips_along_x.at(0)).line_center - upstream_cluster_lines.at(plate_strips_along_x.at(plate_strips_along_x.size() - 1)).line_center;
    initial_values_of_parameters.push_back(track_direction.y() / track_direction.z());
  
    vector<double> errors_in_initial_values_of_parameters;
    errors_in_initial_values_of_parameters.push_back(0.1);
    errors_in_initial_values_of_parameters.push_back(0.1);
    errors_in_initial_values_of_parameters.push_back(0.01);
    errors_in_initial_values_of_parameters.push_back(0.01);
  
    VariableMetricMinimizer the_minimizer;
  
    FunctionMinimum min = the_minimizer.Minimize(the_fcn, initial_values_of_parameters, errors_in_initial_values_of_parameters);
  
    /*
    cout << "Chi Square: " << min.Fval() << endl;
    cout << min.UserState().Value(0) << " " << min.UserState().Error(0) << endl;
    cout << min.UserState().Value(1) << " " << min.UserState().Error(1) << endl;
    cout << min.UserState().Value(2) << " " << min.UserState().Error(2) << endl;
    cout << min.UserState().Value(3) << " " << min.UserState().Error(3) << endl;
    
    cout << min << endl;
    */
    
    // a_vertex.SetCoordinates(min.UserState().Value(2) * midpoint.z() + min.UserState().Value(0), min.UserState().Value(3) * midpoint.z() + min.UserState().Value(1), 0.);
    // cout << "Guessed reco vertex: " << a_vertex << endl;
    
    a_vertex.line_center.SetCoordinates(min.UserState().Value(0), min.UserState().Value(1), midpoint.z());
    a_vertex.line_direction.SetCoordinates(min.UserState().Value(2), min.UserState().Value(3), 1.);
    a_vertex.line_direction = a_vertex.line_direction.Unit();
  }
  
  return good_event;
}
