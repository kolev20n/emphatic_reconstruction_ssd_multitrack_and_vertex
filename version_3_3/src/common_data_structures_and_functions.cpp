#include <map>

#include <TMath.h>
#include <Math/Vector3D.h>
#include <Math/GenVector/Rotation3D.h>
#include <Math/GenVector/RotationX.h>
#include <Math/GenVector/RotationY.h>
#include <Math/GenVector/RotationZ.h>
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameterState.h"
#include "Minuit2/MnPrint.h"
#include "Minuit2/MnMigrad.h"

#include "common_data_structures_and_functions.h"
#include "my_straight_line_fcn.h"

using namespace std;
using namespace ROOT::Math;
using namespace ROOT::Minuit2;

input_root_data_structure::input_root_data_structure()
:plate_number(0), strip_number(0), total_energy_values(0), non_ionization_energy_values(0), contributing_tracks(0),
    track_id(0), particle_code(0),
    parent_track_id(0), creator_process(0), x_vertex(0),
    y_vertex(0), z_vertex(0), px_vertex(0),
    py_vertex(0), pz_vertex(0), ekin_vertex(0),
    visited_plates(0), step_track_id(0), step_volume_name(0),
    step_plate_id(0), step_x_i(0), step_y_i(0), step_z_i(0),
    step_px_i(0), step_py_i(0), step_pz_i(0),
    step_x_f(0), step_y_f(0), step_z_f(0),
    step_px_f(0), step_py_f(0), step_pz_f(0),
    step_initial_energy(0), step_final_energy(0), step_process_name(0)
{
  cout << "Created the input root data structure." << endl;
}

input_root_data_structure::~input_root_data_structure()
{
}

output_root_data_structure::output_root_data_structure()
{
  cout << "Created the output root data structure." << endl;
}

output_root_data_structure::~output_root_data_structure()
{
}

int partial_geometry::get_number_of_plates()
{
  return number_of_plates;
}

partial_geometry::partial_geometry()
:number_of_plates(0), number_of_upstream_plates(0), number_of_midstream_plates(0), number_of_downstream_plates(0)
{
  the_intended_strip_directions.clear();
  the_plate_groups.clear();
  
  cout << "Created an empty geometry." << endl;
}

void partial_geometry::add_plate(int a_plate_id, int a_plate_group, XYZVector an_intended_strip_direction, int a_number_of_strips, double an_efficiency, XYZVector a_position, XYZVector a_rotation, XYZVector a_size)
{
  if (number_of_plates == max_number_of_plates)
  {
    cout << "Max number of ssd plates is being exceeded: " << max_number_of_plates + 1 << endl;
    cout << "Ignoring current plate..." << endl;
    
    return;
  }
  
  the_plates[number_of_plates].plate_id = a_plate_id;
  the_plates[number_of_plates].plate_group = a_plate_group;
  the_plates[number_of_plates].intended_strip_direction = an_intended_strip_direction.Unit();
  the_plates[number_of_plates].number_of_strips = a_number_of_strips;
  the_plates[number_of_plates].efficiency = an_efficiency;
  the_plates[number_of_plates].position = a_position;
  the_plates[number_of_plates].rotation = a_rotation;
  the_plates[number_of_plates].size = a_size;
  number_of_plates++;
}

void partial_geometry::add_target(XYZVector a_position, XYZVector a_rotation, XYZVector a_size, string a_material)
{
  the_target.position = a_position;
  the_target.rotation = a_rotation;
  the_target.size = a_size;
  the_target.material = a_material;
}

ssd_plate partial_geometry::get_plate(int a_plate_number)
{
  if (a_plate_number >= number_of_plates)
  {
    cout << "Requested ssd plate outside the existing number of plates: " << a_plate_number << " from " << number_of_plates << endl;
    cout << "Quitting..." << endl;
    exit(EXIT_FAILURE);
  }
  else
  {
    return the_plates[a_plate_number];
  }
}

void partial_geometry::classify_plates()
{
  const double to_radians = Pi() / 180.;

  // all these are hard coded for the moment
  
  int count_upstream = 0;
  int count_midstream = 0;
  int count_downstream = 0;
  
  RotationX rx;
  RotationY ry;
  RotationZ rz;
  
  Rotation3D total_rotation;
  
  XYZVector trial_vector;

  for (int i = 0; i < number_of_plates; i++)
  {
    if (fabs(the_plates[i].rotation.z()) < 5.) // arbitrary magic number, but allows small misalignments of up to 5 degrees
    {
      x_plates.push_back(i);
    }
    else if (fabs(the_plates[i].rotation.z()) < 95. && fabs(the_plates[i].rotation.z()) > 85.)
    {
      y_plates.push_back(i);
    }
    else if (fabs(the_plates[i].rotation.z()) < 50. && fabs(the_plates[i].rotation.z()) > 40.)
    {
      d_plates.push_back(i);
    }
    
    if (the_plates[i].position.z() < the_target.position.z())
    {
      count_upstream++;
    }
    else if (the_plates[i].position.z() < 16.) // this is a hard coded position of the magnet; change later
    {
      count_midstream++;
    }
    else
    {
      count_downstream++;
    }
    
    rx.SetAngle(the_plates[i].rotation.x() * to_radians);
    ry.SetAngle(the_plates[i].rotation.y() * to_radians);
    rz.SetAngle(the_plates[i].rotation.z() * to_radians);
    
    total_rotation = rz * ry * rx;
    
    the_plates[i].strip_direction.SetCoordinates(0., 1., 0.);
    the_plates[i].strip_direction = total_rotation * the_plates[i].strip_direction;
    
    // hard-coded 2 degrees
    trial_vector.SetCoordinates(1., 0., 0.);

    if (fabs(the_plates[i].strip_direction.Dot(trial_vector)) > 0.99939)
    {
      the_plates[i].plate_type = 1;
      x_plates.push_back(i);
    }
    else
    {
      trial_vector.SetCoordinates(0., 1., 0.);
      
      if (fabs(the_plates[i].strip_direction.Dot(trial_vector)) > 0.99939)
      {
        the_plates[i].plate_type = 2;
        y_plates.push_back(i);
      }
      
      else
      {
        trial_vector.SetCoordinates(1., -1., 0.);
        trial_vector = trial_vector.Unit();
        
        if (fabs(the_plates[i].strip_direction.Dot(trial_vector)) > 0.99939)
        {
          the_plates[i].plate_type = 6;
          d_plates.push_back(i);
        }
        else
        {
          cout << "Uknown ssd plate strip direction: " << i << " " <<  the_plates[i].strip_direction << endl << "Cannot perform reasonable reconstruction. Quitting..." << endl;
          exit(EXIT_FAILURE);
        }
      }
    }
    
    the_plates[i].normal_direction.SetCoordinates(0., 0., 1.);
    the_plates[i].normal_direction = total_rotation * the_plates[i].normal_direction;
  }
  
  number_of_upstream_plates = count_upstream;
  number_of_midstream_plates = count_midstream;
  number_of_downstream_plates = count_downstream;
  
  cout << "<><><><><>" << endl;
  for (int i = 0; i < number_of_plates; i++)
  {
    cout << the_plates[i].strip_direction << " " << the_plates[i].normal_direction << " " << the_plates[i].intended_strip_direction << " " << the_plates[i].plate_type << " " << the_plates[i].plate_group << endl;
  }
  
  // (1, 0, 0) called x
  // (0, 1, 0) called y
  // (1, 1, 0) called d1
  // (1, -1, 0) called d2
  // (x, y, 0) arbitrary, but in the xy plane called a1
  // (x, -y, 0) ditto, called a2
  //
  // consider primary directions x and y, and secondary directions d1 and d2
  // then
  
  vector<int> x_plates;
  vector<int> y_plates;
  vector<int> d1_plates;
  vector<int> d2_plates;

  x_plates.clear();
  y_plates.clear();
  d1_plates.clear();
  d2_plates.clear();
  
  XYZVector v_x(1., 0., 0.);
  XYZVector v_y(0., 1., 0.);
  XYZVector v_d1(1., 1., 0.);
  XYZVector v_d2(1., -1., 0.);
  v_d1 = v_d1.Unit();
  v_d2 = v_d2.Unit();
  
  for (int i = 0; i < number_of_plates; i++)
  {
    if (fabs(the_plates[i].intended_strip_direction.Dot(v_x)) > 0.9999) // arbitrary cut
    {
      x_plates.push_back(i);
      the_plates[i].plate_type = "x";
    }
    else if (fabs(the_plates[i].intended_strip_direction.Dot(v_y)) > 0.9999) // arbitrary cut
    {
      y_plates.push_back(i);
      the_plates[i].plate_type = "y";
    }
    if (fabs(the_plates[i].intended_strip_direction.Dot(v_d1)) > 0.9999) // arbitrary cut
    {
      d1_plates.push_back(i);
      the_plates[i].plate_type = "d1";
    }
    if (fabs(the_plates[i].intended_strip_direction.Dot(v_d2)) > 0.9999) // arbitrary cut
    {
      d2_plates.push_back(i);
      the_plates[i].plate_type = "d2";
    }
    
    // cout << the_plates[i].plate_type << endl;
  }
  
  intended_strip_direction dummy_strip_direction;
  
  dummy_strip_direction.intended_direction_type = "x";
  dummy_strip_direction.plate_index_list = x_plates;
  
  // order is important; so 0 is along x, 1 is along y, 2 is along 110, 3 is along 1-10, etc.
  the_intended_strip_directions.push_back(dummy_strip_direction);
  
  dummy_strip_direction.intended_direction_type = "y";
  dummy_strip_direction.plate_index_list = y_plates;

  the_intended_strip_directions.push_back(dummy_strip_direction);

  dummy_strip_direction.intended_direction_type = "d1";
  dummy_strip_direction.plate_index_list = d1_plates;
  
  the_intended_strip_directions.push_back(dummy_strip_direction);
  
  dummy_strip_direction.intended_direction_type = "d2";
  dummy_strip_direction.plate_index_list = d2_plates;
  
  the_intended_strip_directions.push_back(dummy_strip_direction);
  
  cout << "Intended strip directions: " << endl;
  for (int i = 0; i < the_intended_strip_directions.size(); i++)
  {
    cout << the_intended_strip_directions.at(i).intended_direction_type << " ";
    
    for (int j = 0; j < the_intended_strip_directions.at(i).plate_index_list.size(); j++)
    {
      cout << the_intended_strip_directions.at(i).plate_index_list.at(j) << " ";
    }
    
    cout << endl;
  }
  
  map<int, int> group_ids;
  group_ids.clear();
  
  vector<int> group_id_lookup;
  group_id_lookup.clear();
  
  plate_group dummy_group;

  for (int i = 0; i < number_of_plates; i++)
  {
    cout << group_ids[the_plates[i].plate_group] << endl;
    
    if (group_ids[the_plates[i].plate_group] == 0)
    {
      group_ids[the_plates[i].plate_group]++;
      group_id_lookup.push_back(the_plates[i].plate_group);
    }
  }
  
  for (int i = 0; i < group_id_lookup.size(); i++)
  {
    dummy_group.group_id = group_id_lookup.at(i);
    
    the_plate_groups.push_back(dummy_group);
  }
  
  // this is not optimal...
  for (int i = 0; i < number_of_plates; i++)
  {
    for (int j = 0; j < group_id_lookup.size(); j++)
    {
      if (the_plates[i].plate_group == group_id_lookup.at(j))
      {
        the_plate_groups[j].plate_index_list.push_back(i);
      }
    }
  }
  
  for (int i = 0; i < the_plate_groups.size(); i++)
  {
    // relies that all plates in a group are in the same stream; which should be always true
    //
    if (the_plate_groups[i].plate_index_list[0] < number_of_upstream_plates)
    {
      the_plate_groups[i].group_stream = "upstream";
    }
    else if (the_plate_groups[i].plate_index_list[0] < number_of_upstream_plates + number_of_midstream_plates)
    {
      the_plate_groups[i].group_stream = "midstream";
    }
    else
    {
      the_plate_groups[i].group_stream = "downstream";
    }
    
    if (the_plate_groups[i].plate_index_list.size() == 2)
    {
      if (the_plates[the_plate_groups[i].plate_index_list[0]].plate_type == "x" &&  the_plates[the_plate_groups[i].plate_index_list[1]].plate_type == "y" ||
          the_plates[the_plate_groups[i].plate_index_list[0]].plate_type == "y" &&  the_plates[the_plate_groups[i].plate_index_list[1]].plate_type == "x")
      {
        the_plate_groups[i].group_type = 1;
      }
      else if (the_plates[the_plate_groups[i].plate_index_list[0]].plate_type == "d1" &&  the_plates[the_plate_groups[i].plate_index_list[1]].plate_type == "d2" ||
          the_plates[the_plate_groups[i].plate_index_list[0]].plate_type == "d2" &&  the_plates[the_plate_groups[i].plate_index_list[1]].plate_type == "d1")
      {
        the_plate_groups[i].group_type = 2;
      }
      else
      {
        cout << "Unknown group type of 2 plates. Quitting..." << endl;
        exit(EXIT_FAILURE);
      }
    }
    else if (the_plate_groups[i].plate_index_list.size() == 3)
    {
      int count_x = 0, count_y = 0, count_d1 = 0, count_d2 = 0;
      string temp_string;
      
      for (int j = 0; j < 3; j++)
      {
        temp_string = the_plates[the_plate_groups[i].plate_index_list[j]].plate_type;
        
        cout << temp_string << endl;
        
        if (temp_string == "x")
        {
          count_x++;
        }
        else if (temp_string == "y")
        {
          count_y++;
        }
        else if (temp_string == "d1")
        {
          count_d1++;
        }
        else if (temp_string == "d2")
        {
          count_d2++;
        }
      }
      
      if (count_x == 1 && count_y == 1 && count_d1 == 1)
      {
        the_plate_groups[i].group_type = 3;
      }
      else if (count_x == 1 && count_y == 1 && count_d2 == 1)
      {
        the_plate_groups[i].group_type = 4;
      }
      else if (count_x == 1 && count_d1 == 1 && count_d2 == 1)
      {
        the_plate_groups[i].group_type = 5;
      }
      else if (count_y == 1 && count_d1 == 1 && count_d2 == 1)
      {
        the_plate_groups[i].group_type = 6;
      }
      else
      {
        cout << "Unknown group type of 3 plates. Quitting..." << endl;
        //exit(EXIT_FAILURE);
      }
    }
  }
  
  number_of_upstream_groups = 0;
  number_of_midstream_groups = 0;
  number_of_downstream_groups = 0;
  
  for (int i = 0; i < the_plate_groups.size(); i++)
  {
    if (the_plate_groups.at(i).group_stream == "upstream")
    {
      number_of_upstream_groups++;
    }
    else if (the_plate_groups.at(i).group_stream == "midstream")
    {
      number_of_midstream_groups++;
    }
    else
    {
      number_of_downstream_groups++;
    }
  }
  
  for (int i = 0; i < the_plate_groups.size(); i++)
  {
    cout << "Group " << i << " type " << the_plate_groups.at(i).group_type << ": ";
    
    for (int j = 0; j < the_plate_groups.at(i).plate_index_list.size(); j++)
    {
      cout << the_plate_groups.at(i).plate_index_list.at(j) << " ";
    }
    
    cout << endl;
  }
  
  return;
}

int partial_geometry::get_number_of_upstream_plates()
{
  return number_of_upstream_plates;
}

int partial_geometry::get_number_of_midstream_plates()
{
  return number_of_midstream_plates;
}

int partial_geometry::get_number_of_downstream_plates()
{
  return number_of_downstream_plates;
}

int partial_geometry::get_number_of_groups()
{
  return (number_of_upstream_groups + number_of_midstream_groups + number_of_downstream_groups);
}

int partial_geometry::get_number_of_upstream_groups()
{
  return number_of_upstream_groups;
}

int partial_geometry::get_number_of_midstream_groups()
{
  return number_of_midstream_groups;
}

int partial_geometry::get_number_of_downstream_groups()
{
  return number_of_downstream_groups;
}

vector<int> partial_geometry::get_x_plates()
{
  return x_plates;
}

vector<int> partial_geometry::get_y_plates()
{
  return y_plates;
}

vector<int> partial_geometry::get_d_plates()
{
  return d_plates;
}

XYZVector partial_geometry::get_target_position()
{
  return the_target.position;
}

bool compare_hit_pairs(const hit_pair& first, const hit_pair& second)
{
  return (first.strip_number < second.strip_number);
}

bool compare_map(const pair<string, int> &a, const pair<string, int> &b)
{
  return (a.second > b.second);
}

bool compare_mc_tracks(const sorted_track &a, const sorted_track &b)
{
  return a.track_id < b.track_id;
}

bool compare_z_lines_3d(const line_3d &a, const line_3d &b)
{
  return (a.point.z() < b.point.z());
}

bool geometry_fitters::find_track_3d(std::vector<line_3d> lines, line_3d& result_line, int configuration_key)
{
  XYZVector midpoint;
  XYZVector track_direction;

  vector<double> initial_values_of_parameters;
  vector<double> errors_in_initial_values_of_parameters;

  initial_values_of_parameters.clear();
  errors_in_initial_values_of_parameters.clear();
  
  int number_of_lines = lines.size();

  if (lines.size() < 3)
  {
    cout << "find_track_3d: too few skew lines to fit..." << endl;
    return false;
  }
  
  // they should already be ordered, so this is not needed
  // sort(lines.begin(), lines.end(), compare_z_lines_3d);
  
  // relying that there are no significant angles of rotation about x or y
  // or significant misalignments

  if (configuration_key == 1) // two xy groups in this order
  {
    if (lines.size() < 4)
    {
      cout << "find_track_3d: configuration_key = 1 needs four lines..." << endl;
      return false;
    }
    
    midpoint.SetCoordinates((lines.at(1).point.x() + lines.at(3).point.x()) / 2., (lines.at(0).point.x() + lines.at(2).point.x()) / 2., (lines.at(0).point.z() + lines.at(1).point.z() + lines.at(2).point.z() + lines.at(3).point.z())/ 4.);
    
    initial_values_of_parameters.push_back(midpoint.x());
    initial_values_of_parameters.push_back(midpoint.y());

    track_direction = lines.at(3).point - lines.at(1).point;
    initial_values_of_parameters.push_back(track_direction.x() / track_direction.z());
    track_direction = lines.at(2).point - lines.at(0).point;
    initial_values_of_parameters.push_back(track_direction.y() / track_direction.z());

  
    my_straight_line_fcn the_fcn(lines, midpoint.z());
  
    errors_in_initial_values_of_parameters.push_back(0.1);
    errors_in_initial_values_of_parameters.push_back(0.1);
    errors_in_initial_values_of_parameters.push_back(0.01);
    errors_in_initial_values_of_parameters.push_back(0.01);
  
    VariableMetricMinimizer the_minimizer;
  
    FunctionMinimum min = the_minimizer.Minimize(the_fcn, initial_values_of_parameters, errors_in_initial_values_of_parameters);
  
    result_line.point.SetCoordinates(min.UserState().Value(0), min.UserState().Value(1), midpoint.z());
    result_line.direction.SetCoordinates(min.UserState().Value(2), min.UserState().Value(3), 1.);
    result_line.direction = result_line.direction.Unit();
  }
  
  return true;
}
