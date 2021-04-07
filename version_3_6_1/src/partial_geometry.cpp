#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <utility>

#include <TMath.h>
#include <Math/Vector3D.h>
#include <Math/GenVector/Rotation3D.h>
#include <Math/GenVector/RotationX.h>
#include <Math/GenVector/RotationY.h>
#include <Math/GenVector/RotationZ.h>

#include "common_data_structures_and_functions.h"
#include "global_data_dispatcher.h"
#include "partial_geometry.h"

using namespace std;
using namespace ROOT::Math;

partial_geometry::partial_geometry()
:actual_number_of_plates(0), number_of_tracking_regions(0), number_of_groups(0), primary_direction_1_name(""),
    primary_direction_2_name(""), secondary_direction_1_name(""),
    secondary_direction_2_name("")
{
  cout << "Created an empty geometry." << endl;
}

int partial_geometry::get_number_of_plates()
{
  return actual_number_of_plates;
}

int partial_geometry::get_number_of_plates_in_tracking_region(int a_region)
{
  if (a_region < number_of_tracking_regions)
  {
    return the_tracking_regions[a_region].number_of_plates;
  }
  else
  {
    return -1;
  }
}

int partial_geometry::get_number_of_tracking_regions()
{
  return number_of_tracking_regions;
}

int partial_geometry::get_number_of_groups()
{
  return number_of_groups;
}

int partial_geometry::get_number_of_groups_in_tracking_region(int a_region)
{
  if (a_region < number_of_tracking_regions)
  {
    return the_tracking_regions[a_region].number_of_groups;
  }
  else
  {
    return -1;
  }
}

int partial_geometry::get_number_of_plates_in_group(int a_group)
{
  if (a_group < number_of_groups)
  {
    return the_plate_groups[a_group].number_of_plates_in_group;
  }
  else
  {
    return -1;
  }
}

int partial_geometry::get_number_of_triple_groups()
{
  return number_of_triple_groups;
}

target* partial_geometry::get_target()
{
  return &the_target;
}

XYZVector partial_geometry::get_primary_direction_1()
{
  return primary_direction_1;
}

XYZVector partial_geometry::get_primary_direction_2()
{
  return primary_direction_2;
}

XYZVector partial_geometry::get_secondary_direction_1()
{
  return primary_direction_2;
}

XYZVector partial_geometry::get_secondary_direction_2()
{
  return secondary_direction_2;
}

string partial_geometry::get_primary_direction_1_name()
{
  return primary_direction_1_name;
}

string partial_geometry::get_primary_direction_2_name()
{
  return primary_direction_2_name;
}

string partial_geometry::get_secondary_direction_1_name()
{
  return secondary_direction_1_name;
}

string partial_geometry::get_secondary_direction_2_name()
{
  return secondary_direction_2_name;
}

tracking_region* partial_geometry::get_tracking_region(int a_region)
{
  if (a_region < number_of_tracking_regions)
  {
    return &the_tracking_regions[a_region];
  }
  else
  {
    return nullptr;
  }
}

plate_group* partial_geometry::get_plate_group(int a_group)
{
  if (a_group < number_of_groups)
  {
    return &the_plate_groups[a_group];
  }
  else
  {
    return nullptr;
  }
}

ssd_plate* partial_geometry::get_plate(int a_plate_index)
{
  if (a_plate_index < actual_number_of_plates)
  {
    return &the_plates[a_plate_index];
  }
  else
  {
    return nullptr;
  }
}

ssd_plate* partial_geometry::get_plate_by_id(int a_plate_id)
{
  // implement later
  
  return nullptr;
}

int partial_geometry::get_triple_plate_group_index(int a_group_index)
{
  if (a_group_index < number_of_triple_groups)
  {
    return the_triple_groups[a_group_index];
  }
  else
  {
    cout << "Triple group index requested above the number of triple groups: " << a_group_index << endl;
    return -1;
  }
}

bool partial_geometry::add_plate(ssd_plate a_plate)
{
  if (actual_number_of_plates == max_number_of_plates)
  {
    cout << "Requested plate to add above maxx allowed number. Ignoring..." << endl;
    return false;
  }
  else
  {
    the_plates[actual_number_of_plates] = a_plate;
    actual_number_of_plates++;
    return true;
  }
}

bool partial_geometry::add_target(target a_target)
{
  the_target = a_target;
  
  return true;
}

bool partial_geometry::set_primary_direction_1(XYZVector a_direction)
{
  primary_direction_1 = a_direction;
  primary_direction_1 = primary_direction_1.Unit();
  
  return true;
}

bool partial_geometry::set_primary_direction_2(XYZVector a_direction)
{
  primary_direction_2 = a_direction;
  primary_direction_2 = primary_direction_2.Unit();
  
  return true;
}

bool partial_geometry::set_secondary_direction_1(XYZVector a_direction)
{
  secondary_direction_1 = a_direction;
  secondary_direction_1 = secondary_direction_1.Unit();
  
  return true;
}

bool partial_geometry::set_secondary_direction_2(XYZVector a_direction)
{
  secondary_direction_2 = a_direction;
  secondary_direction_2 = secondary_direction_2.Unit();
  
  return true;
}

bool partial_geometry::set_primary_direction_1_name(string a_name)
{
  primary_direction_1_name = a_name;
  
  return true;
}

bool partial_geometry::set_primary_direction_2_name(string a_name)
{
  primary_direction_2_name = a_name;
  
  return true;
}

bool partial_geometry::set_secondary_direction_1_name(string a_name)
{
  secondary_direction_1_name = a_name;
  
  return true;
}

bool partial_geometry::set_secondary_direction_2_name(string a_name)
{
  secondary_direction_2_name = a_name;
  
  return true;
}

bool partial_geometry::set_tracking_region_names(int a_region, std::string a_name)
{
  if (a_region < max_number_of_tracking_regions)
  {
    tracking_region_names[a_region] = a_name;
    
    return true;
  }
  else
  {
    return false;
  }
}

bool partial_geometry::classify_plates()
{
  bool classification_succeeded = true;
  
  // const double to_radians = Pi() / 180.;
  double dummy_double;
  
  string dummy_string;
  
  RotationX rx;
  RotationY ry;
  RotationZ rz;
  
  Rotation3D total_rotation;
  
  map<int, int> group_index_map;
  group_index_map.clear();
  
  vector<int> group_index_vector;
  group_index_vector.clear();
  
  XYZVector dummy_vector_1, dummy_vector_2;
  
  number_of_triple_groups = 0;
  
  dummy_double = the_global_data_dispatcher->get_reconstruction_options()->maximum_cosangle_between_intended_and_actual_strip_direction;
  
  for (int i = 0; i < max_number_of_tracking_regions; i++)
  {
    the_tracking_regions[i].number_of_groups = 0;
    the_tracking_regions[i].number_of_plates = 0;
  }
  
  for (int i = 0; i < max_number_of_tracking_regions; i++)
  {
    if (tracking_region_names[i] != "NA")
    {
      number_of_tracking_regions++;
    }
  }
  
  for (int i = 0; i < number_of_tracking_regions; i++)
  {
    the_tracking_regions[i].name = tracking_region_names[i];
  }
  
  for (int i = 0; i < actual_number_of_plates; i++)
  {
    if (group_index_map[the_plates[i].plate_group] == 0)
    {
      group_index_map[the_plates[i].plate_group]++;
      group_index_vector.push_back(the_plates[i].plate_group);
    }
    
    the_plates[i].pitch = the_plates[i].size.y() / ((double) the_plates[i].number_of_strips);
  }
  
  for (int i = 0; i < group_index_vector.size(); i++)
  {
    if (group_index_vector[i] != i)
    {
      cout << "Wrong order of plate groups (must be consequtive integers starting at 0)..." << endl;
      classification_succeeded = false;
    }
  }
  
  number_of_groups = group_index_vector.size();
  
  for (int i = 0; i < max_number_of_groups; i++)
  {
    the_plate_groups[i].number_of_plates_in_group = 0;
  }
  
  for (int i = 0; i < actual_number_of_plates; i++)
  {
    rx.SetAngle(the_plates[i].rotation.x() * to_radians);
    ry.SetAngle(the_plates[i].rotation.y() * to_radians);
    rz.SetAngle(the_plates[i].rotation.z() * to_radians);
    
    total_rotation = rz * ry * rx;
    
    // the_plates[i].strip_direction.SetCoordinates(0., 1., 0.);
    the_plates[i].strip_direction = primary_direction_2;
    the_plates[i].strip_direction = total_rotation * the_plates[i].strip_direction;
    
    if (the_plates[i].plate_type == primary_direction_1_name)
    {
      dummy_vector_1 = primary_direction_1.Unit();
      dummy_vector_2 = the_plates[i].strip_direction.Unit();
      
      if (fabs(dummy_vector_1.Dot(dummy_vector_2)) < dummy_double)
      {
        cout << "Plate intended_to_measure and rotation strip directions do not coincide: " << i << endl;
        classification_succeeded = false;
      }
    }
    else if (the_plates[i].plate_type == primary_direction_2_name)
    {
      dummy_vector_1 = primary_direction_2.Unit();
      dummy_vector_2 = the_plates[i].strip_direction.Unit();

      if (fabs(dummy_vector_1.Dot(dummy_vector_2)) < dummy_double)
      {
        cout << "Plate intended_to_measure and rotation strip directions do not coincide: " << i << endl;
        classification_succeeded = false;
      }
    }
    else if (the_plates[i].plate_type == secondary_direction_1_name)
    {
      dummy_vector_1 = secondary_direction_1.Unit();
      dummy_vector_2 = the_plates[i].strip_direction.Unit();
      
      if (fabs(dummy_vector_1.Dot(dummy_vector_2)) < dummy_double)
      {
        cout << "Plate intended_to_measure and rotation strip directions do not coincide: " << i << endl;
        classification_succeeded = false;
      }
    }
    else if (the_plates[i].plate_type == secondary_direction_2_name)
    {
      dummy_vector_1 = secondary_direction_2.Unit();
      dummy_vector_2 = the_plates[i].strip_direction.Unit();
      
      if (fabs(dummy_vector_1.Dot(dummy_vector_2)) < dummy_double)
      {
        cout << "Plate intended_to_measure and rotation strip directions do not coincide: " << i << endl;
        classification_succeeded = false;
      }
    }
    else
    {
      cout << "Unknown intended_to_measure direction..." << endl;
      classification_succeeded = false;
    }
    
    the_plates[i].normal_direction.SetCoordinates(0., 0., 1.);
    the_plates[i].normal_direction = total_rotation * the_plates[i].normal_direction;
    
    the_plates[i].plate_plane_equation.a = the_plates[i].normal_direction.x();
    the_plates[i].plate_plane_equation.b = the_plates[i].normal_direction.y();
    the_plates[i].plate_plane_equation.c = the_plates[i].normal_direction.z();
    the_plates[i].plate_plane_equation.d = -(the_plates[i].normal_direction.x() * the_plates[i].position.x() +
                                            the_plates[i].normal_direction.y() * the_plates[i].position.y() +
                                            the_plates[i].normal_direction.z() * the_plates[i].position.z());
                                            
    if (the_plates[i].tracking_region < number_of_tracking_regions)
    {
      the_tracking_regions[the_plates[i].tracking_region].plate_index_list[the_tracking_regions[the_plates[i].tracking_region].number_of_plates] = i;
      the_tracking_regions[the_plates[i].tracking_region].number_of_plates++;
    }
    else
    {
      cout << "Invalid tracking region: " << the_plates[i].tracking_region << " for plate: " << i << endl;
      classification_succeeded = false;
    }
    
    if (the_plates[i].plate_group < number_of_groups)
    {
      the_plate_groups[the_plates[i].plate_group].plate_index_list[the_plate_groups[the_plates[i].plate_group].number_of_plates_in_group] = i;
      the_plate_groups[the_plates[i].plate_group].number_of_plates_in_group++;
    }
    else
    {
      cout << "Invalid group: " << the_plates[i].plate_group << endl;
      classification_succeeded = false;
    }
  }
  
  for (int i = 0; i < number_of_groups; i++)
  {
    if (the_plate_groups[i].number_of_plates_in_group == 3)
    {
      the_triple_groups[number_of_triple_groups] = i;
      number_of_triple_groups++;
      
      dummy_string  = the_plates[the_plate_groups[i].plate_index_list[0]].plate_type;
      dummy_string += "_";
      dummy_string += the_plates[the_plate_groups[i].plate_index_list[1]].plate_type;
      dummy_string += "_";
      dummy_string += the_plates[the_plate_groups[i].plate_index_list[2]].plate_type;
      the_plate_groups[i].group_type = dummy_string;
      
      if (the_plates[the_plate_groups[i].plate_index_list[0]].tracking_region == the_plates[the_plate_groups[i].plate_index_list[1]].tracking_region &&
          the_plates[the_plate_groups[i].plate_index_list[0]].tracking_region == the_plates[the_plate_groups[i].plate_index_list[2]].tracking_region)
      {
        the_plate_groups[i].tracking_region = the_plates[the_plate_groups[i].plate_index_list[0]].tracking_region;
      }
      else
      {
        cout << "Plates in the same group are in different tracking regions..." << endl;
        classification_succeeded = false;
      }
    }
    else if (the_plate_groups[i].number_of_plates_in_group == 2)
    {
      dummy_string  = the_plates[the_plate_groups[i].plate_index_list[0]].plate_type;
      dummy_string += "_";
      dummy_string += the_plates[the_plate_groups[i].plate_index_list[1]].plate_type;
      the_plate_groups[i].group_type = dummy_string;
      
      if (the_plates[the_plate_groups[i].plate_index_list[0]].tracking_region == the_plates[the_plate_groups[i].plate_index_list[1]].tracking_region)
      {
        the_plate_groups[i].tracking_region = the_plates[the_plate_groups[i].plate_index_list[0]].tracking_region;
      }
      else
      {
        cout << "Plates in the same group are in different tracking regions..." << endl;
        classification_succeeded = false;
      }
    }
    else
    {
      cout << "Number of plates in group different than 2 or 3: group " << i << " with " << the_plate_groups[i].number_of_plates_in_group << " groups." << endl;
      
      classification_succeeded = false;
    }
  }
  
  for (int i = 0; i < number_of_groups; i++)
  {
    the_tracking_regions[the_plate_groups[i].tracking_region].plate_index_list[the_tracking_regions[the_plate_groups[i].tracking_region].number_of_groups] = i;
    the_tracking_regions[the_plate_groups[i].tracking_region].number_of_groups++;
  }
  
  for (int i = 0; i < number_of_groups; i++)
  {
    if (the_plate_groups[i].number_of_plates_in_group == 2)
    {
      // there can be more allowed combinations if arbitrary directions are allowed or a group can be
      // X_D1 or something like that
      if (!(the_plate_groups[i].group_type == "X_Y" || the_plate_groups[i].group_type == "Y_X" ||
            the_plate_groups[i].group_type == "D1_D2" || the_plate_groups[i].group_type == "D2_D1"))
      {
        cout << "A group of plates is found that are not compatible: group " << i << " is type " << the_plate_groups[i].group_type << endl;
        classification_succeeded = false;
      }
    }
    else if (the_plate_groups[i].number_of_plates_in_group == 3)
    {
      // ditto; so a smarter way can be found...
      if (!(the_plate_groups[i].group_type == "X_Y_D1" || the_plate_groups[i].group_type == "Y_X_D1" ||
            the_plate_groups[i].group_type == "X_D1_Y" || the_plate_groups[i].group_type == "Y_D1_X" ||
            the_plate_groups[i].group_type == "D1_X_Y" || the_plate_groups[i].group_type == "D1_Y_X" ||
            the_plate_groups[i].group_type == "X_Y_D2" || the_plate_groups[i].group_type == "Y_X_D2" ||
            the_plate_groups[i].group_type == "X_D2_Y" || the_plate_groups[i].group_type == "Y_D2_X" ||
            the_plate_groups[i].group_type == "D2_X_Y" || the_plate_groups[i].group_type == "D2_Y_X" ||
            the_plate_groups[i].group_type == "D1_D2_X" || the_plate_groups[i].group_type == "D2_D1_X" ||
            the_plate_groups[i].group_type == "D1_X_D2" || the_plate_groups[i].group_type == "D2_X_D1" ||
            the_plate_groups[i].group_type == "X_D1_D2" || the_plate_groups[i].group_type == "X_D2_D1" ||
            the_plate_groups[i].group_type == "D1_D2_Y" || the_plate_groups[i].group_type == "D2_D1_Y" ||
            the_plate_groups[i].group_type == "D1_Y_D2" || the_plate_groups[i].group_type == "D2_Y_D1" ||
            the_plate_groups[i].group_type == "Y_D1_D2" || the_plate_groups[i].group_type == "Y_D2_D1"))
      {
        cout << "A group of plates is found that are not compatible: group " << i << " is type " << the_plate_groups[i].group_type << endl;
        classification_succeeded = false;
      }
    }
  }
   
  for (int i = 0; i < number_of_tracking_regions; i++)
  {
    // can be more cases
    int count_xy = 0;
    int count_dd = 0;
    
    // this part is tricky and may need corrections
    // midstream and downstream can be counter together for this purpose,
    // so correct later if necessary
    //
    for (int j = 0; j < the_tracking_regions[i].number_of_groups; j++)
    {
      if (the_plate_groups[the_tracking_regions[i].group_index_list[j]].group_type == "X_Y" ||
          the_plate_groups[the_tracking_regions[i].group_index_list[j]].group_type == "Y_X" ||
          the_plate_groups[the_tracking_regions[i].group_index_list[j]].group_type == "X_Y_D1" ||
          the_plate_groups[the_tracking_regions[i].group_index_list[j]].group_type == "Y_X_D1" ||
          the_plate_groups[the_tracking_regions[i].group_index_list[j]].group_type == "X_Y_D2" ||
          the_plate_groups[the_tracking_regions[i].group_index_list[j]].group_type == "Y_X_D2" ||
          the_plate_groups[the_tracking_regions[i].group_index_list[j]].group_type == "X_D1_Y" ||
          the_plate_groups[the_tracking_regions[i].group_index_list[j]].group_type == "Y_D1_X" ||
          the_plate_groups[the_tracking_regions[i].group_index_list[j]].group_type == "X_D2_Y" ||
          the_plate_groups[the_tracking_regions[i].group_index_list[j]].group_type == "Y_D2_X" ||
          the_plate_groups[the_tracking_regions[i].group_index_list[j]].group_type == "D1_X_Y" ||
          the_plate_groups[the_tracking_regions[i].group_index_list[j]].group_type == "D1_Y_X" ||
          the_plate_groups[the_tracking_regions[i].group_index_list[j]].group_type == "D2_X_Y" ||
          the_plate_groups[the_tracking_regions[i].group_index_list[j]].group_type == "D2_Y_X")
      {
        count_xy++;
      }
      else if (the_plate_groups[the_tracking_regions[i].group_index_list[j]].group_type == "D1_D2" ||
               the_plate_groups[the_tracking_regions[i].group_index_list[j]].group_type == "D2_D1" ||
               the_plate_groups[the_tracking_regions[i].group_index_list[j]].group_type == "D1_D2_X" ||
               the_plate_groups[the_tracking_regions[i].group_index_list[j]].group_type == "D1_D2_Y" ||
               the_plate_groups[the_tracking_regions[i].group_index_list[j]].group_type == "D1_X_D2" ||
               the_plate_groups[the_tracking_regions[i].group_index_list[j]].group_type == "D1_Y_D2" ||
               the_plate_groups[the_tracking_regions[i].group_index_list[j]].group_type == "X_D1_D2" ||
               the_plate_groups[the_tracking_regions[i].group_index_list[j]].group_type == "Y_D1_D2" ||
               the_plate_groups[the_tracking_regions[i].group_index_list[j]].group_type == "D2_D1_X" ||
               the_plate_groups[the_tracking_regions[i].group_index_list[j]].group_type == "D2_D1_Y" ||
               the_plate_groups[the_tracking_regions[i].group_index_list[j]].group_type == "D2_X_D1" ||
               the_plate_groups[the_tracking_regions[i].group_index_list[j]].group_type == "D2_Y_D1" ||
               the_plate_groups[the_tracking_regions[i].group_index_list[j]].group_type == "X_D2_D1" ||
               the_plate_groups[the_tracking_regions[i].group_index_list[j]].group_type == "Y_D2_D1")
      {
        count_dd++;
      }
    }
    
    if (count_xy < 2 && count_dd < 2)
    {
      cout << "No 2 plate groups of the same type in a region: region " << i << endl;
      classification_succeeded = false;
    }
    
  }
  
  if (PRINT_DEBUG == 3)
  {
    cout << "Geometry directions:" << endl;
    cout << primary_direction_1 << endl;
    cout << primary_direction_1_name << endl;
    cout << primary_direction_2 << endl;
    cout << primary_direction_2_name << endl;
    cout << secondary_direction_1 << endl;
    cout << secondary_direction_1_name << endl;
    cout << secondary_direction_2 << endl;
    cout << secondary_direction_2_name << endl;
  
    cout << "-----------" << endl;
    cout << "-----------" << endl;

    cout << "Tracking regions: " << endl;
    
    for (int i = 0; i < number_of_tracking_regions; i++)
    {
      cout << tracking_region_names[i] << endl;
    }
  
    cout << "-----------" << endl;
    cout << "-----------" << endl;
  
    cout << "Groups: " << number_of_groups << endl;
  
    for (int i = 0; i < number_of_groups; i++)
    {
      cout << the_plate_groups[i].number_of_plates_in_group << " " << the_plate_groups[i].group_type << " " << the_plate_groups[i].tracking_region << endl;
    
      for (int j = 0; j < the_plate_groups[i].number_of_plates_in_group; j++)
      {
        cout << the_plate_groups[i].plate_index_list[j] << endl;
      }
    }
  
    cout << "-----------" << endl;
    cout << "-----------" << endl;
  
    cout << "Plates: " << actual_number_of_plates << endl;
  
    for (int i = 0; i < actual_number_of_plates; i++)
    {
      cout << "-----------" << endl;
      cout << "plate " << i << ":" << endl;
      cout << the_plates[i].plate_id << endl;
      cout << the_plates[i].tracking_region << endl;
      cout << the_plates[i].plate_group << endl;
      cout << the_plates[i].plate_type << endl;
      cout << the_plates[i].efficiency << endl;
      cout << the_plates[i].position << endl;
      cout << the_plates[i].rotation << endl;
      cout << the_plates[i].size << endl;
      cout << the_plates[i].number_of_strips << endl;
      cout << the_plates[i].strip_direction << endl;
      cout << the_plates[i].pitch << endl;
    }
  
    cout << "-----------" << endl;
    cout << "-----------" << endl;
  
    cout << "Target: " << endl;
    cout << the_target.position << endl;
    cout << the_target.rotation << endl;
    cout << the_target.size << endl;
    cout << the_target.material << endl;
  }
  
  return classification_succeeded;
}

bool find_line_plate_intersection_point(line_3d a_track_line, plane_3d plate_plane, XYZVector& intersection_point)
{
  double t;
  double dummy;
  
  dummy = plate_plane.a * a_track_line.direction.x() + plate_plane.b * a_track_line.direction.y() + plate_plane.c * a_track_line.direction.z();
  
  if (dummy == 0) return false;
  
  t = - (plate_plane.d + plate_plane.a * a_track_line.point.x() + plate_plane.b * a_track_line.point.y() + plate_plane.c * a_track_line.point.z()) / dummy;
  
  intersection_point = a_track_line.point + t * a_track_line.direction;
  
  return true;
}

double find_point_to_line_3d_distance(XYZVector a_point, line_3d a_line)
{
  return sqrt(a_line.direction.Cross(a_point - a_line.point).Mag2());
}

bool find_line_to_line_distance(line_3d line_1, line_3d line_2, double& distance)
{
  XYZVector d;
  double t1, t2;
  
  if (line_1.direction.Cross(line_2.direction).Mag2() == 0)
  {
    if (line_1.point == line_2.point)
    {
      distance = 0;
      return true;
    }
    else
    {
      return false;
    }
  }
  
  d = line_1.point - line_2.point;
  
  t1 = (line_1.direction.Cross(line_2.direction)).Dot(d.Cross(line_2.direction)) / sqrt((line_1.direction.Cross(line_2.direction)).Mag2());
  
  t2 = (line_1.direction.Cross(line_2.direction)).Dot(d.Cross(line_1.direction)) / sqrt((line_1.direction.Cross(line_2.direction)).Mag2());
  
  distance = sqrt((line_1.point + t1 * line_1.direction - line_2.point - t2 * line_2.direction).Mag2());
  
  return true;
}

bool find_closest_point(line_3d line_1, line_3d line_2, XYZVector& closest_point)
{
  XYZVector d;
  double t1, t2;
  
  d = line_1.point - line_2.point;
  
  t1 = (line_1.direction.Cross(line_2.direction)).Dot(d.Cross(line_2.direction)) / sqrt((line_1.direction.Cross(line_2.direction)).Mag2());
  
  t2 = (line_1.direction.Cross(line_2.direction)).Dot(d.Cross(line_1.direction)) / sqrt((line_1.direction.Cross(line_2.direction)).Mag2());
  
  closest_point = 0.5 * (line_1.point + t1 * line_1.direction + line_2.point + t2 * line_2.direction);
  
  return true;
}

void track::sort_cluster_lines_by_plate_index()
{
  pair<line_3d, int> a_pair;
  vector<pair<line_3d, int>> the_pairs;
  
  for (int i = 0; i < cluster_lines.size(); i++)
  {
    a_pair.first = cluster_lines.at(i);
    a_pair.second = plate_index_of_cluster_line.at(i);
    the_pairs.push_back(a_pair);
  }
  
  sort(the_pairs.begin(), the_pairs.end(), compare_cluster_line_plate);
  
  cluster_lines.clear();
  plate_index_of_cluster_line.clear();
  
  for (int i = 0; i < the_pairs.size(); i++)
  {
    cluster_lines.push_back(the_pairs.at(i).first);
    plate_index_of_cluster_line.push_back(the_pairs.at(i).second);
  }
}

