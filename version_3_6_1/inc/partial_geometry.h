#ifndef partial_geometry_h
#define partial_geometry_h

#include <vector>
#include <string>

#include <TMath.h>
#include <Math/Vector3D.h>

#include "common_data_structures_and_functions.h"

class line_2d
{
public:
  ROOT::Math::XYZVector point;
  ROOT::Math::XYZVector direction;
  ROOT::Math::XYZVector error_point;
  ROOT::Math::XYZVector error_direction;
};

class plane_3d
{
public:
  // the coefficients in the plane equation ax + by + cz + d = 0
  double a, b, c, d;
};

class ssd_plate
{
public:
  int plate_id; // for reference purposes only
  int tracking_region;
  int plate_group;
  std::string plate_type; // will be compared to the actual orientation
  int number_of_strips;
  double efficiency;
  ROOT::Math::XYZVector position; // position in cm
  ROOT::Math::XYZVector rotation; // rotations in degrees
  ROOT::Math::XYZVector size; // size in cm: along strips, perpendicular to strips, thickness
  ROOT::Math::XYZVector strip_direction; // strip direction unit vector
  ROOT::Math::XYZVector normal_direction; // normal direction unit vector
  // ROOT::Math::XYZVector intended_strip_direction; // remove later; this will be the plate_type
  double pitch;
  plane_3d plate_plane_equation;
};

class plate_group
{
public:
  int number_of_plates_in_group;
  int plate_index_list[max_number_of_plates_in_group]; // maybe have the actual plates?
  
  std::string group_type;
  // xy, yx, xyd1, yxd1, xyd2, yxd2, d1xy, d1yx, d2xy, d2yx, xd1y, yd1x, xd2y, yd2x, etc.
  
  int tracking_region; // 0 for upstream, 1 for midstream, 2 for downstream
};

class tracking_region
{
public:
  std::string name;
  int number_of_groups;
  int group_index_list[max_number_of_groups_in_tracking_region];
  int number_of_plates;
  int plate_index_list[max_number_of_plates_in_tracking_region];
};

class target
{
public:
  ROOT::Math::XYZVector position; // position in cm
  ROOT::Math::XYZVector rotation; // rotations in degrees
  ROOT::Math::XYZVector size; // size in cm
  std::string material; // none, carbon, iron, aluminum
};

class track
{
public:
  std::vector<line_3d> cluster_lines;
  std::vector<int> plate_index_of_cluster_line;
  
  line_3d track_line;
  
  // the following are the position and direction at z = 0;
  // the track may not start or end there;
  // this is irrelevant for the calculation though
  // ROOT::Math::XYZVector position_at_0;
  // ROOT::Math::XYZVector direction_at_0;
  
  double fit_goodness;
  
  // the track shouldn't be extrapolated outside
  // the following z coordinates
  double assumed_start_group;
  double assumed_end_group;
  
  double momentum_magnitude;
  
  void sort_cluster_lines_by_plate_index();
};

class partial_geometry
{
public:
  partial_geometry();

  int get_number_of_plates();
  int get_number_of_plates_in_tracking_region(int a_region);
  int get_number_of_tracking_regions();
  int get_number_of_groups();
  int get_number_of_groups_in_tracking_region(int a_region);
  int get_number_of_plates_in_group(int a_group);
  int get_number_of_triple_groups();
  
  target* get_target();
  
  ROOT::Math::XYZVector get_primary_direction_1();
  ROOT::Math::XYZVector get_primary_direction_2();
  ROOT::Math::XYZVector get_secondary_direction_1();
  ROOT::Math::XYZVector get_secondary_direction_2();
  
  std::string get_primary_direction_1_name();
  std::string get_primary_direction_2_name();
  std::string get_secondary_direction_1_name();
  std::string get_secondary_direction_2_name();

  tracking_region* get_tracking_region(int a_region);
  plate_group* get_plate_group(int a_group);
  ssd_plate* get_plate(int a_plate_index);
  ssd_plate* get_plate_by_id(int a_plate_id); // for reference purposes
  int get_triple_plate_group_index(int a_group_index);
    
  bool add_plate(ssd_plate a_plate);
  bool add_target(target a_target);
  
  bool set_primary_direction_1(ROOT::Math::XYZVector a_direction);
  bool set_primary_direction_2(ROOT::Math::XYZVector a_direction);
  bool set_secondary_direction_1(ROOT::Math::XYZVector a_direction);
  bool set_secondary_direction_2(ROOT::Math::XYZVector a_direction);

  bool set_primary_direction_1_name(std::string a_name);
  bool set_primary_direction_2_name(std::string a_name);
  bool set_secondary_direction_1_name(std::string a_name);
  bool set_secondary_direction_2_name(std::string a_name);
  
  bool set_tracking_region_names(int a_region, std::string a_name);
  
  bool classify_plates();
  
private:
  int actual_number_of_plates;
  ssd_plate the_plates[max_number_of_plates];
  
  int number_of_tracking_regions;
  tracking_region the_tracking_regions[max_number_of_tracking_regions];
  std::string tracking_region_names[max_number_of_tracking_regions];

  int number_of_groups;
  plate_group the_plate_groups[max_number_of_groups];
  
  int number_of_triple_groups;
  int the_triple_groups[max_number_of_groups];

  ROOT::Math::XYZVector primary_direction_1;
  std::string primary_direction_1_name;
  
  ROOT::Math::XYZVector primary_direction_2;
  std::string primary_direction_2_name;
  
  ROOT::Math::XYZVector secondary_direction_1;
  std::string secondary_direction_1_name;
  
  ROOT::Math::XYZVector secondary_direction_2;
  std::string secondary_direction_2_name;

  target the_target;
  

  
};

bool find_line_plate_intersection_point(line_3d a_track_line, plane_3d plate_plane, ROOT::Math::XYZVector& intersection_point);

double find_point_to_line_3d_distance(ROOT::Math::XYZVector a_point, line_3d a_line);

bool find_line_to_line_distance(line_3d line_1, line_3d line_2, double& distance);

bool find_closest_point(line_3d line_1, line_3d line_2, ROOT::Math::XYZVector& closest_point);


#endif
