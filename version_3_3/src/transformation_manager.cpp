#include <iostream>

#include <Math/GenVector/RotationX.h>
#include <Math/GenVector/RotationY.h>
#include <Math/GenVector/RotationZ.h>

#include "transformation_manager.h"
#include "input_output_manager.h"
#include "common_data_structures_and_functions.h"

using namespace std;
using namespace ROOT::Math;

transformation_manager::transformation_manager(partial_geometry* a_geometry)
{
  the_geometry = a_geometry;
}

bool transformation_manager::transform_to_lab_frame(std::vector<double> clusters[], vector<silicon_strip> cluster_lines[])
{
  XYZVector temp_vector;
  double dummy_double;
  
  silicon_strip temp_strip;
  
  RotationX rx;
  RotationY ry;
  RotationZ rz;
  
  Rotation3D total_rotation;
  
  const double to_radians = Pi() / 180.;
  
  for (int i = 0; i < the_geometry->get_number_of_plates(); i++)
  {
    rx.SetAngle(the_geometry->get_plate(i).rotation.x() * to_radians);
    ry.SetAngle(the_geometry->get_plate(i).rotation.y() * to_radians);
    rz.SetAngle(the_geometry->get_plate(i).rotation.z() * to_radians);
    
    total_rotation = rz * ry * rx;
    
    for (int j = 0; j < clusters[i].size(); j++)
    {
      dummy_double = (clusters[i].at(j) - (((double)(the_geometry->get_plate(i).number_of_strips)) - 1.) / 2.) * 60 * um;
      
      temp_vector.SetCoordinates(dummy_double, 0., 0.);
      temp_vector = total_rotation * temp_vector;
      temp_vector = temp_vector + the_geometry->get_plate(i).position;
      temp_strip.line_center = temp_vector;
      
      // this is already done in partial_geometry; use it
      temp_vector.SetCoordinates(0., 1., 0.);
      temp_vector = total_rotation * temp_vector;
      temp_strip.line_direction = temp_vector;
      cluster_lines[i].push_back(temp_strip);
    }
  }
  
  return true;
}

