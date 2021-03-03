#ifndef transformation_manager_h
#define transformation_manager_h

#include <vector>

#include <TMath.h>
#include <Math/Vector3D.h>
#include <Math/GenVector/Rotation3D.h>
#include <Math/GenVector/RotationX.h>
#include <Math/GenVector/RotationY.h>
#include <Math/GenVector/RotationZ.h>

#include "common_data_structures_and_functions.h"

class input_output_manager;

class transformation_manager
{
public:
  transformation_manager(partial_geometry* a_geometry);
  bool transform_to_lab_frame(std::vector<double> clusters[], std::vector<silicon_strip> cluster_lines[]);
  
private:
  partial_geometry* the_geometry;
};

#endif
