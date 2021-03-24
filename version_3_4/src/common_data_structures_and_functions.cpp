#include <map>

#include <TMath.h>
#include <Math/Vector3D.h>
#include <Math/GenVector/Rotation3D.h>
#include <Math/GenVector/RotationX.h>
#include <Math/GenVector/RotationY.h>
#include <Math/GenVector/RotationZ.h>

#include "common_data_structures_and_functions.h"

using namespace std;
using namespace ROOT::Math;

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

bool compare_hit_pairs(const hit_pair& first, const hit_pair& second)
{
  return (first.strip_number < second.strip_number);
}
