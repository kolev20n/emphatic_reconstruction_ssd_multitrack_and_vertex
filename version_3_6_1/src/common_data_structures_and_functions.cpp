#include <map>

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

void event_mc_output_data_structure::print()
{
  cout << "Printing Event MC Data Structure:" << endl;
  cout << "==========================================" << endl;
/*
  cout << "Slected Steps:" << endl;
  cout << "------------------------------------------" << endl;
  
  for (int i = 0; i < the_selected_steps.size(); i++)
  {
    cout << the_selected_steps.at(i).step_index << " "
         << the_selected_steps.at(i).step_track_id << " "
         << the_selected_steps.at(i).step_volume_name << " "
         << the_selected_steps.at(i).step_plate_id  << " "
         << the_selected_steps.at(i).initial_position << " "
         << the_selected_steps.at(i).final_position << " "
         << the_selected_steps.at(i).initial_energy << " "
         << the_selected_steps.at(i).final_energy << " "
         << the_selected_steps.at(i).step_process_name << endl;
  }
*/
  for (int i = 0; i < the_charged_tracks.size(); i++)
  {
    if (the_charged_tracks.at(i).step_plate_index.size() > 0)
    {
      cout << the_charged_tracks.at(i).track_id << " " << the_charged_tracks.at(i).particle_code << " ";
      for (int j = 0; j < the_charged_tracks.at(i).step_plate_index.size(); j++)
      {
        cout << the_charged_tracks.at(i).step_plate_index.at(j) << " "
        << the_charged_tracks.at(i).position.at(j) << " "
        << the_charged_tracks.at(i).energy_lost.at(j) << " ";
      }
      
      cout << endl;
    }
  }
}

void event_reco_output_data_structure::print()
{
  cout << "Printing Event Reco Output Data Structure:" << endl;
  cout << "==========================================" << endl;
  
  for (int j = 0; j < the_global_data_dispatcher->get_partial_geometry()->get_number_of_plates(); j++)
  {
    cout << "Plate: " << j << endl;
    
    for (int k = 0; k < cluster_lines[j].size(); k++)
    {
      cout << cluster_lines[j].at(k).point << " " << cluster_lines[j].at(k).direction << endl;
    }
  }
}

bool compare_hit_pairs(const hit_pair& first, const hit_pair& second)
{
  return (first.strip_number < second.strip_number);
}

bool compare_cluster_line_plate(const pair<line_3d, int> &a, const pair<line_3d, int> &b)
{
  return (a.second < b.second);
}

bool compare_track_options(const track_option& first, const track_option& second)
{
  return (first.distance < second.distance);
}

int temp_point_to_strip_diagonal_plate(ROOT::Math::XYZVector a_point)
{
  double x, y;
  double result;
  
  x = 0.707106781 * a_point.X() + 0.707106781 * a_point.Y();
  y = -0.707106781 * a_point.X() + 0.707106781 * a_point.Y();
    
  return ((int) (x / 0.006 + 320.));
}
