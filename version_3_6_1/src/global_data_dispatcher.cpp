#include "global_data_dispatcher.h"
#include "common_data_structures_and_functions.h"
#include "partial_geometry.h"

using namespace std;
using namespace ROOT::Math;

partial_geometry* global_data_dispatcher::get_partial_geometry()
{
  return the_geometry;
}

void global_data_dispatcher::register_geometry(partial_geometry* a_geometry)
{
  the_geometry = a_geometry;
}

input_root_data_structure* global_data_dispatcher::get_input_root_data_structure()
{
  return the_input_root_data_structure;
}

void global_data_dispatcher::register_input_root_data_structure(input_root_data_structure* an_input_root_data_structure)
{
  the_input_root_data_structure = an_input_root_data_structure;
}

run_output_data_structure* global_data_dispatcher::get_run_output_data_structure()
{
  return the_run_output_data_structure;
}

void global_data_dispatcher::register_run_output_data_structure(run_output_data_structure* a_run_output_data_structure)
{
  the_run_output_data_structure = a_run_output_data_structure;
}

event_mc_output_data_structure* global_data_dispatcher::get_event_mc_output_data_structure()
{
  return the_event_mc_output_data_structure;
}

void global_data_dispatcher::register_event_mc_output_data_structure(event_mc_output_data_structure* a_event_mc_output_data_structure)
{
  the_event_mc_output_data_structure = a_event_mc_output_data_structure;
}

event_reco_output_data_structure* global_data_dispatcher::get_event_reco_output_data_structure()
{
  return the_event_reco_output_data_structure;
}

void global_data_dispatcher::register_event_reco_output_data_structure(event_reco_output_data_structure* a_event_reco_output_data_structure)
{
  the_event_reco_output_data_structure = a_event_reco_output_data_structure;
}

reconstruction_options* global_data_dispatcher::get_reconstruction_options()
{
  return the_reconstruction_options;
}

void global_data_dispatcher::register_reconstruction_options(reconstruction_options* a_reconstruction_options)
{
  the_reconstruction_options = a_reconstruction_options;
}

event_characteristics* global_data_dispatcher::get_event_characteristics()
{
  return the_event_characteristics;
}

void global_data_dispatcher::register_event_characteristics(event_characteristics* an_event_characteristics)
{
  the_event_characteristics = an_event_characteristics;
}

void global_data_dispatcher::record_particle_data()
{
  int dummy_id;
  
  charged_track dummy_track;
  map<int, int> track_id_to_vector_element_map;
  
  double x_sum, y_sum, z_sum, de_sum;
  int counter_steps_in_strip;
  XYZVector dummy_vector;
  bool first_counted = true;
  int previous_step_index;
  
  mc_step dummy_step;
  
  the_event_mc_output_data_structure->the_charged_tracks.clear();
  the_event_mc_output_data_structure->the_selected_steps.clear();
  track_id_to_vector_element_map.clear();
  
  for (int i = 0; i < the_input_root_data_structure->track_id->size(); i++)
  {
    the_run_output_data_structure->particle_codes[the_input_root_data_structure->particle_code->at(i)]++;
    
    if (the_input_root_data_structure->ekin_vertex->at(i) > the_reconstruction_options->cut_min_energy_to_be_high_energy_charged_particle)
    {
      dummy_id = the_input_root_data_structure->particle_code->at(i);
      
      if (abs(dummy_id) == 11 || abs(dummy_id) == 13 || abs(dummy_id) == 211 || abs(dummy_id) == 321 || abs(dummy_id) == 2212 || abs(dummy_id) == 3112 || abs(dummy_id) == 3222 || abs(dummy_id) == 3312 || abs(dummy_id) == 3334)
      {
        dummy_track.track_id = the_input_root_data_structure->track_id->at(i);
        dummy_track.particle_code = dummy_id;
        the_event_mc_output_data_structure->the_charged_tracks.push_back(dummy_track);
        track_id_to_vector_element_map[the_input_root_data_structure->track_id->at(i)] = the_event_mc_output_data_structure->the_charged_tracks.size() - 1;
        
/*
cout << the_input_root_data_structure->track_id->at(i) << " " << dummy_id << " " << the_input_root_data_structure->ekin_vertex->at(i) << endl;
*/
        
/*
if (the_input_root_data_structure->visited_plates->at(i).size() > 0)
{
cout << "Visited plates: ";
for (int j = 0; j < the_input_root_data_structure->visited_plates->at(i).size(); j++)
{
  cout << the_input_root_data_structure->visited_plates->at(i).at(j) << " ";
}
cout << endl;
}
*/
      }
    }
  }

/*
cout << "Track id to vector element map size:" << track_id_to_vector_element_map.size() << endl;
for (map<int, int>::iterator it = track_id_to_vector_element_map.begin(); it != track_id_to_vector_element_map.end(); it++)
{
  cout << it->first << " " << it->second << endl;
}
*/
  
  counter_steps_in_strip = 0;
  x_sum = 0.;
  y_sum = 0.;
  z_sum = 0.;
  de_sum = 0.;
  
  for (int i = 0; i < the_input_root_data_structure->step_track_id->size(); i++)
  {
    if (the_input_root_data_structure->step_volume_name->at(i) == "physical_ssd_strip" && track_id_to_vector_element_map.find(the_input_root_data_structure->step_track_id->at(i)) != track_id_to_vector_element_map.end())
    {
      if (first_counted)
      {
        x_sum  = (the_input_root_data_structure->step_x_f->at(i) + the_input_root_data_structure->step_x_i->at(i)) / 2.;
        y_sum  = (the_input_root_data_structure->step_y_f->at(i) + the_input_root_data_structure->step_y_i->at(i)) / 2.;
        z_sum  = (the_input_root_data_structure->step_z_f->at(i) + the_input_root_data_structure->step_z_i->at(i)) / 2.;
        de_sum =  the_input_root_data_structure->step_final_energy->at(i) - the_input_root_data_structure->step_initial_energy->at(i);
        counter_steps_in_strip++;
        first_counted = false;
        previous_step_index = i;
      }
      else if (the_input_root_data_structure->step_track_id->at(i) == the_input_root_data_structure->step_track_id->at(i - 1) && the_input_root_data_structure->step_plate_id->at(i) == the_input_root_data_structure->step_plate_id->at(i - 1))
      {
        x_sum  += (the_input_root_data_structure->step_x_f->at(i) + the_input_root_data_structure->step_x_i->at(i)) / 2.;
        y_sum  += (the_input_root_data_structure->step_y_f->at(i) + the_input_root_data_structure->step_y_i->at(i)) / 2.;
        z_sum  += (the_input_root_data_structure->step_z_f->at(i) + the_input_root_data_structure->step_z_i->at(i)) / 2.;
        de_sum +=  the_input_root_data_structure->step_initial_energy->at(i) - the_input_root_data_structure->step_final_energy->at(i);
        counter_steps_in_strip++;
        previous_step_index = i;
      }
      else
      {
        dummy_id = track_id_to_vector_element_map[the_input_root_data_structure->step_track_id->at(previous_step_index)];
        
        the_event_mc_output_data_structure->the_charged_tracks.at(dummy_id).step_plate_index.push_back(the_input_root_data_structure->step_plate_id->at(previous_step_index));
        dummy_vector.SetCoordinates(x_sum / ((double) counter_steps_in_strip), y_sum / ((double) counter_steps_in_strip), z_sum / ((double) counter_steps_in_strip));
        the_event_mc_output_data_structure->the_charged_tracks.at(dummy_id).position.push_back(dummy_vector);

        the_event_mc_output_data_structure->the_charged_tracks.at(dummy_id).energy_lost.push_back(de_sum);
        
        counter_steps_in_strip = 0;
        
        x_sum  = (the_input_root_data_structure->step_x_f->at(i) + the_input_root_data_structure->step_x_i->at(i)) / 2.;
        y_sum  = (the_input_root_data_structure->step_y_f->at(i) + the_input_root_data_structure->step_y_i->at(i)) / 2.;
        z_sum  = (the_input_root_data_structure->step_z_f->at(i) + the_input_root_data_structure->step_z_i->at(i)) / 2.;
        de_sum =  the_input_root_data_structure->step_initial_energy->at(i) - the_input_root_data_structure->step_final_energy->at(i);
        counter_steps_in_strip++;
        previous_step_index = i;
      }
    
      dummy_step.step_index = i;
      dummy_step.step_track_id = the_input_root_data_structure->step_track_id->at(i);
      dummy_step.step_volume_name = the_input_root_data_structure->step_volume_name->at(i);
      dummy_step.step_plate_id = the_input_root_data_structure->step_plate_id->at(i);
      dummy_step.initial_position.SetCoordinates(the_input_root_data_structure->step_x_i->at(i), the_input_root_data_structure->step_y_i->at(i), the_input_root_data_structure->step_z_i->at(i));
      dummy_step.final_position.SetCoordinates(the_input_root_data_structure->step_x_f->at(i), the_input_root_data_structure->step_y_f->at(i), the_input_root_data_structure->step_z_f->at(i));
      dummy_step.initial_energy = the_input_root_data_structure->step_initial_energy->at(i);
      dummy_step.final_energy = the_input_root_data_structure->step_final_energy->at(i);
      dummy_step.step_process_name = the_input_root_data_structure->step_process_name->at(i);
      the_event_mc_output_data_structure->the_selected_steps.push_back(dummy_step);
    }
  }
  
  dummy_id = track_id_to_vector_element_map[the_input_root_data_structure->step_track_id->at(previous_step_index)];
  
  /*
  cout << dummy_id << " " << the_input_root_data_structure->step_track_id->at(previous_step_index) << endl;
  cout << counter_steps_in_strip << " " << the_input_root_data_structure->step_track_id->at(previous_step_index) << " " << the_input_root_data_structure->step_plate_id->at(previous_step_index) << " " << x_sum << endl;
  */
  
  the_event_mc_output_data_structure->the_charged_tracks.at(dummy_id).step_plate_index.push_back(the_input_root_data_structure->step_plate_id->at(previous_step_index));
  dummy_vector.SetCoordinates(x_sum / ((double) counter_steps_in_strip), y_sum / ((double) counter_steps_in_strip), z_sum / ((double) counter_steps_in_strip));
  the_event_mc_output_data_structure->the_charged_tracks.at(dummy_id).position.push_back(dummy_vector);
  
  the_event_mc_output_data_structure->the_charged_tracks.at(dummy_id).energy_lost.push_back(de_sum);
}

