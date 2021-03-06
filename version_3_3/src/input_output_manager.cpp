#include <iostream>
#include <string>
#include <sys/stat.h>
#include <sstream>
#include <fstream>

#include <TFile.h>
#include <TTree.h>

#include "input_output_manager.h"
#include "histogram_manager.h"
#include "jsmn.h"
#include "global_data_dispatcher.h"

using namespace std;

input_output_manager::input_output_manager()
:output_message(""), the_input_root_file(0), the_input_ssd_energy_tree(0),
  the_input_particle_truth_tree(0), the_input_steps_truth_tree(0),
  number_of_events_to_process(0), the_reconstruction_options(0)
{
  cout << "Created the_input_output_manager." << endl;
  
  the_reconstruction_options = new reconstruction_options();
  the_global_data_dispatcher->register_reconstruction_options(the_reconstruction_options);
}

input_output_manager::~input_output_manager()
{
}

bool input_output_manager::initialize_files(int argc, char** argv)
{
  int number_of_tokens;
  string current_token;
  
  bool initialization_succeeded = true;
  
  string command_input_file_name = "";
  magnetic_field_file_name = "";
  geometry_file_name = "";
  input_root_file_name = "";
  
  string dummy_string;
  
  struct stat buffer;
  
  ifstream command_input_file_stream;
  
  string line;
  string card;
  string string_card_value;
  int int_card_value;
  
  cout << "File initialization started in the_input_output_manager." << endl;
  
  if (argc == 1)
  {
    output_message += "Too few command line arguments. Usage:\n";
    output_message += argv[0];
    output_message += " command_input_file.txt\n";
    initialization_succeeded = false;
  }
  else if (argc > 2)
  {
    output_message += "Too many command line arguments. Usage:\n";
    output_message += argv[0];
    output_message += " command_input_file.txt\n";
    
    initialization_succeeded = false;
  }
  else
  {
    command_input_file_name = argv[1];
    
    if (stat(command_input_file_name.c_str(), &buffer) != 0)
    {
      output_message += "The command_input_file does not exist: ";
      output_message += command_input_file_name;
      output_message += "\n";
      
      initialization_succeeded = false;
    }
    else
    {
      output_message += "Using command input file name: ";
      output_message += command_input_file_name;
      output_message += "\n";
      
      command_input_file_stream.open(command_input_file_name.c_str(), ifstream::in);
      
      while (getline(command_input_file_stream, line))
      {
        bool empty_line = false;
        
        if (line.empty())
        {
          empty_line = true;
        }
        else
        {
          bool line_of_spaces = true;
          
          for (int i = 0; i < line.length(); i++)
          {
            if (line[i] != ' ')
            {
              line_of_spaces = false;
            }
          }
          
          if (!line_of_spaces) empty_line = false;
        }
        
        if (!empty_line)
        {
          if (line[0] != '#')
          {
            istringstream temp_string_stream(line);
            temp_string_stream >> card;
            if (card == "use_magnetic_field")
            {
              temp_string_stream >> dummy_string;
              
              if (dummy_string == "true")
              {
                the_reconstruction_options->use_magnetic_field = true;
              }
              else if (dummy_string == "false")
              {
                the_reconstruction_options->use_magnetic_field = false;
              }
              else
              {
                output_message += "Uknown option for use_magnetic_field: ";
                output_message += dummy_string;
                output_message += "\n";
                output_message += "Will use default value: false";
                the_reconstruction_options->use_magnetic_field = false;
              }
            }
            else if (card == "magnetic_field_file_name")
            {
              temp_string_stream >> magnetic_field_file_name;
            }
            else if (card == "geometry_file_name")
            {
              temp_string_stream >> geometry_file_name;
            }
            else if (card == "input_root_file_name")
            {
              temp_string_stream >> input_root_file_name;
            }
            else if (card == "input_root_ssd_energy_tree_name")
            {
              temp_string_stream >> input_root_ssd_energy_tree_name;
            }
            else if (card == "input_root_particle_truth_tree_name")
            {
              temp_string_stream >> input_root_particle_truth_tree_name;
            }
            else if (card == "input_root_steps_truth_tree_name")
            {
              temp_string_stream >> input_root_steps_truth_tree_name;
            }
            else if (card == "output_root_file_name")
            {
              temp_string_stream >> output_root_file_name;
            }
            else if (card == "number_of_events_to_process")
            {
              temp_string_stream >> dummy_string;
              
              if (dummy_string == "all")
              {
                number_of_events_to_process = 0;
              }
              else
              {
                number_of_events_to_process = stoi(dummy_string);
              }
            }
            else if (card == "cluster_interpretation_method")
            {
              temp_string_stream >> dummy_string;
              
              if (dummy_string == "weighted")
              {
                the_reconstruction_options->cluster_interpretation_method = dummy_string;
              }
              else
              {
                output_message += "Unknown cluster interpretation method: ";
                output_message += dummy_string;
                output_message += "\n Will use the default weighted method...\n";
                the_reconstruction_options->cluster_interpretation_method = "weighted";
              }
            }
            else
            {
              output_message += "Unknown input item: ";
              output_message += card;
              output_message += "\n";
              output_message += "Ignored...\n";
            }
          }
        }
      }
    }
    
    output_message += "Use magnetic field: ";
    
    if (the_reconstruction_options->use_magnetic_field)
    {
      output_message += "yes.\n";
    }
    else
    {
      output_message += "no.\n";
    }
    
    if (the_reconstruction_options->use_magnetic_field)
    {
      dummy_string = "config/";
      dummy_string += magnetic_field_file_name;
      magnetic_field_file_name = dummy_string;
      output_message += "Magnetic field file: ";
      output_message += magnetic_field_file_name;
      output_message += "\n";
    }
    
    dummy_string = "config/";
    dummy_string += geometry_file_name;
    geometry_file_name = dummy_string;
    output_message += "Geometry file: ";
    output_message += geometry_file_name;
    output_message += "\n";
    
    dummy_string = "../data/";
    dummy_string += input_root_file_name;
    input_root_file_name = dummy_string;
    output_message += "Input root file: ";
    output_message += input_root_file_name;
    output_message += "\n";
    
    std::ostringstream oss;
    oss << number_of_events_to_process;
 
    output_message += "Requested to process ";
    if (number_of_events_to_process == 0)
    {
      output_message += "all";
    }
    else
    {
      output_message += oss.str();
    }
    output_message += " events.\n";
  }
  
  if (initialization_succeeded)
  {
    if (the_reconstruction_options->use_magnetic_field)
    {
      if (stat(magnetic_field_file_name.c_str(), &buffer) != 0)
      {
        output_message += "Magnetic field file does not exist: ";
        output_message += magnetic_field_file_name;
        output_message += "\n";
        initialization_succeeded = false;
      }
    }
    
    if (stat(geometry_file_name.c_str(), &buffer) != 0)
    {
      output_message += "Geometry file does not exist: ";
      output_message += geometry_file_name;
      output_message += "\n";
      initialization_succeeded = false;
    }
    
    if (stat(input_root_file_name.c_str(), &buffer) != 0)
    {
      output_message += "Input root file does not exist: ";
      output_message += input_root_file_name;
      output_message += "\n";
      initialization_succeeded = false;
    }
  }
  
  output_message += "Output histograms will be written to file: ";
  output_message += output_root_file_name;
  output_message +="\n";
  
  cout << output_message;
  
  return initialization_succeeded;
}

bool input_output_manager::initialize_geometry(partial_geometry* the_geometry)
{
  bool initialization_succeeded = true;

  int number_of_tokens;
  int p1_token;
  int p1_size;
  string current_token;
  
  ssd_plate temp_plate;
  target temp_target;
  double temp_x, temp_y, temp_z;
  
  output_message = "";
  
  ifstream geometry_file(geometry_file_name.c_str(), ifstream::in);
  
  string geometry_file_string((std::istreambuf_iterator<char>(geometry_file)), (std::istreambuf_iterator<char>()));
  
  geometry_file.close();
  
  const char* json_string_1 = geometry_file_string.c_str();
  
  jsmn_parser p;
  jsmntok_t t[1024];
  
  jsmn_init(&p);
  number_of_tokens = jsmn_parse(&p, json_string_1, strlen(json_string_1), t, sizeof(t) / sizeof(t[0]));
  
  if (number_of_tokens < 0)
  {
    output_message += "Failed to parse geometry json file: ";
    output_message += geometry_file_name;
    output_message += "\n";
    
    initialization_succeeded = false;
  }
  else
  {
    output_message += "Number of tokens in geometry json file: ";
    std::ostringstream oss;
    oss << number_of_tokens;
    output_message += oss.str();
    output_message += "\n";
  }
  
  if (number_of_tokens < 1 || t[0].type != JSMN_OBJECT)
  {
    output_message += "Object expected as the root of the geometry json file.\n";
    
    initialization_succeeded = false;
  }
  
  int i = 0;
  
  while (i < number_of_tokens - 1)
  {
    i++;
    
    current_token = geometry_file_string.substr(t[i].start, t[i].end - t[i].start);
    
    // cout << i << " " << t[i].type << " " << t[i].size << " " << current_token << endl;
    
    if (current_token == "plate")
    {
      p1_token = i + 1;
      p1_size = t[p1_token].size;
      p1_token++;
    
      for (int j = 0; j < p1_size; j++)
      {
        current_token = geometry_file_string.substr(t[p1_token].start, t[p1_token].end - t[p1_token].start);
     
        if (current_token == "id")
        {
          current_token = geometry_file_string.substr(t[p1_token + 1].start, t[p1_token + 1].end - t[p1_token + 1].start);
          temp_plate.plate_id = stoi(current_token);
          p1_token += 2;
        }
        else if (current_token == "group")
        {
          current_token = geometry_file_string.substr(t[p1_token + 1].start, t[p1_token + 1].end - t[p1_token + 1].start);
          temp_plate.plate_group = stoi(current_token);
          p1_token += 2;
        }
        else if (current_token == "efficiency")
        {
          current_token = geometry_file_string.substr(t[p1_token + 1].start, t[p1_token + 1].end - t[p1_token + 1].start);
          temp_plate.efficiency = stod(current_token);
          p1_token += 2;
        }
        else if (current_token == "intended_strip_direction")
        {
          current_token = geometry_file_string.substr(t[p1_token + 2].start, t[p1_token + 2].end - t[p1_token + 2].start);
          temp_x = stod(current_token);
          current_token = geometry_file_string.substr(t[p1_token + 3].start, t[p1_token + 3].end - t[p1_token + 3].start);
          temp_y = stod(current_token);
          current_token = geometry_file_string.substr(t[p1_token + 4].start, t[p1_token + 4].end - t[p1_token + 4].start);
          temp_z = stod(current_token);
          temp_plate.intended_strip_direction.SetXYZ(temp_x, temp_y, temp_z);
          p1_token += 5;
        }
        else if (current_token == "position")
        {
          current_token = geometry_file_string.substr(t[p1_token + 2].start, t[p1_token + 2].end - t[p1_token + 2].start);
          temp_x = stod(current_token);
          current_token = geometry_file_string.substr(t[p1_token + 3].start, t[p1_token + 3].end - t[p1_token + 3].start);
          temp_y = stod(current_token);
          current_token = geometry_file_string.substr(t[p1_token + 4].start, t[p1_token + 4].end - t[p1_token + 4].start);
          temp_z = stod(current_token);
          temp_plate.position.SetXYZ(temp_x, temp_y, temp_z);
          p1_token += 5;
        }
        else if (current_token == "rotation")
        {
          current_token = geometry_file_string.substr(t[p1_token + 2].start, t[p1_token + 2].end - t[p1_token + 2].start);
          temp_x = stod(current_token);
          current_token = geometry_file_string.substr(t[p1_token + 3].start, t[p1_token + 3].end - t[p1_token + 3].start);
          temp_y = stod(current_token);
          current_token = geometry_file_string.substr(t[p1_token + 4].start, t[p1_token + 4].end - t[p1_token + 4].start);
          temp_z = stod(current_token);
          temp_plate.rotation.SetXYZ(temp_x, temp_y, temp_z);
          p1_token += 5;
        }
        else if (current_token == "size")
        {
          current_token = geometry_file_string.substr(t[p1_token + 2].start, t[p1_token + 2].end - t[p1_token + 2].start);
          temp_x = stod(current_token);
          current_token = geometry_file_string.substr(t[p1_token + 3].start, t[p1_token + 3].end - t[p1_token + 3].start);
          temp_y = stod(current_token);
          temp_plate.size.SetXYZ(temp_x, temp_y, 0.);
          p1_token += 4;
        }
        else if (current_token == "number_of_strips")
        {
          current_token = geometry_file_string.substr(t[p1_token + 1].start, t[p1_token + 1].end - t[p1_token + 1].start);
          temp_plate.number_of_strips = stoi(current_token);
          p1_token += 2;
        }
      }
      
      the_geometry->add_plate(temp_plate.plate_id, temp_plate.plate_group, temp_plate.intended_strip_direction, temp_plate.number_of_strips, temp_plate.efficiency, temp_plate.position, temp_plate.rotation, temp_plate.size);
      i = p1_token - 1;
    }
    else if (current_token == "target")
    {
      p1_token = i + 1;
      p1_size = t[p1_token].size;
      p1_token++;
      
      for (int j = 0; j < p1_size; j++)
      {
        current_token = geometry_file_string.substr(t[p1_token].start, t[p1_token].end - t[p1_token].start);
        
        if (current_token == "position")
        {
          current_token = geometry_file_string.substr(t[p1_token + 2].start, t[p1_token + 2].end - t[p1_token + 2].start);
          temp_x = stod(current_token);
          current_token = geometry_file_string.substr(t[p1_token + 3].start, t[p1_token + 3].end - t[p1_token + 3].start);
          temp_y = stod(current_token);
          current_token = geometry_file_string.substr(t[p1_token + 4].start, t[p1_token + 4].end - t[p1_token + 4].start);
          temp_z = stod(current_token);
          temp_target.position.SetXYZ(temp_x, temp_y, temp_z);
          p1_token += 5;
        }
        if (current_token == "rotation")
        {
          current_token = geometry_file_string.substr(t[p1_token + 2].start, t[p1_token + 2].end - t[p1_token + 2].start);
          temp_x = stod(current_token);
          current_token = geometry_file_string.substr(t[p1_token + 3].start, t[p1_token + 3].end - t[p1_token + 3].start);
          temp_y = stod(current_token);
          current_token = geometry_file_string.substr(t[p1_token + 4].start, t[p1_token + 4].end - t[p1_token + 4].start);
          temp_z = stod(current_token);
          temp_target.rotation.SetXYZ(temp_x, temp_y, temp_z);
          p1_token += 5;
        }
        else if (current_token == "size")
        {
          current_token = geometry_file_string.substr(t[p1_token + 2].start, t[p1_token + 2].end - t[p1_token + 2].start);
          temp_x = stod(current_token);
          current_token = geometry_file_string.substr(t[p1_token + 3].start, t[p1_token + 3].end - t[p1_token + 3].start);
          temp_y = stod(current_token);
          current_token = geometry_file_string.substr(t[p1_token + 4].start, t[p1_token + 4].end - t[p1_token + 4].start);
          temp_z = stod(current_token);
          temp_target.rotation.SetXYZ(temp_x, temp_y, temp_z);
          p1_token += 5;
        }
        else if (current_token == "material")
        {
          temp_target.material = geometry_file_string.substr(t[p1_token + 1].start, t[p1_token + 1].end - t[p1_token + 1].start);
          p1_token += 2;
        }

        the_geometry->add_target(temp_target.position, temp_target.rotation, temp_target.size, temp_target.material);
        i = p1_token - 1;
      }
    }
  }
  
  // create lists of x, y and d plates, and calculate number of upstream, midstream and downstream plates
  the_geometry->classify_plates();
  
  cout << output_message;
  
  return initialization_succeeded;
}

bool input_output_manager::initialize_root_trees(TFile*& an_input_root_file, TTree*& an_input_ssd_energy_tree, TTree*& an_input_particle_truth_tree, TTree*& an_input_steps_truth_tree, input_root_data_structure*& the_input_root_data_structure)
{
  bool initialization_succeeded = true;
  
  output_message = "";
  
  the_input_root_file = new TFile(input_root_file_name.c_str(), "READ");
  the_input_ssd_energy_tree = (TTree*) the_input_root_file->Get(input_root_ssd_energy_tree_name.c_str());
  the_input_particle_truth_tree = (TTree*) the_input_root_file->Get(input_root_particle_truth_tree_name.c_str());
  the_input_steps_truth_tree = (TTree*) the_input_root_file->Get(input_root_steps_truth_tree_name.c_str());

  if (!the_input_ssd_energy_tree)
  {
    output_message += "Root tree not found: ";
    output_message += input_root_ssd_energy_tree_name;
    output_message += "\n";
    initialization_succeeded = false;
  }
  
  if (!the_input_particle_truth_tree)
  {
    output_message += "Root tree not found: ";
    output_message += input_root_particle_truth_tree_name;
    output_message += "\n";
    initialization_succeeded = false;
  }

  if (!the_input_steps_truth_tree)
  {
    output_message += "Root tree not found: ";
    output_message += input_root_steps_truth_tree_name;
    output_message += "\n";
    initialization_succeeded = false;
  }
  
  if (!(the_input_ssd_energy_tree->GetEntries() == the_input_particle_truth_tree->GetEntries() && the_input_ssd_energy_tree->GetEntries() == the_input_steps_truth_tree->GetEntries()))
  {
    output_message += "The root trees have different number of events: ";
    output_message += the_input_ssd_energy_tree->GetEntries();
    output_message += " ";
    output_message += the_input_particle_truth_tree->GetEntries();
    output_message += " ";
    output_message += the_input_steps_truth_tree->GetEntries();
    output_message += "\n Quitting...\n";
    initialization_succeeded = false;
  }
  
  if (initialization_succeeded)
  {
    if (number_of_events_to_process == 0)
    {
      number_of_events_to_process = the_input_ssd_energy_tree->GetEntries();
      output_message += "Will process: ";
      std::ostringstream oss;
      oss << number_of_events_to_process;
      output_message += oss.str();
      output_message += " events.\n";
    }
    else
    {
      if (number_of_events_to_process > the_input_ssd_energy_tree->GetEntries())
      {
        output_message += "Requested number of events to process higher than existing events in the root file. Will process ";
        number_of_events_to_process = the_input_ssd_energy_tree->GetEntries();
        std::ostringstream oss;
        oss << number_of_events_to_process;
        output_message += oss.str();
        output_message += " events.\n";
      }
    }
  }
  
  an_input_root_file = the_input_root_file;
  an_input_ssd_energy_tree = the_input_ssd_energy_tree;
  an_input_particle_truth_tree = the_input_particle_truth_tree;
  an_input_steps_truth_tree = the_input_steps_truth_tree;
  
  the_input_ssd_energy_tree->SetBranchAddress("plate_id", &the_input_root_data_structure->plate_number);
  the_input_ssd_energy_tree->SetBranchAddress("strip_id", &the_input_root_data_structure->strip_number);
  the_input_ssd_energy_tree->SetBranchAddress("ssd_total_energy_deposited", &the_input_root_data_structure->total_energy_values);
  the_input_ssd_energy_tree->SetBranchAddress("ssd_non_ionization_energy_deposited", &the_input_root_data_structure->non_ionization_energy_values);
  the_input_ssd_energy_tree->SetBranchAddress("contributing_tracks", &the_input_root_data_structure->contributing_tracks);
  
  the_input_particle_truth_tree->SetBranchAddress("track_id", &the_input_root_data_structure->track_id);
  the_input_particle_truth_tree->SetBranchAddress("particle_code", &the_input_root_data_structure->particle_code);
  the_input_particle_truth_tree->SetBranchAddress("parent_track_id", &the_input_root_data_structure->parent_track_id);
  the_input_particle_truth_tree->SetBranchAddress("creator_process", &the_input_root_data_structure->creator_process);
  the_input_particle_truth_tree->SetBranchAddress("x_vertex", &the_input_root_data_structure->x_vertex);
  the_input_particle_truth_tree->SetBranchAddress("y_vertex", &the_input_root_data_structure->y_vertex);
  the_input_particle_truth_tree->SetBranchAddress("z_vertex", &the_input_root_data_structure->z_vertex);
  the_input_particle_truth_tree->SetBranchAddress("px_vertex", &the_input_root_data_structure->px_vertex);
  the_input_particle_truth_tree->SetBranchAddress("py_vertex", &the_input_root_data_structure->py_vertex);
  the_input_particle_truth_tree->SetBranchAddress("pz_vertex", &the_input_root_data_structure->pz_vertex);
  the_input_particle_truth_tree->SetBranchAddress("ekin_vertex", &the_input_root_data_structure->ekin_vertex);
  the_input_particle_truth_tree->SetBranchAddress("visited_plates", &the_input_root_data_structure->visited_plates);

  the_input_steps_truth_tree->SetBranchAddress("step_track_id", &the_input_root_data_structure->step_track_id);
  the_input_steps_truth_tree->SetBranchAddress("step_volume_name", &the_input_root_data_structure->step_volume_name);
  the_input_steps_truth_tree->SetBranchAddress("step_plate_id", &the_input_root_data_structure->step_plate_id);
  the_input_steps_truth_tree->SetBranchAddress("step_x_i", &the_input_root_data_structure->step_x_i);
  the_input_steps_truth_tree->SetBranchAddress("step_y_i", &the_input_root_data_structure->step_y_i);
  the_input_steps_truth_tree->SetBranchAddress("step_z_i", &the_input_root_data_structure->step_z_i);
  the_input_steps_truth_tree->SetBranchAddress("step_px_i", &the_input_root_data_structure->step_px_i);
  the_input_steps_truth_tree->SetBranchAddress("step_py_i", &the_input_root_data_structure->step_py_i);
  the_input_steps_truth_tree->SetBranchAddress("step_pz_i", &the_input_root_data_structure->step_pz_i);
  the_input_steps_truth_tree->SetBranchAddress("step_x_f", &the_input_root_data_structure->step_x_f);
  the_input_steps_truth_tree->SetBranchAddress("step_y_f", &the_input_root_data_structure->step_y_f);
  the_input_steps_truth_tree->SetBranchAddress("step_z_f", &the_input_root_data_structure->step_z_f);
  the_input_steps_truth_tree->SetBranchAddress("step_px_f", &the_input_root_data_structure->step_px_f);
  the_input_steps_truth_tree->SetBranchAddress("step_py_f", &the_input_root_data_structure->step_py_f);
  the_input_steps_truth_tree->SetBranchAddress("step_pz_f", &the_input_root_data_structure->step_pz_f);
  the_input_steps_truth_tree->SetBranchAddress("step_initial_energy", &the_input_root_data_structure->step_initial_energy);
  the_input_steps_truth_tree->SetBranchAddress("step_final_energy", &the_input_root_data_structure->step_final_energy);
  the_input_steps_truth_tree->SetBranchAddress("step_process_name", &the_input_root_data_structure->step_process_name);

  cout << output_message;
  
  return initialization_succeeded;
}

int input_output_manager::get_number_of_events_to_process()
{
  return number_of_events_to_process;
}

/*
reconstruction_options input_output_manager::get_reconstruction_options()
{
  return the_reconstruction_options;
}

 
bool input_output_manager::set_reconstruction_option_algorithm(int an_algorithm)
{
  the_reconstruction_options.algorithm_type = an_algorithm;
}
*/

bool input_output_manager::initialize_output(histogram_manager*& a_histogram_manager, output_data_structure* the_output_data_structure)
{
  the_histogram_manager = new histogram_manager(output_root_file_name, the_output_data_structure);
  a_histogram_manager = the_histogram_manager;
  
  the_histogram_manager->book();
}
