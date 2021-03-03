#include "global_data_dispatcher.h"
#include "common_data_structures_and_functions.h"

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

reconstruction_options* global_data_dispatcher::get_reconstruction_options()
{
  return the_reconstruction_options;
}

void global_data_dispatcher::register_reconstruction_options(reconstruction_options* a_reconstruction_options)
{
  the_reconstruction_options = a_reconstruction_options;
}
