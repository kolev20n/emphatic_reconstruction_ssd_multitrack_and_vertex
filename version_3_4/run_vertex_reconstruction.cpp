#include <iostream>

#include "input_output_manager.h"
#include "reconstruction_manager.h"
#include "global_data_dispatcher.h"
#include "partial_geometry.h"

using namespace std;

global_data_dispatcher* the_global_data_dispatcher;

int main(int argc, char **argv)
{
  the_global_data_dispatcher = new global_data_dispatcher;

  input_output_manager* the_input_output_manager = new input_output_manager();

  if (!the_input_output_manager->initialize_files(argc, argv))
  {
    cout << "The initialize_files operation of the_input_output_manager failed. Quitting..." << endl;
    return 101;
  }
  
  partial_geometry* the_geometry = new partial_geometry();
  
  the_global_data_dispatcher->register_geometry(the_geometry);
  
  if (!the_input_output_manager->initialize_geometry(the_geometry))
  {
    cout << "The initialize_geometry operation of the_input_output_manager failed. Quitting..." << endl;
    return 102;
  }
  
  reconstruction_manager the_reconstruction_manager(the_input_output_manager);
  
  the_reconstruction_manager.process_events();
  
  return 0;
}
