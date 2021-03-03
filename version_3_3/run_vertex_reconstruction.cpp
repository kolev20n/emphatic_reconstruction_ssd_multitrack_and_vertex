#include <iostream>

#include "input_output_manager.h"
#include "reconstruction_manager.h"
#include "common_data_structures_and_functions.h"
#include "global_data_dispatcher.h"

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
  
  reconstruction_manager the_reconstruction_manager(the_input_output_manager);
  
  the_reconstruction_manager.process_events();
  
  return 0;
}
