#ifndef global_data_dispatcher_h
#define global_data_dispatcher_h

class partial_geometry;
class input_root_data_structure;
class reconstruction_options;

class global_data_dispatcher
{
public:
  partial_geometry* get_partial_geometry();
  void register_geometry(partial_geometry* a_geometry);
  
  input_root_data_structure* get_input_root_data_structure();
  void register_input_root_data_structure(input_root_data_structure* an_input_root_data_structure);

  reconstruction_options* get_reconstruction_options();
  void register_reconstruction_options(reconstruction_options* a_reconstruction_options);
  
private:
  partial_geometry* the_geometry;
  input_root_data_structure* the_input_root_data_structure;
  reconstruction_options* the_reconstruction_options;
};

#endif
