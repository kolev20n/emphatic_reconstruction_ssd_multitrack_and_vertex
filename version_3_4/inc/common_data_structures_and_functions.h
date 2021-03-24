#ifndef common_data_structures_and_functions_h
#define common_data_structures_and_functions_h

#include <vector>
#include <string>

#include <TMath.h>
#include <Math/Vector3D.h>

const double cm = 1.;     // units used throughout are centimeters
const double mm = 1.e-1;  // millimeter
const double um = 1.e-4; // micrometer

// the actual number of plates will be read from the geometry file
const int max_number_of_plates = 32;
const int max_number_of_plates_in_tracking_region = 32;
const int max_number_of_plates_in_group = 4;
const int max_number_of_groups_in_tracking_region = 16;
const int max_number_of_tracking_regions = 4;
const int max_number_of_groups = 16;

const double pi = 3.141592653589793;
const double to_radians = pi / 180.;

class global_data_dispatcher;

extern global_data_dispatcher* the_global_data_dispatcher;

class reconstruction_options
{
public:
  bool use_magnetic_field;
  std::string cluster_interpretation_method; // default: weighted
  double maximum_cosangle_between_intended_and_actual_strip_direction;
  int max_upstream_plates_with_multiple_clusters;
  int max_clusters_on_upstream_plates_with_multiple_clusters;
  int max_midstream_plus_downstream_plates_with_multiple_clusters;
  int max_midstream_plus_downstream_plates_with_zero_clusters;
  int max_midstream_plus_downstream_plates_with_zero_or_multiple_clusters;
  int max_midstream_plus_downstream_plates_with_less_than_2_clusters;
};

class input_root_data_structure
{
public:
  input_root_data_structure();
  ~input_root_data_structure();
  
  std::vector<int> *plate_number;
  std::vector<int> *strip_number;
  std::vector<double> *total_energy_values;
  std::vector<double> *non_ionization_energy_values;
  std::vector<std::vector<double>> *contributing_tracks;
  
  std::vector<int> *track_id;
  std::vector<int> *particle_code;
  std::vector<int> *parent_track_id;
  std::vector<std::string> *creator_process;
  std::vector<double> *x_vertex;
  std::vector<double> *y_vertex;
  std::vector<double> *z_vertex;
  std::vector<double> *px_vertex;
  std::vector<double> *py_vertex;
  std::vector<double> *pz_vertex;
  std::vector<double> *ekin_vertex;
  std::vector<std::vector<double>> *visited_plates;
  
  std::vector<int> *step_track_id;
  std::vector<std::string> *step_volume_name;
  std::vector<int> *step_plate_id;
  std::vector<double> *step_x_i;
  std::vector<double> *step_y_i;
  std::vector<double> *step_z_i;
  std::vector<double> *step_px_i;
  std::vector<double> *step_py_i;
  std::vector<double> *step_pz_i;
  std::vector<double> *step_x_f;
  std::vector<double> *step_y_f;
  std::vector<double> *step_z_f;
  std::vector<double> *step_px_f;
  std::vector<double> *step_py_f;
  std::vector<double> *step_pz_f;
  std::vector<double> *step_initial_energy;
  std::vector<double> *step_final_energy;
  std::vector<std::string> *step_process_name;
};

class hit_pair
{
public:
  int strip_number;
  double energy_value; // adc_value if experimental hit
};

class track
{
public:
  // these are the silicon strip clusters associated with the tracks;
  // in may be useful to have this information if iterations will be done
  // and for other purposes (e.g. matching studies of mc and reco results)
  std::vector<ROOT::Math::XYZVector> cluster_line_centers;
  std::vector<ROOT::Math::XYZVector> cluster_line_directions;
  
  // the following are the position and direction at z = 0;
  // the track may not start or end there;
  // this is irrelevant for the calculation though
  ROOT::Math::XYZVector position_at_0;
  ROOT::Math::XYZVector direction_at_0;
  
  // the track shouldn't be extrapolated outside
  // the following z coordinates
  double assumed_start_group;
  double assumed_end_group;
  
  double momentum_magnitude;
};

class event_characteristics
{
public:
  int event_class_by_cluster_multiplicity;
  // 10: single track events
  // 11: pure single track but possible bend at the target
  // 12: single track: at most one upstream plate with more than 1 cluster and at most two midstream + downstream plates with more than 1 cluster or with 0 clusters
  // 13: like 12 but one of the exceptions is on a triple group
  //
  // 20: normal multitrack events
  // 21 normal vertex (almost all upstream are 1, almost all mistream + downstream are > 1 (except allowed by cut))
  //
  // 70: abnormal multitrack events (can't be done)
  // 71: multiple tracks upstream (more than 1 plate with multiple clusters)
  //
  // 99 unknown events
  
  int plate_cluster_multiplicity[max_number_of_plates];
  int group_cluster_multiplicity[max_number_of_groups];
};

bool compare_hit_pairs(const hit_pair& first, const hit_pair& second);

#endif
