#ifndef common_data_structures_and_functions_h
#define common_data_structures_and_functions_h

#include <vector>
#include <string>
#include <algorithm>
#include <utility>
#include <map>

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
const int max_number_of_triplet_groups = 4;
const int max_number_of_tracks = 100;

const double pi = 3.141592653589793;
const double to_radians = pi / 180.;

class global_data_dispatcher;

extern global_data_dispatcher* the_global_data_dispatcher;

class line_3d
{
public:
  ROOT::Math::XYZVector point;
  ROOT::Math::XYZVector direction;
  ROOT::Math::XYZVector error_point;
  ROOT::Math::XYZVector error_direction;
};

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
  
  double cut_max_distance_to_accept_in_track;
  double guessed_vertex_weight;
  double cut_max_distance_to_count_to_vertex;
  double cut_min_energy_to_be_high_energy_charged_particle;
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

class event_characteristics
{
public:
  int event_class_by_cluster_multiplicity;
  // 10: single track events
  // 11: pure single track but possible bend at the target
  // 12: single track: at most one upstream plate with more than 1 cluster and at most two midstream + downstream plates with more than 1 cluster or with 0 clusters
  //
  // 20: normal multitrack events
  // 21 normal vertex (all upstream are 1, almost all mistream + downstream are > 1 (except allowed by cut))
  // 22 normal vertex (at most one upstream plate with more than 1 cluster, almost all mistream + downstream are > 1 (except allowed by cut))
  //
  // 70: abnormal multitrack events (can't be done)
  // 71: multiple tracks upstream (more than 1 plate with multiple clusters)
  //
  // 99 unknown events
  
  int event_class_by_reconstructed_tracks;
  // 10: single track events
  // 11: single track no bend
  // 12: single track bend at target
  // 13: single track bend not at target
  
  ROOT::Math::XYZVector the_vertex;
  
  int plate_cluster_multiplicity[max_number_of_plates];
  int group_cluster_multiplicity[max_number_of_groups];
};

class charged_track
{
public:
  int track_id;
  int particle_code;
  int ekin_vertex;
  std::vector<int> step_plate_index;
  std::vector<ROOT::Math::XYZVector> position;
  std::vector<double> energy_lost;
};

class mc_step
{
public:
  int step_index;
  int step_track_id;
  std::string step_volume_name;
  int step_plate_id;
  ROOT::Math::XYZVector initial_position;
  ROOT::Math::XYZVector final_position;
  double initial_energy;
  double final_energy;
  std::string step_process_name;
};

class run_output_data_structure
{
public:
  // geometry
  int number_of_plates;
  
  // algorithm
  int algorithm;
  
  std::map<int, int> particle_codes;
};

class event_mc_output_data_structure
{
public:
  std::vector<charged_track> the_charged_tracks; // high energy charged particles
  std::vector<mc_step> the_selected_steps;
  
  void print();
};

class event_reco_output_data_structure
{
public:
  std::vector<line_3d> cluster_lines[max_number_of_plates];
  
  int event_class_by_cluster_multiplicity;
  // 10: single track events
  // 11: pure single track but possible bend at the target
  // 12: single track: at most one upstream plate with more than 1 cluster and at most two midstream + downstream plates with more than 1 cluster or with 0 clusters
  //
  // 20: normal multitrack events
  // 21 normal vertex (all upstream are 1, almost all mistream + downstream are > 1 (except allowed by cut))
  // 22 normal vertex (at most one upstream plate with more than 1 cluster, almost all mistream + downstream are > 1 (except allowed by cut))
  //
  // 70: abnormal multitrack events (can't be done)
  // 71: multiple tracks upstream (more than 1 plate with multiple clusters)
  //
  // 99 unknown events
  
  void print();
};

class track_option
{
public:
  int track_index;
  int cluster_index_on_triplet_plates[3];
  double distance;
  int strip_on_third_plate;
};

bool compare_hit_pairs(const hit_pair& first, const hit_pair& second);

bool compare_cluster_line_plate(const std::pair<line_3d, int> &a, const std::pair<line_3d, int> &b);

bool compare_track_options(const track_option& first, const track_option& second);

int temp_point_to_strip_diagonal_plate(ROOT::Math::XYZVector a_point);


#endif
