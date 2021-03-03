#ifndef common_data_structures_and_functions_h
#define common_data_structures_and_functions_h

#include <vector>
#include <string>

#include <TMath.h>
#include <Math/Vector3D.h>
#include <Math/Vector2D.h>

class global_data_dispatcher;
class silicon_strip;

const double cm = 1.;     // units used throughout are centimeters
const double mm = 1.e-1;  // millimeter
const double um = 1.e-4; // micrometer

// the actual number of plates will be read from the geometry file
const int max_number_of_plates = 32;

extern global_data_dispatcher* the_global_data_dispatcher;

class line_2d
{
public:
  ROOT::Math::XYVector point;
  ROOT::Math::XYVector direction;
  ROOT::Math::XYVector error_point;
  ROOT::Math::XYVector error_direction;
};

class line_3d
{
public:
  ROOT::Math::XYZVector point;
  ROOT::Math::XYZVector direction;
  ROOT::Math::XYZVector error_point;
  ROOT::Math::XYZVector error_direction;
};

class plane_3d
{
public:
  // the coefficients in the plane equation ax + by + cz + d = 0
  double a, b, c, d;
};

class geometry_fitters
{
public:
  
  static bool find_track_3d(std::vector<line_3d> lines, line_3d& result_line, int configuration_key);
  // 1 two xy groups (as in upstream of real_4 geometry)
  // 2
  // and so on as needed
  
  static bool find_track_2d(std::vector<ROOT::Math::XYVector> points, line_2d& result_line);
  
  static bool find_line_plane_intersection(line_3d line, plane_3d plane, ROOT::Math::XYZVector& intersection_point);

private:
  geometry_fitters();
};

class sorted_track
{
public:
  int track_id;
  int particle_code;
  int parent_track_id;
  double ekin_vertex;
  std::string creator_process;
};

// need to redo line centers and line directions to be this class
// everywhere, so there is one thing to pass around instead of two
class silicon_strip
{
public:
  ROOT::Math::XYZVector line_center;
  ROOT::Math::XYZVector line_direction;
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

// this is the data holder throughout the code;
// all output data should be written here;
// it will be translated in the end to the
// output_root_data_structure, if it is decided
// that it should be recorded
class output_data_structure
{
public:
  // geometry info
  int actual_number_of_plates;
  
  // true info
  //
  // high energy particle > 500 MeV with parent track id = 1
  // have the same vertex (checked on geometry real_3 with
  // 2 to 4 mid- and upstream tracks)
  bool same_vertex_for_high_energy_particles;
  int event_classification;
  ROOT::Math::XYZVector mc_guessed_vertex;
  
  // cluster info
  std::vector<silicon_strip> cluster_lines[max_number_of_plates];
  
  // reco info
  //
  silicon_strip reco_guessed_vertex;
  bool single_cluster_on_all_upstream_plates;
  std::vector<int> clusters_on_midstream_plates;
};

class output_root_data_structure
{
public:
  output_root_data_structure();
  ~output_root_data_structure();
  
  int event_classification;
  //double mc_guessed_vertex[3];
  //double reco_guessed_vertex[3];
  double x_vertex, y_vertex, z_vertex;
  double x_primary_track, y_primary_track;
  std::vector<int> clusters_on_midstream_plates;
  
  std::vector<int> plate_of_cluster;
  std::vector<double> x_of_cluster_center;
  std::vector<double> y_of_cluster_center;
  std::vector<double> z_of_cluster_center;
  std::vector<double> x_of_cluster_direction;
  std::vector<double> y_of_cluster_direction;
  std::vector<double> z_of_cluster_direction;
};

class reconstruction_options
{
public:
  bool use_magnetic_field;
  
  int algorithm_type;
  // 0 no algorithm found
  // 1 single triplet group (normally with diagonal plate)
  std::vector<int> triplet_groups;
  
  bool exclude_multitrack_clusters;
  // more than one track gives the same cluster on a plate; so cannot use a
  // cluster for multiple tracks? Default: yes (do not consider multiple tracks)
  
  int min_plates_per_track_upstream; // default 2
  
  int min_plates_per_track_midstream; // default 2
  
  int min_plates_per_track_downstream; // default 2
  
  int min_plates_per_track_midstream_plus_downstream; // default 4
  
  bool global_track_finding; // default false (so local - each track separately)
  
  std::string cluster_interpretation_method; // default: weighted
};

class ssd_plate
{
public:
  int plate_id;
  int plate_group;
  ROOT::Math::XYZVector intended_strip_direction; // will be checked
  int number_of_strips;
  double efficiency;
  ROOT::Math::XYZVector position; // position in cm
  ROOT::Math::XYZVector rotation; // rotations in degrees
  ROOT::Math::XYZVector size; // size in cm
  ROOT::Math::XYZVector strip_direction; // strip direction unit vector
  ROOT::Math::XYZVector normal_direction; // normal direction unit vector
  std::string plate_type;
};

class target
{
public:
  ROOT::Math::XYZVector position; // position in cm
  ROOT::Math::XYZVector rotation; // rotations in degrees
  ROOT::Math::XYZVector size; // size in cm
  std::string material; // none, carbon, iron, aluminum
};

class intended_strip_direction
{
public:
  std::string intended_direction_type;
  // x, y, d1, d2, a1, a2, ....
  std::vector<int> plate_index_list;
};

class plate_group
{
public:
  int group_id; // don't need to start at 0 in step of 1
  std::vector<int> plate_index_list;
  int group_type;
  // for the moment intended directions should be
  // (1, 0, 0) called x
  // (0, 1, 0) called y
  // (1, 1, 0) called d1
  // (1, -1, 0) called d2
  // (x, y, 0) arbitrary, but in the xy plane called a1
  // (x, -y, 0) ditto, called a2
  //
  // then groups are:
  // 1 x-y group
  // 2 d1-d2 group
  // 3, 4 x-y-d1 or x-y-d2 group resp.
  // 5, 6 d1-d2-x or d1-d2-y group resp.
  // can have many more, but explore later
  //
  std::string group_stream; // upstream, midstream, downstream
};

// for the moment this has only ssd plates and the target;
// it will also have (a reference to?) the magnetic field when it's needed
class partial_geometry
{
public:
  partial_geometry();
  int get_number_of_plates();
  int get_number_of_upstream_plates();
  int get_number_of_midstream_plates();
  int get_number_of_downstream_plates();
  int get_number_of_groups();
  int get_number_of_upstream_groups();
  int get_number_of_midstream_groups();
  int get_number_of_downstream_groups();
  std::vector<int> get_x_plates();
  std::vector<int> get_y_plates();
  std::vector<int> get_d_plates();
  ROOT::Math::XYZVector get_target_position();
  void add_plate(int a_plate_id, int a_plate_group, ROOT::Math::XYZVector an_intended_strip_direction, int a_number_of_strips, double an_efficiency, ROOT::Math::XYZVector a_position, ROOT::Math::XYZVector a_rotation, ROOT::Math::XYZVector a_size);
  void add_target(ROOT::Math::XYZVector a_position, ROOT::Math::XYZVector a_rotation, ROOT::Math::XYZVector size, std::string a_material);
  ssd_plate get_plate(int a_plate_number); // this is the position in the array, not the plate_id,
                                           //which can be different
  void classify_plates();
  
  std::vector<intended_strip_direction> the_intended_strip_directions;
  std::vector<plate_group> the_plate_groups;
  
private:
  int number_of_plates;
  int number_of_upstream_plates;
  int number_of_midstream_plates;
  int number_of_downstream_plates;
  int number_of_upstream_groups;
  int number_of_midstream_groups;
  int number_of_downstream_groups;
  std::vector<int> x_plates; // points along y
  std::vector<int> y_plates; // points along x
  std::vector<int> d_plates; // points along (1, -1, 0)
  ssd_plate the_plates[max_number_of_plates];
  target the_target;
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

class event_class
{
public:
  // normally should have 1 diagonal plate in a triplet group
  int number_of_triplet_groups;
  
  // 1 - single track (with possible disjoint multiple-cluster plates)
  // 2 - normal vertex: all 1 upstream (maybe more on a single plate), at least 2 on more than 1 group
  // 21: all groups the same number except disjoint plates
  // 22: normally decreasing
  // 23: increasing (secondary vertex)
  // 24: messy (ignore for now)
  int general_class;
  
};

bool compare_hit_pairs(const hit_pair& first, const hit_pair& second);

bool compare_map(const std::pair<std::string, int> &a, const std::pair<std::string, int> &b);

bool compare_mc_tracks(const sorted_track &a, const sorted_track &b);

// bool compare_z_lines_3d(const line_3d &a, const line_3d &b);

#endif
