#ifndef track_finding_manager_h
#define track_finding_manager_h

#include <TMath.h>
#include <Math/Vector3D.h>

#include "common_data_structures_and_functions.h"

class histogram_manager;

class track_finding_manager
{
public:
  track_finding_manager(histogram_manager* a_histogram_manager, output_data_structure* an_output_data_structure, partial_geometry* a_geometry);
  bool find_tracks(std::vector<silicon_strip> cluster_lines[], std::vector<track>& the_tracks);
  bool find_preliminary_event_type(std::vector<silicon_strip> cluster_lines[]);
  bool assign_clusters_to_tracks(std::vector<silicon_strip> cluster_lines[], std::vector<track>& the_tracks);
  bool find_guessed_vertex(std::vector<silicon_strip> cluster_lines[], silicon_strip& a_vertex);
  
private:
  histogram_manager* the_histogram_manager;
  output_data_structure* the_output_data_structure;
  partial_geometry* the_geometry;
  std::string preliminary_event_type;
  std::vector<int> group_cluster_multiplicity;
  int potential_vertex_before_group;
};

#endif
