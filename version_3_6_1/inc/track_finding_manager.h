#ifndef track_finding_manager_h
#define track_finding_manager_h

#include <vector>

#include "common_data_structures_and_functions.h"
#include "partial_geometry.h"

class track_finding_manager
{
public:
  track_finding_manager(partial_geometry* a_geometry, int an_algorithm);
  bool find_tracks(std::vector<line_3d> cluster_lines[], std::vector<track>& the_tracks);
  bool classify_event(std::vector<line_3d> cluster_lines[]);
  bool assign_cluster_lines_single_track_upstream(std::vector<line_3d> cluster_lines[]);
  bool assign_cluster_lines_single_track_downstream(std::vector<line_3d> cluster_lines[]);
  bool assign_cluster_lines_multiple_track_downstream(std::vector<line_3d> cluster_lines[], std::vector<track>& downstream_tracks);
  bool find_vertex(std::vector<track>& the_tracks, ROOT::Math::XYZVector the_vertex);
  bool fit_track(track& a_track);
  bool find_guessed_vertex(line_3d upstream_track_line);
  
private:
  partial_geometry* the_geometry;
  
  event_characteristics the_event_characteristics;
  
  int the_algorithm;
  // 1: no field; 1 triplet group; 2 regions
  
  line_3d guessed_vertex;
};

#endif
