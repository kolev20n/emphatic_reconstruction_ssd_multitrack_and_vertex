#ifndef track_finding_manager_h
#define track_finding_manager_h

#include <vector>

#include "common_data_structures_and_functions.h"

class partial_geometry;
class line_3d;

class track_finding_manager
{
public:
  track_finding_manager(partial_geometry* a_geometry, int an_algorithm);
  bool find_tracks(std::vector<line_3d> cluster_lines[], std::vector<track>& the_tracks);
  bool classify_event(std::vector<line_3d> cluster_lines[]);
  bool fit_single_track(std::vector<line_3d> some_cluster_lines[], int first_cluster_line, int last_cluster_line, track& a_track);
  bool assign_cluster_lines_upstream(std::vector<line_3d> cluster_lines[], std::vector<line_3d> assigned_cluster_lines[]);
  
private:
  partial_geometry* the_geometry;
  
  event_characteristics the_event_characteristics;
  
  int the_algorithm;
  // 1: no field; 1 triplet group; 3 regions (but 2 for some purposes - mid and downstream can be combined for most purposes)
};

#endif
