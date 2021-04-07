#ifndef cluster_manager_h
#define cluster_manager_h

#include <vector>
#include <string>

#include "common_data_structures_and_functions.h"

class hit_pair;

class input_output_manager;
class TTree;
class TFile;
class transformation_manager;
class cluster_interpretation_manager;

class cluster_manager
{
public:
  cluster_manager();
  bool cluster(std::vector<hit_pair> hit_pair_list[], std::vector<double> clusters[], int actual_number_of_plates, std::string cluster_interpretation_method);
private:

};

#endif
