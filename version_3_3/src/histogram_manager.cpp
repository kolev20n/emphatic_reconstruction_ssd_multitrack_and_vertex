#include <iostream>
#include <cmath>

#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TTree.h>

#include "histogram_manager.h"

using namespace std;

histogram_manager::histogram_manager(string an_output_file_name, output_data_structure* an_output_data_structure)
:out_file(0), output_tree(0), the_output_root_data_structure(0),
  h101(0)
{
  output_message = "Created the_histogram_manager.\n";
  
  output_root_file_name = an_output_file_name;
  
  the_output_data_structure = an_output_data_structure;
  
  output_message += "Using output root file: ";
  output_message += output_root_file_name;
  output_message += "\n";
  
  cout << output_message;
}

histogram_manager::~histogram_manager()
{
}

bool histogram_manager::book()
{
  out_file = new TFile(output_root_file_name.c_str(), "RECREATE");
  
  h101 = new TH1D("reco_mc_vertex_distance", "reco_mc_vertex_distance", 100, -0.2, 0.8);
  
  the_output_root_data_structure = new output_root_data_structure();

  output_tree = new TTree("ssd_tracks","ssd_tracks");

  output_tree->Branch("event_class", &the_output_root_data_structure->event_classification, "event_class/I");
  output_tree->Branch("plate_of_cluster", &the_output_root_data_structure->plate_of_cluster);
  output_tree->Branch("x_of_cluster_center", &the_output_root_data_structure->x_of_cluster_center);
  output_tree->Branch("y_of_cluster_center", &the_output_root_data_structure->y_of_cluster_center);
  output_tree->Branch("z_of_cluster_center", &the_output_root_data_structure->z_of_cluster_center);
  output_tree->Branch("x_of_cluster_direction", &the_output_root_data_structure->x_of_cluster_direction);
  output_tree->Branch("y_of_cluster_direction", &the_output_root_data_structure->y_of_cluster_direction);
  output_tree->Branch("z_of_cluster_direction", &the_output_root_data_structure->z_of_cluster_direction);
  
  output_tree->Branch("clusters_on_midstream_plates", &the_output_root_data_structure->clusters_on_midstream_plates);
}

bool histogram_manager::fill(input_root_data_structure* the_input_root_data_structure)
{
  double dummy_double_1, dummy_double_2, dummy_double_3, dummy_double_4;
  
  if (the_output_data_structure->single_cluster_on_all_upstream_plates)
  {
    /*
    cout << "--- " << the_output_data_structure->mc_guessed_vertex.x() << " "
  << the_output_data_structure->mc_guessed_vertex.y() << " "
  << the_output_data_structure->mc_guessed_vertex.z() << endl;
  
    cout << the_output_data_structure->reco_guessed_vertex.x() << " "
  << the_output_data_structure->reco_guessed_vertex.y() << " "
  << the_output_data_structure->reco_guessed_vertex.z() << endl;
    */
  }
  
  the_output_root_data_structure->plate_of_cluster.clear();
  the_output_root_data_structure->x_of_cluster_center.clear();
  the_output_root_data_structure->y_of_cluster_center.clear();
  the_output_root_data_structure->z_of_cluster_center.clear();
  the_output_root_data_structure->x_of_cluster_direction.clear();
  the_output_root_data_structure->y_of_cluster_direction.clear();
  the_output_root_data_structure->z_of_cluster_direction.clear();
  
  the_output_root_data_structure->event_classification = the_output_data_structure->event_classification;
  
  for (int i = 0; i < the_output_data_structure->actual_number_of_plates; i++)
  {
    for (int j = 0; j < the_output_data_structure->cluster_lines[i].size(); j++)
    {
      // cout << the_output_data_structure->cluster_lines[i].at(j).line_center << endl;
      
      the_output_root_data_structure->plate_of_cluster.push_back(i);
      the_output_root_data_structure->x_of_cluster_center.push_back(the_output_data_structure->cluster_lines[i].at(j).line_center.x());
      the_output_root_data_structure->y_of_cluster_center.push_back(the_output_data_structure->cluster_lines[i].at(j).line_center.y());
      the_output_root_data_structure->z_of_cluster_center.push_back(the_output_data_structure->cluster_lines[i].at(j).line_center.z());
      the_output_root_data_structure->x_of_cluster_direction.push_back(the_output_data_structure->cluster_lines[i].at(j).line_direction.x());
      the_output_root_data_structure->y_of_cluster_direction.push_back(the_output_data_structure->cluster_lines[i].at(j).line_direction.y());
      the_output_root_data_structure->z_of_cluster_direction.push_back(the_output_data_structure->cluster_lines[i].at(j).line_direction.z());
    }
  }
  
  // cout << the_output_root_data_structure->plate_of_cluster.size() << endl;
  
  output_tree->Fill();
  
  if (the_output_data_structure->event_classification == 210 && the_output_data_structure->same_vertex_for_high_energy_particles)
  {
    //cout << the_output_data_structure->same_vertex_for_high_energy_particles << endl;
    //cout << the_output_data_structure->mc_guessed_vertex << endl;
    //cout << the_output_data_structure->reco_guessed_vertex << endl;
    
    dummy_double_1 = the_output_data_structure->mc_guessed_vertex.x();
    dummy_double_2 = the_output_data_structure->mc_guessed_vertex.y();
    //dummy_double_3 = the_output_data_structure->reco_guessed_vertex.x();
    //dummy_double_4 = the_output_data_structure->reco_guessed_vertex.y();

    h101->Fill(pow((dummy_double_1 - dummy_double_3) * (dummy_double_1 - dummy_double_3) +
                   (dummy_double_2 - dummy_double_4) * (dummy_double_2 - dummy_double_4), 0.5));
  }
  
  the_output_root_data_structure->clusters_on_midstream_plates.clear();
  
  for (int i = 0; i < the_output_data_structure->clusters_on_midstream_plates.size(); i++)
  {
    the_output_root_data_structure->clusters_on_midstream_plates.push_back(the_output_data_structure->clusters_on_midstream_plates.at(i));
  }
  
  return true;
}

bool histogram_manager::save()
{  
  //TFile* out_file = new TFile(output_root_file_name.c_str(), "RECREATE");
  
  h101->Write();
  output_tree->Write();
  
  out_file->Close();
}
