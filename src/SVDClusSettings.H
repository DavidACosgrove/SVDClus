//
// file SmiVSettings.H
//
//  Copyright (C) 2014 AstraZeneca, David Cosgrove
//
//   @@ All Rights Reserved @@
//  This file is part of SVDClus.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the source tree.
//
// This is a class that holds all the settings that the user used to
// start SVDClus with, to make them easier to pass round the place.

#ifndef DAC_SVDCLUS_SETTINGS
#define DAC_SVDCLUS_SETTINGS

#include <iosfwd>
#include <string>
#include <boost/program_options/options_description.hpp>

// **************************************************************************

class SVDClusSettings {

public :

  SVDClusSettings( int argc , char **argv );

  void print_usage( std::ostream &os ) const;

  std::string smi_file() const { return smi_file_; }
  std::string fps_file() const { return fps_file_; }
  std::string usage_text() const { return usage_text_; }

  std::vector<std::string> data_files() const { return data_files_; }

  double gamma() const { return gamma_; }
  double sim_thresh() const { return sim_thresh_; }
  double clus_thresh() const { return clus_thresh_; }
  double tversky_alpha() const { return tversky_alpha_; }
  double tversky_beta() const { return tversky_beta_; }

  int start_num_clus() const { return start_num_clus_; }
  int stop_num_clus() const { return stop_num_clus_ == -1 ? start_num_clus_ : stop_num_clus_; }
  int clus_num_step() const { return clus_num_step_; }

  bool do_svd_clus() const { return do_svd_clus_; }
  bool do_k_means_clus() const { return do_k_means_clus_; }
  bool do_fuzzy_k_means_clus() const { return do_fuzzy_k_means_clus_; }

  bool circular_fps() const { return circular_fps_; }
  bool linear_fps() const { return linear_fps_; }
  float fuzzy_k_means_m() const { return fuzzy_k_means_m_; }

private :

  std::string smi_file_;
  std::string fps_file_;
  std::string usage_text_;
  std::vector<std::string> data_files_;
  double gamma_;
  double sim_thresh_;
  double clus_thresh_;
  double tversky_alpha_ , tversky_beta_;
  int start_num_clus_ , stop_num_clus_ , clus_num_step_;
  bool do_svd_clus_ , do_k_means_clus_ , do_fuzzy_k_means_clus_; // straight away on firing up the program
  bool circular_fps_ , linear_fps_;
  float fuzzy_k_means_m_;

  void build_program_options( boost::program_options::options_description &desc );

};

#endif
