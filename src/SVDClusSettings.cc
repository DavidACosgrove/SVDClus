//
// file SVDClusSettings.cc
//
//  Copyright (C) 2014 AstraZeneca, David Cosgrove
//
//   @@ All Rights Reserved @@
//  This file is part of SVDClus.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the source tree.

#include "SVDClusSettings.H"

#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include <iostream>

using namespace std;
namespace po = boost::program_options;

// *****************************************************************************
SVDClusSettings::SVDClusSettings( int argc , char **argv ) :
  gamma_( 10.0 ) , sim_thresh_( 0.01 ) , clus_thresh_( 0.01 ) ,
  tversky_alpha_( 1.0 ) , tversky_beta_( 1.0 ) , start_num_clus_( -1 ) ,
  stop_num_clus_( -1 ) , clus_num_step_( 1 ) ,
  do_svd_clus_( false ) , do_k_means_clus_( false ) ,
  do_fuzzy_k_means_clus_( false ) ,
  circular_fps_( false ) , linear_fps_( false ) , fuzzy_k_means_m_( 1.05 ) {

  po::options_description desc( "Allowed Options" );
  build_program_options( desc );

  po::variables_map vm;
  po::store( po::parse_command_line( argc , argv , desc ) , vm );
  po::notify( vm );

  if( vm.count( "help" ) ) {
    cout << desc << endl;
    exit( 1 );
  }

  ostringstream oss;
  oss << desc;
  usage_text_ = oss.str();

}

// ***************************************************************************
void SVDClusSettings::print_usage( ostream &os ) const {

  os << usage_text_ << endl;

}

// **************************************************************************
void SVDClusSettings::build_program_options( po::options_description &desc ) {

  desc.add_options()
    ( "help" , "Produce this help text" )
    ( "smiles-file,S" , po::value<string>( &smi_file_ ) , "Input SMILES filename" )
      ( "fingerprint-file,F" , po::value<string>( &fps_file_ ) , "Fingerprint file." )
      ( "data-file,D" , po::value<vector<string> >( &data_files_ ) , "File(s) containing data for molecules." )
      ( "gamma,G" , po::value<double>( &gamma_ ) , "Gamma value for transformation of distances (default 10.0)." )
      ( "similarity-threshold,L" , po::value<double>( &sim_thresh_ ) , "Threshold value for filtering distance matrix (after Gaussian transformation)(default 0.01)." )
      ( "cluster-threshold,C" , po::value<double>( &clus_thresh_ ) , "Threshold for adding molecule to cluster (default 0.01)." )
      ( "start-num-clusters,N" , po::value<int>( &start_num_clus_ ) , "Number of clusters to start with." )
      ( "stop-num-clusters" , po::value<int>( &stop_num_clus_ ) , "Final number of clusters. ")
      ( "cluster-number-step" , po::value<int>( &clus_num_step_ ) , "Step for number of clusters." )
      ( "tversky-alpha,A" , po::value<double>( &tversky_alpha_ ) , "Tversky alpha value (default 1.0).")
      ( "tversky-beta,B" , po::value<double>( &tversky_beta_ ) , "Tversky beta value (default 1.0).")
      ( "do-svd-clusters" , po::value<bool>( &do_svd_clus_ )->zero_tokens() , "Do SVD clustering on program start." )
      ( "do-k-means-clusters" , po::value<bool>( &do_k_means_clus_ )->zero_tokens() , "Do K-Means clustering on program start." )
      ( "do-fuzzy-k-means-clusters" , po::value<bool>( &do_fuzzy_k_means_clus_ )->zero_tokens() , "Do Fuzzy K-Means clustering on program start." )
      ( "circular-fingerprints" , po::value<bool>( &circular_fps_ )->zero_tokens() , "Build circular (ECFP-type) fingerprints." )
      ( "linear-fingerprints" , po::value<bool>( &linear_fps_ )->zero_tokens() , "Build linear (Daylight-like) fingerprints." )
      ( "fuzzy-k-means-m,M" , po::value<float>( &fuzzy_k_means_m_ ) , "M for fuzzy k-means clustering." );

}

