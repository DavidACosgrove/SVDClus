//
// file SVDCluster.H
//
//  Copyright (C) 2014 AstraZeneca, David Cosgrove
//
//   @@ All Rights Reserved @@
//  This file is part of SVDClus.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the source tree.

#include "MoleculeRec.H"
#include "SVDCluster.H"
#include "SVDClusterMember.H"

#include <boost/bind.hpp>
#include <boost/foreach.hpp>

#include <algorithm>
#include <functional>
#include <cmath>

using namespace boost;
using namespace std;

// ****************************************************************************
// sil_scores are the silhouette scores of Rousseeuw or Campello depending on
// whether it's overlapping or non-overlapping clusters.
SVDCluster::SVDCluster( int row_num , double eig_val , int matrix_size ,
                        DMat mat , double cluster_thresh ,
                        const vector<pMolRec> &molecules ,
                        const vector<float> &sil_scores ) :
  thresh_( cluster_thresh ) , eigen_val_( eig_val ) {

  for( int i = 0 ; i < matrix_size ; ++i ) {
    double c = mat->value[row_num][i];
    if( fabs( c ) > thresh_ ) {
      cluster_mems_.push_back( pSVDClusMem( new SVDClusterMember( molecules[i] ,
                                                                  fabs( c ) ,
                                                                  sil_scores[i] ) ) );
    }
  }

  sort_members();

}

// ****************************************************************************
string SVDCluster::member_name( int mem ) const {

  if( mem < 0 || mem >= static_cast<int>( cluster_mems_.size() ) ) {
    return string( "" );
  } else {
    return cluster_mems_[mem]->name();
  }

}

// ****************************************************************************
string SVDCluster::member_smiles( int mem ) const {

  if( mem < 0 || mem >= static_cast<int>( cluster_mems_.size() ) ) {
    return string( "" );
  } else {
    return cluster_mems_[mem]->smiles();
  }

}

// ****************************************************************************
double SVDCluster::member_cont( int mem ) const {

  if( mem < 0 || mem >= static_cast<int>( cluster_mems_.size() ) ) {
    return -1.0;
  } else {
    return cluster_mems_[mem]->cont();
  }

}

// ****************************************************************************
pSVDClusMem SVDCluster::member( int mem ) const {

  if( mem < 0 || mem >= static_cast<int>( cluster_mems_.size() ) ) {
    return pSVDClusMem();
  } else {
    return cluster_mems_[mem];
  }

}

// ****************************************************************************
void SVDCluster::add_member( const pSVDClusMem &new_mem , bool resort ) {

  if( fabs( new_mem->cont() ) > thresh_ ) {
    cluster_mems_.push_back( new_mem );
    if( resort ) {
      sort_members();
    }
  }

}

// ****************************************************************************
// sort cluster into descending order of absolute c values
void SVDCluster::sort_members( bool descending ) {

  if( descending ) {
    sort( cluster_mems_.begin() , cluster_mems_.end() ,
          bind( greater<double>() ,
                bind( &SVDClusterMember::cont , _1 ) ,
                bind( &SVDClusterMember::cont , _2 ) ) );
  } else {
    sort( cluster_mems_.begin() , cluster_mems_.end() ,
          bind( less<double>() ,
                bind( &SVDClusterMember::cont , _1 ) ,
                bind( &SVDClusterMember::cont , _2 ) ) );
  }
}

// ***************************************************************************
// build clusters from the clus vectors
void extract_clusters( const vector<pMolRec> &molecules ,
                       const vector<vector<CLUS_MEM> > &clus ,
                       const vector<float> &sil_scores ,
                       vector<pSVDCluster> &clusters ) {

  BOOST_FOREACH( vector<CLUS_MEM> c , clus ) {
    clusters.push_back( pSVDCluster( new SVDCluster( 0.0 , -1.0 ) ) );
    BOOST_FOREACH( CLUS_MEM m , c ) {
#ifdef NOTYET
      cout << "adding member " << m << " : " << molecules[m.first]->name() << " to cluster " << clusters.size() << endl;
#endif
      clusters.back()->add_member( pSVDClusMem( new SVDClusterMember( molecules[m.first] , m.second ,
                                                                      sil_scores[m.first] ) ), false );
    }
    clusters.back()->sort_members( false );
  }

}
