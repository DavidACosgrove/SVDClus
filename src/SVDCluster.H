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
//
// This is a simple class that holds a single cluster produced by the Spectral Clustering
// method. It's just the eigenvalue for the cluster and the cluster members with their
// contributions to the eigenvector/cluster.  Takes a threshold which selects only the
// important cluster members based on contribution value.

#ifndef SVDCLUSTER_H
#define SVDCLUSTER_H

#include "SVDClusRDKitDefs.H"

#include <string>
#include <vector>

//SVDLIBC
extern "C" {
#include "svdlib.h"
#include "svdutil.h"
}

// ****************************************************************************
class SVDCluster {

public :

  SVDCluster( double eig_val , double cluster_thresh ) :
    thresh_( cluster_thresh ) , eigen_val_( eig_val ) {}

  SVDCluster( int row_num , double eig_val , int matrix_size ,
              DMat mat , double cluster_thresh ,
              const std::vector<pMolRec> &molecules , const std::vector<float> &sil_scores );

  int size() const { return cluster_mems_.size(); }
  const std::vector<pSVDClusMem> &cluster_members() const { return cluster_mems_; }

  std::string member_name( int mem ) const;
  std::string member_smiles( int mem ) const;
  double member_cont( int mem ) const;
  pSVDClusMem member( int mem ) const;

  double eigen_value() const { return eigen_val_; }

  void add_member( const pSVDClusMem &new_mem , bool resort = true );
  void sort_members( bool descending = true );

private :

  double thresh_;
  double eigen_val_;

  std::vector<pSVDClusMem> cluster_mems_;

};

// build clusters from the clus vectors
void extract_clusters( const std::vector<pMolRec> &molecules ,
                       const std::vector<std::vector<CLUS_MEM> > &clus ,
                       const std::vector<float> &sil_scores ,
                       std::vector<pSVDCluster> &clusters );

#endif // SVDCLUSTER_H
