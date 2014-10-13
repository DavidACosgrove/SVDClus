//
// file SVDCluster.cc
//
//  Copyright (C) 2014 AstraZeneca, David Cosgrove
//
//   @@ All Rights Reserved @@
//  This file is part of SVDClus.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the source tree.
//
// This file uses SVDLIBC to find the eigenvectors and eigenvalues of matrix, and uses these
// to create a set of clusters.
// Only components above thresh are used to build the clusters, and if non_overlapping
// is true, a molecule is only put in the cluster where its contribution to the eigenvector
// is highest.

#include "GetRDKitSims.H"
#include "MoleculeRec.H"
#include "SVDClusRDKitDefs.H"
#include "SVDCluster.H"
#include "SVDClusterMember.H"

#include <cmath>
#include <limits>
#include <vector>

#include <boost/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/tuple/tuple.hpp>

//SVDLIBC
extern "C" {
#include "svdlib.h"
#include "svdutil.h"
}

using namespace boost;
using namespace std;

namespace bl = boost::lambda;

// in crisp_silhouette_score.cc
void calc_molecule_cluster_dists( const vector<vector<int> > &clus ,
                                  const vector<pMolRec> &molecules ,
                                  float tversky_alpha ,
                                  float tversky_beta ,
                                  vector<vector<float> > &mol_clus_dists );

// in eponymous file
float crisp_silhouette_score( const vector<vector<int> > &clus ,
                              const vector<vector<float> > dists ,
                              vector<float> &sil_scores );

// in eponymous file
float fuzzy_silhouette_score( const vector<float> &crisp_sil_scores ,
                              const vector<TOP_PAIR> &top_pairs );


// *************************************************************************
SMat create_svd_matrix( int matrix_size , vector<MatrixEl> &matrix ) {

#ifdef NOTYET
  cout << "Matrix elements : " << endl;
  for( int i = 0 , is = matrix.size() ; i < is ; ++i ) {
    cout << matrix[i].get<0>() << " , " << matrix[i].get<1>() << " : " << matrix[i].get<2>() << endl;
  }
#endif

  SMat ret_val = svdNewSMat( matrix_size , matrix_size , matrix.size() );
  for( int i = 0 , n = 0 ; i < matrix_size ; ++i ) {
    // cout << "col " << i << "  n = " << n << endl;
    ret_val->pointr[i] = n;
    while( n < int( matrix.size() ) && matrix[n].get<0>() == i ) {
      ret_val->rowind[n] = matrix[n].get<1>();
      ret_val->value[n] = matrix[n].get<2>();
      ++n;
    }
  }

  ret_val->pointr[ret_val->cols] = ret_val->vals;

  return ret_val;

}

// *************************************************************************
void extract_overlap_clusters( const vector<pMolRec> &molecules ,
                               int rank , DMat mat , double *S ,
                               int matrix_size , double clus_thresh ,
                               const vector<float> &crisp_mol_sil_scores ,
                               vector<pSVDCluster> &clusters ) {

  for( int i = 0 ; i < rank ; ++i ) {
    clusters.push_back( pSVDCluster( new SVDCluster( i , S[i] , matrix_size ,
                                                     mat , clus_thresh ,
                                                     molecules , crisp_mol_sil_scores ) ) );
  }

}

// *************************************************************************
void extract_non_overlap_clusters( const vector<pMolRec> &molecules ,
                                   double *S , double clus_thresh ,
                                   const vector<vector<int> > &clus ,
                                   const vector<float> &crisp_mol_sil_scores ,
                                   const vector<TOP_PAIR> &top_pairs ,
                                   vector<pSVDCluster> &clusters ) {

  for( int i = 0 , is = clus.size() ; i < is ; ++i ) {
    clusters.push_back( pSVDCluster( new SVDCluster( S[i] , clus_thresh ) ) );
    for( int j = 0 , js = clus[i].size() ; j < js ; ++j ) {
      int mem = clus[i][j];
      clusters[i]->add_member( pSVDClusMem( new SVDClusterMember( molecules[mem] ,
                                                                  top_pairs[mem].get<0>() ,
                                                                  crisp_mol_sil_scores[mem] ) ) ,
                               false );
    }
  }

  // sort clusters into descending order of absolute c values
  for( int i = 0 , is = clusters.size() ; i < is ; ++i ) {
    clusters[i]->sort_members();
  }

  // take out any empty clusters
  clusters.erase( remove_if( clusters.begin() , clusters.end() ,
                             bind( equal_to<size_t>() ,
                                   bind( &SVDCluster::size , _1 ) ,
                                   static_cast<size_t>( 0 ) ) ) , clusters.end() );

}

// *************************************************************************
void get_top_pairs( int rank , DMat mat , int matrix_size ,
                    vector<TOP_PAIR> &top_pairs ) {

  top_pairs = vector<TOP_PAIR>( matrix_size ,
                                make_tuple( -numeric_limits<float>::max() , -1 ,
                                            -numeric_limits<float>::max() , -1 ) );

#ifdef NOTYET
  for( int i = 0 ; i < matrix_size ; ++i ) {
    cout << i;
    for( int j = 0 ; j < rank ; ++j ) {
      cout << "," << mat->value[j][i];
    }
    cout << endl;
  }
#endif

  for( int i = 0 ; i < matrix_size ; ++i ) {
    if( fabs( mat->value[0][i] ) > fabs( mat->value[1][i] ) ) {
      top_pairs[i].get<0>() = fabs( mat->value[0][i] );
      top_pairs[i].get<1>() = 0;
      top_pairs[i].get<2>() = fabs( mat->value[1][i] );
      top_pairs[i].get<3>() = 1;
    } else {
      top_pairs[i].get<0>() = fabs( mat->value[1][i] );
      top_pairs[i].get<1>() = 1;
      top_pairs[i].get<2>() = fabs( mat->value[0][i] );
      top_pairs[i].get<3>() = 0;
    }
    for( int j = 2 ; j < rank ; ++j ) {
      if( fabs( mat->value[j][i] ) >= top_pairs[i].get<0>() ) {
        top_pairs[i].get<2>() = top_pairs[i].get<0>();
        top_pairs[i].get<3>() = top_pairs[i].get<1>();
        top_pairs[i].get<0>() = fabs( mat->value[j][i] );
        top_pairs[i].get<1>() = j;
      } else if( fabs( mat->value[j][i] ) > top_pairs[i].get<2>() ) {
        top_pairs[i].get<2>() = fabs( mat->value[j][i] );
        top_pairs[i].get<3>() = j;
      }
    }
  }

#ifdef NOTYET
  for( int i = 0 , is = top_pairs.size() ; i < is ; ++i ) {
    cout << i << " : (" << top_pairs[i].get<0>() << " , " << top_pairs[i].get<1>() << ") "
         << "(" << top_pairs[i].get<2>() << " , " << top_pairs[i].get<3>() << ")" << endl;
  }
#endif

}

// *************************************************************************
void extract_crisp_clusters( vector<TOP_PAIR> &top_pairs , float clus_thresh ,
                             vector<vector<int> > &crisp_clus ) {

  for( int i = 0 , is = top_pairs.size() ; i < is ; ++i ) {
    if( top_pairs[i].get<0>() > clus_thresh ) {
      crisp_clus[top_pairs[i].get<1>()].push_back( i );
    }
  }

#ifdef NOTYET
  cout << "Crisp clusters" << endl;
  for( int i = 0 , is = crisp_clus.size() ; i < is ; ++i ) {
    for_each( crisp_clus[i].begin() , crisp_clus[i].end() , cout << bl::_1 << " " );
    cout << endl;
  }
#endif

}

// *************************************************************************
// get the clusters out of the svd results matrix, returning the silhouette score
float extract_clusters( const vector<pMolRec> &molecules ,
                        float tversky_alpha , float tversky_beta ,
                        int rank , DMat mat , double *S , int matrix_size , double clus_thresh ,
                        bool overlapping_clusters ,
                        vector<pSVDCluster> &clusters ) {

  clusters.clear();

  // to calculate the silhouette scores, both crisp and fuzzy, and also for the
  // non-overlapping clusters, we need the row for which each column has its maximum
  // value.  For the fuzzy score, we need the second highest as well.
  vector<TOP_PAIR> top_pairs;
  get_top_pairs( rank , mat , matrix_size , top_pairs );

  // extract the molecules for the crisp clusters. These are the same as the non-overlapping
  // clusters, and is based on info from top_pairs.
  vector<vector<int> > crisp_clus( rank , vector<int>() );
  extract_crisp_clusters( top_pairs , clus_thresh , crisp_clus );

  // use the crisp clusters to calculated the mean distances from each molecule to each cluster
  vector<vector<float> > mol_clus_dists;
  calc_molecule_cluster_dists( crisp_clus , molecules , tversky_alpha , tversky_beta , mol_clus_dists );

  vector<float> crisp_mol_sil_scores;
  float avg_crisp_sil_score = crisp_silhouette_score( crisp_clus , mol_clus_dists , crisp_mol_sil_scores );

  if( overlapping_clusters ) {
    extract_overlap_clusters( molecules , rank , mat , S , matrix_size ,
                              clus_thresh , crisp_mol_sil_scores , clusters );
    return fuzzy_silhouette_score( crisp_mol_sil_scores , top_pairs );
  } else {
    // we've already extracted the crisp clusters once, no need to do it all again to build the clusters
    // records
    extract_non_overlap_clusters( molecules , S , clus_thresh , crisp_clus , crisp_mol_sil_scores ,
                                  top_pairs , clusters );
    return avg_crisp_sil_score;
  }

}

// ****************************************************************************
void DoSVDCluster( const vector<pMolRec> &molecules ,
                   float tversky_alpha , float tversky_beta ,
                   float gamma , int num_clusters ,
                   double clus_thresh , double sim_thresh ,
                   bool overlapping_clusters ,
                   vector<pSVDCluster> &u_clusters , float &u_sil_score ,
                   vector<pSVDCluster> &v_clusters , float &v_sil_score ) {

  GetRDKitSims grdks( molecules , tversky_alpha , tversky_beta , gamma , sim_thresh );
  vector<MatrixEl> sims = grdks.sims();

  int matrix_size = molecules.size();
  SMat svd_matrix = create_svd_matrix( matrix_size , sims );

  int iterations = 0;
  double las2end[2] = { -1.0e-30 , 1.0e-30 };
  double kappa = 1.0e-6;
  SVDRec svdlib_results = svdLAS2( svd_matrix , num_clusters , iterations , las2end , kappa );

  svdFreeSMat( svd_matrix );

  u_sil_score = extract_clusters( molecules , tversky_alpha , tversky_beta ,
                                  svdlib_results->d , svdlib_results->Ut ,
                                  svdlib_results->S , matrix_size , clus_thresh ,
                                  overlapping_clusters , u_clusters );
  v_sil_score = extract_clusters( molecules , tversky_alpha , tversky_beta ,
                                  svdlib_results->d , svdlib_results->Vt ,
                                  svdlib_results->S , matrix_size , clus_thresh ,
                                  overlapping_clusters , v_clusters );

  svdFreeSVDRec( svdlib_results );

}

