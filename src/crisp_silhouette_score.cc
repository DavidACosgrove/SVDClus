//
// file crisp_silhouette_score.cc
//
//  Copyright (C) 2014 AstraZeneca, David Cosgrove
//
//   @@ All Rights Reserved @@
//  This file is part of SVDClus.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the source tree.
//
// This file does the crisp silhouette score calculation of P. J. Rousseeuw.
// Journal of Computational and Applied Mathematics, 20, 53-65 (1987).
// Campello et al. (Fuzzy Sets and Systems, 157, 2858-2875 (2006) state that
// Hruschka et al have shown that this can be done just on the distances of each
// cluster member to the cluster prototypes rather than needing the full square matrix
// of member distances which is clearly a lot more tractable for large datasets. That's what
// this function is assuming.
// Strictly speaking, this measure is only valid for ratio distances, where doubling the distance
// value means the objects are twice as dissimilar, e.g. Euclidean distances, so it's probably not
// entirely valid for tanimoto type distances. This function doesn't care about niceties like that,
// it assumes you've given it distances that work.

#include <iostream>
#include <limits>
#include <vector>

#include "MoleculeRec.H"
#include "SVDClusRDKitDefs.H"

#include <DataStructs/SparseIntVect.h>
#include <DataStructs/BitOps.h>

#include <boost/foreach.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>

using namespace std;

namespace bl = boost::lambda;

// ****************************************************************************
// for each fingerprint, calculate the mean distance for it to each cluster, using
// a tversky similarity. These can then be used in crisp_silhouette_score.
void calc_molecule_cluster_dists( const vector<vector<int> > &clus ,
                                  const vector<pMolRec> &molecules ,
                                  float tversky_alpha ,
                                  float tversky_beta ,
                                  vector<vector<float> > &mol_clus_dists ) {

  // need to calculate the mean distance from each fingerprint to other members
  // in each cluster.  Start by calculating sum of distances of each fingerprint
  // to clusters.

  // however, not all fingerprints will necessarily be in a cluster. This is certainly true
  // for the spectral clustering. Calculating distances for them is a waste of time.
  vector<char> in_clus( molecules.size() , 0 );
  BOOST_FOREACH( vector<int> c , clus ) {
    BOOST_FOREACH( int cm , c ) {
      in_clus[cm] = 1;
    }
  }

  // dists - rows are fingerprints, columns are clusters. We don't want the distance
  // of fp i to itself in the cluster, but this will be 0.0 so it doesn't matter that
  // we calculate it anyway. It prevents a lot of testing if we do.
  mol_clus_dists = vector<vector<float> >( molecules.size() , vector<float>( clus.size() , 0.0F ) );

  for( int i = 0 , is = molecules.size() ; i < is ; ++i ) {
    if( !in_clus[i] ) {
      continue;
    }
    for( int j = 0 , js = clus.size() ; j < js ; ++j ) {
      for( int k = 0 , ks = clus[j].size() ; k < ks ; ++k ) {
        int mem = clus[j][k];
        float dist = 1.0 - TverskySimilarity<ExplicitBitVect>( *(molecules[i]->get_fingerprint()) ,
                                                               *(molecules[mem]->get_fingerprint()) ,
                                                               tversky_alpha , tversky_beta );
#ifdef NOTYET
        cout << i << " to " << mem << " dist = " << dist << endl;
#endif
        mol_clus_dists[i][j] += dist;
      }
    }
  }

  // now take the means
  for( int i = 0 , is = molecules.size() ; i < is ; ++i ) {
    for( int j = 0 , js = clus.size() ; j < js ; ++j ) {
      mol_clus_dists[i][j] /= float( clus[j].size() );
    }
  }

}

// ****************************************************************************
// each row in clus is a cluster, each row in dists is a molecule.
// there is 1 silhouette score per molecule.  It's a measure of how
// well each molecule fits into a cluster compared with the next nearest
// cluster it might be in.
float crisp_silhouette_score( const vector<vector<int> > &clus ,
                              const vector<vector<float> > dists ,
                              vector<float> &sil_scores ) {

  float css = 0.0;

#ifdef NOTYET
  for( int i = 0 , is = dists.size() ; i < is ; ++i ) {
    cout << i << " : ";
    for_each( dists[i].begin() , dists[i].end() , cout << bl::_1 << " " );
    cout << endl;
  }
#endif

  sil_scores = vector<float>( dists.size() , 0.0F );
  int num_clus_mems = 0;
  for( int i = 0 , is = clus.size() ; i < is ; ++i ) {
#ifdef NOTYET
    cout << "Cluster " << i << " size = " << clus[i].size() << endl;
#endif

    num_clus_mems += clus[i].size();
    if( 1 == clus[i].size() ) {
      continue;
    } else {
      for( int j = 0 , js = clus[i].size() ; j < js ; ++j ) {
        // ai is the distance of cluster member clus[i][j] to it's cluster centre (cluster i)
        // bi is the minimum distance of cluster member to the other cluster centres
        float ai = 0.0 , bi = numeric_limits<float>::max();
        int mem = clus[i][j];
        for( int k = 0 , ks = clus.size() ; k < ks ; ++k ) {
#ifdef NOTYET
          cout << mem << " to cluster " << k << " : " << dists[mem][k] << endl;
#endif
          if( i == k ) {
            ai = dists[mem][k];
          } else {
            if( dists[mem][k] < bi ) {
              bi = dists[mem][k];
            }
          }
        }
        float si = 0.0;
        if( ai != bi  ) {
          si = ( bi - ai ) / std::max( bi , ai );
        }
#ifdef NOTYET
        cout << i << " , " << j << " : mol " << mem << "  ai = " << ai << " bi = " << bi
             << " :: " << si << endl;
#endif
        sil_scores[mem] = si;
        css += si;
      }
    }
  }

  float avg_crisp_sil = css / float( num_clus_mems );
#ifdef NOTYET
  cout << "Crisp Silhouette : " << avg_crisp_sil << endl;
  cout << "individual sil scores : ";
  for_each( sil_scores.begin() , sil_scores.end() , cout << bl::_1 << " " );
  cout << endl;
#endif

  return avg_crisp_sil;
}
