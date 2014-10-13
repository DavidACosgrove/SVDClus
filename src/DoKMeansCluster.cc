//
// File DoKMeansCluster.cc
//
//  Copyright (C) 2014 AstraZeneca, David Cosgrove
//
//   @@ All Rights Reserved @@
//  This file is part of SVDClus.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the source tree.
//
// Takes a vector of MoleculeRec objects and produces k clusters by k-means clustering.
// The nature of the algorithm requires Euclidean distances
// as the cluster centroid must be computed and that won't be a bitstring.

#include "MoleculeRec.H"
#include "SVDClusRDKitDefs.H"
#include "SVDCluster.H"
#include "SVDClusterMember.H"

#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/random/uniform_int.hpp>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/BitOps.h>

#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

using namespace std;
using namespace RDKit;

namespace bl = boost::lambda;

// class for comparing the .first value of the first values in each of two vectors. We're assuming
// neither is empty. I couldn't persuade bind or lambda::bind to do this as
// front() is an overloaded function with const or non-const return types which
// they couldn't resolve unambiguously.
class compare_first_values : public binary_function< vector<CLUS_MEM> , vector<CLUS_MEM> , bool > {
public :
  bool operator()( const vector<CLUS_MEM> &a , const vector<CLUS_MEM> &b ) {
    if( a.front().first < b.front().first ) {
      return true;
    } else {
      return false;
    }
  }
};

// in eponymous file
float crisp_silhouette_score( const vector<vector<int> > &clus ,
                              const vector<vector<float> > dists ,
                              vector<float> &sil_scores );

extern boost::random::mt19937 gen;

// ****************************************************************************
void initialise_centroids( int num_clusters ,
                           const vector<vector<float> > &fps ,
                           vector<vector<float> > &centroids ) {


  // A distribution is a function object.  We generate a random
  // number by calling `dist` with the generator.
  boost::random::uniform_int_distribution<> dist( 0 , fps.size() - 1 );

  vector<unsigned char> selected( fps.size() , 0 );
  for( int i = 0 ; i < num_clusters ; ++i ) {
    while( 1 ) {
      int sel_clus = dist( gen );
      if( !selected[sel_clus] ) {
        centroids.push_back( fps[sel_clus] );
        selected[sel_clus] = 1;
        break;
      }
    }
  }

}

// ****************************************************************************
void rebuild_clusters( const vector<vector<float> > &centroids ,
                       const vector<vector<float> > &fps ,
                       vector<vector<CLUS_MEM> > &clus ) {

  clus = vector<vector<CLUS_MEM> >( centroids.size() , vector<CLUS_MEM>() );

  for( int j = 0 , js = fps.size() ; j < js ; ++j ) {
    float best_dist = numeric_limits<float>::max();
    int best_centroid = -1;
    for( int i = 0 , is = centroids.size() ; i < is ; ++i ) {
      float sq_dist = inner_product( fps[j].begin() , fps[j].end() , centroids[i].begin() ,
                                     0.0 , plus<float>() ,
                                     ( bl::_1 - bl::_2 ) * ( bl::_1 - bl::_2 ) );
#ifdef NOTYET
      cout << "distance from fp " << j << " to centroid " << i << " = " << sq_dist << endl;
#endif
      if( sq_dist < best_dist ) {
        best_dist = sq_dist;
        best_centroid = i;
      }
    }
#ifdef NOTYET
    cout << j << " goes into " << best_centroid << " at dist " << best_dist << endl;
#endif
    clus[best_centroid].push_back( make_pair( j , best_dist ) );
  }

  // remove any empty clusters, though I'm guessing that's not going to happen
  clus.erase( remove_if( clus.begin() , clus.end() ,
                         bl::bind( &vector<CLUS_MEM>::empty , bl::_1 ) ) ,
              clus.end() );

  // sort the clusters to make them easier to check against each other
  // first, by value in the inner vectors
  BOOST_FOREACH( vector<CLUS_MEM> c , clus ) {
    sort( c.begin() , c.end() );
  }

  // then by value of 1st member. It's difficult to see how this shouldn't put
  // identical clusters into identical orders for comparison
  sort( clus.begin() , clus.end() , compare_first_values() );

}

// ****************************************************************************
void recalculate_centroids( const vector<vector<CLUS_MEM> > &clus ,
                            const vector<vector<float> > &fps ,
                            vector<vector<float> > &centroids ) {

  centroids.clear();
  for( int i = 0 , is = clus.size() ; i < is ; ++i ) {
#ifdef NOTYET
    cout << "Cluster " << i << " size : " << clus[i].size() << endl;
#endif
    // assume at least 1 cluster member, which is guaranteed by rebuild_clusters
    centroids.push_back( fps[clus[i].front().first] );
    for( int j = 1 , js = clus[i].size() ; j < js ; ++j ) {
      transform( centroids[i].begin() , centroids[i].end() ,
                 fps[clus[i][j].first].begin() , centroids[i].begin() , plus<float>() );
    }
    float fsize( clus[i].size() );
    transform( centroids[i].begin() , centroids[i].end() ,
               centroids[i].begin() , bl::bind( divides<float>() , bl::_1 , fsize ) );
#ifdef NOTYET
    cout << "New centroid : " << clus[i].size() << " : ";
    for_each( centroids[i].begin() , centroids[i].end() ,
              cout << bl::_1 << " " );
    cout << endl;
#endif
  }

}

// ****************************************************************************
float compute_cluster_sum_of_squares( const vector<vector<float> > &centroids ,
                                      vector<vector<CLUS_MEM> > &clusters ,
                                      const vector<vector<float> > &fps ) {

  float css = 0.0;
  for( int i = 0 , is = clusters.size() ; i < is ; ++i ) {
    for( int j = 0 , js = clusters[i].size() ; j < js ; ++j ) {
      clusters[i][j].second = inner_product( fps[clusters[i][j].first].begin() , fps[clusters[i][j].first].end() ,
                                             centroids[i].begin() , 0.0 , plus<float>() ,
                                             ( bl::_1 - bl::_2 ) * ( bl::_1 - bl::_2 ) );
      css += clusters[i][j].second;
    }
  }

  return css;

}

// ****************************************************************************
void generate_clusters( const vector<vector<float> > &fps ,
                        vector<vector<float> > &centroids ,
                        vector<vector<CLUS_MEM> > &clus ) {

  vector<vector<CLUS_MEM> > curr_clus , new_clus;

  while( 1 ) {
    rebuild_clusters( centroids , fps , new_clus );
#ifdef NOTYET
    cout << "New Clusters" << endl;
    for( int i = 0 , is = new_clus.size() ; i < is ; ++i ) {
      cout << i << " : " << new_clus[i].size() << " :: ";
      for_each( new_clus[i].begin() , new_clus[i].end() , cout << bl::bind( &CLUS_MEM::first , bl::_1 ) << " " );
      cout << endl;
    }
#endif
    if( new_clus == curr_clus ) {
      clus = new_clus;
      return;
    }
    recalculate_centroids( new_clus , fps , centroids );
    curr_clus = new_clus;
#ifdef NOTYET
    curr_css = compute_cluster_sum_of_squares( centroids , curr_clus , fps );
    cout << "Within cluster sum of squares : " << curr_css << endl;
#endif
  }

}

// ****************************************************************************
void calc_mol_to_centroid_scores( const vector<vector<float> > &fps ,
                                  const vector<vector<float> > &centroids ,
                                  vector<vector<float> > &sil_dists ) {

  for( int i = 0 , is = fps.size() ; i < is ; ++i ) {
    for( int j = 0 , js = centroids.size() ; j < js ; ++j ) {
      sil_dists[i][j] = inner_product( fps[i].begin() , fps[i].end() ,
                                       centroids[j].begin() , 0.0 , plus<float>() ,
                                       ( bl::_1 - bl::_2 ) * ( bl::_1 - bl::_2 ) );
#ifdef NOTYET
      cout << "mol " << i << " to centroid " << j << " : " << sil_dists[i][j] << endl;
#endif
    }
  }

}

// ****************************************************************************
void create_sil_clus( const vector<vector<CLUS_MEM> > &clus ,
                      vector<vector<int> > &sil_clus ) {

  BOOST_FOREACH( vector<CLUS_MEM> c , clus ) {
    sil_clus.push_back( vector<int>() );
    transform( c.begin() , c.end() , back_inserter( sil_clus.back() ) ,
               bl::bind( &CLUS_MEM::first , bl::_1 ) );
  }

}

// ****************************************************************************
void DoKMeansCluster( const vector<pMolRec> &molecules ,
                      int num_clusters , int num_iters ,
                      vector<pSVDCluster> &clusters , float &sil_score ) {

  vector<vector<float> > fps;
  get_fingerprints_as_floats( molecules , fps );

  vector<vector<CLUS_MEM> > best_clus;
  vector<vector<float> > best_centroids;
  float best_sil_score = -numeric_limits<float>::max();
  vector<float> mol_sil_scores , best_sil_scores;

  for( int i = 0 ; i < num_iters ; ++i ) {
#ifdef NOTYET
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    cout << "Clustering Iteration " << i << endl;
#endif
    vector<vector<float> > centroids;
    initialise_centroids( num_clusters , fps , centroids );

    vector<vector<CLUS_MEM> > clus;
    generate_clusters( fps , centroids , clus );

    // distances of each molecule to each cluster centroid, for calculation of silhouette score
    vector<vector<float> > sil_dists( molecules.size() , vector<float>( clus.size() , 0.0 ) );
    calc_mol_to_centroid_scores( fps , centroids , sil_dists );

    vector<vector<int> > sil_clus;
    create_sil_clus( clus , sil_clus );

#ifdef NOTYET
    for( int i = 0 , is = sil_clus.size() ; i < is ; ++i ) {
      for_each( sil_clus[i].begin() , sil_clus[i].end() , cout << bl::_1 << " " );
      cout << endl;
    }
    for( int i = 0 , is = sil_dists.size() ; i < is ; ++i ) {
      cout << i << " : " << sil_dists[i].size() << " :: ";
      for_each( sil_dists[i].begin() , sil_dists[i].end() , cout << bl::_1 << " " );
      cout << endl;
    }
#endif

    mol_sil_scores.clear();
    sil_score = crisp_silhouette_score( sil_clus , sil_dists , mol_sil_scores );
#ifdef NOTYET
    cout << "sil_score for iteration " << i << " : " << sil_score << endl;
#endif
    if( sil_score > best_sil_score ) {
#ifdef NOTYET
      cout << "Best sil score " << sil_score << " for iteration " << i << endl;
#endif
      best_sil_score = sil_score;
      best_clus = clus;
      best_sil_scores = mol_sil_scores;
      best_centroids = centroids;
    }
  }

  extract_clusters( molecules , best_clus , best_sil_scores , clusters );
  sort( clusters.begin() , clusters.end() ,
        boost::bind( greater<int>() ,
                     boost::bind( &SVDCluster::size , _1 ) ,
                     boost::bind( &SVDCluster::size , _2 ) ) );

}
