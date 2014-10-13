//
// File DoFuzzyKMeansCluster.cc
//
//  Copyright (C) 2014 AstraZeneca, David Cosgrove
//
//   @@ All Rights Reserved @@
//  This file is part of SVDClus.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the source tree.
//
// Takes a vector of MoleculeRec objects and produces k clusters by fuzzy k-means clustering.
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
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/BitOps.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

using namespace boost;
using namespace std;
using namespace RDKit;

extern boost::random::mt19937 gen;

namespace bl = boost::lambda;
typedef pair<int,float> CLUS_MEM;

// ****************************************************************************
// each molecule's coefficients must add up to 1
void normalise_coefficients( int num_mols ,
                             vector<vector<float> > &clus_coeffs ) {

  for( int i = 0 ; i < num_mols ; ++i ) {
    float sum_of_coeffs = accumulate( clus_coeffs[i].begin() , clus_coeffs[i].end() , 0.0 );
    transform( clus_coeffs[i].begin() , clus_coeffs[i].end() , clus_coeffs[i].begin() ,
               bl::bind( divides<float>() , bl::_1 , sum_of_coeffs ) );
  }

}

// ****************************************************************************
// give each molecule a random contribution to each cluster, between 0.0 and 1.0
void initialise_clusters( int num_mols , int num_clusters ,
                          vector<vector<float> > &clus_coeffs ) {

  clus_coeffs = vector<vector<float> >( num_mols , vector<float>() );

  uniform_real<float> uni_dist( 0 , 1 );
  variate_generator<boost::random::mt19937 &,uniform_real<float> > uni( gen , uni_dist );

  // clus_coeffs[i][j] is the contribution of fingerprint i to cluster j
  for( int i = 0 ; i < num_mols ; ++i ) {
    clus_coeffs[i].reserve( num_clusters );
    for( int j = 0 ; j < num_clusters ; ++j ) {
      clus_coeffs[i].push_back( uni() );
    }
  }

  normalise_coefficients( num_mols , clus_coeffs );

}

// ****************************************************************************
void calculate_fuzzy_centroids( const vector<vector<float> > &fps ,
                                int num_clusters , float m ,
                                const vector<vector<float> > &clus_coeffs ,
                                vector<vector<float> > &centroids ) {

  centroids = vector<vector<float> >( num_clusters ,
                                      vector<float>( fps.front().size() , 0.0 ) );

  // the fuzzy centroid is the weighted mean of all the fps, with the weight being
  // the contribution of the fp to the cluster. The centroid is then normalised by
  // the sum of the weights. There's a parameter m, that the coefficients and other
  // things are raised to the power of, which is some measure of the degree of
  // fuzziness. It can be between 1 and infinity, but most people use 2, which is
  // hard-coded for now.
  for( int i = 0 ; i < num_clusters ; ++i ) {
    for( int j = 0 , js = fps.size() ; j < js ; ++j ) {
      for( int k = 0 , ks = fps[j].size() ; k < ks ; ++k ) {
        centroids[i][k] += fps[j][k] * pow( clus_coeffs[j][i] , m );
      }
    }
  }

  for( int i = 0 ; i < num_clusters ; ++i ) {
    float norm_i = 0.0;
    for( int j = 0 , js = fps.size() ; j < js ; ++j ) {
      norm_i += pow( clus_coeffs[j][i] , m );
    }
#ifdef NOTYET
    cout << i << " : norm_i = " << norm_i << endl;
#endif
    transform( centroids[i].begin() , centroids[i].end() , centroids[i].begin() ,
               bl::bind( divides<float>() , bl::_1 , norm_i ) );
  }

#ifdef NOTYET
  cout << "New centroids" << endl;
  BOOST_FOREACH( vector<float> c , centroids ) {
    for_each( c.begin() , c.end() , cout << bl::_1 << " " );
    cout << endl;
  }
#endif

}

// ****************************************************************************
void calculate_fp_to_centroid_sq_dists( const vector<vector<float> > &fps ,
                                        const vector<vector<float> > &centroids ,
                                        vector<vector<float> > &sq_dists ) {

  sq_dists.clear();
  int num_cents = centroids.size();
  sq_dists.reserve( fps.size() );
  BOOST_FOREACH( const vector<float> fp , fps ) {
    sq_dists.push_back( vector<float>() );
    sq_dists.back().reserve( num_cents );
    BOOST_FOREACH( const vector<float> c , centroids ) {
      float dist = inner_product( fp.begin() , fp.end() , c.begin() ,
                                  0.0 , plus<float>() ,
                                  ( bl::_1 - bl::_2 ) * ( bl::_1 - bl::_2 ) );
      sq_dists.back().push_back( dist );
    }
  }

}

// ****************************************************************************
void update_coefficients( const vector<vector<float> > &fps ,
                          int num_clusters , float m ,
                          const vector<vector<float> > &centroids ,
                          vector<vector<float> > &clus_coeffs ) {

  vector<vector<float> > sq_dists;
  calculate_fp_to_centroid_sq_dists( fps , centroids , sq_dists );

  for( int i = 0 , is = fps.size() ; i < is ; ++i ) {
    for( int j = 0 ; j < num_clusters ; ++j ) {

      clus_coeffs[i][j] = 0.0;
      for( int k = 0 ; k < num_clusters ; ++k ) {
        // raise distances to power 2 / ( m - 1 ) but the distances are already squared
        clus_coeffs[i][j] += pow( sq_dists[i][j] / sq_dists[i][k] , float( 1.0 / ( m - 1.0 ) ) );
      }
      clus_coeffs[i][j] = 1.0 / clus_coeffs[i][j];
    }
  }

#ifdef NOTYET
  int i = 0;
  BOOST_FOREACH( vector<float> c , clus_coeffs ) {
    cout << "coeffs for " << i++ << " : ";
    for_each( c.begin() , c.end() , cout << bl::_1 << " " );
    cout << endl;
  }
#endif

}

// ****************************************************************************
float coeffs_difference( const vector<vector<float> > &coeffs1 ,
                         const vector<vector<float> > &coeffs2 ) {

  float diff = 0.0;

  for( int i = 0 , is = coeffs1.size() ; i < is ; ++i ) {
    diff += inner_product( coeffs1[i].begin() , coeffs1[i].end() ,
                           coeffs2[i].begin() , 0.0 , plus<float>() ,
                           ( bl::_1 - bl::_2 ) * ( bl::_1 - bl::_2 ) );
  }

#ifdef NOTYET
  cout << "coeffs diff = " << diff << endl;
#endif

  return diff;

}

// ****************************************************************************
float calc_j_score( float m , const vector<vector<float> > &fps ,
                    const vector<vector<float> > &centroids ,
                    const vector<vector<float> > &clus_coeffs ) {

  vector<vector<float> > sq_dists;
  calculate_fp_to_centroid_sq_dists( fps , centroids , sq_dists );
  float j_score = 0.0;
  for( int i = 0 , is = fps.size() ; i < is ; ++i ) {
    for( int j = 0 , js = centroids.size() ; j < js ; ++j ) {
      j_score += pow( clus_coeffs[i][j] , m ) * sq_dists[i][j];
    }
  }

  return j_score;

}

// ****************************************************************************
float generate_fuzzy_clusters( const vector<vector<float> > &fps ,
                               int num_clusters , float m ,
                               vector<vector<float> > &centroids ,
                               vector<vector<float> > &clus_coeffs ) {

  for( int i = 0 ; i < 100000 ; ++i ) {
#ifdef NOTYET
    cout << "Clustering ITERATION " << i << endl;
#endif
    calculate_fuzzy_centroids( fps , num_clusters , m , clus_coeffs , centroids );
    vector<vector<float> > new_coeffs = clus_coeffs;
    update_coefficients( fps , num_clusters , m , centroids , new_coeffs );
    if( coeffs_difference( clus_coeffs , new_coeffs ) < 1.0e-6 ) {
      break;
    }
    clus_coeffs = new_coeffs;
  }

  return calc_j_score( m , fps , centroids , clus_coeffs );

}

// *************************************************************************
void get_top_pairs( const vector<vector<float> > &coeffs ,
                    vector<TOP_PAIR> &top_pairs ) {

  top_pairs = vector<TOP_PAIR>( coeffs.size() ,
                                make_tuple( -numeric_limits<float>::max() , -1 ,
                                            -numeric_limits<float>::max() , -1 ) );

  int rank = coeffs.front().size();

#ifdef NOTYET
  cout << "coeffs.size() : " << coeffs.size() << " by rank = " << rank << endl;
#endif

  for( int i = 0 , is = coeffs.size() ; i < is ; ++i ) {
    if( coeffs[i][0] > coeffs[i][1] ) {
      top_pairs[i].get<0>() = coeffs[i][0];
      top_pairs[i].get<1>() = 0;
      top_pairs[i].get<2>() = coeffs[i][1];
      top_pairs[i].get<3>() = 1;
    } else {
      top_pairs[i].get<0>() = coeffs[i][1];
      top_pairs[i].get<1>() = 1;
      top_pairs[i].get<2>() = coeffs[i][0];
      top_pairs[i].get<3>() = 0;
    }
    for( int j = 2 ; j < rank ; ++j ) {
      if( coeffs[i][j] >= top_pairs[i].get<0>() ) {
        top_pairs[i].get<2>() = top_pairs[i].get<0>();
        top_pairs[i].get<3>() = top_pairs[i].get<1>();
        top_pairs[i].get<0>() = coeffs[i][j];
        top_pairs[i].get<1>() = j;
      } else if( coeffs[i][j] > top_pairs[i].get<2>() ) {
        top_pairs[i].get<2>() = coeffs[i][j];
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

// ****************************************************************************
// for each fingerprint, calculate the mean distance for it to each cluster, for
// use calculating the crisp silhouette score
void calc_fp_cluster_dists( const vector<vector<int> > &clus ,
                            const vector<vector<float> > &fps ,
                            vector<vector<float> > &fp_clus_dists ) {

  fp_clus_dists = vector<vector<float> >( fps.size() ,
                                          vector<float>( clus.size() , 0.0 ) );

  for( int i = 0 , is = fps.size() ; i < is ; ++i ) {
    for( int j = 0 , js = clus.size() ; j < js ; ++j ) {
      for( int k = 0 , ks = clus[j].size() ; k < ks ; ++k ) {
        int mem = clus[j][k];
        float dist = inner_product( fps[i].begin() , fps[i].end() , fps[mem].begin() ,
                                    0.0 , plus<float>() ,
                                    ( bl::_1 - bl::_2 ) * ( bl::_1 - bl::_2 ) );
        fp_clus_dists[i][j] += dist;
      }
    }
  }

  // now take the means
  for( int i = 0 , is = fps.size() ; i < is ; ++i ) {
    for( int j = 0 , js = clus.size() ; j < js ; ++j ) {
      fp_clus_dists[i][j] /= float( clus[j].size() );
    }
  }

}

// ****************************************************************************
float calculate_fuzzy_sil_score( const vector<vector<float> > &fps ,
                                 const vector<vector<float> > &coeffs ,
                                 double clus_thresh ) {

  // this is in DoSVDCluster.cc, for historical reasons.
  void extract_crisp_clusters( vector<TOP_PAIR> &top_pairs , float clus_thresh ,
                               vector<vector<int> > &crisp_clus );

  // in eponymous file
  float crisp_silhouette_score( const vector<vector<int> > &clus ,
                                const vector<vector<float> > dists ,
                                vector<float> &sil_scores );
  // likewise
  float fuzzy_silhouette_score( const vector<float> &crisp_sil_scores ,
                                const vector<TOP_PAIR> &top_pairs );
  // to calculate the silhouette scores, both crisp and fuzzy, and also for the
  // we need the row/molecule for which each column has its maximum
  // value.  For the fuzzy score, we need the second highest as well.
  vector<TOP_PAIR> top_pairs;
  get_top_pairs( coeffs , top_pairs );

  // we need the crisp clusters for the the crisp silhouette score
  vector<vector<int> > crisp_clus( coeffs.front().size() , vector<int>() );
  extract_crisp_clusters( top_pairs , clus_thresh , crisp_clus );

  // use the crisp clusters to calculated the mean distances from each molecule to each cluster
  vector<vector<float> > mol_clus_dists;
  calc_fp_cluster_dists( crisp_clus , fps , mol_clus_dists );

  vector<float> crisp_mol_sil_scores;
  crisp_silhouette_score( crisp_clus , mol_clus_dists , crisp_mol_sil_scores );
  return fuzzy_silhouette_score( crisp_mol_sil_scores , top_pairs );

}

// ****************************************************************************
float extract_fuzzy_k_means_clusters( const vector<pMolRec> &molecules ,
                                      const vector<vector<float> > &fps ,
                                      const vector<vector<float> > &coeffs ,
                                      double clus_thresh ,
                                      vector<pSVDCluster> &clusters ) {

#ifdef NOTYET
  int i = 0;
  BOOST_FOREACH( vector<float> c , coeffs ) {
    cout << "coeffs for " << i++ << " : ";
    for_each( c.begin() , c.end() , cout << bl::_1 << " " );
    cout << endl;
  }
#endif

  int num_clusters = coeffs.front().size();
  clusters = vector<pSVDCluster>();
  for( int i = 0 ; i < num_clusters ; i++ ) {
    clusters.push_back( pSVDCluster( new SVDCluster( 0.0 , clus_thresh ) ) );
  }

  for( int j = 0 , js = molecules.size() ; j < js ; ++j ) {
    for( int i = 0 , is = num_clusters ; i < is ; ++i ) {
      if( coeffs[j][i] > clus_thresh ) {
#ifdef NOTYET
        cout << "adding molecule " << molecules[j]->name() << " to cluster " << i << endl;
#endif
        clusters[i]->add_member( pSVDClusMem( new SVDClusterMember( molecules[j] ,
                                                                    coeffs[j][i] , 0.0 ) ) ) ;
#ifdef NOTYET
        for( int i = 0 , is = clusters.size() ; i < is ; ++i ) {
          cout << i << " :";
          for( int j = 0 , js = clusters[i]->size() ; j < js ; ++j ) {
            cout << " " << clusters[i]->member_name( j );
          }
          cout << endl;
        }
        cout << "finished adding" << endl;
#endif
      }
    }
  }

#ifdef NOTYET
  for( int i = 0 , is = clusters.size() ; i < is ; ++i ) {
    cout << i << " :";
    for( int j = 0 , js = clusters[i]->size() ; j < js ; ++j ) {
      cout << " " << clusters[i]->member_name( j );
    }
    cout << endl;
  }
#endif

  return calculate_fuzzy_sil_score( fps , coeffs , clus_thresh );

}

// ****************************************************************************
void DoFuzzyKMeansCluster( const vector<pMolRec> &molecules ,
                           int num_clusters , int num_iters ,
                           double clus_thresh , float m ,
                           vector<pSVDCluster> &clusters , float &sil_score ) {

  if( molecules.empty() ) {
    return;
  }

  vector<vector<float> > fps;
  get_fingerprints_as_floats( molecules , fps );

#ifdef NOTYET
  // sneath's data
  fps.clear();
  fps = vector<vector<float> >( 16 , vector<float>( 2 , 0.0 ) );
  fps[0][0] = 0; fps[0][1] = 4;
  fps[1][0] = 0; fps[1][1] = 3;
  fps[2][0] = 1; fps[2][1] = 5;
  fps[3][0] = 2; fps[3][1] = 4;
  fps[4][0] = 3; fps[4][1] = 3;
  fps[5][0] = 2; fps[5][1] = 2;
  fps[6][0] = 2; fps[6][1] = 1;
  fps[7][0] = 1; fps[7][1] = 0;
  fps[8][0] = 5; fps[8][1] = 5;
  fps[9][0] = 6; fps[9][1] = 5;
  fps[10][0] = 7; fps[10][1] = 6;
  fps[11][0] = 5; fps[11][1] = 3;
  fps[12][0] = 7; fps[12][1] = 3;
  fps[13][0] = 6; fps[13][1] = 2;
  fps[14][0] = 6; fps[14][1] = 1;
  fps[15][0] = 8; fps[15][1] = 1;
  num_clusters = 3;
#endif

#ifdef NOTYET
  cout << "Input fps " << endl;
  BOOST_FOREACH( vector<float> cd , fps ) {
    for_each( cd.begin() , cd.end() , cout << bl::_1 << " " );
    cout << endl;
  }
#endif

  vector<vector<float> > best_clus_coeffs;
  vector<vector<float> > best_centroids;
  float best_j_score = numeric_limits<float>::max();

  int num_bests = 0;
#ifdef NOTYET
  cout << "Number of clusters : " << num_clusters << endl;
#endif

  for( int i = 0 ; i < num_iters ; ++i ) {

#ifdef NOTYET
    cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << endl;
    cout << "Clustering Iteration " << i << endl;
#endif

    vector<vector<float> > clus_coeffs; // coefficient of each molecule in each cluster
    vector<vector<float> > centroids;
    // clus_coeffs comes back with rows as molecules, columns as clusters, so that
    // clus_coeffs[i][j] is the contribution of molecule i to cluster j.
    initialise_clusters( fps.size() , num_clusters , clus_coeffs );
    float this_j = generate_fuzzy_clusters( fps , num_clusters , m , centroids , clus_coeffs );

    cout << "j score = " << this_j << " best score " << best_j_score
         << "  discrim : " << fabs( this_j - best_j_score ) << endl;
    if( fabs( this_j - best_j_score ) < 1.0e-3 ) {
      // assume if we get the same score three times it's a minimum
      best_clus_coeffs = clus_coeffs;
      best_centroids = centroids;
      ++num_bests;
      if( num_bests == 3 ) {
        cout << "Breaking at iteration " << i << " num_bests = " << num_bests << endl;
        break;
      }
    }
    if( this_j < best_j_score ) {
      best_j_score = this_j;
      best_clus_coeffs = clus_coeffs;
      best_centroids = centroids;
      num_bests = 1;
    }

  }

#ifdef NOTYET
  cout << "Best j score : " << best_j_score << endl;
  int i = 0;
  BOOST_FOREACH( vector<float> c , best_clus_coeffs ) {
    cout << "coeffs for " << i++ << " : ";
    for_each( c.begin() , c.end() , cout << bl::_1 << " " );
    cout << endl;
  }
#endif

  sil_score = extract_fuzzy_k_means_clusters( molecules , fps , best_clus_coeffs ,
                                              clus_thresh , clusters );

}
