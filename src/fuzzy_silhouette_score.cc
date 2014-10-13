//
// file fuzzy_silhouette_score.cc
//
//  Copyright (C) 2014 AstraZeneca, David Cosgrove
//
//   @@ All Rights Reserved @@
//  This file is part of SVDClus.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the source tree.
//
// Computes the fuzzy silhouette score of a set of molecules given the corresponding crisp
// score and the top 2 contributions each molecule makes to a cluster.
// Taken from Campanello et al, Fuzzy Sets and Systems, 157, 2858-2875 (2006)
// Using a hard-coded value of 1.0 for alpha, as they suggest.

#include <iostream>
#include <vector>

#include <boost/tuple/tuple.hpp>

using namespace boost;
using namespace std;

// to hold the two highest values of each molecule's contribution to the fuzzy clusters
typedef tuple<float,int,float,int> TOP_PAIR;

// ****************************************************************************
float fuzzy_silhouette_score( const vector<float> &crisp_sil_scores ,
                              const vector<TOP_PAIR> &top_pairs ) {

  float fs = 0.0F , normaliser = 0.0F;

  for( int i = 0 , is = crisp_sil_scores.size() ; i < is ; ++i ) {
#ifdef NOTYET
    cout << i << " : " << crisp_sil_scores[i] << " : " << top_pairs[i].get<0>() << " and "
         << top_pairs[i].get<2>() << endl;
#endif
    fs += ( top_pairs[i].get<0>() - top_pairs[i].get<2>() ) * crisp_sil_scores[i];
    normaliser += ( top_pairs[i].get<0>() - top_pairs[i].get<2>() );
  }

  cout << "Fuzzy silhouette score : " << fs / normaliser << endl;

  return fs / normaliser;

}
