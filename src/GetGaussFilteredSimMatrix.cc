//
// file GetGaussFilteredSimMatrix.cc
//
//  Copyright (C) 2014 AstraZeneca, David Cosgrove
//
//   @@ All Rights Reserved @@
//  This file is part of SVDClus.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the source tree.
//
// Takes a set of RDKit ROMols and computes a tversky distance matrix,
// applies the gaussian tranformation function, filters by threshold,
// and returns those matrix elements that pass the threshold.
// If gamma <= -0.5 (which is the default) no filtering is done.
// If tversky_alpha and tversky_beta are left at default values of 1.0,
// it's a tanimoto distance.

#include "SVDClusRDKitDefs.H"
#include "boost_tuples_and_bind.H"
#include "MoleculeRec.H"

#include <cmath>
#include <iostream>
#include <vector>

#include <boost/bind.hpp>

#include <DataStructs/BitOps.h>

using namespace boost;
using namespace std;

// ****************************************************************************
void GetGaussFilteredSimMatrix( const vector<pMolRec> &molecules ,
                                vector<MatrixEl> &dists ,
                                double sim_thresh ,
                                double gamma = -1.0 ,
                                double tversky_alpha = 1.0 ,
                                double tversky_beta = 1.0 ) {

  for( int i = 0 , is = static_cast<int>( molecules.size() ) - 1 ; i < is ; ++i ) {
    for( int j = i + 1 , js = molecules.size() ; j < js ; ++j ) {
      double sim = TverskySimilarity<ExplicitBitVect>( *(molecules[i]->get_fingerprint()) ,
                                                       *(molecules[j]->get_fingerprint()) ,
                                                tversky_alpha , tversky_beta );
      if( gamma > -0.5 ) {
        sim = exp( -1.0 * gamma * ( sim - 1.0 ) * ( sim - 1.0 ) );
      }
      if( sim > sim_thresh ) {
        dists.push_back( make_tuple( i , j , sim ) );
        dists.push_back( make_tuple( j , i , sim ) );
      }
    }
  }

  // dists must be sorted into ascending order of first value, so that it can be fed into the
  // SVD functions correctly. Using the magic spell in boost_tuples_and_bind
  sort( dists.begin() , dists.end() ,
        bind( less<int>() ,
              bind( &getTupleElement<0,MatrixEl> , _1 ) ,
              bind( &getTupleElement<0,MatrixEl> , _2 ) ) );

#ifdef NOTYET
  for( int i = 0 , is = dists.size() ; i < is ; ++i ) {
    cout << dists[i].get<0>() << " to " << dists[i].get<1>() << " : " << dists[i].get<2>() << endl;
  }
#endif

}
