//
// file GetRDKitSims.cc
//
//  Copyright (C) 2014 AstraZeneca, David Cosgrove
//
//   @@ All Rights Reserved @@
//  This file is part of SVDClus.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the source tree.
//

#include "GetRDKitSims.H"
#include "MoleculeRec.H"
#include "boost_tuples_and_bind.H"

#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/tuple/tuple.hpp>

#include <DataStructs/BitOps.h>

#include <cmath>

using namespace boost;
using namespace std;

// ****************************************************************************
GetRDKitSims::GetRDKitSims( const vector<pMolRec> &molecules ,
                            double tversky_alpha , double tversky_beta ,
                            double gamma , double sim_thresh ) :
  molecules_( molecules ) , tversky_alpha_( tversky_alpha ) , tversky_beta_( tversky_beta ) ,
  gamma_( gamma ) , sim_thresh_( sim_thresh ) {

  build_sim_matrix();

}

// ****************************************************************************
void GetRDKitSims::build_sim_matrix() {

  sims_.clear();

  for( int i = 0 , is = static_cast<int>( molecules_.size() ) - 1 ; i < is ; ++i ) {
    if( !molecules_[i]->get_fingerprint() ) {
      continue;
    }
    for( int j = i + 1 , js = molecules_.size() ; j < js ; ++j ) {
      if( !molecules_[j]->get_fingerprint() ) {
        continue;
      }
      double sim = TverskySimilarity<ExplicitBitVect>( *(molecules_[i]->get_fingerprint()) ,
                                                       *(molecules_[j]->get_fingerprint()) ,
                                                       tversky_alpha_ , tversky_beta_ );
#ifdef NOTYET
      cout << i << " , " << j << " -> " << sim << endl;
#endif
      sim = exp( -1.0 * gamma_ * ( sim - 1.0 ) * ( sim - 1.0 ) );
      if( sim > sim_thresh_ ) {
        sims_.push_back( make_tuple( i , j , sim ) );
      }

      if( tversky_alpha_ == tversky_beta_ ) {
        sims_.push_back( make_tuple( j , i , sim ) );
      } else {
        sim = TverskySimilarity<ExplicitBitVect>( *(molecules_[j]->get_fingerprint()) ,
                                                  *(molecules_[i]->get_fingerprint()) ,
                                                  tversky_alpha_ , tversky_beta_ );
#ifdef NOTYET
        cout << j << " , " << i << " -> " << sim << endl;
#endif
        sims_.push_back( make_tuple( j , i , sim ) );
      }
    }
  }

  // dists must be sorted into ascending order of first value, so that it can be fed into the
  // SVD functions correctly. Using the magic spell in boost_tuples_and_bind
  sort( sims_.begin() , sims_.end() ,
        bind( less<int>() ,
              bind( &getTupleElement<0,MatrixEl> , _1 ) ,
              bind( &getTupleElement<0,MatrixEl> , _2 ) ) );

}
