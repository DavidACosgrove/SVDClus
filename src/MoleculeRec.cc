//
// file MoleculeRec.cc
//
//  Copyright (C) 2014 AstraZeneca, David Cosgrove
//
//   @@ All Rights Reserved @@
//  This file is part of SVDClus.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the source tree.

#include "MoleculeRec.H"

#include <boost/bind.hpp>
#include <boost/foreach.hpp>

using namespace boost;
using namespace std;

// ****************************************************************************
bool MoleculeRec::data( const string &data_name , float &data_val ) const {

  vector<pair<string,float> >::const_iterator p = find_if( data_.begin() , data_.end() ,
                                                           bind( equal_to<string>() ,
                                                                 bind( &pair<string,float>::first , _1 ) ,
                                                                 data_name ) );
  if( p == data_.end() ) {
    return false;
  } else {
    data_val = p->second;
    return true;
  }

}

// ****************************************************************************
bool MoleculeRec::add_data( const string &data_name , float data_val , bool overwrite ) {

  vector<pair<string,float> >::iterator p = find_if( data_.begin() , data_.end() ,
                                                     bind( equal_to<string>() ,
                                                           bind( &pair<string,float>::first , _1 ) ,
                                                           data_name ) );

  if( p == data_.end() ) {
    data_.push_back( make_pair( data_name , data_val ) );
    return true;
  } else {
    if( overwrite ) {
      p->second = data_val;
    } else {
      return false;
    }
  }

  return false;

}

// ****************************************************************************
// return the contents of each molecule's fingerprint as a vector of floats. Some
// clustering algorithms such as k-means don't work with the bitstrings
void get_fingerprints_as_floats( const vector<pMolRec> &molecules ,
                                 vector<vector<float> > &fps ) {

  BOOST_FOREACH( pMolRec mol , molecules ) {

    pRD_FP fp = mol->get_fingerprint();
    if( !fp ) {
      continue;
    }

    fps.push_back( vector<float>( fp->getNumBits() ) );
    for( unsigned int i = 0 ; i < fp->getNumBits() ; ++i ) {
      if( (*fp)[i] ) {
        fps.back()[i] = 1.0;
      } else {
        fps.back()[i] = 0.0;
      }
    }

  }

}
