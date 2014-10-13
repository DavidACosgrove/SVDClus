//
// file MoleculeRec.H
//
//  Copyright (C) 2014 AstraZeneca, David Cosgrove
//
//   @@ All Rights Reserved @@
//  This file is part of SVDClus.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the source tree.
//
// A simple class for carrying stuff about molecules about.

#ifndef MOLECULEREC_H
#define MOLECULEREC_H

#include "SVDClusRDKitDefs.H"

#include <string>
#include <vector>

// ****************************************************************************

class MoleculeRec {

public :

  MoleculeRec() {}
  MoleculeRec( const std::string &smi , const std::string &nam ) :
    smiles_( smi ) , name_( nam ) {}

  std::string smiles() const { return smiles_; }
  std::string name() const { return name_; }

  const std::vector<std::pair<std::string,float> > &data() const { return data_; }
  // returns true if data of given name found, false otherwise
  bool data( const std::string &data_name , float &data_val ) const;
  // add new data, returning true if successful, false if not. Unsuccessful in this
  // case means overwrite was false and the datum of given name already existed.
  bool add_data( const std::string &data_name , float data_val , bool overwrite = false );

  void set_fingerprint( pRD_FP fp ) {
    fp_ = fp;
  }
  pRD_FP get_fingerprint() { return fp_; }

private :

  std::string smiles_;
  std::string name_;
  std::vector<std::pair<std::string,float> > data_;
  pRD_FP fp_; // fingerprint

};

// return the contents of each molecule's fingerprint as a vector of floats. Some
// clustering algorithms such as k-means don't work with the bitstrings
void get_fingerprints_as_floats( const std::vector<pMolRec> &molecules ,
                                 std::vector<std::vector<float> > &fps );

#endif // MOLECULEREC_H
