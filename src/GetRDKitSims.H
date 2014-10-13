//
// file GetRDKitSims.H
//
//  Copyright (C) 2014 AstraZeneca, David Cosgrove
//
//   @@ All Rights Reserved @@
//  This file is part of SVDClus.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the source tree.
//
// This class takes a set of SMILES strings and returns a set of similarities,
// filtered by a Gaussian function and then a threshold applied. It uses
// RDKit's Morgan fingerprints and a Tversky similarity to calculate the initial distances.

#ifndef GETRDKITDISTS_H
#define GETRDKITDISTS_H

#include <string>
#include <vector>

#include "SVDClusRDKitDefs.H" // sundry definitions including MatrixEl.

// ****************************************************************************

class GetRDKitSims {

public :

  GetRDKitSims( const std::vector<pMolRec> &molecules ,
                double tversky_alpha = 1.0 , double tversky_beta = 1.0 ,
                double gamma = 10.0 , double sim_thresh = 0.01 );

  std::vector<MatrixEl> sims() const { return sims_; }

private :

  std::vector<pMolRec> molecules_;

  double tversky_alpha_ , tversky_beta_; // defaults to 1.0, 1.0 i.e. tanimoto sim
  double gamma_; // for the gaussian transformation
  double sim_thresh_; // for filtering transformed similarities

  std::vector<MatrixEl> sims_;

  void build_sim_matrix();

};

#endif // GETRDKITDISTS_H
