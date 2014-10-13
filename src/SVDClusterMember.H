//
// file SVDClusterMember.H
//
//  Copyright (C) 2014 AstraZeneca, David Cosgrove
//
//   @@ All Rights Reserved @@
//  This file is part of SVDClus.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the source tree.
//
// This class just holds relevant info for individual cluster members.

#ifndef SVDCLUSTERMEMBER_H
#define SVDCLUSTERMEMBER_H

#include "SVDClusRDKitDefs.H"
#include "MoleculeRec.H"

#include <string>

// ****************************************************************************
class SVDClusterMember {

public :

  SVDClusterMember( const pMolRec &mol , double c , float s = 0.0F ) :
    mol_( mol ) , cont_( c ) , sil_score_( s ) {}

  std::string name() const { return mol_->name(); }
  std::string smiles() const { return mol_->smiles(); }
  double cont() const { return cont_; }
  float sil_score() const { return sil_score_; }
  pMolRec molecule() { return mol_; }

private :

  pMolRec mol_;
  double cont_; // contribution to eigenvector that defines cluster
  float  sil_score_; // silhouette score (crisp of fuzzy, depending) of member in cluster

};

#endif // SVDCLUSTERMEMBER_H
