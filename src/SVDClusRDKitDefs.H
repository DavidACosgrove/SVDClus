//
// file SVDClusRDKitDefs.H
//
//  Copyright (C) 2014 AstraZeneca, David Cosgrove
//
//   @@ All Rights Reserved @@
//  This file is part of SVDClus.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the source tree.
//
// This is a set of definitions used here and there in RDKitSVDClus

#ifndef SVDCLUSRDKITDEFS_H
#define SVDCLUSRDKITDEFS_H

#include <boost/cstdint.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/tuple/tuple.hpp>

#include <DataStructs/ExplicitBitVect.h>

#include <QVariant>

typedef boost::shared_ptr<ExplicitBitVect> pRD_FP;

// column number, row number, matrix value, which is the order
// the svd matrix wants them in
typedef boost::tuples::tuple<int,int,double> MatrixEl;

class SVDCluster;
typedef boost::shared_ptr<SVDCluster> pSVDCluster;

class SVDClusterMember;
typedef boost::shared_ptr<SVDClusterMember> pSVDClusMem;

class MoleculeRec;
typedef boost::shared_ptr<MoleculeRec> pMolRec;

// so we can put a pSVDClusMem and pMolRec into a QVariant
Q_DECLARE_METATYPE( pSVDClusMem );
Q_DECLARE_METATYPE( pMolRec );

// for K-Means clustering, this is a cluster member. The int is a molecule number,
// the float is the distance the molecule is from the cluster centroid
typedef std::pair<int,float> CLUS_MEM;

typedef enum { CIRCULAR_FPS , LINEAR_FPS , USER_FPS , NO_FPS } FP_TYPE;
typedef enum { HIGHER_BETTER , LOWER_BETTER } COLOUR_DATA_SENSE;

// to hold the two highest values of each column in the coefficients/contributions matrix
// of a fuzzy clustering
typedef boost::tuple<float,int,float,int> TOP_PAIR;

#endif // SVDCLUSRDKITDEFS_H
