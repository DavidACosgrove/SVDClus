//
// file svdclus_main.cc
//
//  Copyright (C) 2014 AstraZeneca, David Cosgrove
//
//   @@ All Rights Reserved @@
//  This file is part of SVDClus.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the source tree.
//
// This is the main function for the program for spectral clustering, using RDKit
// fingerprints and 2D molecule display.

#include "SVDClusRDKit.H"

#include <RDGeneral/versions.h>

#include <QApplication>

#include <boost/random/mersenne_twister.hpp>

#include <iostream>

// global random number generator, so we only have 1 on the go
boost::random::mt19937 gen;

// ****************************************************************************
int main( int argc , char **argv ) {

  QApplication a( argc , argv );

  std::cout << "\nSVDClus\n\nCopyright (C) 2014 AstraZeneca\n\n"
	    << "Built using RDKit version " << RDKit::rdkitVersion << "\n"
	    << "and Qt version " << QT_VERSION_STR << ".\n"
	    << "Running with Qt version " << qVersion() << "\n\n";

  SVDClusRDKit *w = new SVDClusRDKit( argc , argv );

  w->setWindowTitle( "Spectral Clustering with RDKit" );
  w->setGeometry( 50 , 50 , 700 , 650 );
  w->show();

  return a.exec();

}
