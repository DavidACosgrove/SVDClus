//
// file FuzzyKMeansClustersDialog.cc
//
//  Copyright (C) 2014 AstraZeneca, David Cosgrove
//
//   @@ All Rights Reserved @@
//  This file is part of SVDClus.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the source tree.
//

#include "FuzzyKMeansClustersDialog.H"
#include "SVDClusSettings.H"

#include <QFormLayout>
#include <QLabel>
#include <QLayout>
#include <QLineEdit>

#include <iostream>

using namespace std;

// ****************************************************************************
FuzzyKMeansClustersDialog::FuzzyKMeansClustersDialog( SVDClusSettings *initial_settings ,
                                                      QWidget *parent , Qt::WindowFlags f ) :
  KMeansClustersDialog( initial_settings , parent , f ) {

  build_widget( initial_settings );

}

// ****************************************************************************
void FuzzyKMeansClustersDialog::get_settings( int &start_num_clus , int &stop_num_clus ,
                                              int &num_clus_step , int &num_iters ,
                                              float &m ) const {

  KMeansClustersDialog::get_settings( start_num_clus , stop_num_clus ,
                                      num_clus_step , num_iters );

  m = m_->text().toFloat();

}

// ****************************************************************************
void FuzzyKMeansClustersDialog::build_widget( SVDClusSettings *initial_settings ) {

  if( !main_vbox_ ) {
    main_vbox_ = new QVBoxLayout;
    setLayout( main_vbox_ );
  }

  if( !main_form_ ) {
    main_form_ = new QFormLayout;
    main_vbox_->addLayout( main_form_ );
  }

  m_ = new QLineEdit( QString( "%1" ).arg( initial_settings->fuzzy_k_means_m() ) );
  main_form_->addRow( "M (fuzziness)" , m_ );

  setWindowTitle( "Fuzzy K-Means Clusters" );

}
