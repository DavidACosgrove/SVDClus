//
// file SVDClustersDialog.cc
//
//  Copyright (C) 2014 AstraZeneca, David Cosgrove
//
//   @@ All Rights Reserved @@
//  This file is part of SVDClus.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the source tree.

#include "SVDClusSettings.H"
#include "SVDClustersDialog.H"

#include <QCheckBox>
#include <QDoubleValidator>
#include <QFormLayout>
#include <QFrame>
#include <QLabel>
#include <QLayout>
#include <QLineEdit>
#include <QPushButton>
#include <QString>

#include <limits>

using namespace std;

// ****************************************************************************
SVDClustersDialog::SVDClustersDialog( SVDClusSettings *initial_settings ,
                                      QWidget *parent , Qt::WindowFlags f ) :
  BuildClustersDialog( initial_settings , parent , f ) {

  build_widget( initial_settings );

}

// ****************************************************************************
void SVDClustersDialog::get_settings( int &start_num_clus , int &stop_num_clus ,
                                      int &num_clus_step , double &tv_alpha , double &tv_beta ,
                                      double &gamma ,
                                      double &sim_thresh , double &clus_thresh ,
                                      bool &overlapping_clusters ) const {

  BuildClustersDialog::get_settings( start_num_clus , stop_num_clus , num_clus_step );

  tv_alpha = tv_alpha_->text().toDouble();
  tv_beta = tv_beta_->text().toDouble();
  gamma = gamma_->text().toDouble();
  sim_thresh = sim_thresh_->text().toDouble();
  clus_thresh = clus_thresh_->text().toDouble();
  overlapping_clusters = overlap_clusters_->isChecked();

}

// ****************************************************************************
void SVDClustersDialog::build_widget( SVDClusSettings *initial_settings ) {

  if( !main_vbox_ ) {
    main_vbox_ = new QVBoxLayout;
    setLayout( main_vbox_ );
  }

  if( !main_form_ ) {
    main_form_ = new QFormLayout;
    main_vbox_->addLayout( main_form_ );
  }

#ifdef NOTYET
  cout << "SVDClustersDialog::build_widget : " << endl;
#endif

  // still to do - make use of these validators!
  QDoubleValidator *dval = new QDoubleValidator( 0.0 , 1.0 , 3 , this );

  tv_alpha_ = new QLineEdit( QString( "%1" ).arg( initial_settings->tversky_alpha() ) );
  tv_alpha_->setValidator( dval );
  main_form_->addRow( "Tversky Alpha" , tv_alpha_ );

  tv_beta_ = new QLineEdit( QString( "%1" ).arg( initial_settings->tversky_beta() ) );
  tv_beta_->setValidator( dval );
  main_form_->addRow( "Tversky Beta" , tv_beta_ );

  gamma_ = new QLineEdit( QString( "%1" ).arg( initial_settings->gamma() ) );
  main_form_->addRow( "Gaussian filter (gamma)" , gamma_ );

  sim_thresh_ = new QLineEdit( QString( "%1" ).arg( initial_settings->sim_thresh() ) );
  main_form_->addRow( "Similarity threshold" , sim_thresh_ );

  clus_thresh_ = new QLineEdit( QString( "%1" ).arg( initial_settings->clus_thresh() ) );
  main_form_->addRow( "Cluster threshold" , clus_thresh_ );

  overlap_clusters_ = new QCheckBox;
  overlap_clusters_->setChecked( true );
  main_form_->addRow( "Overlapping clusters" , overlap_clusters_ );

  setWindowTitle( "SVD Clusters" );

}
