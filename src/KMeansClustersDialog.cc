//
// file KMeansClustersDialog.cc
//
//  Copyright (C) 2014 AstraZeneca, David Cosgrove
//
//   @@ All Rights Reserved @@
//  This file is part of SVDClus.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the source tree.
//

#include "KMeansClustersDialog.H"
#include "SVDClusSettings.H"

#include <QFormLayout>
#include <QFrame>
#include <QIntValidator>
#include <QLabel>
#include <QLayout>
#include <QLineEdit>
#include <QPushButton>

#include <iostream>
#include <limits>

using namespace std;

// ****************************************************************************
KMeansClustersDialog::KMeansClustersDialog( SVDClusSettings *initial_settings ,
                                            QWidget *parent , Qt::WindowFlags f ) :
  BuildClustersDialog( initial_settings , parent , f ) {

  build_widget( initial_settings );

}

// ****************************************************************************
void KMeansClustersDialog::get_settings( int &start_num_clus , int &stop_num_clus ,
                                         int &num_clus_step ,  int &num_iters ) const {

  BuildClustersDialog::get_settings( start_num_clus , stop_num_clus , num_clus_step );

  num_iters = num_iters_->text().toInt();

}

// ****************************************************************************
void KMeansClustersDialog::build_widget( SVDClusSettings *initial_settings ) {

  if( !main_vbox_ ) {
    main_vbox_ = new QVBoxLayout;
    setLayout( main_vbox_ );
  }

  if( !main_form_ ) {
    main_form_ = new QFormLayout;
    main_vbox_->addLayout( main_form_ );
  }

#ifdef NOTYET
  cout << "KMeansClustersDialog::build_widget() : " << endl;
#endif

  // still to do - make use of this validator
  QIntValidator *ncval = new QIntValidator( 0 , std::numeric_limits<int>::max() , this );
  num_iters_ = new QLineEdit( QString( "%1" ).arg( 10 ) );
  num_iters_->setValidator( ncval );
  main_form_->addRow( "Num Iterations" , num_iters_ );

  setWindowTitle( "K-Means Clusters" );

}
