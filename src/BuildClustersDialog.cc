//
// file BuildClustersDialog.cc
//
//  Copyright (C) 2014 AstraZeneca, David Cosgrove
//
//   @@ All Rights Reserved @@
//  This file is part of SVDClus.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the source tree.
//

#include "BuildClustersDialog.H"
#include "SVDClusSettings.H"

#include <iostream>

#include <QFormLayout>
#include <QIntValidator>
#include <QLabel>
#include <QLayout>
#include <QLineEdit>
#include <QPushButton>

using namespace std;

// ****************************************************************************
BuildClustersDialog::BuildClustersDialog( SVDClusSettings *initial_settings ,
                                          QWidget *parent , Qt::WindowFlags f ) :
  QDialog( parent , f ) , main_form_( 0 ) , main_vbox_( 0 ) {

  build_widget( initial_settings );

}

// ****************************************************************************
void BuildClustersDialog::get_settings( int &start_num_clus , int &stop_num_clus ,
                                        int &num_clus_step ) const {

  start_num_clus = start_num_clus_->text().toInt();
  stop_num_clus = stop_num_clus_->text().toInt();
  num_clus_step = num_clus_step_->text().toInt();

}

// ****************************************************************************
void BuildClustersDialog::build_widget( SVDClusSettings *initial_settings ) {

  if( !main_vbox_ ) {
    main_vbox_ = new QVBoxLayout;
    setLayout( main_vbox_ );
  }

  if( !main_form_ ) {
    main_form_ = new QFormLayout;
    main_vbox_->addLayout( main_form_ );
  }

#ifdef NOTYET
  cout << "BuildClustersDialog::build_widget() : " << endl;
#endif

  QIntValidator *ncval = new QIntValidator( 0 , numeric_limits<int>::max() , this );

  start_num_clus_ = new QLineEdit();
  if( initial_settings->start_num_clus() > 0 ) {
    start_num_clus_->setText( QString( "%1" ).arg( initial_settings->start_num_clus() ) );
  }
  start_num_clus_->setValidator( ncval );
  main_form_->addRow( "First Num. Clusters" , start_num_clus_ );

  stop_num_clus_ = new QLineEdit();
  if( initial_settings->stop_num_clus() > 0 ) {
    stop_num_clus_->setText( QString( "%1" ).arg( initial_settings->stop_num_clus() ) );
  } else if( initial_settings->start_num_clus() > 0 ) {
    stop_num_clus_->setText( QString( "%1" ).arg( initial_settings->start_num_clus() ) );
  }
  stop_num_clus_->setValidator( ncval );
  main_form_->addRow( "Last Num. Clusters" , stop_num_clus_ );

  num_clus_step_ = new QLineEdit();
  if( initial_settings->clus_num_step() > 0 ) {
     num_clus_step_->setText( QString( "%1" ).arg( initial_settings->clus_num_step() ) );
  }
  num_clus_step_->setValidator( ncval );
  main_form_->addRow( "Num. Clusters Step" , num_clus_step_ );

  main_vbox_->addWidget( build_action_box() );

}

// ****************************************************************************
QWidget *BuildClustersDialog::build_action_box() {

  QFrame *action_frame = new QFrame;
  action_frame->setFrameStyle( QFrame::Box );

  QHBoxLayout *hlayout = new QHBoxLayout;

  QPushButton *button = new QPushButton( "Ok" );
  hlayout->addWidget( button );
  button->setDefault( true );
  connect( button , SIGNAL( clicked() ) , this , SLOT( accept() ) );

  button = new QPushButton( "Cancel" );
  hlayout->addWidget( button );
  button->setDefault( false );
  button->setAutoDefault( false );
  connect( button , SIGNAL( clicked() ) , this , SLOT( reject() ) );

  action_frame->setLayout( hlayout );

  return action_frame;

}
