//
// file ColourClusterMolsDialog.cc
//
//  Copyright (C) 2014 AstraZeneca, David Cosgrove
//
//   @@ All Rights Reserved @@
//  This file is part of SVDClus.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the source tree.
//

#include "ColourClusterMolsDialog.H"

#include <QComboBox>
#include <QFrame>
#include <QLabel>
#include <QLayout>
#include <QLineEdit>
#include <QPushButton>

#include <boost/foreach.hpp>

using namespace boost;
using namespace std;

// ****************************************************************************
ColourClusterMolsDialog::ColourClusterMolsDialog( const vector<string> &names ,
                                                  float init_cutoff ,
                                                  QWidget *parent ,
                                                  Qt::WindowFlags f ) :
  QDialog( parent , f ) {

  build_widget( names , init_cutoff );

}

// ****************************************************************************
void ColourClusterMolsDialog::build_widget( const vector<string> &names ,
                                            float init_cutoff ) {

  QVBoxLayout *vbox = new QVBoxLayout;
  QGridLayout *grid = new QGridLayout;

  grid->addWidget( new QLabel( "Data for Colouring" ) , 0 , 0 );
  data_names_ = new QComboBox;
  BOOST_FOREACH( string dn , names ) {
    data_names_->addItem( dn.c_str() );
  }
  grid->addWidget( data_names_ , 0 , 1 );

  grid->addWidget( new QLabel( "Active cutoff" ) , 1 , 0 );
  act_cutoff_ = new QLineEdit;
  act_cutoff_->setText( QString( "%1" ).arg( init_cutoff ) );
  grid->addWidget( act_cutoff_ , 1 , 1 );

  grid->addWidget( new QLabel( "Data direction" ) , 2 , 0 );
  data_sense_ = new QComboBox;
  data_sense_->addItem( "Higher Better" );
  data_sense_->addItem( "Lower Better" );
  grid->addWidget( data_sense_ , 2 , 1 );

  vbox->addLayout( grid );
  vbox->addWidget( build_action_box() );
  setLayout( vbox );

}

// ****************************************************************************
QWidget *ColourClusterMolsDialog::build_action_box() {

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

// ****************************************************************************
void ColourClusterMolsDialog::get_settings( int &name_num ,
                                            float &act_cutoff ,
                                            COLOUR_DATA_SENSE &ds ) const {

  name_num = data_names_->currentIndex();
  act_cutoff = act_cutoff_->text().toFloat();

  switch( data_sense_->currentIndex() ) {
  case 0 :
    ds = HIGHER_BETTER;
    break;
  case 1 :
    ds = LOWER_BETTER;
    break;
  }

}
