//
// file QTHelpViewer.H
//
//  Copyright (C) 2014 AstraZeneca, David Cosgrove
//
//   @@ All Rights Reserved @@
//  This file is part of SVDClus.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the source tree.
//
// This is a general widget that displays html text, intended to be program help
// text.

#include <cstdlib>
#include <fstream>
#include <iostream>

#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QDockWidget>
#include <QFileInfo>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QMessageBox>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QTreeWidget>
#include <QtWidgets/QTreeWidgetItem>
#include <QString>
#include <QtWidgets/QTextBrowser>

#include "QTHelpViewer.H"

// icons for the buttons - taken from the QT helpviewer example
#include "back_xpm.h"
#include "forward_xpm.h"
#include "home_xpm.h"
#include "panel-arrow-up_xpm.h"
#include "panel-arrow-down_xpm.h"

using namespace std;

// **************************************************************************
QTHelpViewer::QTHelpViewer( const string &prog_path , const string &env_name ,
			    const string &init_help_file_guess ,
			    QWidget *parent , Qt::WindowFlags flags ) :
  QMainWindow( parent , flags ) ,
  prog_path_( prog_path ) , env_name_( env_name ) ,
  help_file_guess_( init_help_file_guess ) {

  build_widget();
  read_help_text();
  build_index();

  text_view_->home();

  connect( text_view_ , SIGNAL( anchorClicked( const QUrl & ) ) ,
	   this , SLOT( slotLinkClicked( const QUrl & ) ) );
}

// **************************************************************************
void QTHelpViewer::build_widget() {

  QDockWidget *dock = new QDockWidget( "Index" , this );
  index_ = new QTreeWidget( dock );
  dock->setWidget( index_ );
  addDockWidget( Qt::LeftDockWidgetArea , dock );
  // don't want the dock to be closable
  dock->setFeatures( QDockWidget::DockWidgetMovable | QDockWidget::DockWidgetFloatable );

  text_view_ = new QTextBrowser();
  setCentralWidget( text_view_ );

  build_toolbar();

}

// **************************************************************************
void QTHelpViewer::build_toolbar() {

  QAction *back_act = new QAction( QIcon( QPixmap( back ) ) , "Back" , this );
  connect( back_act , SIGNAL( triggered() ) , text_view_ , SLOT( backward() ) );

  QAction *forward_act = new QAction( QIcon( QPixmap( forward ) ) , "Forward" ,
				      this );
  connect( forward_act , SIGNAL( triggered() ) , text_view_ , SLOT( forward() ) );

  QAction *home_act = new QAction( QIcon( QPixmap( home ) ) , "Home" , this );
  connect( home_act , SIGNAL( triggered() ) , text_view_ , SLOT( home() ) );

  QToolBar *toolbar = addToolBar( "Help ToolBar" );
  toolbar->addAction( back_act );
  toolbar->addAction( forward_act );
  toolbar->addAction( home_act );

  toolbar->addSeparator();
  toolbar->addWidget( new QLabel( "Search : " ) );
  text_search_ = new QLineEdit( this );
  toolbar->addWidget( text_search_ );
  connect( text_search_ , SIGNAL( returnPressed() ) ,
	   this , SLOT( slot_search_down() ) );

  QAction *search_down = new QAction( QIcon( QPixmap( panel_arrow_down_xpm ) ) ,
				      "Search Down" , this );
  connect( search_down , SIGNAL( triggered() ) ,
	   this , SLOT( slot_search_down() ) );
  toolbar->addAction( search_down );

  QAction *search_up = new QAction( QIcon( QPixmap( panel_arrow_up_xpm ) ) ,
				    "Search Up" , this );
  connect( search_up , SIGNAL( triggered() ) , this , SLOT( slot_search_up() ) );
  toolbar->addAction( search_up );

}

// **************************************************************************
void QTHelpViewer::read_help_text() {

  // first try to find help_text
  string help_file;
  if( !prog_path_.empty() ) {
    size_t slash_pos = prog_path_.rfind( '/' );
    if( string::npos != slash_pos ) {
      help_file = prog_path_.substr( 0 , slash_pos + 1 ) +
	help_file_guess_;
    }
  }
  
  if( help_file.empty() ) {
    // try env_name_.
    char *help_path = getenv( env_name_.c_str() );
    if( !help_path ) {
      return;
    }
    help_file = help_path;
  } else {
    QFileInfo qfi( help_file.c_str() );
    if( !qfi.exists() ) {
      // try env_name_.
      char *help_path = getenv( env_name_.c_str() );
      if( !help_path ) {
	return;
      }
      help_file = help_path;
    }
  }
  QFileInfo qfi( help_file.c_str() );
  source_file_ = qfi.absoluteFilePath();
  text_view_->setSource( QUrl( QUrl::fromLocalFile( source_file_ ) ) );

}

// ***************************************************************************
void QTHelpViewer::slotLinkClicked( const QUrl &new_link ) {

  QString url( new_link.toString() );
  if( !url.startsWith( "http://" ) )
    text_view_->setSource( url );
  else {
    QMessageBox::information( this , "Help Information." ,
			      QString( "This weedy little browser can't follow\n"\
				       "URLs properly, so can't connect you to\n"\
				       "%1\nfor which it is abjectly sorry." ).
			      arg( url ) );
    return;
  }

}

// ***************************************************************************
void QTHelpViewer::slot_search_up() {

  QString query( text_search_->text() );
  if( !text_view_->find( query , QTextDocument::FindBackward ) ) {
    QApplication::beep();
    text_view_->find( query , 0 );
  }

}

// ***************************************************************************
void QTHelpViewer::slot_search_down() {

  QString query( text_search_->text() );
  if( !text_view_->find( query , 0 ) ) {
    QApplication::beep();
    text_view_->find( query , QTextDocument::FindBackward );    
  }

}

// ***************************************************************************
void QTHelpViewer::slotIndexItemSelected( QTreeWidgetItem *sel_item ) {

  map<QString,QString>::iterator p = anchor_map_.find( sel_item->text( 0 ) );
  if( p == anchor_map_.end() ) {
    QApplication::beep();
  } else {
    text_view_->setSource( p->second );
  }

}

// **************************************************************************
void QTHelpViewer::build_index() {

  QRegExp start_header( "<[Hh][2-4]" );
  QRegExp stop_header( "</H" );
  QRegExp start_anchor( "<[Aa] name=\"" );
  QRegExp stop_anchor( "\">" );
  QRegExp stop_label( "</[Aa]>" );

  ifstream ifs( source_file_.toLocal8Bit().constData() );
  if( !ifs )
    return; // we'll already know if it can't be found, probably

  QString file_text , qs_next_word , space( " " );
  string  next_word;

  QStringList header;
  header << "Index";
  index_->setHeaderLabels( header );

  connect( index_ ,
	   SIGNAL( currentItemChanged( QTreeWidgetItem * , QTreeWidgetItem * ) ) ,
	   this , SLOT( slotIndexItemSelected( QTreeWidgetItem * ) ) );

  QTreeWidgetItem *level2_item = 0 , *level3_item = 0 , *level4_item = 0;

  while( 1 ) {
    while( ifs >> next_word ) {
      qs_next_word = next_word.c_str();
      if( qs_next_word.contains( start_header ) ) {
	file_text = qs_next_word;
	while( ifs >> next_word ) {
	  qs_next_word = next_word.c_str();
	  file_text += space + qs_next_word;
	  if( qs_next_word.contains( stop_header ) )
	    break;
	}
	break;
      }
    }
    if( file_text.isEmpty() )
      break;

    int level , next_start , next_stop;
    QString anchor , label;

    level = file_text[2].digitValue();

    next_start = file_text.indexOf( start_anchor ) + 9;
    next_stop = file_text.indexOf( stop_anchor );
    anchor = file_text.mid( next_start , next_stop - next_start );

    next_start = next_stop + 2;
    next_stop = file_text.indexOf( stop_label );
    label = file_text.mid( next_start , next_stop - next_start );

    switch( level ) {
      case 2 :
	if( !level2_item )
	  level2_item = new QTreeWidgetItem( index_ );
	else
	  level2_item = new QTreeWidgetItem( index_ , level2_item );
	level2_item->setText( 0 , label );
	level3_item = level4_item = 0;
	break;
      case 3 :
	if( !level3_item )
	  level3_item = new QTreeWidgetItem( level2_item );
	else
	  level3_item = new QTreeWidgetItem( level2_item , level3_item );
	level3_item->setText( 0 , label );
	break;
      case 4 :
	if( !level3_item )
	  level3_item = new QTreeWidgetItem( level2_item );
	if( !level4_item )
	  level4_item = new QTreeWidgetItem( level3_item );
	else
	  level4_item = new QTreeWidgetItem( level3_item , level4_item );
	level4_item->setText( 0 , label );
	break;
    }

    anchor.prepend( "#" );
    anchor_map_.insert( make_pair( label , anchor ) );
    file_text.truncate( 0 );
  }

  QTreeWidgetItemIterator it( index_ );
  while( *it ) {
    (*it)->setExpanded( true );
    ++it;
  }

}
