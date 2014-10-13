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
// text. The Help text is expected to be either help_index.html in the same
// directory as the executable, if that can be deduced, or pointed to by the
// environment variable env_name, both of which are passed into the constructor.

#ifndef DAC_QT_HELP_VIEWER__
#define DAC_QT_HELP_VIEWER__

#include <map>
#include <string>

#include <QMainWindow>
#include <QString>

// ***************************************************************************

class QLineEdit;
class QTreeWidget;
class QTreeWidgetItem;
class QTextBrowser;
class QUrl;

class QTHelpViewer : public QMainWindow {

  Q_OBJECT

private slots :

  void slotLinkClicked( const QUrl & );
  void slot_search_up();
  void slot_search_down();
  void slotIndexItemSelected( QTreeWidgetItem *sel_item );

protected :

  std::string prog_path_ , env_name_; // where the help text might be
  std::string help_file_guess_;
  QString source_file_; // where the help text was found to be

  QTreeWidget *index_;
  QTextBrowser *text_view_;
  QLineEdit *text_search_;

  std::map<QString , QString> anchor_map_;

  void build_widget();
  void read_help_text();
  void build_index();
  void build_toolbar();

public :

  QTHelpViewer( const std::string &prog_path , const std::string &env_name ,
		const std::string &init_help_file_guess = std::string( "" ) ,
		QWidget *parent = 0 , Qt::WindowFlags flags = 0 );

};

#endif
