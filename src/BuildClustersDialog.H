//
// file BuildClustersDialog.H
//
//  Copyright (C) 2014 AstraZeneca, David Cosgrove
//
//   @@ All Rights Reserved @@
//  This file is part of SVDClus.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the source tree.
//
// This virtual base class has some of the clustering info that both k-means
// and svd clustering needs.

#ifndef BUILDCLUSTERSDIALOG_H
#define BUILDCLUSTERSDIALOG_H

#include <QDialog>

class SVDClusSettings;
class QFormLayout;
class QVBoxLayout;
class QLineEdit;

// ****************************************************************************

class BuildClustersDialog : public QDialog {

public :

  BuildClustersDialog( SVDClusSettings *initial_settings ,
                       QWidget *parent = 0 , Qt::WindowFlags f = 0 );

  void get_settings( int &start_num_clus , int &stop_num_clus ,
                     int &num_clus_step ) const;

protected :

  QFormLayout *main_form_;
  QVBoxLayout *main_vbox_;

  virtual void build_widget( SVDClusSettings *initial_settings );

private :

  QLineEdit *start_num_clus_ , *stop_num_clus_ , *num_clus_step_;

  QWidget *build_action_box();

};

#endif // BUILDCLUSTERSDIALOG_H
