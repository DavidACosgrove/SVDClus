//
// file KMeansClustersDialog.H
//
//  Copyright (C) 2014 AstraZeneca, David Cosgrove
//
//   @@ All Rights Reserved @@
//  This file is part of SVDClus.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the source tree.
//
// This class puts up a dialog that collects information for k means clustering.

#ifndef KMEANSCLUSTERSDIALOG_H
#define KMEANSCLUSTERSDIALOG_H

#include "BuildClustersDialog.H"

class SVDClusSettings;
class QLineEdit;

// ****************************************************************************

class KMeansClustersDialog : public BuildClustersDialog {

public :

  KMeansClustersDialog( SVDClusSettings *initial_settings , QWidget *parent = 0 ,
                        Qt::WindowFlags f = 0 );

  void get_settings( int &start_num_clus , int &stop_num_clus ,
                     int &num_clus_step , int &num_iters ) const;

protected :

  void build_widget( SVDClusSettings *initial_settings );

private :

  QLineEdit *num_iters_;

};

#endif // KMEANSCLUSTERSDIALOG_H
