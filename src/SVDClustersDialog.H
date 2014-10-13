//
// file SVDClustersDialog.H
//
//  Copyright (C) 2014 AstraZeneca, David Cosgrove
//
//   @@ All Rights Reserved @@
//  This file is part of SVDClus.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the source tree.
//
// This class pops up a modal dialog to get settings for a spectral clustering.

#ifndef SVDCLUSTERSDIALOG_H
#define SVDCLUSTERSDIALOG_H

#include "BuildClustersDialog.H"

class SVDClusSettings;
class QCheckBox;
class QLineEdit;

// ****************************************************************************

class SVDClustersDialog : public BuildClustersDialog {

public :

  SVDClustersDialog( SVDClusSettings *initial_settings , QWidget *parent ,
                     Qt::WindowFlags f = 0 );

  void get_settings( int &start_num_clus , int &stop_num_clus ,
                     int &num_clus_step , double &tv_alpha , double &tv_beta ,
                     double &gamma , double &sim_thresh , double &clus_thresh ,
                     bool &overlapping_clusters ) const;

private :

  QLineEdit *tv_alpha_ , *tv_beta_;
  QLineEdit *gamma_;
  QLineEdit *sim_thresh_ , *clus_thresh_;
  QCheckBox *overlap_clusters_;

  void build_widget( SVDClusSettings *initial_settings );

};

#endif // SVDCLUSTERSDIALOG_H
