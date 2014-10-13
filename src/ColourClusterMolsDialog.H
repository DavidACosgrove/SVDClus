//
// file ColourClusterMolsDialog.H
//
//  Copyright (C) 2014 AstraZeneca, David Cosgrove
//
//   @@ All Rights Reserved @@
//  This file is part of SVDClus.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the source tree.
//
// This class defines a dialog used to pick a data set and cutoff for colouring
// clusters in a cluster view.

#ifndef COLOURCLUSTERMOLSDIALOG_H
#define COLOURCLUSTERMOLSDIALOG_H

#include <string>
#include <vector>

#include <QDialog>

#include "SVDClusRDKitDefs.H"

class QComboBox;
class QLineEdit;

// ****************************************************************************
class ColourClusterMolsDialog : public QDialog {

public :

  ColourClusterMolsDialog( const std::vector<std::string> &names ,
                           float init_cutoff , QWidget *parent = 0 ,
                           Qt::WindowFlags f = 0 );

  void get_settings( int &name_num , float &act_cutoff ,
                     COLOUR_DATA_SENSE &ds ) const;

private :

  QComboBox *data_names_;
  QComboBox *data_sense_;
  QLineEdit *act_cutoff_;

  void build_widget( const std::vector<std::string> &names , float init_cutoff );
  QWidget *build_action_box();

};

#endif // COLOURCLUSTERMOLSDIALOG_H
