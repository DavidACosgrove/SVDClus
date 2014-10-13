//
// file MoleculeTableView.H
//
//  Copyright (C) 2014 AstraZeneca, David Cosgrove
//
//   @@ All Rights Reserved @@
//  This file is part of SVDClus.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the source tree.

// This class, derived from QTableView assists in the display of the molecules and associated data

#ifndef MOLECULETABLEVIEW_H
#define MOLECULETABLEVIEW_H

#include <QCloseEvent>
#include <QSize>
#include <QTableView>

class QEvent;
class QHelpEvent;

// ****************************************************************************

class MoleculeTableView : public QTableView {

public :

  QSize minimumSizeHint() const {
    return QSize( 150 , 150 );
  }

  QSize sizeHint() const {
    return QSize( 200 , 200 );
  }

protected :

  // over-ride the closeEvent handler so it can't be closed
  void closeEvent( QCloseEvent *e ) {
    e->ignore();
  }

};


#endif // MOLECULETABLEVIEW_H
