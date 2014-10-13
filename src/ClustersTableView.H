//
// file ClustersTableView.H
//
//  Copyright (C) 2014 AstraZeneca, David Cosgrove
//
//   @@ All Rights Reserved @@
//  This file is part of SVDClus.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the source tree.

// This class, derived from QTableView assists in the display of the clusters in svdclus.

#ifndef CLUSTERSTABLEVIEW_H
#define CLUSTERSTABLEVIEW_H

#include <QSize>
#include <QTableView>

class QAction;
class QEvent;
class QHelpEvent;

// ****************************************************************************

class ClustersTableView : public QTableView {

public :

  QSize minimumSizeHint() const;
  QSize sizeHint() const;

};

#endif // CLUSTERSTABLEVIEW_H
