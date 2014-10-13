//
// file ClustersTableView.cc
//
//  Copyright (C) 2014 AstraZeneca, David Cosgrove
//
//   @@ All Rights Reserved @@
//  This file is part of SVDClus.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the source tree.

#include "ClustersTableView.H"

#include <QAction>
#include <QContextMenuEvent>
#include <QMenu>

#include <iostream>

using namespace std;

// ****************************************************************************
QSize ClustersTableView::minimumSizeHint() const {
  return QSize( 150 , 150 );
}

// ****************************************************************************
QSize ClustersTableView::sizeHint() const {
  return QSize( 200 , 200 );
}
