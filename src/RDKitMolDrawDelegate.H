//
// file RDKitMolDrawDelegate.H
//
//  Copyright (C) 2014 AstraZeneca, David Cosgrove
//
//   @@ All Rights Reserved @@
//  This file is part of SVDClus.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the source tree.
//
// This class is derived from QAbstractItemDelegate and draws a cluster member,
// which is a SMILES string and associated data. It's not particularly well-named;
// it would be better named RDKitClusMemDrawDelegate or some such but it is as it
// is for now.

#ifndef RDKITMOLDRAWDELEGATE_H
#define RDKITMOLDRAWDELEGATE_H

#include "SVDClusRDKitDefs.H"

#include <QStyledItemDelegate>
#include <QSize>

// ****************************************************************************

class RDKitMolDrawDelegate : public QStyledItemDelegate {

public :

  RDKitMolDrawDelegate( QObject *parent = 0 ) : QStyledItemDelegate( parent ) {}

  void paint( QPainter *qp , const QStyleOptionViewItem &option ,
              const QModelIndex &index ) const;
  QSize sizeHint( const QStyleOptionViewItem &option ,
                  const QModelIndex &index ) const;

private :

  void draw_cluster_member( QPainter &qp , const QStyleOptionViewItem &option ,
                            const QModelIndex &index , pSVDClusMem &mem ) const;
  void draw_molecule_record( QPainter &qp , const QStyleOptionViewItem &option ,
                             const QModelIndex &index , pMolRec &mol ) const;
  void draw_string( QPainter &qp , const QStyleOptionViewItem &option ,
                    const QModelIndex &index , QString lab ) const;

  // at the moment, to show selected/not selected
  void set_background( const QModelIndex &index ,
                       const QStyleOptionViewItem &option ,
                       QPainter &qp ) const;

};

#endif // RDKITMOLDRAWDELEGATE_H
