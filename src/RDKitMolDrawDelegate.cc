//
// file RDKitMolDrawDelegate.cc
//
//  Copyright (C) 2014 AstraZeneca, David Cosgrove
//
//   @@ All Rights Reserved @@
//  This file is part of SVDClus.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the source tree.
//

#include <iostream>

#include <SVDClusterMember.H>
#include <SVDClusRDKitDefs.H>
#include <RDKitMolDrawDelegate.H>

#include <QApplication>
#include <QModelIndex>
#include <QPainter>
#include <QStyleOptionViewItem>

#include <boost/lexical_cast.hpp>

using namespace std;

// in RDKitMolToQPainter.cc
void smiles_to_qpainter( const string &smiles , const QString &label ,
                         int width , int height , bool bow , QPainter &qp );

// ****************************************************************************
void RDKitMolDrawDelegate::paint( QPainter *qp , const QStyleOptionViewItem &option ,
                                  const QModelIndex &index ) const {

  if( index.model()->data( index , Qt::DisplayRole).canConvert<pSVDClusMem>() ) {
    pSVDClusMem mem = index.model()->data( index , Qt::DisplayRole ).value<pSVDClusMem>();
    draw_cluster_member( *qp , option , index , mem );
  } else if( index.model()->data( index , Qt::DisplayRole).canConvert<pMolRec>() ) {
    pMolRec mol = index.model()->data( index , Qt::DisplayRole ).value<pMolRec>();
    draw_molecule_record( *qp , option , index , mol );
  } else if( index.model()->data( index , Qt::DisplayRole ).canConvert<QString>() ) {
    draw_string( *qp , option , index , index.model()->data( index , Qt::DisplayRole ).value<QString>() );
  } else {
    qp->fillRect( option.rect , QColor( "White" ) );
    return;
  }

}

// ****************************************************************************
void RDKitMolDrawDelegate::draw_cluster_member( QPainter &qp , const QStyleOptionViewItem &option ,
                                                const QModelIndex &index ,
                                                pSVDClusMem &mem ) const {

  set_background( index , option , qp );

  // if there's nothing to draw, or the thing's so small there's no point that's it.
  if( !mem || option.rect.width() < 20 || option.rect.height() < 20 ) {
    return;
  }

  qp.save();
  qp.resetMatrix();
  qp.translate( option.rect.x() , option.rect.y() );
  QString label = QString( "%1\nCont=%2, S=%3" ).arg( mem->name().c_str() ).arg( QString::number( mem->cont() ) )
      .arg( QString::number( mem->sil_score() ) );
  smiles_to_qpainter( mem->smiles() , label , option.rect.width() , option.rect.height() , false , qp );

  qp.restore();

}

// ****************************************************************************
void RDKitMolDrawDelegate::draw_molecule_record( QPainter &qp , const QStyleOptionViewItem &option ,
                                                 const QModelIndex &index , pMolRec &mol ) const {

  set_background( index , option , qp );

  // if there's nothing to draw, or the thing's so small there's no point that's it.
  if( !mol || option.rect.width() < 20 || option.rect.height() < 20 ) {
    return;
  }

  qp.save();
  qp.resetMatrix();
  qp.translate( option.rect.x() , option.rect.y() );
  smiles_to_qpainter( mol->smiles() , QString() , option.rect.width() , option.rect.height() , false , qp );

  qp.restore();

}

// ****************************************************************************
void RDKitMolDrawDelegate::draw_string( QPainter &qp , const QStyleOptionViewItem &option ,
                                        const QModelIndex &index , QString lab ) const {

  qp.drawText( option.rect , Qt::AlignCenter , lab );

}

// ****************************************************************************
QSize RDKitMolDrawDelegate::sizeHint( const QStyleOptionViewItem &option ,
                                      const QModelIndex &index ) const {

  return QSize( 200 , 230 );

}

// ****************************************************************************
void RDKitMolDrawDelegate::set_background( const QModelIndex &index ,
                                           const QStyleOptionViewItem &option ,
                                           QPainter &qp ) const {

  if( option.state & QStyle::State_Selected ) {
    qp.setBackground( QApplication::palette().highlight().color() );
    qp.fillRect( option.rect , QApplication::palette().highlight().color() );
  } else {

    QVariant bgv = index.data( Qt::BackgroundRole );
    QBrush bg = bgv.value<QBrush>();
    qp.setBackground( bg.color() );
    qp.fillRect( option.rect , bg.color() );

  }

}
