//
// file MoleculeTableModel.cc
//
//  Copyright (C) 2014 AstraZeneca, David Cosgrove
//
//   @@ All Rights Reserved @@
//  This file is part of SVDClus.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the source tree.
//

#include "MoleculeRec.H"
#include "MoleculeTableModel.H"

#include <iostream>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>

#include <QBrush>

using namespace boost;
using namespace std;

// ****************************************************************************
MoleculeTableModel::MoleculeTableModel( QObject *parent ) :
  QAbstractTableModel( parent ) , col_names_( vector<QString>( 1 , QString( "SMILES" ) ) ) {

}

// ****************************************************************************
void MoleculeTableModel::add_molecules( const vector<pMolRec> &new_mols ) {

  if( new_mols.empty() ) {
    return;
  }

  rebuild_column_names();

  beginInsertRows( QModelIndex() , molecules_.size() , new_mols.size() - 1 );
  molecules_.insert( molecules_.end() , new_mols.begin() , new_mols.end() );
  endInsertRows();

}

// ****************************************************************************
int MoleculeTableModel::rowCount( const QModelIndex &parent ) const {

  return molecules_.size();

}

// ****************************************************************************
int MoleculeTableModel::columnCount( const QModelIndex &parent ) const {

  return col_names_.size();

}

// ****************************************************************************
QVariant MoleculeTableModel::data( const QModelIndex &index , int role ) const {

#ifdef NOTYET
  cout << "data for " << index.row() << " , " << index.column() << endl;
#endif

  if( index.column() >= static_cast<int>( col_names_.size() ) ||
      index.row() >= static_cast<int>( molecules_.size() ) ) {
    return QVariant();
  }

  if( role == Qt::DisplayRole ) {
    float data_val;
    if( 0 == index.column() ) {
      QVariant ret_val;
      ret_val.setValue( molecules_[index.row()] );
      return ret_val;
      return QVariant( QString( molecules_[index.row()]->smiles().c_str() ) );
    } else {
      if( molecules_[index.row()]->data( col_names_[index.column()].toLocal8Bit().data() ,
                                         data_val ) ) {
        return QVariant( data_val );
      } else {
        return QVariant(); // this molecule doesn't have any data of that name
      }
    }
  } else if( role == Qt::BackgroundRole ) {
    return QBrush( QColor( "White" ) );
  }

  return QVariant();

}

// ****************************************************************************
QVariant MoleculeTableModel::headerData( int section , Qt::Orientation orientation ,
                                         int role ) const {

  if( role != Qt::DisplayRole ) {
    return QVariant();
  }

  if( orientation == Qt::Horizontal ) {
#ifdef NOTYET
    cout << "Horizontal header for " << section << endl;
#endif
    if( section >= 0 && section < static_cast<int>( col_names_.size() ) ) {
#ifdef NOTYET
      cout << "returning " << col_names_[section].toStdString() << endl;
#endif
      return QVariant( col_names_[section] );
    } else {
#ifdef NOTYET
      cout << "returning empty wotsit" << endl;
#endif
      return QVariant();
    }
  } else {
    if( section >= 0 && section < static_cast<int>( molecules_.size() ) ) {
      return QVariant( molecules_[section]->name().c_str() );
    } else {
      return QVariant();
    }
  }

  return QVariant();

}

// ****************************************************************************
// take the new set of molecules and add to the col_names_ any new column names.
void MoleculeTableModel::rebuild_column_names() {

#ifdef NOTYET
  cout << "rebuilding column names : " << columnCount() << " and " << rowCount() << endl;
#endif
  typedef pair<string,float> DATA_PAIR;
  BOOST_FOREACH( pMolRec mol , molecules_ ) {
    BOOST_FOREACH( DATA_PAIR dp , mol->data() ) {
      QString dp_name( dp.first.c_str() );
      if( col_names_.end() == find( col_names_.begin() , col_names_.end() , dp_name ) ) {
#ifdef NOTYET
        cout << "New column " << dp.first << " at " << col_names_.size() - 1 << endl;
#endif
        beginInsertColumns( QModelIndex() , col_names_.size() , col_names_.size() + 1 );
        col_names_.push_back( dp_name );
        endInsertColumns();
      }
    }
  }

#ifdef NOTYET
  cout << "Number of columns now " << columnCount() << endl;
#endif

}

// ****************************************************************************
// count the number of molecules with fingerprints
int MoleculeTableModel::count_fingerprints() const {

  return count_if( molecules_.begin() , molecules_.end() ,
                   bind( &MoleculeRec::get_fingerprint , _1 ) );

}
