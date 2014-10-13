//
// file MoleculeTableModel.H
//
//  Copyright (C) 2014 AstraZeneca, David Cosgrove
//
//   @@ All Rights Reserved @@
//  This file is part of SVDClus.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the source tree.
//
// This class holds a bunch of MoleculeRecs and presents them in a manner acceptable to
// a QTableView or some such.

#ifndef MOLECULETABLEMODEL_H
#define MOLECULETABLEMODEL_H

#include "SVDClusRDKitDefs.H"

#include <vector>

#include <QAbstractTableModel>
#include <QString>

// ****************************************************************************

class MoleculeTableModel : public QAbstractTableModel {

public :

  MoleculeTableModel( QObject *parent = 0 );

  // functions that must be implemented
  int rowCount( const QModelIndex &parent = QModelIndex() ) const;
  int columnCount( const QModelIndex &parent = QModelIndex() ) const;
  QVariant data( const QModelIndex &index , int role = Qt::DisplayRole ) const;

  // functions that are better if implemented
  QVariant headerData( int section , Qt::Orientation orientation , int role = Qt::DisplayRole ) const;

  // functions that are specific to MoleculeTableModel
  void add_molecules( const std::vector<pMolRec> &new_mols );
  const std::vector<pMolRec> &molecules() const { return molecules_; }

  // go through the molecules make a new set of column names
  void rebuild_column_names();

  // count the number of molecules with fingerprints
  int count_fingerprints() const;

private :

  std::vector<QString> col_names_; // in the order they come from the molecule records
  std::vector<pMolRec> molecules_;

};

#endif // MOLECULETABLEMODEL_H
