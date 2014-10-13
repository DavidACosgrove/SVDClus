//
// file ClustersTableModel.H
//
//  Copyright (C) 2014 AstraZeneca, David Cosgrove
//
//   @@ All Rights Reserved @@
//  This file is part of SVDClus.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the source tree.
//
// This class, derived from QAbstractTableModel, is the model that holds a set
// of clusters. It's a ragged-edged table: each cluster is in its own row,
// so each row is of a different length. The width of the table is thus the size
// of the largest cluster, which won't necessarily be the first one.

#ifndef CLUSTERSTABLEMODEL_H
#define CLUSTERSTABLEMODEL_H

#include "SVDClusRDKitDefs.H"

#include <vector>

#include <QAbstractTableModel>

// ****************************************************************************

class ClustersTableModel : public QAbstractTableModel {

  Q_OBJECT

  typedef enum { SORT_BY_SIZE , SORT_BY_EIGVAL } SortMode;

public :

  ClustersTableModel( QObject *parent = 0 );

  // functions that must be implemented
  int rowCount( const QModelIndex &parent = QModelIndex() ) const;
  int columnCount( const QModelIndex &parent = QModelIndex() ) const;
  QVariant data( const QModelIndex &index , int role = Qt::DisplayRole ) const;

  // functions that are better if implemented
  QVariant headerData( int section , Qt::Orientation orientation , int role = Qt::DisplayRole ) const;

  // this sort only works on the cluster size or eigen value, both in descending order
  void sort( int column , Qt::SortOrder order = Qt::DescendingOrder );

  // functions particular to this Model.
  // overlapping denotes whether molecules can be in more than 1 cluster. True if they can.
  void set_clusters( std::vector<pSVDCluster> new_clus , bool overlapping );

  void set_background_data_details( const std::string &data_name ,
                                    float data_cutoff ,
                                    COLOUR_DATA_SENSE ds );
  void set_cluster_colouring( bool new_val );

  // return a vector of all the names for data held by the cluster molecules
  std::vector<std::string> get_data_names() const;

  // do the qci score, only relevant for non-overlapping clusters with activity data,
  // after a cluster colouring has occurred (so we can distinguish between actives and inactives).
  // It's a measure of how well the actives and inactives are separated by the clustering.
  float calc_qci_score();

  // change the sorting mode and resort.  We'll only be sorting on the one column, and always
  // in descending order for now.
  void sort_by_size() { sort_mode_ = SORT_BY_SIZE; sort( 0 , Qt::DescendingOrder ); }
  void sort_by_eigval() { sort_mode_ = SORT_BY_EIGVAL; sort( 0 , Qt::DescendingOrder ); }

private :

  mutable int num_cols_; // because it's not necessarily a rectangular table, this must be calculated, so cache it.
  // SVDCluster has 1 cluster, all members held as SMILES strings, names and other data.
  std::vector<pSVDCluster> clusters_;
  bool overlapping_clusters_;
  COLOUR_DATA_SENSE colour_data_sense_;

  // data for colouring the background. The string is the name of the data in the molecule
  // records. The max, min and cutoff values are used to grade the colour. Max and min
  // are deduced from the data, cutoff is supplied by user.
  bool cluster_colouring_;
  std::string background_data_name_;
  float background_data_max_ , background_data_min_ , background_data_cutoff_;

  SortMode sort_mode_;

  // return data for a tooltip. Assumes that index is valid.
  QString build_tooltip( const QModelIndex &index ) const;
  // this one sets a background colour depending on value of data at index. Again, assumes
  // index is valid.
  QBrush background_brush( const QModelIndex &index ) const;
  QColor build_green( float val , float bottom_val ) const;
  QColor build_red( float val , float top_val ) const;

signals :

  void text_for_output( QString text );

};

#endif // CLUSTERSTABLEMODEL_H
