//
// file ClusterWindow.H
//
//  Copyright (C) 2014 AstraZeneca, David Cosgrove
//
//   @@ All Rights Reserved @@
//  This file is part of SVDClus.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the source tree.
//
// This Widget holds a set of clusters as a ClustersTableModel and
// ClustersTableView. It is intended to be fed into a QMdiArea so that multiple
// clusters can be viewed at the same time.

#ifndef CLUSTERWINDOW_H
#define CLUSTERWINDOW_H

#include "SVDClusRDKitDefs.H"

#include <vector>

#include <QWidget>

class ClustersTableModel;
class ClustersTableView;
class RDKitMolDrawDelegate;
class SVDClusRDKit;

class QAction;
class QContextMenuEvent;
class QString;
class QTextEdit;

// ****************************************************************************

class ClusterWindow : public QWidget {

  Q_OBJECT

public :

  ClusterWindow( const std::vector<pSVDCluster> &clus ,
                 bool overlapping_clusters ,
                 RDKitMolDrawDelegate *mdd , const QString &label ,
                 QWidget *parent = 0 , Qt::WindowFlags f = 0 );

  void connect_selection( SVDClusRDKit *parent_wid );
  void disconnect_selection( SVDClusRDKit *parent_wid );
  void select_molecules( const std::vector<std::string> &sel_mol_names );
  bool is_this_my_model( const ClustersTableModel *mod ) const;
  std::vector<std::string> selected_molecules() const;

public slots :

  void slot_text_to_show( QString text );

protected :

  void contextMenuEvent( QContextMenuEvent *event );

private :

  ClustersTableModel *clusters_model_;
  ClustersTableView *clusters_view_;
  QTextEdit *text_window_;
  QAction *write_clusters_ , *write_eigen_vecs_;
  QAction *colour_mols_ , *dont_colour_mols_ , *heat_map_;
  QAction *sort_by_size_ , *sort_by_eigen_val_;

  std::vector<pSVDCluster> clusters_;
  static QString last_dir_;

  void build_widget( RDKitMolDrawDelegate *mdd , const QString &label );
  void build_actions();

private slots :

  void slot_write_clusters();
  void slot_write_eigenvectors();
  void slot_colour_mols();
  void slot_dont_colour_mols();
  void slot_heat_map();
  void slot_sort_clusters_by_size();
  void slot_sort_clusters_by_eigval();

};

#endif // CLUSTERWINDOW_H
