//
// file SVDClusRDKit.H
//
//  Copyright (C) 2014 AstraZeneca, David Cosgrove
//
//   @@ All Rights Reserved @@
//  This file is part of SVDClus.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the source tree.
//
// This is the main Qt window for svdclus.

#include <iosfwd>
#include <string>
#include <vector>

#include <QMainWindow>
#include <QString>

#include "SVDClusRDKitDefs.H"

#ifndef SVDCLUSRDKIT_H
#define SVDCLUSRDKIT_H

class ClusterWindow;
class Chronograph;
class FuzzyKMeansClustersDialog;
class KMeansClustersDialog;
class SVDClustersDialog;
class ClustersTableModel;
class ClustersTableView;
class MoleculeTableModel;
class MoleculeTableView;
class RDKit2DMolDisplay;
class RDKitMolDrawDelegate;
class SVDCluster;
class SVDClusSettings;
class QTHelpViewer; // one of my own devising, not a Qt one

class QAction;
class QItemSelection;
class QMdiArea;
class QMenu;
class QSignalMapper;
class QSlider;
class QString;
class QTableView;

// ****************************************************************************
class SVDClusRDKit : public QMainWindow {

  Q_OBJECT

public :

  SVDClusRDKit( int argc = 0 , char **argv = 0);

public slots :

  void slot_cluster_selection_changed( const QItemSelection &sel_items , const QItemSelection &desel_items );

private :

  MoleculeTableModel *mol_table_;
  MoleculeTableView *mol_table_view_;

  SVDClusSettings *settings_;
  SVDClustersDialog *svd_clusters_dialog_;
  KMeansClustersDialog *kmeans_clusters_dialog_;
  FuzzyKMeansClustersDialog *fuzzy_kmeans_clusters_dialog_;

  RDKitMolDrawDelegate *mol_draw_del_;
  QString last_dir_;

  QAction *file_read_smiles_ , *file_read_data_ , *file_quit_;
  QAction *build_svd_clusters_ , *build_k_means_clusters_ , *build_fuzzy_k_means_clusters_;
  QAction *cascade_windows_ , *tile_windows_ , *separator_act_;
  QAction *circular_fps_ , *linear_fps_ , *user_fps_;
  QAction *help_about_ , *help_index_;
  QTHelpViewer *help_viewer_;

  QMdiArea *mdi_area_;
  QMenu *windows_menu_;
  QSignalMapper *windows_mapper_;

  FP_TYPE fp_type_;

  void build_widget();
  void build_actions();
  void build_menubar();
  void build_windows_menu();

  void parse_args( int argc , char **argv );

  void read_smiles_file( const std::string &smi_file );
  void read_data_file( const std::string &data_file );

  void do_svd_clustering( double tv_alpha , double tv_beta , int num_clus_start ,
                          int num_clus_stop , int clus_num_step ,
                          double gamma , double sim_thresh , double clus_thresh ,
                          bool overlapping_clusters );
  void do_k_means_clustering( int start_num_clus , int stop_num_clus ,
                              int clus_num_step , int num_iters );
  void do_fuzzy_k_means_clustering( int start_num_clus , int stop_num_clus ,
                                    int clus_num_step , int num_iters , float m );

  pMolRec get_molecule( const std::string &mol_name );

  // send summary information to the cluster window for display in the text widget.
  void report_clus_statistics( ClusterWindow *clus_win , const std::vector<pSVDCluster> &clusters ,
                               float sil_score , bool overlapping_clusters , Chronograph &chrono );

  void build_fingerprints();
  void build_circular_fingerprints();
  void build_linear_fingerprints();
  // read a file of user-generated fingerprints, assumed to be a name followed by
  // an stream of 0 and 1, possibly with whitespace inbetween. The name can't contain
  // whitespace.
  void read_user_fingerprints( const QString &fp_file );
  bool check_fingerprints_for_clustering();

  // decode fp_type_ into a string
  QString fingerprint_label() const;

  // connect and disconnect the selection signals of the ClusterWindows from
  // this object, so that when passing selections round we don't get in a loop
  void disconnect_cluster_window_selections();
  void connect_cluster_window_selections();

  void select_molecules_in_mol_table( const std::vector<std::string> &sel_mol_names );

private slots :

  void slot_read_smiles();
  void slot_read_data_file();
  void slot_build_svd_clusters();
  void slot_build_k_means_clusters();
  void slot_build_fuzzy_k_means_clusters();
  void slot_quit();
  void slot_update_windows_menu();
  void slot_set_active_sub_window( QWidget *window );
  void slot_circular_fps();
  void slot_linear_fps();
  void slot_user_fps();
  void slot_show_help_index();
  void slot_show_help_about();
  void slot_mol_selection_changed( const QItemSelection &sel_items , const QItemSelection &desel_items );

};

#endif // SVDCLUSRDKIT_H
