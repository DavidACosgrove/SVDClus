//
// file SVDClusRDKit.cc
//
//  Copyright (C) 2014 AstraZeneca, David Cosgrove
//
//   @@ All Rights Reserved @@
//  This file is part of SVDClus.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the source tree.
//
// This code uses RDKit Morgan (i.e. circular) fingerprints to generate
// a distance matrix that is fed into SVDLIBC to do a spectral clustering,
// the results being rendered into a Qt Widget using RDKit's rendering
// engine.
// Options now also include using RDKit path fingerprints, and K-Means
// and Fuzzy K-Means clustering for comparison with spectral clustering.

#include "chrono.h"
#include "ClustersTableModel.H"
#include "ClustersTableView.H"
#include "ClusterWindow.H"
#include "FuzzyKMeansClustersDialog.H"
#include "KMeansClustersDialog.H"
#include "MoleculeRec.H"
#include "MoleculeTableModel.H"
#include "MoleculeTableView.H"
#include "RDKitMolDrawDelegate.H"
#include "SVDClusSettings.H"
#include "SVDClusRDKitDefs.H"
#include "SVDCluster.H"
#include "SVDClusterMember.H"
#include "SVDClustersDialog.H"
#include "SVDClusRDKit.H"
#include "QTHelpViewer.H"

#include <fstream>
#include <iostream>
#include <numeric>
#include <sstream>
#include <vector>

#include <boost/bind.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

#include <GraphMol/RDKitBase.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include <QAction>
#include <QApplication>
#include <QFileDialog>
#include <QFileInfo>
#include <QLayout>
#include <QMdiArea>
#include <QMdiSubWindow>
#include <QMenu>
#include <QMenuBar>
#include <QMessageBox>
#include <QSignalMapper>
#include <QSlider>
#include <QStatusBar>
#include <QTableView>

using namespace boost;
using namespace std;

// in eponymous file
void DoSVDCluster( const vector<pMolRec> &molecules ,
                   float tversky_alpha , float tversky_beta ,
                   float gamma , int num_clusters ,
                   double clus_thresh , double sim_thresh ,
                   bool overlapping_clusters ,
                   vector<pSVDCluster> &u_clusters , float &u_sil_score ,
                   vector<pSVDCluster> &v_clusters , float &v_sil_score );

// in eponymous file
void DoKMeansCluster( const vector<pMolRec> &molecules ,
                      int num_clusters , int num_iters ,
                      vector<pSVDCluster> &clusters , float &sil_score );
void DoKMeansCluster2( const vector<pMolRec> &molecules ,
                       int num_clusters , int num_iters ,
                       vector<pSVDCluster> &clusters , float &sil_score );

// in eponymous file
void DoFuzzyKMeansCluster( const vector<pMolRec> &molecules ,
                           int num_clusters , int num_iters , double clus_thresh ,
                           float m ,
                           vector<pSVDCluster> &clusters , float &sil_score );

// in file ClusterWindow.cc
void write_clusters( ostream &os , const vector<pSVDCluster> &clusters );

namespace RDKit {
  extern const char *rdkitVersion;
}

// *************************************************************************
SVDClusRDKit::SVDClusRDKit( int argc , char **argv ) :
  svd_clusters_dialog_( 0 ) , kmeans_clusters_dialog_( 0 ) ,
  fuzzy_kmeans_clusters_dialog_( 0 ) , last_dir_( QString( "." ) ) ,
  fp_type_( NO_FPS ) {

  build_widget();
  build_actions();
  build_menubar();

  parse_args( argc , argv );

  help_viewer_ = new QTHelpViewer( argv[0] , "SVDCLUS_HELP_TEXT" ,
                                   "svdclus_help_index.html" , this );

}

// ***************************************************************************
void SVDClusRDKit::slot_cluster_selection_changed( const QItemSelection &sel_items ,
                                                   const QItemSelection &desel_items ) {

  // get the list of selected molecules. sel_items and desel_items just give
  // ones that changed, we want all of them.
  const ClustersTableModel *clus_model = 0;
  if( !sel_items.indexes().isEmpty() ) {
    clus_model = qobject_cast<const ClustersTableModel *>( sel_items.indexes().front().model() );
  } else if( !desel_items.indexes().isEmpty() ) {
    clus_model = qobject_cast<const ClustersTableModel *>( desel_items.indexes().front().model() );
  }

  vector<string> sel_mol_names;
  if( clus_model ) {
    QList<QMdiSubWindow *> windows = mdi_area_->subWindowList();

    for( int i = 0 , is = windows.size() ; i < is ; ++i ) {
      QWidget *child = windows.at( i )->widget();
      ClusterWindow *cw = qobject_cast<ClusterWindow *>( child );
      if( cw && cw->is_this_my_model( clus_model ) ) {
        sel_mol_names = cw->selected_molecules();
        break;
      }
    }
    disconnect_cluster_window_selections();
    for( int i = 0 , is = windows.size() ; i < is ; ++i ) {
      QWidget *child = windows.at( i )->widget();
      ClusterWindow *cw = qobject_cast<ClusterWindow *>( child );
      // do the originating cluster window as well, as for overlapping clusters
      // the same molecule may occur more than once, and we want to select all
      // of them
      if( cw ) {
        cw->select_molecules( sel_mol_names );
      }
    }
    connect_cluster_window_selections();
  }

  select_molecules_in_mol_table( sel_mol_names );

}

// *************************************************************************
void SVDClusRDKit::build_widget() {

  mol_draw_del_ = new RDKitMolDrawDelegate;

  mdi_area_ = new QMdiArea;
  setCentralWidget( mdi_area_ );

  mol_table_ = new MoleculeTableModel;
  mol_table_view_ = new MoleculeTableView;
  mol_table_view_->setWindowTitle( "Molecules" );
  mol_table_view_->setModel( mol_table_ );
  mol_table_view_->setItemDelegate( mol_draw_del_ );
  connect( mol_table_view_->selectionModel() ,
           &QItemSelectionModel::selectionChanged ,
           this ,
           &SVDClusRDKit::slot_mol_selection_changed );

  mdi_area_->addSubWindow( mol_table_view_ );
  mol_table_view_->show();

}

// *************************************************************************
void SVDClusRDKit::build_actions() {

  file_read_smiles_ = new QAction( "Read SMILES" , this );
  connect( file_read_smiles_ , SIGNAL( triggered() ) , this , SLOT( slot_read_smiles() ) );
  file_read_data_ = new QAction( "Read Data" , this );
  connect( file_read_data_ , SIGNAL( triggered() ) , this , SLOT( slot_read_data_file() ) );

  file_quit_ = new QAction( "Quit" , this );
  file_quit_->setShortcut( QString( "Ctrl+Q"  ) );
  connect( file_quit_ , SIGNAL( triggered() ) , this , SLOT( slot_quit() ) );

  build_svd_clusters_ = new QAction( "Build SVD" , this );
  connect( build_svd_clusters_ , SIGNAL( triggered() ) ,
           this , SLOT( slot_build_svd_clusters() ) );

  build_k_means_clusters_ = new QAction( "Build K-Means" , this );
  connect( build_k_means_clusters_ , SIGNAL( triggered() ) ,
           this , SLOT( slot_build_k_means_clusters() ) );

  build_fuzzy_k_means_clusters_ = new QAction( "Build Fuzzy K-Means" , this );
  connect( build_fuzzy_k_means_clusters_ , SIGNAL( triggered() ) ,
           this , SLOT( slot_build_fuzzy_k_means_clusters() ) );

  circular_fps_ = new QAction( "Circular (ECFP-like)" , this );
  connect( circular_fps_ , SIGNAL( triggered() ) ,
           this , SLOT( slot_circular_fps() ) );
  linear_fps_ = new QAction( "Linear (Daylight-like)" , this );
  connect( linear_fps_ , SIGNAL( triggered() ) ,
           this , SLOT( slot_linear_fps() ) );
  user_fps_ = new QAction( "From File" , this );
  connect( user_fps_ , SIGNAL( triggered() ) ,
           this , SLOT( slot_user_fps() ) );

  help_index_ = new QAction( "Using SVDClus" , this );
  connect( help_index_ , SIGNAL( triggered() ) ,
           this , SLOT( slot_show_help_index() ) );
  help_about_ = new QAction( "About SVDClus" , this );
  connect( help_about_ , SIGNAL( triggered() ) ,
           this , SLOT( slot_show_help_about() ) );

}

// *************************************************************************
void SVDClusRDKit::build_menubar() {

  QMenu *menu = menuBar()->addMenu( "File" );

  menu->addAction( file_read_smiles_ );
  menu->addAction( file_read_data_ );
  menu->addSeparator();
  menu->addAction( file_quit_ );

  menu = menuBar()->addMenu( "Fingerprints" );
  menu->addAction( circular_fps_ );
  menu->addAction( linear_fps_ );
  menu->addAction( user_fps_ );

  menu = menuBar()->addMenu( "Clusters" );
  menu->addAction( build_svd_clusters_ );
  menu->addAction( build_k_means_clusters_ );
  menu->addAction( build_fuzzy_k_means_clusters_ );

  build_windows_menu();

  menu = menuBar()->addMenu( "Help" );
  menu->addAction( help_index_ );
  menu->addSeparator();
  menu->addAction( help_about_ );

}

// *************************************************************************
void SVDClusRDKit::build_windows_menu() {

  windows_menu_ = menuBar()->addMenu( "Windows" );
  cascade_windows_ = new QAction( "Cascade" , this );
  tile_windows_ = new QAction( "Tile" , this );
  separator_act_ = new QAction( this );
  separator_act_->setSeparator( true );

  connect( windows_menu_ , SIGNAL( aboutToShow() ) ,
           this , SLOT( slot_update_windows_menu() ) );
  connect( cascade_windows_ , SIGNAL( triggered() ) ,
           mdi_area_ , SLOT( cascadeSubWindows() ) );
  connect( tile_windows_ , SIGNAL( triggered() ) ,
           mdi_area_ , SLOT( tileSubWindows() ) );

  windows_mapper_ = new QSignalMapper( this );
  connect( windows_mapper_ , SIGNAL( mapped( QWidget * ) ) ,
           this , SLOT( slot_set_active_sub_window( QWidget * ) ) );

}

// *************************************************************************
void SVDClusRDKit::parse_args( int argc , char **argv ) {

  settings_ = new SVDClusSettings( argc , argv );


  if( !settings_->smi_file().empty() ) {
    read_smiles_file( settings_->smi_file() );
  }

  if( !settings_->fps_file().empty() ) {
    read_user_fingerprints( settings_->fps_file().c_str() );
  }

  if( settings_->circular_fps() ) {
    if( mol_table_->count_fingerprints() ) {
      cerr << "Warning - already have fingerprints. Previous ones will be over-written." << endl;
    }
    build_circular_fingerprints();
  }

  if( settings_->linear_fps() ) {
    if( mol_table_->count_fingerprints() ) {
      cerr << "Warning - already have fingerprints. Previous ones will be over-written." << endl;
    }
    build_linear_fingerprints();
  }

  vector<string> data_files = settings_->data_files();
  if( !data_files.empty() ) {
    BOOST_FOREACH( string df , data_files ) {
      read_data_file( df );
    }
  }

  if( settings_->do_svd_clus() ) {
    if( !mol_table_->count_fingerprints() ) {
      cerr << "Error - can't do SVD clustering, no fingerprints." << endl;
    } else {
      do_svd_clustering( settings_->tversky_alpha() , settings_->tversky_beta() ,
                         settings_->start_num_clus() , settings_->stop_num_clus() , settings_->clus_num_step() ,
                         settings_->gamma() , settings_->sim_thresh() ,
                         settings_->clus_thresh() , true );
    }
  }

  if( settings_->do_k_means_clus() ) {
    if( !mol_table_->count_fingerprints() ) {
      cerr << "Error - can't do K-Means clustering, no fingerprints." << endl;
    } else {
      do_k_means_clustering( settings_->start_num_clus() , settings_->stop_num_clus() ,
                             settings_->clus_num_step() , 10 );
    }
  }

  if( settings_->do_fuzzy_k_means_clus() ) {
    if( !mol_table_->count_fingerprints() ) {
      cerr << "Error - can't do Fuzzy K-Means clustering, no fingerprints." << endl;
    } else {
      do_fuzzy_k_means_clustering( settings_->start_num_clus() , settings_->stop_num_clus() ,
                                   settings_->clus_num_step() , 2 , settings_->fuzzy_k_means_m() );
    }
  }

}

// *************************************************************************
void SVDClusRDKit::read_smiles_file( const string &smi_file ) {

  ifstream ifs( smi_file.c_str() );

  vector<pMolRec> new_mols;
  while( 1 ) {
    string next_line , smi , name;
    getline( ifs , next_line );
    if( ifs.eof() || !ifs.good() ) {
      break;
    }

    istringstream iss( next_line );
    iss >> smi;
    iss >> name;
    if( iss.fail() ) {
      name = string( "Str" ) + boost::lexical_cast<string>( mol_table_->rowCount() + 1 );
    }
    new_mols.push_back( pMolRec( new MoleculeRec( smi , name ) ) );
  }

  mol_table_->add_molecules( new_mols );
  mol_table_view_->resizeColumnsToContents();
  mol_table_view_->resizeRowsToContents();
  mdi_area_->tileSubWindows();
  QString msg = QString( "Now have %1 SMILES strings" ).arg( mol_table_->rowCount() );
  statusBar()->showMessage( msg , 0 );

}

// *************************************************************************
void SVDClusRDKit::read_data_file( const string &data_file ) {

  ifstream ifs( data_file.c_str() );
  if( !ifs.good() ) {
    QMessageBox::warning( this , "File not found." ,
                          QString( "Couldn't open %1 for reading." ).arg( data_file.c_str() ) );
    return;
  }

#ifdef NOTYET
  cout << "Reading file " << data_file << endl;
#endif
  // first line is assumed to be headers.
  string first_line;
  getline( ifs , first_line );
  vector<string> headers;
  boost::algorithm::trim( first_line );
  if( first_line.empty() ) {
    return; // but it's a duff file
  }
  boost::algorithm::split( headers , first_line , boost::algorithm::is_any_of( " \t" ) ,
                           boost::algorithm::token_compress_on );
#ifdef NOTYET
  cout << "File headers :";
  BOOST_FOREACH( string header , headers ) {
    cout << " " << header;
  }
  cout << endl;
#endif

  while( 1 ) {
    string next_line;
    getline( ifs , next_line );
    if( ifs.eof() || !ifs.good() ) {
      break;
    }
    boost::algorithm::trim( next_line );
    if( next_line.empty() ) {
      continue;
    }

    list<string> splits;
    boost::algorithm::split( splits , next_line , boost::algorithm::is_any_of( " \t" ) ,
                             boost::algorithm::token_compress_on );
    string mol_name = splits.front();
    splits.pop_front();
    pMolRec mol = get_molecule( mol_name );

    if( !mol ) {
#ifdef NOTYET
      cout << data_file << " has data for " << mol_name << " not currently in molecule dataset." << endl;
#endif
      continue;
    }
    int i = 1;
    BOOST_FOREACH( string dp , splits ) {
      float data_val;
      try {
        data_val = lexical_cast<float>( dp );
      } catch( bad_lexical_cast &e ) {
#ifdef NOTYET
        cout << "Bad data value " << dp << " for " << headers[i] << " of molecule " << mol_name << endl;
#endif
        continue;
      }
      if( !mol->add_data( headers[i] , data_val ) ) {
#ifdef NOTYET
        cout << data_file << " has value " << data_val << " for " << headers[i] << " for "
             << mol_name << " which already has datum called that." << endl;
#endif
      }
      ++i;
    }
  }

  mol_table_->rebuild_column_names();

}

// *************************************************************************
void SVDClusRDKit::do_svd_clustering( double tv_alpha , double tv_beta , int start_num_clus ,
                                      int stop_num_clus , int num_clus_step ,
                                      double gamma , double sim_thresh ,
                                      double clus_thresh , bool overlapping_clusters ) {

  if( start_num_clus < 0 || stop_num_clus < 0 ) {
    QMessageBox::warning( this , "Bad cluster number" , "Number of clusters not specified." );
    return;
  }

  QApplication::setOverrideCursor( Qt::WaitCursor );
  for( int dims = start_num_clus ; dims <= stop_num_clus ; dims += num_clus_step ) {

    Chronograph chrono2;

    std::vector<pSVDCluster> u_clusters , v_clusters;
    float u_sil_score , v_sil_score;

    chrono2.start();
    DoSVDCluster( mol_table_->molecules() , tv_alpha , tv_beta , gamma , dims ,
                  clus_thresh , sim_thresh , overlapping_clusters ,
                  u_clusters , u_sil_score , v_clusters , v_sil_score );
    chrono2.stop();

#ifdef NOTYET
    write_clusters( cout , u_clusters );
#endif

    QString label = QString( "U Clusters : Num. Clusters = %1, FPs = %8, Alpha = %2 , Beta = %3 , Gamma = %4 , Clus. Thresh = %5 , Sim. Thresh = %6 , Overlapping = %7" )
        .arg( dims ).arg( tv_alpha ).arg( tv_beta ).arg( gamma ).arg( clus_thresh ).arg( sim_thresh ).arg( overlapping_clusters ).arg( fingerprint_label() );

    ClusterWindow *new_win = new ClusterWindow( u_clusters , overlapping_clusters , mol_draw_del_ , label );
    new_win->connect_selection( this );

    mdi_area_->addSubWindow( new_win );
    new_win->show();
    report_clus_statistics( new_win , u_clusters , u_sil_score , overlapping_clusters , chrono2 );

    if( tv_alpha != tv_beta ) {
      // the v clusters will be different, so show those, too
      QString label = QString( "V Clusters : Num. Clusters = %1, FPs = %8, Alpha = %2 , Beta = %3 , Gamma = %4 , Clus. Thresh = %5 , Sim. Thresh = %6 , Overlapping = %7" )
          .arg( dims ).arg( tv_alpha ).arg( tv_beta ).arg( gamma ).arg( clus_thresh ).arg( clus_thresh ).arg( overlapping_clusters ).arg( fingerprint_label() );
      ClusterWindow *new_win = new ClusterWindow( v_clusters , overlapping_clusters , mol_draw_del_ , label );
      new_win->connect_selection( this );
      mdi_area_->addSubWindow( new_win );
      new_win->show();
      report_clus_statistics( new_win , v_clusters , v_sil_score , overlapping_clusters , chrono2 );
    }

    cout << "Clustering with " << dims << " clusters took " << chrono2.elapsed() << " seconds." << endl;

    mdi_area_->tileSubWindows();

  }

  QApplication::restoreOverrideCursor();

}

// *************************************************************************
void SVDClusRDKit::do_k_means_clustering( int start_num_clus , int stop_num_clus ,
                                          int clus_num_step , int num_iters ) {

  if( start_num_clus < 0 || stop_num_clus < 0 ) {
    QMessageBox::warning( this , "Bad cluster number" , "Number of clusters not specified." );
    return;
  }

  QApplication::setOverrideCursor( Qt::WaitCursor );

  for( int num_clus = start_num_clus ; num_clus <= stop_num_clus ; num_clus += clus_num_step ) {

    cout << "clusters : " << num_clus << " of " << start_num_clus << " to " << stop_num_clus << endl;
    vector<pSVDCluster> clusters;
    float sil_score;

    Chronograph chrono1;
    chrono1.start();
    DoKMeansCluster( mol_table_->molecules() , num_clus , num_iters , clusters , sil_score );
    chrono1.stop();

    QString fp_lab = fingerprint_label();
    QString label( QString( "K-Means Clustering. FPs = %2, Num. clusters = %1." ).arg( clusters.size() )
                   .arg( fingerprint_label() ) );
    ClusterWindow *new_win = new ClusterWindow( clusters , false , mol_draw_del_ , label );
    new_win->connect_selection( this );
    mdi_area_->addSubWindow( new_win );
    new_win->show();
    mdi_area_->tileSubWindows();

    report_clus_statistics( new_win , clusters , sil_score , false , chrono1 );

  }

  QApplication::restoreOverrideCursor();

}

// *************************************************************************
void SVDClusRDKit::do_fuzzy_k_means_clustering( int start_num_clus , int stop_num_clus ,
                                                int clus_num_step , int num_iters ,
                                                float m ) {

  if( start_num_clus < 0 || stop_num_clus < 0 ) {
    QMessageBox::warning( this , "Bad cluster number" , "Number of clusters not specified." );
    return;
  }

  QApplication::setOverrideCursor( Qt::WaitCursor );

  for( int num_clus = start_num_clus ; num_clus <= stop_num_clus ; num_clus += clus_num_step ) {

    vector<pSVDCluster> clusters;
    float sil_score;

    Chronograph chrono1;
    chrono1.start();
    DoFuzzyKMeansCluster( mol_table_->molecules() , num_clus , num_iters , 1.0e-6 ,
                          m , clusters , sil_score );
    chrono1.stop();

    QString fp_lab = fingerprint_label();
    QString label( QString( "Fuzzy K-Means Clustering. FPs = %2, Num. clusters = %1." ).arg( clusters.size() )
                   .arg( fingerprint_label() ) );
    ClusterWindow *new_win = new ClusterWindow( clusters , false , mol_draw_del_ , label );
    new_win->connect_selection( this );
    mdi_area_->addSubWindow( new_win );
    new_win->show();
    mdi_area_->tileSubWindows();

    report_clus_statistics( new_win , clusters , sil_score , true , chrono1 );

  }

  QApplication::restoreOverrideCursor();

}

// *************************************************************************
pMolRec SVDClusRDKit::get_molecule( const string &mol_name ) {

  if( !mol_table_->rowCount() ) {
    return pMolRec();
  }

  static vector<pMolRec>::const_iterator p = mol_table_->molecules().begin();

#ifdef NOTYET
  cout << "looking for " << mol_name << " in " << mol_table_->rowCount() << " mols."
       << " starting at " << (*p)->name() << endl;
#endif

  // do a circular search from where we left off last time, just in case
  // the molecules are in the same order in the data file as in the
  // SMILES file, for example.
  vector<pMolRec>::const_iterator p_start = p;
  do {
    if( mol_name == (*p)->name() ) {
      return *p;
    }
    ++p;
    if( p == mol_table_->molecules().end() ) {
      p = mol_table_->molecules().begin();
    }
  } while( p != p_start );

  return pMolRec();

}

// *************************************************************************
// send summary information to the cluster window for display in the text widget.
void SVDClusRDKit::report_clus_statistics( ClusterWindow *clus_win , const vector<pSVDCluster> &clusters ,
                                           float sil_score , bool overlapping_clusters , Chronograph &chrono ) {

  int num_clustered = 0;

  if( overlapping_clusters ) {
    // calculate number of unique molecules in clusters (molecules can be in more than 1 cluster)
    vector<string> mol_names;
    BOOST_FOREACH( pSVDCluster clus , clusters ) {
      const vector<pSVDClusMem> &mems = clus->cluster_members();
      transform( mems.begin() , mems.end() , back_inserter( mol_names ) ,
                 bind( &SVDClusterMember::name , _1 ) );
    }
    sort( mol_names.begin() , mol_names.end() );
    mol_names.erase( unique( mol_names.begin() , mol_names.end() ) , mol_names.end() );
    num_clustered = mol_names.size();
  } else {
    accumulate( clusters.begin() , clusters.end() , 0 ,
                bind( plus<int>() , _1 , bind( &SVDCluster::size , _2 ) ) );
  }

  QString label = QString( "Clustered %1 molecules into %2 clusters.  " ).arg( mol_table_->count_fingerprints() ).arg( clusters.size() );
  // The number of molecules can be more or less than number of molecules in clusters
  // if they are overlapping clusters, where a molecule can be in more than one
  // cluster. Also, the number of molecules in non-overlapping clusters might be
  // fewer than the total number for SVD clustering as some may not have a
  // contribution to any cluster that is above the threshold.
  if( num_clustered < static_cast<int>( mol_table_->rowCount()) ) {
    label += QString( "%2 of them ended up in clusters.\n" ).arg( num_clustered );
  } else if( num_clustered > static_cast<int>( mol_table_->rowCount() ) ) {
    label += QString( "Each molecule on average in %1 clusters.\n" )
        .arg( float( num_clustered ) / static_cast<float>( mol_table_->rowCount() ) );
  } else {
    label += QString( "\n" );
  }

  if( overlapping_clusters ) {
    label += QString( "Fuzzy silhouette score = %1.\n" ).arg( sil_score );
  } else {
    label += QString( "Crisp silhouette score = %1.\n" ).arg( sil_score );
  }
  label += QString( "Time to cluster = %1s.\n" ).arg( chrono.elapsed() );
  clus_win->slot_text_to_show( label );

}

// *************************************************************************
void SVDClusRDKit::build_fingerprints() {

  switch( fp_type_ ) {
  case CIRCULAR_FPS : build_circular_fingerprints(); break;
  case LINEAR_FPS : build_linear_fingerprints(); break;
  case USER_FPS : slot_user_fps(); break;
  case NO_FPS : break; // to stop the compiler whinging
  }

}

// *************************************************************************
void SVDClusRDKit::build_circular_fingerprints() {

#ifdef NOTYET
  cout << "Building circular fingerprints" << endl;
#endif
  BOOST_FOREACH( pMolRec mol , mol_table_->molecules() ) {
    RDKit::ROMol *newmol = RDKit::SmilesToMol( mol->smiles() );
    mol->set_fingerprint( pRD_FP( RDKit::MorganFingerprints::getFingerprintAsBitVect( *newmol , 3 , 2048 ) ) );
    delete newmol;
  }

  fp_type_ = CIRCULAR_FPS;

}

// *************************************************************************
void SVDClusRDKit::build_linear_fingerprints() {

#ifdef NOTYET
  cout << "building linear fingerprints" << endl;
#endif
  BOOST_FOREACH( pMolRec mol , mol_table_->molecules() ) {
    RDKit::ROMol *newmol = RDKit::SmilesToMol( mol->smiles() );
    mol->set_fingerprint( pRD_FP( RDKit::RDKFingerprintMol( *newmol ) ) );
    delete newmol;
  }

  fp_type_ = LINEAR_FPS;

}

// *************************************************************************
// read a file of user-generated fingerprints, assumed to be a name followed by
// an stream of 0 and 1, possibly with whitespace inbetween. The name can't contain
// whitespace.
void SVDClusRDKit::read_user_fingerprints( const QString &fp_file ) {

#ifdef NOTYET
  cout << "reading fingerprints from " << fp_file.toStdString() << endl;
#endif
  ifstream ifs( fp_file.toLocal8Bit().data() );
  while( 1 ) {
    string next_line;
    getline( ifs , next_line );
    if( ifs.eof() || !ifs.good() ) {
      break;
    }
    boost::trim( next_line );
    if( next_line.empty() ) {
      continue;
    }
    string::iterator i = next_line.begin();
    string mol_name;
    while( !isspace( *i ) && i != next_line.end() ) {
      mol_name += *i;
      ++i;
    }

    // find the right molecule
    pMolRec mol = get_molecule( mol_name );
    if( mol ) {
      vector<char> fp_bits;
      for( ; i != next_line.end() ; ++i ) {
        if( *i == '0' || *i == '1' ) {
          fp_bits.push_back( *i );
        } else if( isspace( *i ) ) {
          continue;
        } else {
          break;
        }
      }

#ifdef NOTYET
      cout << mol_name << " : " << fp_bits.size() << endl;
#endif
      pRD_FP new_fp( new ExplicitBitVect( fp_bits.size() ) );
      for( int i = 0 , is = fp_bits.size() ; i < is ; ++i ) {
        if( '1' == fp_bits[i] ) {
          new_fp->setBit( i );
        }
      }
#ifdef NOTYET
      cout << "Num on bits : " << new_fp->getNumOnBits() << " and off " << new_fp->getNumOffBits() << endl;
      cout << "Count of input bits : " << count( fp_bits.begin() , fp_bits.end() , '1' ) << endl;
#endif
      mol->set_fingerprint( new_fp );
    } else {
      cout << "Error : molecule " << mol_name << " not found in dataset. Skipping." << endl;
    }

  }

  fp_type_ = USER_FPS;

}

// *************************************************************************
// make sure there are some fingerprints for clustering
bool SVDClusRDKit::check_fingerprints_for_clustering() {

  if( NO_FPS == fp_type_ ) {
    QMessageBox::warning( this , "No fingerprints" , "No fingerprints yet. Can't continue." );
    return false;
  }

  int num_fps = mol_table_->count_fingerprints();
#ifdef NOTYET
  cout << "num fps = " << num_fps << endl;
#endif
  if( !num_fps ) {
    QMessageBox::warning( this , "No fingerprints" , "No fingerprints yet. Can't continue." );
    return false;
  }
  int num_mols = mol_table_->rowCount();
  if( num_fps < num_mols &&
      QMessageBox::Ok != QMessageBox::question( this , "Missing fingerprints" ,
                                                QString( "Only %1 of the %2 molecules have fingerprints. Proceed anyway?" ).arg( num_fps ).arg( num_mols) ) ) {
    return false;
  }

  return true;

}

// *************************************************************************
QString SVDClusRDKit::fingerprint_label() const {

  switch( fp_type_ ) {
  case CIRCULAR_FPS : return QString( "Circular" );
  case LINEAR_FPS : return QString( "Linear" );
  case USER_FPS : return QString( "User" );
  default : return QString( "Unknown" ); // to appease the compiler as much as anything.
  }

}

// *************************************************************************
void SVDClusRDKit::disconnect_cluster_window_selections() {

  QList<QMdiSubWindow *> windows = mdi_area_->subWindowList();

  for( int i = 0 , is = windows.size() ; i < is ; ++i ) {
    QWidget *child = windows.at( i )->widget();
    ClusterWindow *cw = qobject_cast<ClusterWindow *>( child );
    if( cw ) {
      cw->disconnect_selection( this );
    }
  }

}

// *************************************************************************
void SVDClusRDKit::connect_cluster_window_selections() {

  QList<QMdiSubWindow *> windows = mdi_area_->subWindowList();

  for( int i = 0 , is = windows.size() ; i < is ; ++i ) {
    QWidget *child = windows.at( i )->widget();
    ClusterWindow *cw = qobject_cast<ClusterWindow *>( child );
    if( cw ) {
      cw->connect_selection( this );
    }
  }

}

// *************************************************************************
void SVDClusRDKit::select_molecules_in_mol_table( const vector<string> &sel_mol_names ) {

  disconnect( mol_table_view_->selectionModel() , &QItemSelectionModel::selectionChanged ,
              this , &SVDClusRDKit::slot_mol_selection_changed );

  QItemSelection qis;
  for( int i = 0 , is = mol_table_->rowCount() ; i < is ; ++i ) {
    string row_name = string( mol_table_->headerData( i , Qt::Vertical ).toString().toLocal8Bit().data() );
    if( sel_mol_names.end() != std::find( sel_mol_names.begin() , sel_mol_names.end() , row_name ) ) {
      QModelIndex ind = mol_table_->index( i , 0 );
      qis.merge( QItemSelection( ind , ind ) , QItemSelectionModel::Select );
    }
  }

  mol_table_view_->selectionModel()->clearSelection();
  mol_table_view_->selectionModel()->select( qis , QItemSelectionModel::Select );
  if( !qis.indexes().isEmpty() ) {
    mol_table_view_->scrollTo( qis.indexes().front() );
  }

  connect( mol_table_view_->selectionModel() , &QItemSelectionModel::selectionChanged ,
           this , &SVDClusRDKit::slot_mol_selection_changed );

}

// *************************************************************************
void SVDClusRDKit::slot_read_smiles() {

  QString smi_file = QFileDialog::getOpenFileName( this , "SMILES file" , last_dir_ , "SMILES (*.smi)" );
  if( smi_file.isEmpty() ) {
    return;
  }

  read_smiles_file( smi_file.toLocal8Bit().data() );

  QFileInfo fi( smi_file );
  last_dir_ = fi.absolutePath();

}

// *************************************************************************
void SVDClusRDKit::slot_read_data_file() {

  QString data_file = QFileDialog::getOpenFileName( this , "Data file" , last_dir_ , "Any file (*.*)" );
  if( data_file.isEmpty() ) {
    return;
  }

  read_data_file( data_file.toLocal8Bit().data() );

  QFileInfo fi( data_file );
  last_dir_ = fi.absolutePath();

}

// *************************************************************************
void SVDClusRDKit::slot_build_svd_clusters() {

  if( !check_fingerprints_for_clustering() ) {
    return;
  }

  if( !svd_clusters_dialog_ ) {
    svd_clusters_dialog_ = new SVDClustersDialog( settings_ , this );
  }

  if( QDialog::Accepted != svd_clusters_dialog_->exec() ) {
    return;
  }

  double tv_alpha , tv_beta , gamma , sim_thresh , clus_thresh;
  int start_num_clus , stop_num_clus , num_clus_step;
  bool overlapping_clusters;
  svd_clusters_dialog_->get_settings( start_num_clus , stop_num_clus , num_clus_step ,
                                      tv_alpha , tv_beta , gamma ,
                                      sim_thresh , clus_thresh , overlapping_clusters );

  do_svd_clustering( tv_alpha , tv_beta , start_num_clus , stop_num_clus , num_clus_step ,
                     gamma , sim_thresh , clus_thresh , overlapping_clusters );

}

// *************************************************************************
void SVDClusRDKit::slot_build_k_means_clusters() {

  if( !check_fingerprints_for_clustering() ) {
    return;
  }

  if( !kmeans_clusters_dialog_ ) {
    kmeans_clusters_dialog_ = new KMeansClustersDialog( settings_ , this );
  }

  if( QDialog::Accepted != kmeans_clusters_dialog_->exec() ) {
    return;
  }

  int start_num_clus , stop_num_clus , num_clus_step , num_iters;
  kmeans_clusters_dialog_->get_settings( start_num_clus , stop_num_clus , num_clus_step , num_iters );
  do_k_means_clustering( start_num_clus , stop_num_clus , num_clus_step , num_iters );

}

// *************************************************************************
void SVDClusRDKit::slot_build_fuzzy_k_means_clusters() {

  if( !check_fingerprints_for_clustering() ) {
    return;
  }

  if( !fuzzy_kmeans_clusters_dialog_ ) {
    fuzzy_kmeans_clusters_dialog_ = new FuzzyKMeansClustersDialog( settings_ , this );
  }

  if( QDialog::Accepted != fuzzy_kmeans_clusters_dialog_->exec() ) {
    return;
  }

  int start_num_clus , stop_num_clus , num_clus_step , num_iters;
  float m;
  fuzzy_kmeans_clusters_dialog_->get_settings( start_num_clus , stop_num_clus , num_clus_step , num_iters , m );
  do_fuzzy_k_means_clustering( start_num_clus , stop_num_clus , num_clus_step , num_iters , m );

}

// *************************************************************************
void SVDClusRDKit::slot_quit() {

  exit( 0 );

}

// *************************************************************************
void SVDClusRDKit::slot_update_windows_menu() {

  windows_menu_->clear();
  windows_menu_->addAction( tile_windows_ );
  windows_menu_->addAction( cascade_windows_ );
  windows_menu_->addSeparator();

  QList<QMdiSubWindow *> windows = mdi_area_->subWindowList();
  separator_act_->setVisible( !windows.isEmpty() );

  for( int i = 0 , is = windows.size() ; i < is ; ++i ) {
    QWidget *child = windows.at( i )->widget();
    QAction *action = windows_menu_->addAction( child->windowTitle() );
    connect( action , SIGNAL( triggered() ) , windows_mapper_ , SLOT( map() ) );
    windows_mapper_->setMapping( action , windows.at( i ) );
  }

}

// ***************************************************************************
void SVDClusRDKit::slot_set_active_sub_window( QWidget *window ) {

  if( !window ) {
    return;
  }

  mdi_area_->setActiveSubWindow( qobject_cast<QMdiSubWindow *>( window ) );

}

// ***************************************************************************
void SVDClusRDKit::slot_circular_fps() {

  if( fp_type_ != CIRCULAR_FPS ) {
    fp_type_ = CIRCULAR_FPS;
    build_fingerprints();
  }

}

// ***************************************************************************
void SVDClusRDKit::slot_linear_fps() {

  if( fp_type_ != LINEAR_FPS ) {
    fp_type_ = LINEAR_FPS;
    build_fingerprints();
  }

}

// ***************************************************************************
void SVDClusRDKit::slot_user_fps() {

  fp_type_ = USER_FPS;
  QString fp_file = QFileDialog::getOpenFileName( this , "Fingerprint file" , last_dir_ , "Any file (*.*)" );
  if( fp_file.isEmpty() ) {
    return;
  }

  read_user_fingerprints( fp_file );

  QFileInfo fi( fp_file );
  last_dir_ = fi.absolutePath();

}

// ***************************************************************************
void SVDClusRDKit::slot_show_help_about() {

  QString msg = QString( "SVDClus\n\n\nCopyright (C) 2014 AstraZeneca\n\nBuilt using RDKit version %3 and Qt version %1.\nRunning with Qt version %2\n\n" ).arg( QT_VERSION_STR ).arg( qVersion() ).arg( RDKit::rdkitVersion );

  msg += settings_->usage_text().c_str();

  QMessageBox::information( this , "About SVDClus" , msg ,
                            QMessageBox::Ok | QMessageBox::Default );

}

// ***************************************************************************
void SVDClusRDKit::slot_show_help_index() {

  help_viewer_->show();

}

// ***************************************************************************
void SVDClusRDKit::slot_mol_selection_changed( const QItemSelection &sel_items ,
                                               const QItemSelection &desel_items ) {

#ifdef NOTYET
  cout << "slot_mol_selection_changed" << endl;
#endif

  // disconnect the cluster windows so we don't go round in infinite circles
  disconnect_cluster_window_selections();

  // get the names of all the molecules selected in the mol_table_view_
  QModelIndexList sel_inds = mol_table_view_->selectionModel()->selectedIndexes();
  vector<string> sel_mol_names;
  for( int i = 0 , is = sel_inds.count() ; i < is ; ++i ) {
    if( mol_table_->data( sel_inds[i] , Qt::DisplayRole ).canConvert<pMolRec>() ) {
      pMolRec mol = mol_table_->data( sel_inds[i] , Qt::DisplayRole ).value<pMolRec>();
      sel_mol_names.push_back( mol->name().c_str() );
#ifdef NOTYET
      cout << i << " : " << sel_mol_names.back() << endl;
#endif
    }
  }

  // pass the new selections to the ClusterViews
  QList<QMdiSubWindow *> windows = mdi_area_->subWindowList();
  for( int i = 0 , is = windows.size() ; i < is ; ++i ) {
    QWidget *child = windows.at( i )->widget();
    ClusterWindow *cw = qobject_cast<ClusterWindow *>( child );
    if( cw ) {
      cw->select_molecules( sel_mol_names );
    }
  }

  // and reconnect the cluster windows
  connect_cluster_window_selections();

}
