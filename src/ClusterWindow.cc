//
// file ClusterWindow.cc
//
//  Copyright (C) 2014 AstraZeneca, David Cosgrove
//
//   @@ All Rights Reserved @@
//  This file is part of SVDClus.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the source tree.
//

#include "stddefs.H"
#include "ClustersTableModel.H"
#include "ClustersTableView.H"
#include "ClusterWindow.H"
#include "ColourClusterMolsDialog.H"
#include "RDKitMolDrawDelegate.H"
#include "SVDCluster.H"
#include "SVDClusterMember.H"
#include "SVDClusRDKit.H"

#include <QAction>
#include <QContextMenuEvent>
#include <QFileDialog>
#include <QLayout>
#include <QMenu>
#include <QMessageBox>
#include <QItemSelectionModel>
#include <QSplitter>
#include <QTextEdit>

#include <fstream>

using namespace std;

QString ClusterWindow::last_dir_ = QString( "." );

// *************************************************************************
void write_clusters( ostream &os , const vector<pSVDCluster> &clusters ) {

  os << "Molecule name : Cluster size : Cluster Members" << endl;
  for( int i = 0 , is = clusters.size() ; i < is ; ++i ) {
    const vector<pSVDClusMem> &this_clus = clusters[i]->cluster_members();
    if( this_clus.empty() ) {
      continue;
    }
    os << this_clus.front()->name() << " : "
        << this_clus.size() << "(" << this_clus.size() << ") : ";
    for( int j = 0 , js = this_clus.size() ; j < js ; ++j ) {
      os << this_clus[j]->name() << " ";
    }
    os << endl;
  }

}

// ****************************************************************************
void write_eigenvectors( ostream &os , const vector<pSVDCluster> &clusters ) {

  for( int i = 0 , is = clusters.size() ; i < is ; ++i ) {
    os << clusters[i]->eigen_value();
    const vector<pSVDClusMem> &this_clus = clusters[i]->cluster_members();
    if( this_clus.empty() ) {
      continue;
    }
    for( int j = 0 , js = this_clus.size() ; j < js ; ++j ) {
      os << " " << this_clus[j]->cont();
    }
    os << endl;
  }

}


// ****************************************************************************
ClusterWindow::ClusterWindow( const vector<pSVDCluster> &clus ,
                              bool overlapping_clusters ,
                              RDKitMolDrawDelegate *mdd ,
                              const QString &label ,
                              QWidget *parent ,
                              Qt::WindowFlags f ) :
  QWidget( parent , f ) , clusters_( clus ) {

  build_widget( mdd , label );

  clusters_model_->set_clusters( clusters_ , overlapping_clusters );
  clusters_view_->resizeColumnsToContents();
  clusters_view_->resizeRowsToContents();

  build_actions();

}

// *************************************************************************
void ClusterWindow::connect_selection( SVDClusRDKit *parent_wid ) {

#ifdef NOTYET
  cout << "connecting " << windowTitle().toStdString() << endl;
#endif
  connect( clusters_view_->selectionModel() , &QItemSelectionModel::selectionChanged ,
           parent_wid , &SVDClusRDKit::slot_cluster_selection_changed );

}

// *************************************************************************
void ClusterWindow::disconnect_selection( SVDClusRDKit *parent_wid ) {

#ifdef NOTYET
  cout << "disconnecting " << windowTitle().toStdString() << endl;
#endif
  disconnect( clusters_view_->selectionModel() , &QItemSelectionModel::selectionChanged ,
              parent_wid , &SVDClusRDKit::slot_cluster_selection_changed );

}

// *************************************************************************
void ClusterWindow::select_molecules(const vector<string> &sel_mol_names) {

#ifdef NOTYET
  cout << "ClusterWindow::select_molecules : ";
  copy( sel_mol_names.begin() , sel_mol_names.end() , stringOut );
  cout << " :: " << windowTitle().toStdString() << endl;
#endif

  QItemSelection qis;
  for( int i = 0 , is = clusters_model_->rowCount() ; i < is ; ++i ) {
    for( int j = 0 , js = clusters_model_->columnCount() ; j < js ; ++j ) {
      QModelIndex ind = clusters_model_->index( i , j );
      QVariant var = clusters_model_->data( ind , Qt::DisplayRole );
      if( var.isValid() && var.canConvert<pSVDClusMem>() ) {
        if( sel_mol_names.end() != std::find( sel_mol_names.begin() , sel_mol_names.end() ,
                                              var.value<pSVDClusMem>()->name() ) ) {
          qis.merge( QItemSelection( ind , ind ) , QItemSelectionModel::Select );
        }
      }
    }
  }

  clusters_view_->selectionModel()->clearSelection();
  clusters_view_->selectionModel()->select( qis , QItemSelectionModel::Select );
  if( !qis.indexes().isEmpty() ) {
    clusters_view_->scrollTo( qis.indexes().front() );
  }

#ifdef NOTYET
  cout << "leaving ClusterWindow::select_molecules : " << windowTitle().toStdString() << endl;
#endif

}

// *************************************************************************
bool ClusterWindow::is_this_my_model( const ClustersTableModel *mod ) const {

  return( mod == clusters_model_ );

}

// *************************************************************************
vector<string> ClusterWindow::selected_molecules() const {

  vector<string> sel_mol_names;
  QModelIndexList sel_inds = clusters_view_->selectionModel()->selectedIndexes();
  for( int i = 0 , is = sel_inds.size() ; i < is ; ++i ) {
    QVariant var = clusters_model_->data( sel_inds[i] , Qt::DisplayRole );
    if( var.isValid() && var.canConvert<pSVDClusMem>() ) {
      sel_mol_names.push_back( var.value<pSVDClusMem>()->name() );
    }
  }

  return sel_mol_names;

}

// *************************************************************************
void ClusterWindow::slot_text_to_show( QString text ) {

  text_window_->append( text );

}

// ****************************************************************************
void ClusterWindow::build_widget( RDKitMolDrawDelegate *mdd ,
                                  const QString &label ) {

  clusters_model_ = new ClustersTableModel;
  clusters_view_ = new ClustersTableView;
  clusters_view_->setModel( clusters_model_ );
  clusters_view_->setItemDelegate( mdd );

  QVBoxLayout *vbox = new QVBoxLayout;
  QSplitter *splitter = new QSplitter( Qt::Vertical );
  splitter->addWidget( clusters_view_ );

  text_window_ = new QTextEdit;
  text_window_->setReadOnly( true );
  splitter->addWidget( text_window_ );
  splitter->setStretchFactor( 0 , 1 );

  vbox->addWidget( splitter );
  setLayout( vbox );

  // set sensible areas to view and text window
  QList<int> win_sizes;
  win_sizes.push_back( 50 );
  win_sizes.push_back( 1 );
  splitter->setSizes( win_sizes );

  clusters_model_->set_cluster_colouring( false );

  setWindowTitle( label );

  connect( clusters_model_ , SIGNAL( text_for_output( QString ) ) ,
           this , SLOT( slot_text_to_show( QString) ) );

}

// ****************************************************************************
void ClusterWindow::build_actions() {

  write_clusters_ = new QAction( "Write clusters" , this );
  connect( write_clusters_ , SIGNAL( triggered() ) ,
           this , SLOT( slot_write_clusters() ) );

  write_eigen_vecs_ = new QAction( "Write eigenvectors" , this );
  connect( write_eigen_vecs_ , SIGNAL( triggered() ) ,
           this , SLOT( slot_write_eigenvectors() ) );

  colour_mols_ = new QAction( "Colour Mols" , this );
  connect( colour_mols_ , SIGNAL( triggered() ) ,
           this , SLOT( slot_colour_mols() ) );

  dont_colour_mols_ = new QAction( "Don't Colour Mols" , this );
  connect( dont_colour_mols_ , SIGNAL( triggered() ) ,
           this , SLOT( slot_dont_colour_mols() ) );

  heat_map_ = new QAction( "Heat Map" , this );
  connect( heat_map_ , SIGNAL( triggered() ) ,
           this , SLOT( slot_heat_map() ) );
  heat_map_->setCheckable( true );
  heat_map_->setChecked( false );

  sort_by_eigen_val_ = new QAction( tr( "Sort By Eigenvalue" ) , this );
  connect( sort_by_eigen_val_ , SIGNAL( triggered() ) ,
           this , SLOT( slot_sort_clusters_by_eigval() ) );

  sort_by_size_ = new QAction( tr( "Sort By Size" ) , this );
  connect( sort_by_size_ , SIGNAL( triggered() ) ,
           this , SLOT( slot_sort_clusters_by_size() ) );

}

// ****************************************************************************
void ClusterWindow::contextMenuEvent( QContextMenuEvent *e ) {

  QMenu menu( this );
  menu.addAction( write_clusters_ );
  menu.addAction( write_eigen_vecs_ );
  menu.addAction( colour_mols_ );
  menu.addAction( dont_colour_mols_ );
  menu.addAction( heat_map_ );
  menu.addAction( sort_by_eigen_val_ );
  menu.addAction( sort_by_size_ );

  menu.exec( e->globalPos() );

}

// *************************************************************************
void ClusterWindow::slot_write_clusters() {

  QString clus_file = QFileDialog::getSaveFileName( this , "Clusters File" ,
                                                    last_dir_ , "samples (*.samples);;All (*)" );
  if( clus_file.isEmpty() ) {
    return;
  }

  ofstream ofs( clus_file.toLocal8Bit().data() );
  if( !ofs.good() ) {
    return;
  }

  write_clusters( ofs , clusters_ );

  QFileInfo fi( clus_file );
  last_dir_ = fi.absolutePath();

}

// *************************************************************************
void ClusterWindow::slot_write_eigenvectors() {

  QString ev_file = QFileDialog::getSaveFileName( this , "Eigenvectors File" ,
                                                  last_dir_ , "All (*)" );
  if( ev_file.isEmpty() ) {
    return;
  }

  ofstream ofs( ev_file.toLocal8Bit().data() );
  if( !ofs.good() ) {
    return;
  }

  write_eigenvectors( ofs , clusters_ );

  QFileInfo fi( ev_file );
  last_dir_ = fi.absolutePath();

}

// *************************************************************************
void ClusterWindow::slot_colour_mols() {

  vector<string> names = clusters_model_->get_data_names();
  if( names.empty() ) {
    QMessageBox::warning( this , "No Data" , "No data, so nothing to colour by." );
    return;
  }

  ColourClusterMolsDialog ccmd( names , 300.0 , this );
  if( QDialog::Rejected == ccmd.exec() ) {
    return;
  }

  int name_num;
  float act_cutoff;
  COLOUR_DATA_SENSE ds;
  ccmd.get_settings( name_num , act_cutoff , ds );
  clusters_model_->set_background_data_details( names[name_num] , act_cutoff , ds );
  clusters_model_->set_cluster_colouring( true );

  cout << "QCI score " << clusters_model_->calc_qci_score() << endl;
  text_window_->append( QString( "QCI score = %1\n" ).arg( clusters_model_->calc_qci_score() ) );

}

// *************************************************************************
void ClusterWindow::slot_dont_colour_mols() {

  clusters_model_->set_cluster_colouring( false );

}

// *************************************************************************
void ClusterWindow::slot_heat_map() {

  int cw , rh;

  if( heat_map_->isChecked() ) {
    cw = width() / clusters_model_->columnCount();
    cw = cw < 5 ? 5 : cw;
    rh = height() / clusters_model_->rowCount();
    rh = rh < 5 ? 5 : rh;
  } else {
    QSize s = clusters_view_->sizeHint();
    cw = s.width();
    rh = s.height();
  }

  for( int i = 0 , is = clusters_model_->columnCount() ; i < is ; ++i ) {
    clusters_view_->setColumnWidth( i , cw );
  }

  for( int i = 0 , is = clusters_model_->rowCount() ; i < is ; ++i ) {
    clusters_view_->setRowHeight( i , rh );
  }

}

// *************************************************************************
void ClusterWindow::slot_sort_clusters_by_eigval() {

  clusters_model_->sort_by_eigval();

}

// *************************************************************************
void ClusterWindow::slot_sort_clusters_by_size() {

  clusters_model_->sort_by_size();

}
