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

#include "ClustersTableModel.H"
#include "SVDCluster.H"
#include "SVDClusterMember.H"

#include <QBrush>
#include <QColor>
#include <QVariant>

#include <boost/bind.hpp>
#include <boost/foreach.hpp>

#include <algorithm>
#include <set>

using namespace boost;
using namespace std;

// ****************************************************************************
ClustersTableModel::ClustersTableModel( QObject *parent ) :
  QAbstractTableModel( parent ) , num_cols_( 0 ) , overlapping_clusters_( false ) ,
  colour_data_sense_( HIGHER_BETTER ) , cluster_colouring_( false ) ,
  sort_mode_( SORT_BY_EIGVAL ) {

}

// ****************************************************************************
int ClustersTableModel::rowCount( const QModelIndex &parent ) const {

  return clusters_.size();

}

// ****************************************************************************
int ClustersTableModel::columnCount( const QModelIndex &parent ) const {

  num_cols_ = 0; // not strictly necessary as it's initialised in the constructor and table is read-only at present
  BOOST_FOREACH( pSVDCluster clus , clusters_ ) {
    if( clus->size() > num_cols_ ) {
      num_cols_ = clus->size();
    }
  }

  return num_cols_;

}

// ****************************************************************************
QVariant ClustersTableModel::data( const QModelIndex &index , int role ) const {

  if( !index.isValid() ) {
    return QVariant();
  }

  if( index.row() >= static_cast<int>( clusters_.size() ) ||
      index.column() >= clusters_[index.row()]->size() ) {
    return QVariant();
  }

  if( role == Qt::DisplayRole ) {
    QVariant ret_val;
    ret_val.setValue( clusters_[index.row()]->member( index.column() ) );
    return ret_val;
  } else if( role == Qt::ToolTipRole ) {
    QString ret = build_tooltip( index );
    return QVariant( ret );
  } else if( role == Qt::BackgroundRole ) {
    return background_brush( index );
  }

  return QVariant();

}

// ****************************************************************************
QVariant ClustersTableModel::headerData( int section , Qt::Orientation orientation ,
                                         int role ) const {

  if( role != Qt::DisplayRole ) {
    return QVariant();
  }

  if( orientation == Qt::Horizontal ) {
    return QVariant( section + 1 );
  } else {
    // get the eigenvalue of the cluster for this row and add it to the row number
    // it might be an empty row, though, if the contents of the model have changed to a smaller number
    // of clusters
    if( section >= static_cast<int>( clusters_.size() ) || !clusters_[section] ) {
      return QVariant( section + 1 );
    }
    QString ret_val = QString( "%1\nEig. Val. : %2\nClus. Size: %3" )
        .arg( section + 1 )
        .arg( QString::number( clusters_[section]->eigen_value() ) )
        .arg( QString::number( clusters_[section]->size() ) );
    return ret_val;
  }

}

// ****************************************************************************
// this sort only works on the cluster size or eigen value, both in descending order
void ClustersTableModel::sort( int column , Qt::SortOrder order ) {

  if( SORT_BY_SIZE == sort_mode_ ) {
    std::sort( clusters_.begin() , clusters_.end() ,
               bind( greater<int>() ,
                     bind( &SVDCluster::size , _1 ) ,
                     bind( &SVDCluster::size , _2 ) ) );
  } else if( SORT_BY_EIGVAL == sort_mode_ ) {
    std::sort( clusters_.begin() , clusters_.end() ,
               bind( greater<double>() ,
                     bind( &SVDCluster::eigen_value , _1 ) ,
                     bind( &SVDCluster::eigen_value , _2 ) ) );
  }

  emit layoutChanged();

}

// ****************************************************************************
// overlapping denotes whether molecules can be in more than 1 cluster. True if they can.
void ClustersTableModel::set_clusters( std::vector<pSVDCluster> new_clus ,
                                       bool overlapping ) {

  beginRemoveRows( QModelIndex() , 0 , rowCount() );
  beginRemoveColumns( QModelIndex() , 0 , columnCount() );
  clusters_.clear();
  endRemoveColumns();
  endRemoveRows();

  num_cols_ = -1;
  BOOST_FOREACH( pSVDCluster clus , new_clus ) {
    if( clus->size() > num_cols_ ) {
      num_cols_ = clus->size();
    }
  }

  beginInsertRows( QModelIndex() , 0 , static_cast<int>( new_clus.size() ) - 1 );
  beginInsertColumns( QModelIndex() , 0 , num_cols_  - 1 );
  clusters_ = new_clus;
  endInsertColumns();
  endInsertRows();

  overlapping_clusters_ = overlapping;

}

// ****************************************************************************
void ClustersTableModel::set_background_data_details( const string &data_name ,
                                                      float data_cutoff ,
                                                      COLOUR_DATA_SENSE ds ) {

  beginResetModel();

  background_data_name_ = data_name;
  background_data_cutoff_ = data_cutoff;
  colour_data_sense_ = ds;

  background_data_max_ = -numeric_limits<float>::max();
  background_data_min_ = numeric_limits<float>::max();
  BOOST_FOREACH( pSVDCluster clus , clusters_ ) {
    BOOST_FOREACH( const pSVDClusMem clus_mem , clus->cluster_members() ) {
      float val;
      if( clus_mem->molecule()->data( background_data_name_ , val ) ) {
        if( val > background_data_max_ ) {
          background_data_max_ = val;
        }
        if( val < background_data_min_ ) {
          background_data_min_ = val;
        }
      }
    }
  }

#ifdef NOTYET
  cout << "Min : " << background_data_min_ << " Max : " << background_data_max_
       << "  Cutoff : " << background_data_cutoff_ << endl;
#endif

  endResetModel();

}

// ****************************************************************************
void ClustersTableModel::set_cluster_colouring( bool new_val ) {

  beginResetModel();
  cluster_colouring_ = new_val;
  endResetModel();

}

// ****************************************************************************
// return a vector of all the names for data held by the cluster molecules
vector<string> ClustersTableModel::get_data_names() const {

  set<string> data_names;
  typedef pair<string,float> PAIR_DATA;
  BOOST_FOREACH( pSVDCluster clus , clusters_ ) {
    BOOST_FOREACH( const pSVDClusMem clus_mem , clus->cluster_members() ) {
      const vector<pair<string,float> > &mol_data = clus_mem->molecule()->data();
      BOOST_FOREACH( const PAIR_DATA pd , mol_data ) {
        data_names.insert( pd.first );
      }
    }
  }

  return vector<string>( data_names.begin() , data_names.end() );

}

// ****************************************************************************
// do the qci score, only relevant for non-overlapping clusters with activity data.
// It's a measure of how well the actives and inactives are separated by the clustering.
// Varin, T., Bureau, R., Mueller, C. & Willett, P. J. of Mol. Graph. 28 (2), 187-195 (2009)
float ClustersTableModel::calc_qci_score() {

  if( overlapping_clusters_ || !cluster_colouring_ ) {
    return -1.0;
  }

  // qci = p / ( p + q + r + s )
  // p = number of actives in active clusters (not singletons)
  // q = number of inactives in active clusters
  // r = actives in inactive clusters
  // s = number of singleton actives
  // An active cluster is one where the proportion of actives is greater than
  // the proportion of actives in the whole set.
  vector<int> actives_in_clusters( clusters_.size() , 0 );
  int num_acts = 0 , num_tot = 0 , s = 0;

  for( int i = 0 , is = clusters_.size() ; i < is ; ++i ) {
    for( int j = 0 , js = clusters_[i]->size() ; j < js ; ++j ) {
      pMolRec mol = clusters_[i]->member( j )->molecule();
      float val;
      if( mol->data( background_data_name_ , val ) ) {
        if( LOWER_BETTER == colour_data_sense_ ) {
          if( val < background_data_cutoff_ ) {
            actives_in_clusters[i]++;
          }
        } else {
          if( val > background_data_cutoff_ ) {
            actives_in_clusters[i]++;
          }
        }
      }
    }
    if( 1 == clusters_[i]->size() && actives_in_clusters[i] ) {
      ++s;
    }
    num_acts += actives_in_clusters[i];
    num_tot += clusters_[i]->size();
  }

  float gp_act = float( num_acts ) / float( num_tot );
  int p = 0 , q = 0 , r = 0;
  for( int i = 0 , is = clusters_.size() ; i < is ; ++i ) {
    if( float( actives_in_clusters[i] ) / float( clusters_[i]->size() ) > gp_act ) {
      if( clusters_[i]->size()> 1 ) {
        // it's an active cluster
        p += actives_in_clusters[i];
      }
      q += clusters_[i]->size() - actives_in_clusters[i];
    } else {
      // it's an inactive cluster
      r += actives_in_clusters[i];
    }
  }

  QString label = QString( "%1 molecules were in clusters, of which %2 were active.\n" ).arg( num_tot ).arg( num_acts );
  emit text_for_output( label );

#ifdef NOTYET
  cout << num_tot << " molecules were in clusters, of which " << num_acts << " were active." << endl;
  cout << "p = " << p << "  q = " << q << "  r = " << r << "  s = " << s << endl;
#endif
  return float( p ) / float( p + q + r + s );

}

// ****************************************************************************
// assume that index is valid.
QString ClustersTableModel::build_tooltip( const QModelIndex &index ) const {

  QString ret_val;

  const pSVDClusMem &mem = clusters_[index.row()]->member( index.column() );

  ret_val = QString( "Name : %1" ).arg( mem->name().c_str() );
  const pMolRec &mol = mem->molecule();
  const vector<pair<string,float> > &mol_data = mol->data();
  typedef pair<string,float> DATA_PAIR;
  BOOST_FOREACH( DATA_PAIR dp , mol_data ) {
    ret_val += QString( "\n%1 : %2" ).arg( dp.first.c_str() ).arg( dp.second );
  }

#ifdef NOTYET
  // put colour values to tooltip, for debugging the colour scaling
  float val;
  if( mol->data( background_data_name_ , val ) ) {
    QColor retval;
    if( LOWER_BETTER == colour_data_sense_ ) {
      cout << "lower better" << endl;
      if( val < background_data_cutoff_ ) {
        // build a shade of green
        retval = build_green( val , background_data_min_ );
      } else {
        // or a shade of red
        retval = build_red( val , background_data_max_ );
      }
    } else {
      cout << "higher better" << endl;
      if( val > background_data_cutoff_ ) {
        // build a shade of green
        retval = build_green( val , background_data_max_ );
      } else {
        // or a shade of red
        retval = build_red( val , background_data_min_ );
      }
    }
    ret_val += QString( "\n%1 , %2 , %3" ).arg( retval.red() ).arg( retval.green() ).arg( retval.blue() );
  }
#endif

  return ret_val;
}

// ****************************************************************************
// assume that index is valid.
QBrush ClustersTableModel::background_brush( const QModelIndex &index ) const {

  if( !cluster_colouring_ ) {
    return QBrush( QColor( "White" ) );
  }

  const pSVDClusMem &mem = clusters_[index.row()]->member( index.column() );
  const pMolRec &mol = mem->molecule();

  float val;
  if( mol->data( background_data_name_ , val ) ) {
    if( val == background_data_cutoff_ ) {
      return QBrush( QColor( "White" ) ); // it's on the cutoff, it must be white, though it might look odd
    }
    if( LOWER_BETTER == colour_data_sense_ ) {
      if( val < background_data_cutoff_ ) {
        // build a shade of green
        return QBrush( build_green( val , background_data_min_ ) );
      } else {
        // or a shade of red
        return QBrush( build_red( val , background_data_max_ ) );
      }
    } else {
      if( val > background_data_cutoff_ ) {
        // build a shade of green
        return QBrush( build_green( val , background_data_max_ ) );
      } else {
        // or a shade of red
        return QBrush( build_red( val , background_data_min_ ) );
      }
    }
  } else {
    return QBrush( QColor( "Yellow" ) ); // to show there's no data
  }

  return QBrush( QColor( "White" ) );

}

// ****************************************************************************
QColor ClustersTableModel::build_green( float val , float bottom_val ) const {

  float red;
  if( background_data_cutoff_ != bottom_val ) {
    red = 255.0 * ( val - bottom_val ) / ( background_data_cutoff_ - bottom_val );
  } else {
    if( val == bottom_val ) {
      red = 0.0; // 0 for numerator wins
    } else {
      red = 255.0; // 0 for denominator wins
    }
  }
  red = red < 0.0 ? 0.0 : red;
  red = red > 255.0 ? 255.0 : red;
  float blue = red;
  float green = 255.0;

#ifdef NOTYET
  cout << "building green for " << val << " : " << red << " , " << green << " , " << blue
       << " data cutoff = " << background_data_cutoff_ << " and bottom val = " << bottom_val << endl;
#endif

  return QColor( int( red ) , int( green ) , int( blue ) );

}

// ****************************************************************************
QColor ClustersTableModel::build_red( float val , float top_val  ) const {

  float green;
  if( top_val != background_data_cutoff_ ) {
    green = 255.0 - 255.0 * ( val - background_data_cutoff_ ) / ( top_val - background_data_cutoff_ );
  } else {
    if( val == background_data_cutoff_ ) {
      green = 0.0; // 0 for numerator wins
    } else {
      green = 255.0; // 0 for denominator wins
    }
  }
  float blue = green;
  float red = 255.0;

#ifdef NOTYET
  cout << "building red for " << val << " : " << red << " , " << green << " , " << blue
       << " data cutoff = " << background_data_cutoff_ << " and top val = " << top_val << endl;
#endif

  return QColor( int( red ) , int( green ) , int( blue ) );

}

