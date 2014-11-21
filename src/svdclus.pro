SOURCES += SVDClustersDialog.cc SVDClusRDKit.cc \
  ClustersTableModel.cc RDKitMolDrawDelegate.cc SVDClusSettings.cc \
  DoSVDCluster.cc RDKitMolToQPainter.cc SVDCluster.cc \
  GetGaussFilteredSimMatrix.cc \
  GetRDKitSims.cc svdclus_main.cc \
    ClusterWindow.cc \
    MoleculeRec.cc \
    ColourClusterMolsDialog.cc \
    DoKMeansCluster.cc \
    KMeansClustersDialog.cc \
    crisp_silhouette_score.cc \
    fuzzy_silhouette_score.cc \
    BuildClustersDialog.cc \
    MoleculeTableModel.cc \
    QTHelpViewer.cc \
    DoFuzzyKMeansCluster.cc \
    FuzzyKMeansClustersDialog.cc \
    ClustersTableView.cc

HEADERS += SVDClustersDialog.H SVDClusSettings.H \
ClustersTableModel.H RDKitMolDrawDelegate.H SVDCluster.H \
ClustersTableView.H SVDClusRDKitDefs.H SVDClusterMember.H \
GetRDKitSims.H SVDClusRDKit.H \
    MoleculeRec.H \
    ClusterWindow.H \
    ColourClusterMolsDialog.H \
    KMeansClustersDialog.H \
    BuildClustersDialog.H \
    MoleculeTableModel.H \
    MoleculeTableView.H \
    QTHelpViewer.H \
    FuzzyKMeansClustersDialog.H

TARGET = svdclus

CONFIG += qt debug_and_release
CONFIG(debug, debug|release) {
     TARGET = $$join(TARGET,,,_debug)
}

QT += widgets

INCLUDEPATH += ${RDBASE}/Code ${BOOST_ROOT}/include ${SVDLIBC_HOME}

RD_STATIC_LIBS = -lSmilesParse_static -lFingerprints_static \
-lSubstructMatch_static \
-lDepictor_static -lSubgraphs_static -lGraphMol_static \
-lRDGeometryLib_static -lDataStructs_static -lRDGeneral_static

LIBS += -L${RDBASE}/lib $${RD_STATIC_LIBS}

LIBS += -L${SVDLIBC_HOME} -lsvd

LIBS += ${BOOST_ROOT}/lib/libboost_program_options.a ${BOOST_ROOT}/lib/libboost_regex.a




