/*
 * ntuple_pixelclusters.h
 *
 *  Created on: 1 November 2021
 *      Author: Alexandre De Moor
 */

#ifndef DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_PIXEL_CLUSTERS_
#define DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_PIXEL_CLUSTERS_

#include "ntuple_content.h"

// For pixel clusters and topology
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"

// Geometry
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "Geometry/CommonDetUnit/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

class ntuple_pixelclusters: public ntuple_content{
public:

    ntuple_pixelclusters(double jetR = 0.4);
    ~ntuple_pixelclusters();

    void getInput(const edm::ParameterSet& iConfig);
    void initBranches(TTree* );
    void readEvent(const edm::Event& iEvent);
    void readSetup(const edm::EventSetup& iSetup);
    void checkEventSetup(const edm::EventSetup & iSetup);

    //use either of these functions

    bool fillBranches(const pat::Jet &, const size_t& jetidx, const  edm::View<pat::Jet> * coll=0);


    void setCandidatesToken(const edm::EDGetTokenT<edm::View<pat::PackedCandidate> > & t){
        CandidateToken=t;
     }

private:

    // seed candidates
    static constexpr size_t max_gtracks=4;
    
    unsigned int n_layers;
    
    float r004[max_gtracks];
    float r006[max_gtracks];
    float r008[max_gtracks];
    float r010[max_gtracks];
    float r016[max_gtracks];
    float rvar[max_gtracks];
    float rvwt[max_gtracks];

    edm::EDGetTokenT<edm::View<reco::Jet> > m_jets;
    edm::EDGetTokenT<reco::VertexCollection> m_vertices;
    edm::EDGetTokenT<edmNew::DetSetVector<SiPixelCluster> > m_pixelhit; //To header !
    edm::Handle<edmNew::DetSetVector<SiPixelCluster> > collectionHandle; //To header !
    edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> m_geomToken; //To header !
    edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> m_topoToken; //To header !

    //tokens to be defined from main analyzer                                                                                                                                                              
    edm::EDGetTokenT<edm::View<pat::PackedCandidate> > CandidateToken;
    // Open Geometry                                                                                                                                                                                       
    const TrackerGeometry* theTracker;
    // Retrieve tracker topology from geometry                                                                                                                                                             
    const TrackerTopology* tTopo;
};

#endif /* DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_PIXEL_CLUSTERS_ */
