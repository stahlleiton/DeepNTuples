/*
 * ntuple_pixelclusters.h
 *
 *  Created on: 1 November 2021
 *      Author: Alexandre De Moor
 */

#ifndef DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_PIXEL_CLUSTERS_
#define DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_PIXEL_CLUSTERS_

#include "ntuple_content.h"

// TagInfo                                                                                                                                                                                                 
#include "DataFormats/BTauReco/interface/PixelClusterTagInfo.h"

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

private:

    // seed candidates
    static constexpr size_t max_gtracks=10;
    
    unsigned int n_layers;

    unsigned int r004[max_gtracks];
    unsigned int r006[max_gtracks];
    unsigned int r008[max_gtracks];
    unsigned int r010[max_gtracks];
    unsigned int r016[max_gtracks];
    unsigned int rvar[max_gtracks];
    unsigned int rvwt[max_gtracks];

};

#endif /* DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_PIXEL_CLUSTERS_ */
