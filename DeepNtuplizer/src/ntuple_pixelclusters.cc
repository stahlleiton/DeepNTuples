/*
 * ntuple_pixelclusters.cc
 *
 *  Created on: 22 December 2021
 *      Author: Alexandre De Moor

 */


#include "../interface/ntuple_pixelclusters.h"

// system include files
#include <memory>


#include "DataFormats/Common/interface/CMS_CLASS_VERSION.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/ESGetToken.h"

// TagInfo
#include "DataFormats/BTauReco/interface/PixelClusterTagInfo.h"

// For vertices
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// For jet
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

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

struct PixelClusterProperties {
    float x = 0;
    float y = 0;
    float z = 0;
    int charge = 0;
    unsigned int layer = 0;
  };

struct PixelClusterData {
    std::vector<int8_t> r004;
    std::vector<int8_t> r006;
    std::vector<int8_t> r008;
    std::vector<int8_t> r010;
    std::vector<int8_t> r016;
    std::vector<int8_t> rvar;
    std::vector<unsigned int> rvwt;
    PixelClusterData(unsigned int l = 4) {
      r004 = std::vector<int8_t>(l, 0);
      r006 = std::vector<int8_t>(l, 0);
      r008 = std::vector<int8_t>(l, 0);
      r010 = std::vector<int8_t>(l, 0);
      r016 = std::vector<int8_t>(l, 0);
      rvar = std::vector<int8_t>(l, 0);
      rvwt = std::vector<unsigned int>(l, 0);
    }
    CMS_CLASS_VERSION(3)
  };

ntuple_pixelclusters::ntuple_pixelclusters(double jetR):ntuple_content(jetR){}

ntuple_pixelclusters::~ntuple_pixelclusters(){}

void ntuple_pixelclusters::getInput(const edm::ParameterSet& iConfig){
}

void ntuple_pixelclusters::initBranches(TTree* tree){
    
    addBranch(tree,"n_layers",&n_layers, "n_layers/I");
    
    addBranch(tree,"r004", &r004, "r004[n_layers]/F");
    addBranch(tree,"r006", &r006, "r006[n_layers]/F");
    addBranch(tree,"r008", &r008, "r008[n_layers]/F");
    addBranch(tree,"r010", &r010, "r010[n_layers]/F");
    addBranch(tree,"r016", &r016, "r016[n_layers]/F");
    
    addBranch(tree,"rvar", &rvar, "rvar[n_layers]/F");
    addBranch(tree,"rvwt", &rvwt, "rvwt[n_layers]/F");
}

void ntuple_pixelclusters::readEvent(const edm::Event& iEvent){
    iEvent.getByToken(m_pixelhit, collectionHandle);
}

void ntuple_pixelclusters::readSetup(const edm::EventSetup& iSetup){
    // Open Geometry
    iSetup.get<TrackerDigiGeometryRecord>().get(geom);
    // Retrieve tracker topology from geometry 
    iSetup.get<TrackerTopologyRcd>().get(tTopoH);
}

void ntuple_pixelclusters::checkEventSetup(const edm::EventSetup & iSetup) {
    // Open Geometry                                                                                                                                                                                       
    //theTracker = &iSetup.getData(m_geomToken);
    // Retrieve tracker topology from geometry
    //tTopo = &iSetup.getData(m_topoToken);
}

bool ntuple_pixelclusters::fillBranches(const pat::Jet & jet, const size_t& jetidx, const  edm::View<pat::Jet> * coll){

    const bool m_addFPIX = true;
    const int m_minADC = -1;
    unsigned int m_nLayers = 4;
    float m_hadronMass = 12.0;
    
    // pv info
    const reco::Vertex &pv = vertices()->at(0);
    GlobalPoint v3(pv.x(),pv.y(),pv.z());
    
    // Open Pixel Cluster collection
    const edmNew::DetSetVector<SiPixelCluster>& collectionClusters(*collectionHandle);
    const TrackerGeometry& theTracker(*geom);
    const TrackerTopology* tTopo = tTopoH.product();
    
    std::vector<reco::PixelClusterProperties> clusters;
    
    // Get vector of detunit ids, and fill a vector of PixelClusterProperties in the loop
    for (auto const& detUnit : collectionClusters) {
        if (detUnit.empty())
            continue;
        DetId detId = DetId(detUnit.detId());  // Get the Detid object for pixel detector selection
        if (detId.det() != DetId::Tracker)
            continue;
        if (!(detId.subdetId() == PixelSubdetector::PixelBarrel ||
            (m_addFPIX && detId.subdetId() == PixelSubdetector::PixelEndcap)))
            continue;
        unsigned int layer = tTopo->layer(detId);  // The layer index is in range 1-4 or 1-3
        if (layer == 0 || layer > m_nLayers)
            continue;
        
        // Get the geom-detector
        const auto* geomDet = theTracker.idToDet(detId);
	const auto* topol = &geomDet->topology();
        
        for (auto const& clUnit : detUnit) {
            if (m_minADC > 0 and clUnit.charge() < m_minADC)
                continue;  // skip cluster if below threshold
            // Get global position of the cluster
            LocalPoint lp = topol->localPosition(MeasurementPoint(clUnit.x(), clUnit.y()));
            GlobalPoint gp = geomDet->surface().toGlobal(lp);
            // Fill PixelClusterProperties vector for matching
            reco::PixelClusterProperties cp = {gp.x(), gp.y(), gp.z(), clUnit.charge(), layer};
            clusters.push_back(cp);
        }
    }
    
    float cR = m_hadronMass * 2. / (jet.pt());

    reco::PixelClusterData data(m_nLayers);
    size_t counter= 0;
    //if(counter >= 4){continue;}
    for (auto const& cluster : clusters) {
        GlobalPoint c3(cluster.x, cluster.y, cluster.z);  // Get cluster 3D position
        float dR = reco::deltaR(c3 - v3, jet.momentum());
        // Match pixel clusters to jets and fill Data struct
        if (cluster.layer >= 1 && cluster.layer <= m_nLayers) {
            int idx(cluster.layer - 1);
            if (dR < 0.16) {
                data.r016[idx]++;
                if (dR < 0.10) {
                    data.r010[idx]++;
                    if (dR < 0.08) {
                        data.r008[idx]++;
                        if (dR < 0.06) {
                            data.r006[idx]++;
                            if (dR < 0.04){
                                data.r004[idx]++;
                            }
                        }
                    }
                }
            }
            if (dR < cR) {
                data.rvar[idx]++;
                data.rvwt[idx] += cluster.charge;
            }
        }
    }
    
    for (long unsigned int ind=0; ind<data.r016.size(); ind++){
    
        r004[counter] = catchInfsAndBound(data.r004[ind],-1.0,-100.0,100.0);
        r006[counter] = catchInfsAndBound(data.r006[ind],-1.0,-100.0,100.0);
        r008[counter] = catchInfsAndBound(data.r008[ind],-1.0,-100.0,100.0);
        r010[counter] = catchInfsAndBound(data.r010[ind],-1.0,-100.0,100.0);
        r016[counter] = catchInfsAndBound(data.r016[ind],-1.0,-100.0,100.0);
        
        rvar[counter] = catchInfsAndBound(data.rvar[ind],-1.0,-100.0,100.0);
        rvwt[counter] = catchInfsAndBound(data.rvwt[ind],-1.0,-100.0,100.0);
        
        counter++;
    }
    
    n_layers = counter;

    return true;
}
