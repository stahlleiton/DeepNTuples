/*
 * ntuple_pixelclusters.cc
 *
 *  Created on: 16 February 2022
 *      Author: Alexandre De Moor & Joshuha Thomas-Wilsker

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


ntuple_pixelclusters::ntuple_pixelclusters(double jetR):ntuple_content(jetR){}

ntuple_pixelclusters::~ntuple_pixelclusters(){}

void ntuple_pixelclusters::getInput(const edm::ParameterSet& iConfig){
}

void ntuple_pixelclusters::initBranches(TTree* tree){
    
    addBranch(tree,"n_layers",&n_layers, "n_layers/I");
    
    addBranch(tree,"r004", &r004, "r004[n_layers]/I");
    addBranch(tree,"r006", &r006, "r006[n_layers]/I");
    addBranch(tree,"r008", &r008, "r008[n_layers]/I");
    addBranch(tree,"r010", &r010, "r010[n_layers]/I");
    addBranch(tree,"r016", &r016, "r016[n_layers]/I");
    
    addBranch(tree,"rvar", &rvar, "rvar[n_layers]/I");
    addBranch(tree,"rvwt", &rvwt, "rvwt[n_layers]/I");
}

void ntuple_pixelclusters::readEvent(const edm::Event& iEvent){
}

void ntuple_pixelclusters::readSetup(const edm::EventSetup& iSetup){
}

void ntuple_pixelclusters::checkEventSetup(const edm::EventSetup & iSetup) {
}

bool ntuple_pixelclusters::fillBranches(const pat::Jet & jet, const size_t& jetidx, const  edm::View<pat::Jet> * coll){

  if(jet.hasTagInfo("pixelCluster")){
    const reco::PixelClusterTagInfo *b = static_cast<const reco::PixelClusterTagInfo*>( jet.tagInfo("pixelCluster") );
    const reco::PixelClusterData data = b->data();
    //std::cout << "Yes " << jet.pt() << std::endl;                                                                                                                                                       

    unsigned int counter1 = 0;

    for (size_t ind=0; ind < data.r016.size(); ind++){

      r004[counter1] = catchInfsAndBound((int)data.r004[ind],-1.0,-100.0,100.0);
      r006[counter1] = catchInfsAndBound((int)data.r006[ind],-1.0,-100.0,100.0);
      r008[counter1] = catchInfsAndBound((int)data.r008[ind],-1.0,-100.0,100.0);
      r010[counter1] = catchInfsAndBound((int)data.r010[ind],-1.0,-100.0,100.0);
      r016[counter1] = catchInfsAndBound((int)data.r016[ind],-1.0,-100.0,100.0);

      rvar[counter1] = catchInfsAndBound((int)data.rvar[ind],-1.0,-100.0,100.0);
      rvwt[counter1] = catchInfsAndBound((int)data.rvwt[ind],-1.0,-100.0,100.0);

      counter1++;
    }

    n_layers = counter1;
  }

  else{
    //std::cout << "No " << jet.pt() << std::endl;                                                                                                                                                        

    unsigned int counter1 = 0;

    for (size_t ind=0; ind < 4; ind++){

      r004[counter1] = 0; //catchInfsAndBound((int)data.r004[ind],-1.0,-100.0,100.0);                                                                                                                     
      r006[counter1] = 0; //catchInfsAndBound((int)data.r006[ind],-1.0,-100.0,100.0);                                                                                                                     
      r008[counter1] = 0; //catchInfsAndBound((int)data.r008[ind],-1.0,-100.0,100.0);                                                                                                                     
      r010[counter1] = 0; //catchInfsAndBound((int)data.r010[ind],-1.0,-100.0,100.0);                                                                                                                     
      r016[counter1] = 0; //catchInfsAndBound((int)data.r016[ind],-1.0,-100.0,100.0);                                                                                                                     

      rvar[counter1] = 0; //catchInfsAndBound((int)data.rvar[ind],-1.0,-100.0,100.0);                                                                                                                     
      rvwt[counter1] = 0; //catchInfsAndBound((int)data.rvwt[ind],-1.0,-100.0,100.0);                                                                                                                     

      counter1++;
    }

    n_layers = counter1;
  }

  return true;
}
