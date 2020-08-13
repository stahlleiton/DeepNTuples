/*
 * ntuple_GraphB.cc
 *
 *  Created on: 23 June 2017
 *      Author: Seth Moortgat

 */


#include "../interface/ntuple_GraphB.h"

#include "DataFormats/GeometrySurface/interface/Line.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "../interface/sorting_modules.h"
#include "CondFormats/BTauObjects/interface/TrackProbabilityCalibration.h"
#include "CondFormats/DataRecord/interface/BTagTrackProbability2DRcd.h"
#include "CondFormats/DataRecord/interface/BTagTrackProbability3DRcd.h"
#include "FWCore/Framework/interface/EventSetupRecord.h"
#include "FWCore/Framework/interface/EventSetupRecordImplementation.h"
#include "FWCore/Framework/interface/EventSetupRecordKey.h"

ntuple_GraphB::ntuple_GraphB(double jetR):ntuple_content(jetR){}

ntuple_GraphB::~ntuple_GraphB(){}

void ntuple_GraphB::getInput(const edm::ParameterSet& iConfig){}

void ntuple_GraphB::initBranches(TTree* tree){
    
    addBranch(tree,"n_gtracks",&n_gtracks, "n_gtracks/i");
    addBranch(tree,"nGtracks",&nGtracks, "nGtracks/f");

    addBranch(tree,"gtrack_pt",&gtrack_pt, "gtrack_pt[n_gtracks]/f");
    addBranch(tree,"gtrack_eta",&gtrack_eta, "gtrack_eta[n_gtracks]/f");
    addBranch(tree,"gtrack_phi",&gtrack_phi, "gtrack_phi[n_gtracks]/f");
    addBranch(tree,"gtrack_mass",&gtrack_mass, "gtrack_mass[n_gtracks]/f");
    
    addBranch(tree,"gtrack_dz", &gtrack_dz, "gtrack_dz[n_gtracks]/f");
    addBranch(tree,"gtrack_dxy", &gtrack_dxy, "gtrack_dxy[n_gtracks]/f");
    addBranch(tree,"gtrack_3D_ip", &gtrack_3D_ip, "gtrack_3D_ip[n_gtracks]/f");
    addBranch(tree,"gtrack_3D_sip", &gtrack_3D_sip, "gtrack_3D_sip[n_gtracks]/f");
    addBranch(tree,"gtrack_2D_ip", &gtrack_2D_ip, "gtrack_2D_ip[n_gtracks]/f");
    addBranch(tree,"gtrack_2D_sip", &gtrack_2D_sip, "gtrack_2D_sip[n_gtracks]/f");
    addBranch(tree,"gtrack_dR", &gtrack_dR, "gtrack_dR[n_gtracks]/f");
    
    addBranch(tree,"gtrack_3D_TrackProbability", &gtrack_3D_TrackProbability, "gtrack_3D_TrackProbability[n_gtracks]/f");
    addBranch(tree,"gtrack_2D_TrackProbability", &gtrack_2D_TrackProbability, "gtrack_2D_TrackProbability[n_gtracks]/f");
    
    addBranch(tree,"gtrack_chi2reduced",&gtrack_chi2reduced, "gtrack_chi2reduced[n_gtracks]/f");
    addBranch(tree,"gtrack_nPixelHits",&gtrack_nPixelHits, "gtrack_nPixelHits[n_gtracks]/f");
    addBranch(tree,"gtrack_nHits",&gtrack_nHits, "gtrack_nHits[n_gtracks]/f");
    addBranch(tree,"gtrack_jetAxisDistance",&gtrack_jetAxisDistance, "gtrack_jetAxisDistance[n_gtracks]/f");
    addBranch(tree,"gtrack_jetAxisDlength",&gtrack_jetAxisDlength, "gtrack_jetAxisDlength[n_gtracks]/f");
    addBranch(tree,"gtrack_PCAtrackFromPV",&gtrack_PCAtrackFromPV, "gtrack_PCAtrackFromPV[n_gtracks]/f");
    addBranch(tree,"gtrack_dotProdTrack",&gtrack_dotProdTrack, "gtrack_dotProdTrack[n_gtracks]/f");
    addBranch(tree,"gtrack_dotProdTrack2D",&gtrack_dotProdTrack2D, "gtrack_dotProdTrack2D[n_gtracks]/f");
    

}


void ntuple_GraphB::readEvent(const edm::Event& iEvent)
{
    iEvent.getByToken(CandidateToken, tracks);
}


void ntuple_GraphB::readSetup(const edm::EventSetup& iSetup){
    
    //this part was to be in checkEventSetup, but idk how to call it
    using namespace edm;
    using namespace edm::eventsetup;

    const EventSetupRecord & re2D= iSetup.get<BTagTrackProbability2DRcd>();
    const EventSetupRecord & re3D= iSetup.get<BTagTrackProbability3DRcd>();
    unsigned long long cacheId2D= re2D.cacheIdentifier();
    unsigned long long cacheId3D= re3D.cacheIdentifier();

    if(cacheId2D!=m_calibrationCacheId2D || cacheId3D!=m_calibrationCacheId3D  )  //Calibration changed
    {
        //iSetup.get<BTagTrackProbabilityRcd>().get(calib);
        ESHandle<TrackProbabilityCalibration> calib2DHandle;
        iSetup.get<BTagTrackProbability2DRcd>().get(calib2DHandle);
        ESHandle<TrackProbabilityCalibration> calib3DHandle;
        iSetup.get<BTagTrackProbability3DRcd>().get(calib3DHandle);

        const TrackProbabilityCalibration *  ca2D= calib2DHandle.product();
        const TrackProbabilityCalibration *  ca3D= calib3DHandle.product();

        m_probabilityEstimator.reset(new HistogramProbabilityEstimator(ca3D,ca2D));

    }

    m_calibrationCacheId3D=cacheId3D;
    m_calibrationCacheId2D=cacheId2D;
    
    //readEvent only line
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);

}

void ntuple_GraphB::checkEventSetup(const edm::EventSetup & iSetup) {
  
  using namespace edm;
  using namespace edm::eventsetup;

   const EventSetupRecord & re2D= iSetup.get<BTagTrackProbability2DRcd>();
   const EventSetupRecord & re3D= iSetup.get<BTagTrackProbability3DRcd>();
   unsigned long long cacheId2D= re2D.cacheIdentifier();
   unsigned long long cacheId3D= re3D.cacheIdentifier();

   if(cacheId2D!=m_calibrationCacheId2D || cacheId3D!=m_calibrationCacheId3D  )  //Calibration changed
   {
     //iSetup.get<BTagTrackProbabilityRcd>().get(calib);
     ESHandle<TrackProbabilityCalibration> calib2DHandle;
     iSetup.get<BTagTrackProbability2DRcd>().get(calib2DHandle);
     ESHandle<TrackProbabilityCalibration> calib3DHandle;
     iSetup.get<BTagTrackProbability3DRcd>().get(calib3DHandle);

     const TrackProbabilityCalibration *  ca2D= calib2DHandle.product();
     const TrackProbabilityCalibration *  ca3D= calib3DHandle.product();

     m_probabilityEstimator.reset(new HistogramProbabilityEstimator(ca3D,ca2D));

   }
   
   m_calibrationCacheId3D=cacheId3D;
   m_calibrationCacheId2D=cacheId2D;
   
}


bool ntuple_GraphB::fillBranches(const pat::Jet & jet, const size_t& jetidx, const  edm::View<pat::Jet> * coll){

    // pv info
    const reco::Vertex &pv = vertices()->at(0);
    GlobalPoint pvp(pv.x(),pv.y(),pv.z());

    
    std::vector<reco::TransientTrack> selectedTracks;
    std::vector<float> masses;
    

   for(size_t k = 0; k<tracks->size(); ++k) {
        if((*tracks)[k].bestTrack() != 0 &&  (*tracks)[k].pt()>0.5 && std::fabs(pvp.z()-builder->build(tracks->ptrAt(k)).track().vz())<0.5) {
            selectedTracks.push_back(builder->build(tracks->ptrAt(k)));
            masses.push_back(tracks->ptrAt(k)->mass());
        }
    }
    
    double jet_radius = jetR();
    GlobalVector direction(jet.px(), jet.py(), jet.pz());
    std::vector<sorting::sortingClass<size_t> > sorted_tracks;
    for(std::vector<reco::TransientTrack>::const_iterator it = selectedTracks.begin(); it != selectedTracks.end(); it++){  
      float angular_distance=reco::deltaR(jet,it->track());
      sorted_tracks.push_back(sorting::sortingClass<size_t>
			      (it-selectedTracks.begin(), -angular_distance,
			       -1, -1));
    }
    std::sort(sorted_tracks.begin(),sorted_tracks.end(),sorting::sortingClass<size_t>::compareByABCInv);
    std::vector<size_t> sorted_track_indices;
    sorted_track_indices=sorting::invertSortingVector(sorted_tracks);
    //n_gtracks = std::min(sorted_tracks.size(),max_gtracks);
    //nGtracks = n_gtracks;
 
    size_t counter= 0;
    for(std::vector<reco::TransientTrack>::const_iterator it = selectedTracks.begin(); it != selectedTracks.end(); it++){  
      //is the track in the jet cone?
      float angular_distance=reco::deltaR(jet,it->track());
      //std::sqrt(std::pow(jet.eta()-it->track().eta(),2) + std::pow(jet.phi()-it->track().phi(),2) );
      bool hasNeighbour = false;
      bool include = false;
      if (angular_distance>jet_radius) {
	if(std::fabs(pvp.z()-it->track().vz())>0.1) continue;
	for(std::vector<reco::TransientTrack>::const_iterator tt = selectedTracks.begin();tt!=selectedTracks.end(); ++tt ) {
	  TwoTrackMinimumDistance dist;
	  float near_angular_distance=reco::deltaR(jet,tt->track());
	  if(*tt==*it) continue;
	  if(near_angular_distance<jet_radius){
	    std::pair<bool,Measurement1D> ip = IPTools::absoluteImpactParameter3D(*tt, pv);
	    if(ip.second.significance() < 1.0) continue;
	    if(dist.calculate(tt->impactPointState(),it->impactPointState())) {
	      float distance = dist.distance();
	      if(distance < 0.02){
		hasNeighbour = true;
	      }
	    }
	  }
	}
	if(hasNeighbour){
	  int matched_jets = 0;
	  for (std::size_t jet_n = 0; jet_n < coll->size(); jet_n++) {
	    const auto& test_jet = coll->at(jet_n);
	    if(test_jet.pt() < 5.0){continue;}
	    float new_angular_distance=reco::deltaR(test_jet,it->track());
	    if(new_angular_distance < jet_radius){
	      matched_jets = matched_jets + 1;
	    }
	  }
	  if(matched_jets != 0){continue;}
	  include = true;
	}
      }
      if(angular_distance<jet_radius){
	include = true;
      }
      if(include){
	std::pair<bool,Measurement1D> ip = IPTools::signedImpactParameter3D(*it, direction, pv);        
        std::pair<bool,Measurement1D> ip2d = IPTools::signedTransverseImpactParameter(*it, direction, pv);
        std::pair<double, Measurement1D> jet_dist =IPTools::jetTrackDistance(*it, direction, pv);                   
        TrajectoryStateOnSurface closest = IPTools::closestApproachToJet(it->impactPointState(),pv, direction,it->field());
	if(counter >= max_gtracks){continue;}	
	gtrack_dR[counter] = catchInfsAndBound( reco::deltaR(jet,it->track()), -1.0,0.0,10.0);
	float length=999;
        if (closest.isValid()) length=(closest.globalPosition() - pvp).mag();
	gtrack_jetAxisDlength[counter] = catchInfsAndBound(length,-1.0,-100.0,100.0);
        gtrack_jetAxisDistance[counter] = catchInfsAndBound(jet_dist.second.value(),-1.0,-100.0,100.0);
	gtrack_chi2reduced[counter] = catchInfsAndBound(it->track().normalizedChi2(),-1.0,-100.0,300.0);        
	gtrack_3D_ip[counter] = catchInfsAndBound(ip.second.value(),0.0,-500.0,500.0); 
	gtrack_3D_sip[counter] = catchInfsAndBound(ip.second.significance(),0.0,-500.0,500.0);
	gtrack_2D_ip[counter] = catchInfsAndBound(ip2d.second.value(),0.0,-500.0,500.0); 
	gtrack_2D_sip[counter] = catchInfsAndBound(ip2d.second.significance(),0.0,-500.0,500.0);
        gtrack_pt[counter] = catchInfsAndBound(it->track().pt(),0.0,-100.0,1000.0);
        gtrack_eta[counter] = catchInfsAndBound(it->track().eta(),0.0,-2.5,2.5);
        gtrack_phi[counter] = catchInfsAndBound(it->track().phi(),0.0,-5,5);
	gtrack_mass[counter] = catchInfsAndBound(masses[it-selectedTracks.begin()],0.0,-1.0,500.0);
	gtrack_dxy[counter] = catchInfsAndBound(it->track().dxy(pv.position()),0.0,-100.0,100.0);
	gtrack_dz[counter] = catchInfsAndBound(it->track().dz(pv.position()),0.0,-100.0,100.0);
	gtrack_PCAtrackFromPV[counter] =  catchInfsAndBound( (it->impactPointState().globalPosition()-pvp).mag() , -1.0,-100.0,100.0);               
	gtrack_dotProdTrack[counter] = catchInfsAndBound( (it->impactPointState().globalPosition()-pvp).unit().dot(it->impactPointState().globalDirection().unit()), -3.0,-3.0,3.0);
	GlobalVector trackDir2D(it->impactPointState().globalDirection().x(),it->impactPointState().globalDirection().y(),0.);
	GlobalPoint pvp2d(pv.x(),pv.y(),0.0);
	GlobalPoint trackPos2D(it->impactPointState().globalPosition().x(),it->impactPointState().globalPosition().y(),0.0);
	gtrack_dotProdTrack2D[counter] = catchInfsAndBound( (trackPos2D-pvp2d).unit().dot(trackDir2D.unit()), -3.0,-3.0,3.0);
	std::pair<bool,double> probability;
	//probability with 3D ip
	probability = m_probabilityEstimator->probability(0, 0,ip.second.significance(),it->track(),jet,pv);
	double prob3D=(probability.first ? probability.second : -1.);
	gtrack_3D_TrackProbability[counter] = catchInfsAndBound(prob3D,-2.0,-100.0,100.0);
	//probability with 2D ip
	probability = m_probabilityEstimator->probability(0, 1,ip2d.second.significance(),it->track(),jet,pv);
	double prob2D=(probability.first ? probability.second : -1.);
	gtrack_2D_TrackProbability[counter] = catchInfsAndBound(prob2D,-2.0,-100.0,100.0);

        gtrack_nPixelHits[counter] = catchInfsAndBound(it->hitPattern().numberOfValidPixelHits(),-1.0,0.0,100.0);
        gtrack_nHits[counter] = catchInfsAndBound(it->hitPattern().numberOfValidHits(),-1.0,0.0,100.0);
	counter++;
      }
    }
    n_gtracks = counter;
    nGtracks = n_gtracks;
 
    masses.clear();

    return true;
}

