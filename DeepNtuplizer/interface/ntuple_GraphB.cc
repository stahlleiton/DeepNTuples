/*
 * ntuple_GraphB.cc
 *
 *  Created on: 23 June 2017
 *      Author: Seth Moortgat

 */


#include "../interface/ntuple_GraphB.h"

#include "DataFormats/GeometrySurface/interface/Line.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
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
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TVector3.h"

class TrackInfoBuilder{
public:
    TrackInfoBuilder(edm::ESHandle<TransientTrackBuilder> & build):
        builder(build),
        trackMomentum_(0),
        trackEta_(0),
        trackEtaRel_(0),
        trackPtRel_(0),
        trackPPar_(0),
        trackDeltaR_(0),
        trackPtRatio_(0),
        trackPParRatio_(0),
        trackSip2dVal_(0),
        trackSip2dSig_(0),
        trackSip3dVal_(0),
        trackSip3dSig_(0),
        trackJetDistVal_(0),
        trackJetDistSig_(0),
	trackImpactPointState_(0)
{


}

    void buildTrackInfo(const pat::PackedCandidate* PackedCandidate_ ,const math::XYZVector&  jetDir, GlobalVector refjetdirection, const reco::Vertex & pv){
			TVector3 jetDir3(jetDir.x(),jetDir.y(),jetDir.z());
			if(!PackedCandidate_->hasTrackDetails()) {
				TVector3 trackMom3(
					PackedCandidate_->momentum().x(),
					PackedCandidate_->momentum().y(),
					PackedCandidate_->momentum().z()
					);
				trackMomentum_=PackedCandidate_->p();
				trackEta_= PackedCandidate_->eta();
				trackEtaRel_=reco::btau::etaRel(jetDir, PackedCandidate_->momentum());
				trackPtRel_=trackMom3.Perp(jetDir3);
				trackPPar_=jetDir.Dot(PackedCandidate_->momentum());
				trackDeltaR_=reco::deltaR(PackedCandidate_->momentum(), jetDir);
				trackPtRatio_=trackMom3.Perp(jetDir3) / PackedCandidate_->p();
				trackPParRatio_=jetDir.Dot(PackedCandidate_->momentum()) / PackedCandidate_->p();
				trackSip2dVal_=0.;
				trackSip2dSig_=0.;
				trackSip3dVal_=0.;
				trackSip3dSig_=0.;
				trackJetDistVal_=0.;
				trackJetDistSig_=0.;
				//trackImpactPointState_=PackedCandidate_->impactPointState();
				return;
			}

        const reco::Track & PseudoTrack =  PackedCandidate_->pseudoTrack();

        reco::TransientTrack transientTrack;
        transientTrack=builder->build(PseudoTrack);
        Measurement1D meas_ip2d=IPTools::signedTransverseImpactParameter(transientTrack, refjetdirection, pv).second;
        Measurement1D meas_ip3d=IPTools::signedImpactParameter3D(transientTrack, refjetdirection, pv).second;
        Measurement1D jetdist=IPTools::jetTrackDistance(transientTrack, refjetdirection, pv).second;
        math::XYZVector trackMom = PseudoTrack.momentum();
        double trackMag = std::sqrt(trackMom.Mag2());
        TVector3 trackMom3(trackMom.x(),trackMom.y(),trackMom.z());


        trackMomentum_=std::sqrt(trackMom.Mag2());
        trackEta_= trackMom.Eta();
        trackEtaRel_=reco::btau::etaRel(jetDir, trackMom);
        trackPtRel_=trackMom3.Perp(jetDir3);
        trackPPar_=jetDir.Dot(trackMom);
        trackDeltaR_=reco::deltaR(trackMom, jetDir);
        trackPtRatio_=trackMom3.Perp(jetDir3) / trackMag;
        trackPParRatio_=jetDir.Dot(trackMom) / trackMag;
        trackSip2dVal_=(meas_ip2d.value());

        trackSip2dSig_=(meas_ip2d.significance());
        trackSip3dVal_=(meas_ip3d.value());


        trackSip3dSig_=meas_ip3d.significance();
        trackJetDistVal_= jetdist.value();
        trackJetDistSig_= jetdist.significance();
	trackImpactPointState_=transientTrack.impactPointState();

    }

    const float& getTrackDeltaR() const {return trackDeltaR_;}
    const float& getTrackEta() const {return trackEta_;}
    const float& getTrackEtaRel() const {return trackEtaRel_;}
    const float& getTrackJetDistSig() const {return trackJetDistSig_;}
    const float& getTrackJetDistVal() const {return trackJetDistVal_;}
    const float& getTrackMomentum() const {return trackMomentum_;}
    const float& getTrackPPar() const {return trackPPar_;}
    const float& getTrackPParRatio() const {return trackPParRatio_;}
    const float& getTrackPtRatio() const {return trackPtRatio_;}
    const float& getTrackPtRel() const {return trackPtRel_;}
    const float& getTrackSip2dSig() const {return trackSip2dSig_;}
    const float& getTrackSip2dVal() const {return trackSip2dVal_;}
    const float& getTrackSip3dSig() const {return trackSip3dSig_;}
    const float& getTrackSip3dVal() const {return trackSip3dVal_;}
    const TrajectoryStateOnSurface& getImpactPointState() const {return trackImpactPointState_;}

private:

    edm::ESHandle<TransientTrackBuilder>& builder;

    float trackMomentum_;
    float trackEta_;
    float trackEtaRel_;
    float trackPtRel_;
    float trackPPar_;
    float trackDeltaR_;
    float trackPtRatio_;
    float trackPParRatio_;
    float trackSip2dVal_;
    float trackSip2dSig_;
    float trackSip3dVal_;
    float trackSip3dSig_;

    float trackJetDistVal_;
    float trackJetDistSig_;
    TrajectoryStateOnSurface trackImpactPointState_;

};

ntuple_GraphB::ntuple_GraphB(double jetR):ntuple_content(jetR){}

ntuple_GraphB::~ntuple_GraphB(){}

void ntuple_GraphB::getInput(const edm::ParameterSet& iConfig){
  min_candidate_pt_ = (iConfig.getParameter<double>("minCandidatePt"));
}

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
    addBranch(tree,"gtrack_dist_neigh", &gtrack_dist_neigh, "gtrack_dist_neigh[n_gtracks]/f");
    
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

void ntuple_GraphB::readEvent(const edm::Event& iEvent){
    iEvent.getByToken(CandidateToken, tracks);
    n_Npfcand_=0;
    n_Cpfcand_=0;
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
	iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);	

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

    double jet_radius = jetR();

    math::XYZVector jetDir = jet.momentum().Unit();
    GlobalVector direction(jet.px(), jet.py(), jet.pz());
    std::vector<sorting::sortingClass<size_t> > sortedcharged;
    std::vector<size_t> sortedchargedindices;

    const float jet_uncorr_pt=jet.correctedJet("Uncorrected").pt();
    //const float jet_uncorr_e=jet.correctedJet("Uncorrected").energy();

    for(size_t k = 0; k<tracks->size(); ++k) {
      if((*tracks)[k].bestTrack() != 0 &&  (*tracks)[k].pt()>0.5 && std::fabs(pvp.z()-builder->build(tracks->ptrAt(k)).track().vz())<0.5) {
	selectedTracks.push_back(builder->build(tracks->ptrAt(k)));
	masses.push_back(tracks->ptrAt(k)->mass());
        }
    }
        
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


    TrackInfoBuilder trackinfo(builder);
    //create collection of cpf/npf for the seeding
    for (unsigned int i = 0; i <  jet.numberOfDaughters(); i++){
        const pat::PackedCandidate* PackedCandidate = dynamic_cast<const pat::PackedCandidate*>(jet.daughter(i));
        if(PackedCandidate){
	  if(PackedCandidate->pt() < min_candidate_pt_) continue; 
            if(PackedCandidate->charge()!=0){
                trackinfo.buildTrackInfo(PackedCandidate,jetDir,direction,pv);
                sortedcharged.push_back(sorting::sortingClass<size_t>
                (i, trackinfo.getTrackSip2dSig(),
                        -mindrsvpfcand(PackedCandidate), PackedCandidate->pt()/jet_uncorr_pt));
            }
        }
    }
    std::sort(sortedcharged.begin(),sortedcharged.end(),sorting::sortingClass<size_t>::compareByABCInv);
    long unsigned int max_cpf = 30;
    n_Cpfcand_ = std::min(sortedcharged.size(),max_cpf);
    sortedchargedindices = sorting::invertSortingVector(sortedcharged);
 
    size_t counter= 0;

    for(std::vector<reco::TransientTrack>::const_iterator it = selectedTracks.begin(); it != selectedTracks.end(); it++){  
      //is the track in the jet cone?
      float angular_distance=reco::deltaR(jet,it->track());
      //std::sqrt(std::pow(jet.eta()-it->track().eta(),2) + std::pow(jet.phi()-it->track().phi(),2) );
      bool hasNeighbour = false;
      bool include = false;
      float dist_part = -0.01;

      //matching the track with the jet radius
      if(/*(angular_distance < 1.50*jet_radius) &&*/ (angular_distance > jet_radius)){
	if((std::fabs(pvp.z() - it->track().vz()) > 0.1)) {continue;}  	
	   //matching the track with a cpf_seed of not

	  for (unsigned int i = 0; i <  jet.numberOfDaughters(); i++){
	    const pat::PackedCandidate* PackedCandidate = dynamic_cast<const pat::PackedCandidate*>(jet.daughter(i));
	    if(PackedCandidate->hasTrackDetails()){
	      if(PackedCandidate->pt() < min_candidate_pt_) {continue;} 
	      if(PackedCandidate->charge()!=0){
		const reco::Track & PseudoTrack =  PackedCandidate->pseudoTrack();
		reco::TransientTrack transientTrack;
		transientTrack=builder->build(PseudoTrack);
	      
		TwoTrackMinimumDistance dist;
		std::pair<bool, Measurement1D> ip=IPTools::absoluteImpactParameter3D(transientTrack, pv);
		float near_angular_dist = reco::deltaR(jet,transientTrack.track());
		//std::pair<double, Measurement1D> jet_dist =IPTools::jetTrackDistance(transientTrack, direction, pv);
		//float length = 999;
		//TrajectoryStateOnSurface closest = IPTools::closestApproachToJet(transientTrack.impactPointState(), pv, direction, transientTrack.field());
		//if(closest.isValid()){
		// length = (closest.globalPosition() - pvp).mag();
		//}
		if (transientTrack == *it) {continue;}
		if (near_angular_dist < jet_radius){
		  /*if (!(ip.first && ip.second.value() >= 0.0 && ip.second.significance() >= 1.0 && ip.second.value() <= 9999. 
			&& ip.second.significance() <= 9999. && transientTrack.track().normalizedChi2() < 5. && std::fabs(transientTrack.track().dxy(pv.position())) < 2 
			&& std::fabs(transientTrack.track().dz(pv.position())) < 17 && jet_dist.second.value() < 0.07 && length < 5.)){continue;}*/
		  if(ip.second.significance() < 1.0) {continue;}
		  if(dist.calculate(transientTrack.impactPointState(),it->impactPointState())){
		    float distance = dist.distance();
		    if(distance < 0.02){
		      hasNeighbour = true;
		      dist_part = distance;
		    }
		  }
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
      if(angular_distance<=jet_radius){
	include = true;
	}
      
      if(include){

	std::pair<bool,Measurement1D> ip = IPTools::signedImpactParameter3D(*it, direction, pv);        
        std::pair<bool,Measurement1D> ip2d = IPTools::signedTransverseImpactParameter(*it, direction, pv);
        std::pair<double, Measurement1D> jet_dist =IPTools::jetTrackDistance(*it, direction, pv);                   
        TrajectoryStateOnSurface closest = IPTools::closestApproachToJet(it->impactPointState(),pv, direction,it->field());
	if(counter >= max_gtracks){continue;}	
	gtrack_dR[counter] = catchInfsAndBound(reco::deltaR(jet,it->track()), -1.0,0.0,10.0);
	gtrack_dist_neigh[counter] = catchInfsAndBound(dist_part, -10.0,-5.0,100.0);
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

float ntuple_GraphB::mindrsvpfcand(const pat::PackedCandidate* pfcand) {

  float mindr_ = jetradius_;
  for (unsigned int i=0; i<secVertices()->size(); ++i) {
    if(!pfcand) continue;
    //if(!svs.at(i)) continue;
    float tempdr_ = reco::deltaR(secVertices()->at(i),*pfcand);
    if (tempdr_<mindr_) { mindr_ = tempdr_; }

  }
  return mindr_;
}
