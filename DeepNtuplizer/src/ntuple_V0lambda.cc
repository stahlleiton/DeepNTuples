/*
 * ntuple_V0lambda.cc
 *
 *  Created on: 3rd August 2023
 *      Author: Alexandre De Moor
 */


#include "../interface/ntuple_V0lambda.h"
// for ivf
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "TLorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "TVector3.h"

const reco::Vertex * ntuple_V0lambda::spvp_;

ntuple_V0lambda::ntuple_V0lambda(std::string prefix, double jetR):ntuple_content(jetR),v0lambda_num_(0){
    prefix_ = prefix;
}
ntuple_V0lambda::~ntuple_V0lambda(){}


void ntuple_V0lambda::getInput(const edm::ParameterSet& iConfig){

}

void ntuple_V0lambda::initBranches(TTree* tree){
  // SV candidates
  addBranch(tree,(prefix_+"n_v0lambda").c_str()           ,&v0lambda_num_         ,(prefix_+"v0lambda_num_/I").c_str()     );
  addBranch(tree,(prefix_+"nv0lambda").c_str()            ,&nv0lambda_            ,(prefix_+"nv0lambda_/F").c_str()         );
  addBranch(tree,(prefix_+"v0lambda_pt").c_str()          ,&v0lambda_pt_          ,(prefix_+"v0lambda_pt_["+prefix_+"v0lambda_num_]/F").c_str()        );
  addBranch(tree,(prefix_+"v0lambda_eta").c_str()         ,&v0lambda_eta_         ,(prefix_+"v0lambda_eta_["+prefix_+"v0lambda_num_]/F").c_str()        );
  addBranch(tree,(prefix_+"v0lambda_phi").c_str()         ,&v0lambda_phi_         ,(prefix_+"v0lambda_phi_["+prefix_+"v0lambda_num_]/F").c_str()        );
  addBranch(tree,(prefix_+"v0lambda_e").c_str()           ,&v0lambda_e_           ,(prefix_+"v0lambda_e_["+prefix_+"v0lambda_num_]/F").c_str()        );
  addBranch(tree,(prefix_+"v0lambda_etarel").c_str()      ,&v0lambda_etarel_      ,(prefix_+"v0lambda_etarel_["+prefix_+"v0lambda_num_]/F").c_str()         );
  addBranch(tree,(prefix_+"v0lambda_phirel").c_str()      ,&v0lambda_phirel_      ,(prefix_+"v0lambda_phirel_["+prefix_+"v0lambda_num_]/F").c_str()         );
  addBranch(tree,(prefix_+"v0lambda_deltaR").c_str()      ,&v0lambda_deltaR_      ,(prefix_+"v0lambda_deltaR_["+prefix_+"v0lambda_num_]/F").c_str()         );
  addBranch(tree,(prefix_+"v0lambda_mass").c_str()        ,&v0lambda_mass_        ,(prefix_+"v0lambda_mass_["+prefix_+"v0lambda_num_]/F").c_str()        );
  addBranch(tree,(prefix_+"v0lambda_ntracks").c_str()     ,&v0lambda_ntracks_     ,(prefix_+"v0lambda_ntracks_["+prefix_+"v0lambda_num_]/F").c_str()     );
  addBranch(tree,(prefix_+"v0lambda_chi2").c_str()        ,&v0lambda_chi2_        ,(prefix_+"v0lambda_chi2_["+prefix_+"v0lambda_num_]/F").c_str()        );
  addBranch(tree,(prefix_+"v0lambda_ndf").c_str()         ,&v0lambda_ndf_         ,(prefix_+"v0lambda_ndf_["+prefix_+"v0lambda_num_]/F").c_str()         );
  addBranch(tree,(prefix_+"v0lambda_normchi2").c_str()    ,&v0lambda_normchi2_    ,(prefix_+"v0lambda_normchi2_["+prefix_+"v0lambda_num_]/F").c_str()     );
  addBranch(tree,(prefix_+"v0lambda_dxy").c_str()         ,&v0lambda_dxy_         ,(prefix_+"v0lambda_dxy_["+prefix_+"v0lambda_num_]/F").c_str()         );
  addBranch(tree,(prefix_+"v0lambda_dxyerr").c_str()      ,&v0lambda_dxyerr_      ,(prefix_+"v0lambda_dxyerr_["+prefix_+"v0lambda_num_]/F").c_str()      );
  addBranch(tree,(prefix_+"v0lambda_dxysig").c_str()      ,&v0lambda_dxysig_      ,(prefix_+"v0lambda_dxysig_["+prefix_+"v0lambda_num_]/F").c_str()      );
  addBranch(tree,(prefix_+"v0lambda_d3d").c_str()         ,&v0lambda_d3d_         ,(prefix_+"v0lambda_d3d_["+prefix_+"v0lambda_num_]/F").c_str()         );
  addBranch(tree,(prefix_+"v0lambda_d3derr").c_str()      ,&v0lambda_d3derr_      ,(prefix_+"v0lambda_d3err_["+prefix_+"v0lambda_num_]/F").c_str()       );
  addBranch(tree,(prefix_+"v0lambda_d3dsig").c_str()      ,&v0lambda_d3dsig_      ,(prefix_+"v0lambda_d3dsig_["+prefix_+"v0lambda_num_]/F").c_str()      );
  addBranch(tree,(prefix_+"v0lambda_costhetasvpv").c_str(),&v0lambda_costhetasvpv_,(prefix_+"v0lambda_costhetasvpv_["+prefix_+"v0lambda_num_]/F").c_str());
  addBranch(tree,(prefix_+"v0lambda_enratio").c_str()     ,&v0lambda_enratio_     ,(prefix_+"v0lambda_enratio_["+prefix_+"v0lambda_num_]/F").c_str());

}

void ntuple_V0lambda::readSetup(const edm::EventSetup& iSetup){

}

void ntuple_V0lambda::readEvent(const edm::Event& iEvent){


}


bool ntuple_V0lambda::compareDxyDxyErr(const reco::VertexCompositePtrCandidate &sva,const reco::VertexCompositePtrCandidate &svb){
    reco::Vertex pv = *spvp_;
    float adxy = ntuple_V0lambda::vertexDxy(sva,pv).value();
    float bdxy = ntuple_V0lambda::vertexDxy(svb,pv).value();
    float aerr = ntuple_V0lambda::vertexDxy(sva,pv).error();
    float berr = ntuple_V0lambda::vertexDxy(svb,pv).error();

    float asig = ntuple_V0lambda::catchInfs(adxy/aerr,0.);
    float bsig = ntuple_V0lambda::catchInfs(bdxy/berr,0.);
    return bsig<asig;
}

bool ntuple_V0lambda::fillBranches(const pat::Jet & jet, const size_t& jetidx, const  edm::View<pat::Jet> * coll){

    const float jet_uncorr_e = jet.correctedJet("Uncorrected").energy();
    const reco::Vertex & pv = vertices()->at(0);
    GlobalVector jetRefTrackDir(jet.px(),jet.py(),jet.pz());

    v0lambda_num_ = 0;

    reco::VertexCompositePtrCandidateCollection cpvtx=*V0lambda();
    spvp_ =   & vertices()->at(0);
    std::sort(cpvtx.begin(),cpvtx.end(),ntuple_V0lambda::compareDxyDxyErr);

    float etasign=1;
    etasign++; //avoid unused warning
    if(jet.eta()<0)etasign=-1;

    double jet_radius = jetR();
    if (jet_radius<0){
      // subjets: use maxDR(subjet, pfcand)
      for (unsigned idau=0; idau<jet.numberOfDaughters(); ++idau){
        double dR = reco::deltaR(*jet.daughter(idau), jet);
        if (dR>jet_radius)
          jet_radius = dR;
      }
    }

    for (const reco::VertexCompositePtrCandidate &sv : cpvtx) {

        if (reco::deltaR(sv,jet)>jet_radius) { continue; }
        if((int)max_sv>v0lambda_num_){

	  v0lambda_pt_[v0lambda_num_]           = sv.pt();
	  v0lambda_eta_[v0lambda_num_]          = sv.eta();
	  v0lambda_phi_[v0lambda_num_]          = sv.phi();
	  v0lambda_etarel_[v0lambda_num_]       = catchInfsAndBound(fabs(sv.eta()-jet.eta())-0.5,0,-2,0);
	  v0lambda_phirel_[v0lambda_num_]       = catchInfsAndBound(fabs(reco::deltaPhi(sv.phi(),jet.phi()))-0.5,0,-2,0);
	  v0lambda_deltaR_[v0lambda_num_]       = catchInfsAndBound(fabs(reco::deltaR(sv,jet))-0.5,0,-2,0);
	  v0lambda_mass_[v0lambda_num_]         = sv.mass();
	  v0lambda_ntracks_[v0lambda_num_]      = sv.numberOfDaughters();
	  v0lambda_chi2_[v0lambda_num_]         = sv.vertexChi2();
	  v0lambda_ndf_[v0lambda_num_]          = sv.vertexNdof();
	  v0lambda_normchi2_[v0lambda_num_]     = catchInfsAndBound(v0lambda_chi2_[v0lambda_num_]/v0lambda_ndf_[v0lambda_num_],1000,-1000,1000);
	  v0lambda_dxy_[v0lambda_num_]          = vertexDxy(sv,pv).value();
	  v0lambda_dxyerr_[v0lambda_num_]       = catchInfsAndBound(vertexDxy(sv,pv).error()-2,0,-2,0);
	  v0lambda_dxysig_[v0lambda_num_]       = catchInfsAndBound(v0lambda_dxy_[v0lambda_num_]/vertexDxy(sv,pv).error() ,0,-1,800);
	  v0lambda_d3d_[v0lambda_num_]          = vertexD3d(sv,pv).value();
	  v0lambda_d3derr_[v0lambda_num_]       = catchInfsAndBound(vertexD3d(sv,pv).error()-2,0,-2,0);
	  v0lambda_d3dsig_[v0lambda_num_]       = catchInfsAndBound(vertexD3d(sv,pv).value()/vertexD3d(sv,pv).error(),0,-1,800);
	  v0lambda_costhetasvpv_[v0lambda_num_] = vertexDdotP(sv,pv); // the pointing angle (i.e. the angle between the sum of the momentum                                                                                                                          
	  v0lambda_enratio_[v0lambda_num_]      = sv.energy()/jet_uncorr_e;
	  v0lambda_e_[v0lambda_num_]            = sv.energy();

            v0lambda_num_++;
        }
    } // end of looping over the secondary vertices
    nv0lambda_=v0lambda_num_;

    return true;
}

///helpers seldomly touched
Measurement1D ntuple_V0lambda::vertexDxy(const reco::VertexCompositePtrCandidate &svcand, const reco::Vertex &pv)  {
    VertexDistanceXY dist;
    reco::Vertex::CovarianceMatrix csv; svcand.fillVertexCovariance(csv);
    reco::Vertex svtx(svcand.vertex(), csv);
    return dist.distance(svtx, pv);
}

Measurement1D ntuple_V0lambda::vertexD3d(const reco::VertexCompositePtrCandidate &svcand, const reco::Vertex &pv)  {
    VertexDistance3D dist;
    reco::Vertex::CovarianceMatrix csv; svcand.fillVertexCovariance(csv);
    reco::Vertex svtx(svcand.vertex(), csv);
    return dist.distance(svtx, pv);
}

float ntuple_V0lambda::vertexDdotP(const reco::VertexCompositePtrCandidate &sv, const reco::Vertex &pv)  {
    reco::Candidate::Vector p = sv.momentum();
    reco::Candidate::Vector d(sv.vx() - pv.x(), sv.vy() - pv.y(), sv.vz() - pv.z());
    return p.Unit().Dot(d.Unit());
}
