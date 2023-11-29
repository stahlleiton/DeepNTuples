/*
 * ntuple_V0lambda.h
 *
 *  Created on: 3rd August 2023
 *      Author: Alexandre De Moor
 */

#ifndef DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_V0LAMBDA_H_
#define DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_V0LAMBDA_H_

#include "ntuple_content.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

class ntuple_V0lambda: public ntuple_content{
public:

    ntuple_V0lambda(std::string prefix = "", double jetR = 0.4);
    ~ntuple_V0lambda();

    void getInput(const edm::ParameterSet& iConfig);
    void initBranches(TTree* );
    void readEvent(const edm::Event& iEvent);
    void readSetup(const edm::EventSetup& iSetup);

    //use either of these functions

    bool fillBranches(const pat::Jet &, const size_t& jetidx, const  edm::View<pat::Jet> * coll=0);



private:

    // SV candidates
    int   v0lambda_num_;
    float nv0lambda_;
    std::string prefix_;

    edm::ESHandle<TransientTrackBuilder> builder;

    static constexpr size_t max_sv=10;

    float v0lambda_pt_[max_sv];
    float v0lambda_eta_[max_sv];
    float v0lambda_phi_[max_sv];
    float v0lambda_e_[max_sv];
    float v0lambda_etarel_[max_sv];
    float v0lambda_phirel_[max_sv];
    float v0lambda_deltaR_[max_sv];
    float v0lambda_mass_[max_sv];                                                                                                                                                                       
    float v0lambda_ntracks_[max_sv];
    float v0lambda_chi2_[max_sv];
    float v0lambda_ndf_[max_sv];
    float v0lambda_normchi2_[max_sv];
    float v0lambda_dxy_[max_sv];
    float v0lambda_dxyerr_[max_sv];
    float v0lambda_dxysig_[max_sv];
    float v0lambda_d3d_[max_sv];
    float v0lambda_d3derr_[max_sv];
    float v0lambda_d3dsig_[max_sv];
    float v0lambda_costhetasvpv_[max_sv];
    float v0lambda_enratio_[max_sv];

    static const reco::Vertex * spvp_;

    static bool compareDxyDxyErr(const reco::VertexCompositePtrCandidate &sva,const reco::VertexCompositePtrCandidate &svb);

    //helper functions:
    static Measurement1D vertexDxy(const reco::VertexCompositePtrCandidate &svcand, const reco::Vertex &pv)  ;
    static Measurement1D vertexD3d(const reco::VertexCompositePtrCandidate &svcand, const reco::Vertex &pv)  ;
    static float vertexDdotP(const reco::VertexCompositePtrCandidate &sv, const reco::Vertex &pv)  ;

};

#endif /* DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_V0LAMBDA_H_ */
