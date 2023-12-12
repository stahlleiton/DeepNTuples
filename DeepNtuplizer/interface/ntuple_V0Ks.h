/*
 * ntuple_V0Ks.h
 *
 *  Created on: 3rd August 2023
 *      Author: Alexandre De Moor
 */

#ifndef DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_V0Ks_H_
#define DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_V0Ks_H_

#include "ntuple_content.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

class ntuple_V0Ks: public ntuple_content{
public:

    ntuple_V0Ks(std::string prefix = "", double jetR = 0.4);
    ~ntuple_V0Ks();

    void getInput(const edm::ParameterSet& iConfig);
    void initBranches(TTree* );
    void readEvent(const edm::Event& iEvent);
    void readSetup(const edm::EventSetup& iSetup);

    //use either of these functions

    bool fillBranches(const pat::Jet &, const size_t& jetidx, const  edm::View<pat::Jet> * coll=0);



private:

    // SV candidates
    int   v0ks_num_;
    float nv0ks_;
    std::string prefix_;

    edm::ESHandle<TransientTrackBuilder> builder;

    static constexpr size_t max_sv=10;

    float v0ks_pt_[max_sv];
    float v0ks_px_[max_sv];
    float v0ks_py_[max_sv];
    float v0ks_pz_[max_sv];
    float v0ks_eta_[max_sv];
    float v0ks_phi_[max_sv];
    float v0ks_e_[max_sv];
    float v0ks_etarel_[max_sv];
    float v0ks_phirel_[max_sv];
    float v0ks_deltaR_[max_sv];
    float v0ks_mass_[max_sv];
    float v0ks_ntracks_[max_sv];
    float v0ks_chi2_[max_sv];
    float v0ks_ndf_[max_sv];
    float v0ks_normchi2_[max_sv];
    float v0ks_dxy_[max_sv];
    float v0ks_dxyerr_[max_sv];
    float v0ks_dxysig_[max_sv];
    float v0ks_d3d_[max_sv];
    float v0ks_d3derr_[max_sv];
    float v0ks_d3dsig_[max_sv];
    float v0ks_costhetasvpv_[max_sv];
    float v0ks_enratio_[max_sv];

    static const reco::Vertex * spvp_;

    static bool compareDxyDxyErr(const reco::VertexCompositePtrCandidate &sva,const reco::VertexCompositePtrCandidate &svb);

    //helper functions:
    static Measurement1D vertexDxy(const reco::VertexCompositePtrCandidate &svcand, const reco::Vertex &pv)  ;
    static Measurement1D vertexD3d(const reco::VertexCompositePtrCandidate &svcand, const reco::Vertex &pv)  ;
    static float vertexDdotP(const reco::VertexCompositePtrCandidate &sv, const reco::Vertex &pv)  ;

};

#endif /* DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_V0KS_H_ */
