/*
 * ntuple_GraphB.h
 *
 *  Created on: 23 June 2017
 *      Author: Seth Moortgat
 */

#ifndef DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_GRAPHB_H_
#define DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_GRAPHB_H_

#include "ntuple_content.h"
#include "neighbourTrackVars.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "RecoBTag/TrackProbability/interface/HistogramProbabilityEstimator.h"

class HistogramProbabilityEstimator;
#include <typeinfo>

class ntuple_GraphB: public ntuple_content{
public:

    ntuple_GraphB(double jetR = 0.4);
    ~ntuple_GraphB();

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
    float jetradius_;
    float min_candidate_pt_ = -1;

    unsigned int n_Cpfcand_;
    float nCpfcand_;
    unsigned int n_Npfcand_;
    float nNpfcand_;

    static constexpr size_t max_gtracks=100;
    
    unsigned int n_gtracks=0;
    unsigned int n_match10=0;
    unsigned int n_match20=0;
    unsigned int n_match30=0;
    unsigned int n_match40=0;
    unsigned int n_match50=0;

    unsigned int n_gen_match10=0;
    unsigned int n_gen_match20=0;
    unsigned int n_gen_match30=0;
    unsigned int n_gen_match40=0;
    unsigned int n_gen_match50=0;

    float nGtracks=0;
    float n_daughters=0;
    float n_gen_daughters=0;
    float n_gen_constit=0;
    
    float gtrack_isMatched10[max_gtracks];
    float gtrack_isMatched20[max_gtracks];
    float gtrack_isMatched30[max_gtracks];
    float gtrack_isMatched40[max_gtracks];
    float gtrack_isMatched50[max_gtracks];
    float gtrack_pt[max_gtracks];
    float gtrack_eta[max_gtracks];
    float gtrack_phi[max_gtracks];
    float gtrack_mass[max_gtracks];
    float gtrack_dR[max_gtracks];
    float gtrack_dist_neigh[max_gtracks];
    float gtrack_dz[max_gtracks];
    float gtrack_dxy[max_gtracks];
    float gtrack_3D_ip[max_gtracks];
    float gtrack_3D_sip[max_gtracks];
    float gtrack_2D_ip[max_gtracks];
    float gtrack_2D_sip[max_gtracks];
    float gtrack_3D_TrackProbability[max_gtracks];
    float gtrack_2D_TrackProbability[max_gtracks];
    float gtrack_chi2reduced[max_gtracks];
    float gtrack_nPixelHits[max_gtracks];
    float gtrack_nHits[max_gtracks];
    float gtrack_jetAxisDistance[max_gtracks];
    float gtrack_jetAxisDlength[max_gtracks];
    float gtrack_dotProdTrack[max_gtracks];
    float gtrack_dotProdTrack2D[max_gtracks];
    float gtrack_PCAtrackFromPV[max_gtracks];

    
    // IVF cut parameters (HARDCODED?? OR CONFIGURABLE IN PYTHON CONFIG)
    float min3DIPValue=0.005;
    float min3DIPSignificance=1.2;
    int max3DIPValue=9999.;
    int max3DIPSignificance=9999.;
    

    //tokens to be defined from main analyzer
    edm::EDGetTokenT<edm::View<pat::PackedCandidate> > CandidateToken;

    //helper:
    edm::Handle<edm::View<pat::PackedCandidate> > tracks;
    
    // builder    
    edm::ESHandle<TransientTrackBuilder> builder;
        
    std::auto_ptr<HistogramProbabilityEstimator> m_probabilityEstimator;
    bool m_computeProbabilities=1;
    unsigned long long  m_calibrationCacheId2D; 
    unsigned long long m_calibrationCacheId3D;
    float mindrsvpfcand(const pat::PackedCandidate* pfcand);

};



#endif /* DEEPNTUPLES_DEEPNTUPLIZER_INTERFACE_NTUPLE_GRAPHB_H_ */
