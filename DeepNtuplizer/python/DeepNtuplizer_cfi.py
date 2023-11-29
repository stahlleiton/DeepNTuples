import FWCore.ParameterSet.Config as cms

deepntuplizer = cms.EDAnalyzer('DeepNtuplizer',
                               vertices   = cms.InputTag("offlineSlimmedPrimaryVertices"),
                               secVertices = cms.InputTag("slimmedSecondaryVertices"),
                               V0_ks = cms.InputTag("slimmedKshortVertices"),
                               V0_lambda = cms.InputTag("slimmedLambdaVertices"),
                               jets       = cms.InputTag("slimmedJetsPuppi"),
                               losttracks = cms.InputTag("lostTracks"),
                               packed   = cms.InputTag("packedGenParticles"),
                               jetR       = cms.double(0.4),
                               runFatJet = cms.bool(False),
                               eta = cms.bool(True),
                               puppi = cms.bool(True),
                               runDeepVertex = cms.bool(False),
                               pupInfo = cms.InputTag("slimmedAddPileupInfo"),
                               rhoInfo = cms.InputTag("fixedGridRhoFastjetAll"),	
                               SVs  = cms.InputTag("slimmedSecondaryVertices"),
                               LooseSVs = cms.InputTag("inclusiveCandidateSecondaryVertices"),
                               genJetMatchWithNu = cms.InputTag("patGenJetMatchWithNu"),
                               genJetMatchRecluster = cms.InputTag("patGenJetMatchRecluster"),
                               genJetMatchAllowDuplicates = cms.InputTag("patGenJetMatchAllowDuplicates"),
                               pruned = cms.InputTag("prunedGenParticles"),
                               fatjets = cms.InputTag('slimmedJetsAK8'),
                               muons = cms.InputTag("slimmedMuons"),
                               electrons = cms.InputTag("slimmedElectrons"),
                               jetPtMin     = cms.double(10.0),
                               jetPtMax     = cms.double(2000),
                               jetAbsEtaMin = cms.double(0.0),
                               jetAbsEtaMax = cms.double(5.0),
                               gluonReduction = cms.double(0.0),
                               tagInfoName = cms.string('deepNN'),
                               tagInfoFName = cms.string('pfBoostedDoubleSVAK8'),
                               bDiscriminators = cms.vstring(),
                               qgtagger        = cms.string("QGTagger"),
                               candidates      = cms.InputTag("packedPFCandidates"),
                               minCandidatePt  = cms.double(0.95),
                               useHerwigCompatible=cms.bool(False),
                               isHerwig=cms.bool(False),
                               useOffsets=cms.bool(True),
                               applySelection=cms.bool(True)
)
