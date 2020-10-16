import FWCore.ParameterSet.Config as cms

MuonPogTree = cms.EDAnalyzer("MuonPogTreeProducer",
                             TrigResultsTag = cms.untracked.InputTag("TriggerResults::HLT"),
                             TrigSummaryTag = cms.untracked.InputTag("hltTriggerSummaryAOD::HLT"),
                             TriggerObjectTag = cms.untracked.InputTag("none"),

                             TrigFilterCut = cms.untracked.string("all"),
                             TrigPathCut   = cms.untracked.string("all"),

                             MuonTag          = cms.untracked.InputTag("muons"),
                             TrackTag         = cms.untracked.InputTag("generalTracks"),
                             simMuonTag       = cms.untracked.InputTag("muonSimClassifier"),
                             
                             PrimaryVertexTag = cms.untracked.InputTag("offlinePrimaryVertices"),
                             secondaryKsVertexTag = cms.untracked.InputTag("none"),
                             secondaryVertexTag = cms.untracked.InputTag("none"),
                             BeamSpotTag      = cms.untracked.InputTag("offlineBeamSpot"),
                             
                             PFMetTag         = cms.untracked.InputTag("pfMet"), 
                             PFChMetTag       = cms.untracked.InputTag("pfChMet"), 
                             CaloMetTag       = cms.untracked.InputTag("caloMet"),
                             ScalersTag = cms.untracked.InputTag("scalersRawToDigi"),
                             l1MuonsTag = cms.untracked.InputTag("gmtStage2Digis:Muon:"),

                             GenTag = cms.untracked.InputTag("genParticles"), # pruned
                             PileUpInfoTag = cms.untracked.InputTag("addPileupInfo"),
                             GenInfoTag = cms.untracked.InputTag("generator"),
                             MinMuPtCut = cms.untracked.double(0.),
                             MinNMuCut  = cms.untracked.int32(0),
                             miniAODRun = cms.untracked.bool(False),
                             fillKsVertices = cms.untracked.bool(False),
                             fillPhiVertices = cms.untracked.bool(False)
                             )
