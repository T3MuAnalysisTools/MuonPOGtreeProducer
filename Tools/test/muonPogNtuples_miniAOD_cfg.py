import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import subprocess

import sys

options = VarParsing.VarParsing()

options.register('globalTag',
                 '102X_dataRun2_Prompt_v16', #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Global Tag")

options.register('nEvents',
                 10, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "Maximum number of processed events")

options.register('eosInputFolder',
                 '/store/relval/CMSSW_8_0_3/RelValZMM_13/MINIAODSIM/PU25ns_80X_mcRun2_asymptotic_2016_v3_gs71xNewGtHcalCust-v1/00000', #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "EOS folder with input files")

options.register('ntupleName',
                 './muonPOGNtuple_miniAOD_singleMuon_charmonium.root', #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Folder and name ame for output ntuple")

options.register('runOnMC',
                 False, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Run on DATA or MC")

options.register('hltPathFilter',
                 "all", #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Filter on paths (now only accepts all or IsoMu20)")

options.register('minMuPt',
                 0.0, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 "Skim the ntuple selecting only STA || TRK || GLB muons with pT > of this value")

options.register('minNMu',
                 1, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "number of TRK or GLB muons with pT > minMuPt to pass the skim")


options.parseArguments()

if options.hltPathFilter == "all" :
    pathCut   = "all"
    filterCut = "all"
elif options.hltPathFilter == "IsoMu20" :
    pathCut   = "HLT_IsoMu20_v"
    filterCut = "hltL3crIsoL1sMu16L1f0L2f10QL3f20QL3trkIsoFiltered0p09"
else :
    print "[" + sys.argv[0] + "]:", "hltPathFilter=", options.hltPathFilter, "is not a valid parameter!"
    sys.exit(100)

process = cms.Process("NTUPLES")
#process = cms.Process("Ntuples")

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.nEvents))

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = cms.string(options.globalTag)

process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(),
        secondaryFileNames = cms.untracked.vstring()
)

#files = subprocess.check_output([ "/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select", "ls", options.eosInputFolder ])
#process.source.fileNames = [ options.eosInputFolder+"/"+f for f in files.split() ]    

process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
#process.load("Geometry.CommonDetUnit.globalTrackingGeometry_cfi")
#process.load("RecoMuon.DetLayers.muonDetLayerGeometry_cfi")
process.load("MuonPOGtreeProducer.Tools.TrackCollectionProducer_cfi")
process.TrackCollection.PFCandidateTag = cms.untracked.InputTag("packedPFCandidates")
process.TrackCollection.LostTrackTag = cms.untracked.InputTag("lostTracks")

from MuonPOGtreeProducer.Tools.MuonPogNtuples_cff import appendMuonPogNtuple, customiseHlt, customiseMuonCuts
    
appendMuonPogNtuple(process,options.runOnMC,"HLT",options.ntupleName)
customiseHlt(process,pathCut,filterCut)
customiseMuonCuts(process,options.minMuPt,options.minNMu)

process.AOutput.replace(process.muonPogNtuple, process.TrackCollection + process.muonPogNtuple)

process.goodOfflinePrimaryVertices.src = cms.InputTag("offlineSlimmedPrimaryVertices")
process.MuonPogTree.MuonTag = cms.untracked.InputTag("slimmedMuons")
process.MuonPogTree.TrackTag = cms.untracked.InputTag("TrackCollection:pfTracks:NTUPLES")
process.MuonPogTree.PrimaryVertexTag = cms.untracked.InputTag("offlineSlimmedPrimaryVertices")
process.MuonPogTree.secondaryKsVertexTag = cms.untracked.InputTag("slimmedKshortVertices")
process.MuonPogTree.secondaryVertexTag = cms.untracked.InputTag("slimmedSecondaryVertices")
process.MuonPogTree.GenTag = cms.untracked.InputTag("prunedGenParticles")
process.MuonPogTree.TrigResultsTag = cms.untracked.InputTag("none")
process.MuonPogTree.TrigSummaryTag = cms.untracked.InputTag("none")
process.MuonPogTree.PFMetTag = cms.untracked.InputTag("none")
process.MuonPogTree.PFChMetTag = cms.untracked.InputTag("none")
process.MuonPogTree.CaloMetTag = cms.untracked.InputTag("none")
process.MuonPogTree.TriggerObjectTag = cms.untracked.InputTag("slimmedPatTrigger")
process.MuonPogTree.miniAODRun = cms.untracked.bool(True)
process.MuonPogTree.doKsVertices = cms.untracked.bool(True)
process.MuonPogTree.doPhiVertices = cms.untracked.bool(True)

#process.source.fileNames = ['/store/mc/RunIIAutumn18MiniAOD/BuToJpsiK_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/90000/10311389-8843-0C4E-BB57-CD2003467E55.root']
process.source.fileNames = ['file:/tmp/bjoshi/4FC506A2-247C-234F-A377-33939B80303C.root'] # double muon low mass
#process.source.fileNames = ['file:/tmp/bjoshi/0731BA23-725E-8540-BA91-D6B06D8D7826.root'] # single muon

#process.source.fileNames = ['file:./Skim.root']
