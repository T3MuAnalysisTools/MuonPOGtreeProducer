import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
import subprocess
import sys

options = VarParsing.VarParsing()

options.register('globalTag',
                 '102X_dataRun2_Prompt_v16', #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Global Tag")

options.register('nEvents',
                 -1, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "Maximum number of processed events")

options.register('ntupleName',
                 './muonPOGNtuple_singleMuon.root', #default value
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
                 1.2, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 "Skim the ntuple saving only STA || TRK || GLB muons with pT > of this value")

options.register('minNMu',
                 0, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "number of TRK or GLB muons with pT > minMuPt to pass the skim")

options.parseArguments()

if options.hltPathFilter == "all" :
    pathCut   = "all"
    filterCut = "all"
elif options.hltPathFilter == "IsoMu20" :
    pathCut   = "HLT_IsoMu20_v"
    if options.runOnMC :
        filterCut = "hltL3crIsoL1sMu16L1f0L2f10QL3f20QL3trkIsoFiltered0p09"
    else :
        filterCut = "hltL3crIsoL1sMu18L1f0L2f10QL3f20QL3trkIsoFiltered0p09"
        
else :
    print "[" + sys.argv[0] + "]:", "hltPathFilter=", options.hltPathFilter, "is not a valid parameter!"
    sys.exit(100)

process = cms.Process("NTUPLES")

process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.nEvents))

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration.Geometry.GeometrySimDB_cff')
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.GlobalTag.globaltag = cms.string(options.globalTag)

process.source = cms.Source("PoolSource",
                            
        fileNames = cms.untracked.vstring(),
        secondaryFileNames = cms.untracked.vstring()

)

#files = subprocess.check_output([ "/afs/cern.ch/project/eos/installation/0.3.15/bin/eos.select", "ls", options.eosInputFolder ])
#process.source.fileNames = [ options.eosInputFolder+"/"+f for f in files.split() ]  
#process.source.fileNames = ['/store/mc/RunIIAutumn18DRPremix/BdToPiPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/102X_upgrade2018_realistic_v15-v2/1120000/2EB6F3EE-5988-8D42-9723-E6C882DEDD00.root']

process.source.fileNames = ['/store/data/Run2018D/SingleMuon/AOD/22Jan2019-v2/1110000/9EF1066D-C552-474C-AE57-AA53A6DD6218.root'] # single muon
from MuonPOGtreeProducer.Tools.MuonPogNtuples_cff import appendMuonPogNtuple, customiseHlt, customiseMuonCuts
    
appendMuonPogNtuple(process,options.runOnMC,"HLT",options.ntupleName)

customiseHlt(process,pathCut,filterCut)
customiseMuonCuts(process,options.minMuPt,options.minNMu)

process.MuonPogTree.miniAODRun = cms.untracked.bool(False)
process.MuonPogTree.doKsVertices = cms.untracked.bool(True)
process.MuonPogTree.doPhiVertices = cms.untracked.bool(True)
process.MuonPogTree.TriggerObjectTag = cms.untracked.InputTag("none")
process.MuonPogTree.TrackTag = cms.untracked.InputTag("generalTracks")
process.MuonPogTree.secondaryVertexTag = cms.untracked.InputTag("inclusiveCandidateSecondaryVertices")
process.MuonPogTree.secondaryKsVertexTag = cms.untracked.InputTag("generalV0Candidates:Kshort:RECO")
