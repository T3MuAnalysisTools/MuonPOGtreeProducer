import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

import subprocess

import sys

options = VarParsing.VarParsing()

options.register('globalTag',
                 '80X_mcRun2_asymptotic_2016_v3', #default value
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
                 './muonPOGNtuple_miniAOD_tracks.root', #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "Folder and name ame for output ntuple")

options.register('runOnMC',
                 False, #default value
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.bool,
                 "Run on DATA or MC")

options.parseArguments()


process = cms.Process("NTUPLES")

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
process.source.fileNames = ['file:/tmp/bjoshi/4FC506A2-247C-234F-A377-33939B80303C.root']

process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("MuonPOGtreeProducer.Tools.TrackCollectionProducer_cfi")

process.output = cms.OutputModule( "PoolOutputModule", outputCommands = cms.untracked.vstring('keep *'), fileName = cms.untracked.string("Skim.root"))
process.TrackCollection.PFCandidateTag = cms.untracked.InputTag("packedPFCandidates")
process.TrackCollection.LostTrackTag = cms.untracked.InputTag("lostTracks")
process.p = cms.Path(process.TrackCollection)
process.outputPath = cms.EndPath(process.output)
process.schedule = cms.Schedule(process.p, process.outputPath)
#process.source.fileNames = ['/store/mc/RunIIAutumn18MiniAOD/BuToJpsiK_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/90000/10311389-8843-0C4E-BB57-CD2003467E55.root']
