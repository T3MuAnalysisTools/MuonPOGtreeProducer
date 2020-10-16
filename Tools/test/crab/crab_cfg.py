from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.workArea = 'Production_10_06_2020' 
config.General.requestName = ''

config.section_('JobType')
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = '../muonPogNtuples_miniAOD_cfg.py'
#config.JobType.outputFiles = ['muonNTuple.root']
config.JobType.pyCfgParams = ['globalTag=102X_upgrade2018_realistic_v20',
                              'ntupleName=muonPOGNtuple.root',
                              'nEvents=-1',
                              'runOnMC=True',
                              'hltPathFilter=all',
                              'minMuPt=1.2',
                              'minNMu=1'
               ]

config.JobType.allowUndistributedCMSSW = True  # To fix cmssw releases

config.section_('Data')
config.Data.inputDataset = '/DsEtaMuNu_EtaMuMuGamma/wangjian-CRAB3_RunIIAutumn18DR_AODSIM-c5ffe086db9ac96902dde93cd44a5aa0/USER'
config.Data.outputDatasetTag = 'Production2018_19_05_2020_DsEtaMuNu_EtaMuMuGamma__wangjian-CRAB3_RunIIAutumn18DR_AODSIM-c5ffe086db9ac96902dde93cd44a5aa0'


config.Data.lumiMask = 'Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt'
#config.Data.runRange = '273299'

config.Data.splitting    = 'FileBased'
config.Data.inputDBS  = 'global'
config.Data.totalUnits = -1
config.Data.unitsPerJob  = 5  # Since files based, 5 files per job
config.Data.publication = True
config.Data.outLFNDirBase  = '/store/user/bjoshi'


#config.Data.ignoreLocality = True

config.section_('Site')
config.Site.storageSite = 'T2_US_Florida'
#config.Site.whitelist = ['T2_US_Wisconsin','T2_US_Purdue','T1_US_FNAL']
