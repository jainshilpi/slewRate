from CRABClient.UserUtilities import config
import sys

config = config()


#**************************submit function***********************                                                      
from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
'''
from httplib import HTTPException
def submit(config):
        try:
                crabCommand('submit', config = config)
        except HTTPException as hte:
                print("Failed submitting task: %s" % (hte.headers))
        except ClientException as cle:
                print("Failed submitting task: %s" % (cle))
'''
#****************************************************************                                




# Common configuration

config.General.workArea     = 'crab_projects_ntuples'
config.General.transferLogs = False
config.JobType.pluginName   = 'Analysis' # PrivateMC
#    config.JobType.inputFiles   = ['Summer16_23Sep2016V4_MC_L2Relative_AK8PFchs.txt', 'Summer16_23Sep2016V4_MC_L3Absolute_AK8PFchs.txt', 'Summer16_23Sep2016V4_MC.db']
config.JobType.sendExternalFolder = True
config.Data.inputDBS        = 'global'    
#config.Data.splitting       = 'LumiBased' # EventBased, FileBased, LumiBased (1 lumi ~= 300 events)
config.Data.splitting       = 'FileBased' # EventBased, FileBased, LumiBased (1 lumi ~= 300 events)
config.Data.totalUnits      = -1
config.Data.allowNonValidInputDataset      = True
config.Data.publication     = False
config.Site.storageSite     = 'T2_CH_CERN'
config.JobType.allowUndistributedCMSSW = True
#config.Data.useParent                   = True
config.Data.useParent                   = False
#config.Data.inputBlocks = '/DoubleElectron_FlatPt-1To300/Run3Winter20DRMiniAOD-FlatPU0to80_110X_mcRun3_2021_realistic_v6-v3/GEN-SIM-RAW#d5be51ad-1f40-4565-a8cd-fe0de08338f3'
# dataset dependent configuration

config.JobType.psetName     = 'ConfFile_cfg.py'
config.General.requestName = 'EGamma0_Run2023C'
config.Data.unitsPerJob    = 5

### to initiate the data transfer - 11th Dec, 2023
config.Data.inputDataset   = '/EGamma0/Run2023C-EcalUncalZElectron-PromptReco-v1/ALCARECO'
 
#config.Data.outLFNDirBase  = '/store/user/shilpi' ### first version here
config.Data.outLFNDirBase  = '/store/user/shilpi/slewRate'
#submit(config)

