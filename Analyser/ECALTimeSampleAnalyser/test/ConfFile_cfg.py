import FWCore.ParameterSet.Config as cms
#from Configuration.StandardSequences.Eras import eras
#process = cms.Process("EcalTimeSample", eras.Run2_2017)

from Configuration.Eras.Era_Run3_cff import Run3
process = cms.Process('Analyse',Run3)

process.load("FWCore.MessageService.MessageLogger_cfi")
#######################
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
#process.load('Configuration.StandardSequences.RawToDigi_cff')
#process.load('Configuration.StandardSequences.L1Reco_cff')
#process.load('Configuration.StandardSequences.Reconstruction_cff')
#process.load('Configuration.StandardSequences.RecoSim_cff')
#process.load('CommonTools.ParticleFlow.EITopPAG_cff')
#process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
#process.load('Configuration.StandardSequences.PATMC_cff')
#process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load("RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi")

#############################

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '106X_mc2017_realistic_v8', '') ###UL 2017 BH 
#process.GlobalTag = GlobalTag(process.GlobalTag, '105X_upgrade2018_realistic_IdealEcalIC_v4', '') ###UL 2017 BH 
#process.GlobalTag = GlobalTag(process.GlobalTag, '110X_mcRun3_2021_realistic_v6', '')  
#process.GlobalTag = GlobalTag(process.GlobalTag, '113X_mcRun3_2021_realistic_v10', '')  
#process.GlobalTag = GlobalTag(process.GlobalTag, '122X_mcRun3_2021_realistic_v9', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, '130X_dataRun3_PromptAnalysis_v1', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, '130X_dataRun3_Prompt_v4', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '133X_dataRun3_PromptAnalysis_v1', '')  

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5000) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )


process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
                                    #'/store/user/shilpi/V21_AODSIM_BEAM1ON_PU/BeamHalo_2017_Beam1ON/BeamHalo_AODSIM_V21_AODSIM_BEAM1ON_PU/200215_230732/0000/BH_3_1.root'
                                    #'/store/user/shilpi/V21_AODSIM_BEAM1ON_PU/BeamHalo_2017_Beam1ON/BeamHalo_AODSIM_V21_AODSIM_BEAM1ON_PU/200215_230732/0000/BH_3_907.root'
                                    #'/store/mc/RunIIWinter19PFCalibDR/DoubleElectron_FlatPt-1To300/AODSIM/2018ConditionsFlatPU0to70ECALGT_105X_upgrade2018_realistic_IdealEcalIC_v4-v1/40000/0038F7DC-41D2-7A4B-BE93-5F20220400C5.root'
                                    #'/store/mc/Run3Winter20DRMiniAOD/DoubleElectron_FlatPt-1To300/AODSIM/FlatPU0to80_110X_mcRun3_2021_realistic_v6-v3/10000/030B2A09-E305-7A48-AEE6-2BB967D727B9.root'
                                    #'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/030B2A09-E305-7A48-AEE6-2BB967D727B9.root'
                                    #'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/045196a3-6768-412f-bef9-96e07364fee6_alcareco.root'
                                    #'file:/eos/cms/store/user/shilpi/08c6dee4-08af-42f4-86aa-92857e4d2cfc.root'
                                    'file:/eos/cms/store/user/shilpi/49725c63-390d-4d17-bb29-ad75b3e2e658.root'
                                )
                            
                            #secondaryFileNames = cms.untracked.vstring(
                            #)


                        )


process.InterimOutput = cms.OutputModule("PoolOutputModule",
                                         fileName = cms.untracked.string('myOutputFile.root'), 
                                          SelectEvents = cms.untracked.PSet(
                                              SelectEvents = cms.vstring("p")
                                              ),
                                         outputCommands = cms.untracked.vstring('keep *')
                                         
                                     )

                                                        


process.TFileService = cms.Service("TFileService", fileName = cms.string('timeSampleTree.root'))

###https://cmssdt.cern.ch/lxr/source/RecoEcal/EgammaClusterProducers/python/ecalDigiSelector_cfi.py
process.timeSample = cms.EDAnalyzer('ECALTimeSampleAnalyser',
                                    #EBdigiCollection = cms.InputTag("selectDigi","selectedEcalEBDigiCollection"),
                                    #EEdigiCollection = cms.InputTag("selectDigi","selectedEcalEEDigiCollection"),
                                    #ebRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
                                    #eeRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
                                    
                                    EBdigiCollection = cms.InputTag("ecalDigis","ebDigis", "RECO"),
                                    EEdigiCollection = cms.InputTag("ecalDigis","eeDigis", "RECO"),
                                    ebRecHitCollection = cms.InputTag("ecalRecHit", "EcalRecHitsEB"),
                                    eeRecHitCollection = cms.InputTag("ecalRecHit", "EcalRecHitsEE"),

                                    ebRecHitWeightCollection = cms.InputTag("ecalRecHitWeight", "EcalRecHitsEBWeight"),
                                    eeRecHitWeightCollection = cms.InputTag("ecalRecHitWeight", "EcalRecHitsEEWeight"),
                                    rhoLabel   = cms.InputTag("fixedGridRhoFastjetAll"),
                                    VtxLabel             = cms.InputTag("offlinePrimaryVertices"),
                                    #electronSrc          = cms.InputTag("selectedPatElectrons"),
                                    #photonSrc            = cms.InputTag("selectedPatPhotons")
                                    #recoEleSrc            = cms.InputTag("gedGsfElectrons") ,
                                    recoEleSrc            = cms.InputTag("electronRecalibSCAssociator") ,
                                    genParticleSrc        = cms.InputTag("genParticles")
                                    #recoPhotonSrc            = cms.InputTag("gedPhotons") ,
                                )


#"ecalMultiFitUncalibRecHit"   "EcalUncalibRecHitsEB"
#"ecalWeightUncalibRecHit"   "EcalUncalibRecHitsEBWeight"
# "ecalRecHitWeight"          "EcalRecHitsEBWeight"
 
process.load('RecoLuminosity.LumiProducer.bunchSpacingProducer_cfi')
#process.load('Configuration.StandardSequences.Reconstruction_cff')
#process.load("RecoEcal.EgammaClusterProducers.reducedRecHitsSequence_cff")
#process.load("RecoLocalCalo.EcalRecProducers.ecalLocalRecoSequence_cff")
process.load("RecoLocalCalo.Configuration.ecalLocalRecoSequence_cff")
#process.load("RecoEcal.EgammaClusterProducers.particleFlowSuperClusterECAL_cfi")
process.load("RecoEcal.Configuration.RecoEcal_cff")
process.load("Calibration.EcalCalibAlgos.electronRecalibSCAssociator_cfi")
process.load("RecoParticleFlow.PFClusterProducer.particleFlowClusterECAL_cfi")
process.load("RecoParticleFlow.PFClusterProducer.particleFlowClusterECALUncorrected_cfi")
process.load("RecoParticleFlow.PFClusterProducer.particleFlowClusterPS_cfi")
process.load("RecoParticleFlow.PFClusterProducer.particleFlowRecHitECAL_cfi")
process.load("RecoParticleFlow.PFClusterProducer.particleFlowRecHitPS_cfi")

###electron-SC matching needs to be run since I am re-running the SC collection and hence when ele-->SC, it doesnot recognize. So I recreate another electron collection. SC algo needed to be run because it didnt exist and needed for ele-->SC. 
###https://cmssdt.cern.ch/lxr/source/Calibration/EcalCalibAlgos/src/ElectronRecalibSuperClusterAssociator.cc
### https://cmssdt.cern.ch/lxr/source/Calibration/EcalCalibAlgos/python/electronRecalibSCAssociator_cfi.py


###Need to run the weight method
#ECAL reconstruction
#from RecoLocalCalo.EcalRecProducers.ecalWeightUncalibRecHit_cfi import *
#from RecoLocalCalo.EcalRecProducers.ecalFixedAlphaBetaFitUncalibRecHit_cfi import *
#from RecoLocalCalo.EcalRecProducers.ecalRecHit_cff import *

process.load('RecoLocalCalo.EcalRecProducers.ecalWeightUncalibRecHit_cfi')

process.ecalWeightUncalibRecHit.EBhitCollection = cms.string('EcalUncalibRecHitsEBWeight')
process.ecalWeightUncalibRecHit.EEhitCollection = cms.string('EcalUncalibRecHitsEEWeight')


# get rechits e.g. from the weights
#process.load("CalibCalorimetry.EcalLaserCorrection.ecalLaserCorrectionService_cfi")
#process.load("RecoLocalCalo.EcalRecProducers.ecalRecHit_cfi")

process.ecalRecHitWeight = process.ecalRecHit.clone()

process.ecalRecHitWeight.cpu.EBuncalibRecHitCollection = cms.InputTag("ecalWeightUncalibRecHit", "EcalUncalibRecHitsEBWeight")#
process.ecalRecHitWeight.cpu.EEuncalibRecHitCollection = cms.InputTag("ecalWeightUncalibRecHit", "EcalUncalibRecHitsEEWeight")

process.ecalRecHitWeight.cpu.EBrechitCollection = cms.string('EcalRecHitsEBWeight')
process.ecalRecHitWeight.cpu.EErechitCollection = cms.string('EcalRecHitsEEWeight')


process.p = cms.Path(
    process.bunchSpacingProducer *
    process.ecalLocalRecoSequence *
    process.ecalWeightUncalibRecHit*
    process.ecalRecHitWeight * 
    process.ecalClusters *
    process.particleFlowRecHitPS *
    process.particleFlowRecHitECAL *
    process.particleFlowClusterECALUncorrected *
    process.particleFlowClusterPS *
    process.particleFlowClusterECAL *
    process.electronRecalibSCAssociator * 
    process.timeSample
)

#process.e = cms.EndPath(process.InterimOutput)  


#print(process.dumpPython())
