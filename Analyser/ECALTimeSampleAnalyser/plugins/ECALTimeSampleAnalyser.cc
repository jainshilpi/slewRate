// -*- C++ -*-
//
// Package:    Analyser/ECALTimeSampleAnalyser
// Class:      ECALTimeSampleAnalyser
//
/**\class ECALTimeSampleAnalyser ECALTimeSampleAnalyser.cc Analyser/ECALTimeSampleAnalyser/plugins/ECALTimeSampleAnalyser.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Shilpi Jain
//         Created:  Fri, 12 May 2023 19:29:28 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalRawData/interface/EcalRawDataCollections.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "CondFormats/EcalObjects/interface/EcalGainRatios.h"
#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
#include "CondFormats/DataRecord/interface/EcalGainRatiosRcd.h"
#include "CondFormats/DataRecord/interface/EcalPedestalsRcd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "CondFormats/EcalObjects/interface/EcalPFRecHitThresholds.h"
#include "CondFormats/DataRecord/interface/EcalPFRecHitThresholdsRcd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
//#include "RecoEgamma/EgammaTools/interface/HGCalEgammaIDHelper.h"

#include "RecoLocalCalo/EcalRecAlgos/interface/EigenMatrixTypes.h"

#include "CondFormats/DataRecord/interface/EcalPulseShapesRcd.h"
#include "CondFormats/EcalObjects/interface/EcalPulseShapes.h"

#include <TTree.h>
#include <TLorentzVector.h>

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


using reco::TrackCollection;

class ECALTimeSampleAnalyser : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit ECALTimeSampleAnalyser(const edm::ParameterSet&);
      ~ECALTimeSampleAnalyser();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      static const int NMAXSAMPLES = 10;
      static const int NMAXCRYS = 100;
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
  //virtual std::vector<std::array<double, NMAXSAMPLES> > getTimeSamplesAroundEle(std::vector<DetId> v_id, edm::Handle<EBDigiCollection> pEBDigi, edm::Handle<EEDigiCollection> pEEDigi);
  //virtual std::vector<std::vector<double> > getTimeSamplesAroundEle(std::vector<DetId> v_id, edm::Handle<EBDigiCollection> pEBDigi, edm::Handle<EEDigiCollection> pEEDigi, const EcalRecHitCollection* EBRecHits, const EcalRecHitCollection* EERecHits,  const EcalPFRecHitThresholds* thresholds, std::vector<double> &hitsEnergy, std::vector<double> &hitsThr);
  virtual std::vector<std::vector<double> > getTimeSamplesAroundEle(std::vector<DetId> v_id, edm::Handle<EBDigiCollection> pEBDigi, edm::Handle<EEDigiCollection> pEEDigi, const EcalRecHitCollection* EBRecHits, const EcalRecHitCollection* EERecHits,  const EcalRecHitCollection* EBRecHitsWeight, const EcalRecHitCollection* EERecHitsWeight, const EBUncalibratedRecHitCollection* EBUncalibRecHit, const EEUncalibratedRecHitCollection* EEUncalibRecHit, std::vector<double> &hitsEnergy, std::vector<double> &hitsEnergyWeight, std::vector<double> &hitsAmp,  std::vector<double> &hitsAmpError, std::vector<std::vector<double>> &hitsOOAmp, std::vector<std::vector<double>> &hitsPUSubSamples);

  virtual std::vector<reco::GenParticle>::const_iterator  getGenMatch(std::vector<std::vector<reco::GenParticle>::const_iterator> genLep, reco::GsfElectron gsfele, double &dRmin);
  
  
      // ----------member data ---------------------------
      edm::ESGetToken<CaloTopology, CaloTopologyRecord> caloTopoToken_;
      edm::ESHandle<CaloTopology> caloTopology_;
      edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
      edm::EDGetTokenT<EBDigiCollection> ebDigiToken_;
      edm::EDGetTokenT<EEDigiCollection> eeDigiToken_;
      const std::string ebdigiCollection_;
  //edm::EDGetTokenT<reco::PhotonCollection>         recophotonCollection_;
      edm::EDGetTokenT<reco::GsfElectronCollection>         recoelectronCollection_;
  
      edm::EDGetTokenT<EcalRecHitCollection>           ebRecHitCollection_;
      edm::EDGetTokenT<EcalRecHitCollection>           eeRecHitCollection_;

      edm::EDGetTokenT<EcalRecHitCollection>           ebRecHitCollectionWeight_;
      edm::EDGetTokenT<EcalRecHitCollection>           eeRecHitCollectionWeight_;
  
  edm::EDGetTokenT<EBUncalibratedRecHitCollection> EBuncalibrechitCollection_;
  edm::EDGetTokenT<EEUncalibratedRecHitCollection> EEuncalibrechitCollection_;
  
  //edm::EDGetTokenT<EcalRecHitCollection>           esRecHitCollection_;
  
      edm::ESHandle<EcalPedestals> peds;
      edm::ESGetToken<EcalPedestals, EcalPedestalsRcd> pedsToken_;
      edm::ESHandle<EcalGainRatios> gains;
      edm::ESGetToken<EcalGainRatios, EcalGainRatiosRcd> gainsToken_;
  edm::EDGetTokenT<reco::VertexCollection>         vtxLabel_;
  edm::EDGetTokenT<double>                         rhoLabel_;

   edm::ESHandle<EcalPulseShapes> pulseshapes;
   edm::ESGetToken<EcalPulseShapes, EcalPulseShapesRcd> pulseShapesToken_;
  
  //std::vector<double> m_EBPulseShapeTemplate, m_EEPulseShapeTemplate;
   std::vector<double> ebPulseShape_;
  std::vector<double> eePulseShape_;

  std::vector<int> activeBX;
  
  //const std::string ebdigiProducer_;

  TTree   *treeEle;
  TTree   *treePho;

  Int_t       run_;
  Long64_t    event_;
  Int_t       lumis_;
  std::vector<double> eleEta_;
  std::vector<double> elePhi_;
  std::vector<double> elePt_;
  std::vector<double> eleE_;
  std::vector<double> eleSCEta_;
  std::vector<double> eleSCPhi_;
  std::vector<float>  eleHoverE_;
  std::vector<float>  eleSigmaIEtaIEtaFull5x5_;
  std::vector<float>  eleSigmaIPhiIPhiFull5x5_;
  std::vector<float>  eleE1x5Full5x5_;
  std::vector<float>  eleE2x5Full5x5_;
  std::vector<float>  eleE5x5Full5x5_;
  std::vector<float>  eleR9Full5x5_;
  std::vector<float>  eledEtaAtVtx_;
  std::vector<float>  eledPhiAtVtx_;
  std::vector<float>  eledEtaAtCalo_;
  std::vector<float>  eleEoverP_;
  std::vector<float>  eleEoverPout_;
  std::vector<float>  eleEoverPInv_;
  std::vector<float>  eleBrem_;
    
  
  std::vector<float>  elePFChIso_;
  std::vector<float>  elePFPhoIso_;
  std::vector<float>  elePFNeuIso_;
  std::vector<float>  elePFPUIso_;

  
  std::vector<double> e5x5_; 
  std::vector<int> nCrys_;
  int  nsamples_;
  //double genPt_, genEta_, genPhi_, gendR_, genE_, genStatus_;
  //std::vector<std::array<double, NMAXSAMPLES> > hitsAmplitudes_;
  //std::vector<std::array<double, 10> > hitsAmplitudes_;
  //double hitsAmplitudes_[NMAXCRYS][NMAXSAMPLES];
  std::vector<std::vector<std::vector<double>>> hitsAmplitudes_;
  //std::vector<std::map<int,std::vector<double>>> hitsAmplitudes_;
  std::vector<std::vector<double>> hitsEnergy_;
  std::vector<std::vector<double>> hitsEnergyWeight_;
  std::vector<std::vector<double>> hitsMFAmp_;
  std::vector<std::vector<double>> hitsMFAmpErr_;
  std::vector<std::vector<std::vector<double>>> hitsMFOOTAmp_;
  std::vector<std::vector<std::vector<double>>> hitsAmplitudesPUSub_;
  
  Int_t       nVtx_;
  float       rho_;
  
  //std::vector<double> hitsThr_;

  //edm::EDGetTokenT<std::vector<reco::GenParticle> > genParticlesCollection_;



};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
ECALTimeSampleAnalyser::ECALTimeSampleAnalyser(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed
  ebDigiToken_ = consumes<EBDigiCollection>(iConfig.getParameter<edm::InputTag>("EBdigiCollection"));
  eeDigiToken_ = consumes<EEDigiCollection>(iConfig.getParameter<edm::InputTag>("EEdigiCollection"));
  //recophotonCollection_       = consumes<reco::PhotonCollection>        (iConfig.getParameter<edm::InputTag>("recoPhotonSrc"));
  recoelectronCollection_       = consumes<reco::GsfElectronCollection>        (iConfig.getParameter<edm::InputTag>("recoEleSrc"));
  ebRecHitCollection_ = consumes<EcalRecHitCollection>          (iConfig.getParameter<edm::InputTag>("ebRecHitCollection"));
  eeRecHitCollection_ = consumes<EcalRecHitCollection>          (iConfig.getParameter<edm::InputTag>("eeRecHitCollection"));

  EBuncalibrechitCollection_ = consumes<EBUncalibratedRecHitCollection>          (iConfig.getParameter<edm::InputTag>("EBuncalibrechitCollection"));
  EEuncalibrechitCollection_ = consumes<EBUncalibratedRecHitCollection>          (iConfig.getParameter<edm::InputTag>("EEuncalibrechitCollection"));

  
  ebRecHitCollectionWeight_ = consumes<EcalRecHitCollection>          (iConfig.getParameter<edm::InputTag>("ebRecHitWeightCollection"));
  eeRecHitCollectionWeight_ = consumes<EcalRecHitCollection>          (iConfig.getParameter<edm::InputTag>("eeRecHitWeightCollection"));

  rhoLabel_                  = consumes<double>(iConfig.getParameter<edm::InputTag>("rhoLabel"));
  vtxLabel_                  = consumes<reco::VertexCollection>        (iConfig.getParameter<edm::InputTag>("VtxLabel"));
  
  //esRecHitCollection_ = consumes<EcalRecHitCollection>          (iConfig.getParameter<edm::InputTag>("esRecHitCollection"));
  //m_EBPulseShapeTemplate = iConfig.getParameter<std::vector<double> >("EBPulseShapeTemplate");
  //m_EEPulseShapeTemplate = iConfig.getParameter<std::vector<double> >("EEPulseShapeTemplate");
   pulseShapesToken_ = esConsumes<EcalPulseShapes, EcalPulseShapesRcd>();

  caloTopoToken_ = esConsumes();
  pedsToken_ = esConsumes<EcalPedestals, EcalPedestalsRcd>();
  gainsToken_ = esConsumes<EcalGainRatios, EcalGainRatiosRcd>();

  //activeBX.resize(10);
  activeBX.push_back(-5);
  activeBX.push_back(-4);
  activeBX.push_back(-3);
  activeBX.push_back(-2);
  activeBX.push_back(-1);
  activeBX.push_back(0);
  activeBX.push_back(1);
  activeBX.push_back(2);
  activeBX.push_back(3);
  activeBX.push_back(4);
  //genParticlesCollection_   = consumes<std::vector<reco::GenParticle> >    (iConfig.getParameter<edm::InputTag>("genParticleSrc"));
  
}


ECALTimeSampleAnalyser::~ECALTimeSampleAnalyser()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ECALTimeSampleAnalyser::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  ///https://cmssdt.cern.ch/lxr/source/CalibCalorimetry/EcalLaserAnalyzer/plugins/EcalTestPulseAnalyzer.cc
  /// https://cmssdt.cern.ch/lxr/source/RecoEcal/EgammaClusterProducers/src/EcalDigiSelector.cc ---> selected digis
  ////seems to save 3x3 around the max hit only: https://cmssdt.cern.ch/lxr/source/RecoEcal/EgammaClusterProducers/src/EcalDigiSelector.cc#0164
  ///pedestal subtraction ---> 
  ///https://cmssdt.cern.ch/lxr/source/RecoLocalCalo/EcalRecAlgos/src/EcalUncalibRecHitMultiFitAlgo.cc#0064

  ///https://cmssdt.cern.ch/lxr/source/Validation/EcalRecHits/python/ecalRecHitsValidation_cfi.py

  ///ECAL Pulse shape: https://cmssdt.cern.ch/lxr/source/CondTools/Ecal/src/EcalPulseShapesHandler.cc#0018

  ///https://cmssdt.cern.ch/lxr/source/RecoLocalCalo/EcalRecProducers/plugins/EcalUncalibRecHitWorkerMultiFit.cc#0440
  
   using namespace edm;
   run_    = iEvent.id().run();
   event_  = iEvent.id().event();
   lumis_  = iEvent.luminosityBlock();
   
   eleEta_.clear();
   elePhi_.clear();

   eleSCEta_.clear();
   eleSCPhi_.clear();
   
   elePt_.clear();
   eleE_.clear();

   e5x5_.clear();
   nCrys_.clear();

   eleHoverE_                  .clear();
   eleEoverP_                  .clear();
   eleEoverPout_               .clear();
   eleEoverPInv_               .clear();
   eleBrem_                    .clear();
   eledEtaAtVtx_               .clear();
   eledPhiAtVtx_               .clear();
   eledEtaAtCalo_              .clear();
   eleSigmaIEtaIEtaFull5x5_    .clear();
   eleSigmaIPhiIPhiFull5x5_    .clear();
   elePFChIso_                 .clear();
   elePFPhoIso_                .clear();
   elePFNeuIso_                .clear();
   elePFPUIso_                 .clear();
   eleE1x5Full5x5_             .clear();
   eleE2x5Full5x5_             .clear();
   eleE5x5Full5x5_             .clear();
   eleR9Full5x5_               .clear();
   
   
   // Iterate over each element of the outer vector
   for (auto& outerVec : hitsAmplitudes_) {
     // Iterate over each element of the middle vector
     for (auto& middleVec : outerVec) {
       // Clear each inner vector
       middleVec.clear();
     }
     // Clear the middle vector
     outerVec.clear();
   }
   // Clear the outer vector
   hitsAmplitudes_.clear();


   // Iterate over each element of the outer vector
   for (auto& outerVec : hitsMFOOTAmp_) {
     // Iterate over each element of the middle vector
     for (auto& middleVec : outerVec) {
       // Clear each inner vector
       middleVec.clear();
     }
     // Clear the middle vector
     outerVec.clear();
   }
   // Clear the outer vector
   hitsMFOOTAmp_.clear();

   // Iterate over each element of the outer vector
   for (auto& outerVec : hitsAmplitudesPUSub_) {
     // Iterate over each element of the middle vector
     for (auto& middleVec : outerVec) {
       // Clear each inner vector
       middleVec.clear();
     }
     // Clear the middle vector
     outerVec.clear();
   }
   // Clear the outer vector
   hitsAmplitudesPUSub_.clear();


   
   
 
   
   hitsEnergy_.clear();
   // Iterate through each inner vector and clear them
   for (auto& innerVec : hitsEnergy_) {
     innerVec.clear();
   }
   
   hitsEnergyWeight_.clear();
   // Iterate through each inner vector and clear them
   for (auto& innerVec : hitsEnergyWeight_) {
     innerVec.clear();
   }

   hitsMFAmp_.clear();
   // Iterate through each inner vector and clear them
   for (auto& innerVec : hitsMFAmp_) {
     innerVec.clear();
   }

   hitsMFAmpErr_.clear();
   // Iterate through each inner vector and clear them
   for (auto& innerVec : hitsMFAmpErr_) {
     innerVec.clear();
   }

   
   
   //hitsThr_.clear();

   //hitsAmplitudes_.clear();
   nsamples_ = NMAXSAMPLES;

   //noZS::EcalClusterLazyTools lazyToolnoZS(iEvent, iSetup, ebRecHitCollection_, eeRecHitCollection_, esRecHitCollection_);
   
   caloTopology_ = iSetup.getHandle(caloTopoToken_);
   gains = iSetup.getHandle(gainsToken_);
   peds = iSetup.getHandle(pedsToken_);


   ///
   edm::Handle<double> rhoHandle;
   iEvent.getByToken(rhoLabel_, rhoHandle);
  
   rho_    = *(rhoHandle.product());
   
   edm::Handle<reco::VertexCollection> vtxHandle;
   iEvent.getByToken(vtxLabel_, vtxHandle);
   
   nVtx_     = -1;
   if (vtxHandle.isValid()) {
     nVtx_     = 0;
     
     for (std::vector<reco::Vertex>::const_iterator v = vtxHandle->begin(); v != vtxHandle->end(); ++v) {

       /*
       if (nVtx_ == 0) {
	 vtx_     = v->x();
	 vty_     = v->y();
	 vtz_     = v->z();
	 
	 isPVGood_ = false;
	 if (!v->isFake() && v->ndof() > 4. && fabs(v->z()) <= 24. && fabs(v->position().rho()) <= 2.) isPVGood_ = true;
       }
       */
       //if (!v->isFake() && v->ndof() > 4. && fabs(v->z()) <= 24. && fabs(v->position().rho()) <= 2.) nGoodVtx_++;
       nVtx_++;
     }
   } else
     edm::LogWarning("ggNtuplizer") << "Primary vertices info not unavailable";


   ///ECAL pulse shapes
   pulseshapes = iSetup.getHandle(pulseShapesToken_);

   ///ECAL multifit uncalib rechits
   const EBUncalibratedRecHitCollection *EBUncalibRecHit = nullptr;
   Handle<EBUncalibratedRecHitCollection> EcalUncalibRecHitEB;
   iEvent.getByToken(EBuncalibrechitCollection_, EcalUncalibRecHitEB);
   if (EcalUncalibRecHitEB.isValid()) {
     EBUncalibRecHit = EcalUncalibRecHitEB.product();
   } else {
     std::cout<<"EBUncalibRecHit collection now found!!!!"<<std::endl;
   }
 
   const EEUncalibratedRecHitCollection *EEUncalibRecHit = nullptr;
   Handle<EEUncalibratedRecHitCollection> EcalUncalibRecHitEE;
   iEvent.getByToken(EEuncalibrechitCollection_, EcalUncalibRecHitEE);
   if (EcalUncalibRecHitEE.isValid()) {
     EEUncalibRecHit = EcalUncalibRecHitEE.product();
   } else {
     std::cout<<"EEUncalibRecHit collection now found!!!!"<<std::endl;
   }


   
   //// nominal ECAL rechits from multifit method
   edm::Handle<EcalRecHitCollection> barrelRecHitsHandle;
   edm::Handle<EcalRecHitCollection> endcapRecHitsHandle;
   //edm::Handle<EcalRecHitCollection> esRecHitsHandle;
   iEvent.getByToken(ebRecHitCollection_,barrelRecHitsHandle);
   iEvent.getByToken(eeRecHitCollection_,endcapRecHitsHandle);
   //iEvent.getByToken(esRecHitCollection_,esRecHitsHandle);

   const EcalRecHitCollection* EBRecHits = nullptr;
   const EcalRecHitCollection* EERecHits = nullptr;
   //const EcalRecHitCollection* ESRecHits = nullptr;
  
   if ( !barrelRecHitsHandle.isValid() ){
     LogDebug("") << "Error! EB rechits can't get product!" << std::endl;
   } else{
     EBRecHits = barrelRecHitsHandle.product();
   }

   if ( !endcapRecHitsHandle.isValid() ){
     LogDebug("") << "Error! EE rechits can't get product!" << std::endl;
   } else{
     EERecHits = endcapRecHitsHandle.product();
   }




   //// ECAL rechits from weight method
   edm::Handle<EcalRecHitCollection> barrelRecHitsWeightHandle;
   edm::Handle<EcalRecHitCollection> endcapRecHitsWeightHandle;
   iEvent.getByToken(ebRecHitCollectionWeight_,barrelRecHitsWeightHandle);
   iEvent.getByToken(eeRecHitCollectionWeight_,endcapRecHitsWeightHandle);

   const EcalRecHitCollection* EBRecHitsWeight = nullptr;
   const EcalRecHitCollection* EERecHitsWeight = nullptr;
  
   if ( !barrelRecHitsWeightHandle.isValid() ){
     LogDebug("") << "Error! EB weight method rechits can't get product!" << std::endl;
   } else{
     EBRecHitsWeight = barrelRecHitsWeightHandle.product();
   }

   if ( !endcapRecHitsWeightHandle.isValid() ){
     LogDebug("") << "Error! EE weight rechits can't get product!" << std::endl;
   } else{
     EERecHitsWeight = endcapRecHitsWeightHandle.product();
   }
   
   
   /*if ( !esRecHitsHandle.isValid() ){
     LogDebug("") << "Error! ES rechits can't get product!" << std::endl;
   } else{
     ESRecHits = esRecHitsHandle.product();
   }
   */

   // retrieving crystal data from Event
   edm::Handle<EBDigiCollection> pEBDigi;
   edm::Handle<EEDigiCollection> pEEDigi;

   
   iEvent.getByToken(ebDigiToken_, pEBDigi);
   iEvent.getByToken(eeDigiToken_, pEEDigi);
   
   if (!pEBDigi.isValid()) {
     std::cout<<"Error! can't get the product retrieving EB crystal data, i.e. EBDigiCollection " <<std::endl;
   } 
   
   if (!pEEDigi.isValid()) {
     std::cout<<"Error! can't get the product retrieving EE crystal data, i.e. EEDigiCollection " <<std::endl;
   }

   

   /*//https://github.com/swagata87/OldLocalCovMiniAOD/blob/main/plugins/OldLocalCovMiniAOD.cc
   edm::ESHandle<EcalPFRecHitThresholds> pThresholds;
   iSetup.get<EcalPFRecHitThresholdsRcd>().get(pThresholds);
   const EcalPFRecHitThresholds* thresholds = pThresholds.product();


   /// gen collection
   edm::Handle<std::vector<reco::GenParticle> > genParticlesHandle;
   iEvent.getByToken(genParticlesCollection_, genParticlesHandle);

   if (!genParticlesHandle.isValid()) {
     edm::LogWarning("ggNtuplizer") << "no reco::GenParticles in event";
     return;
   }
   
     std::vector<std::vector<reco::GenParticle>::const_iterator> genLep;

     
     for (std::vector<reco::GenParticle>::const_iterator ip = genParticlesHandle->begin(); ip != genParticlesHandle->end(); ++ip) {
     

     

     Int_t status = ip->status();
     bool photonOrLepton =
       (ip->pdgId() == 22 && (ip->isPromptFinalState() || ip->isLastCopy() || status == 1)) ||
       (status == 1 && abs(ip->pdgId()) == 11 && (ip->isPromptFinalState() || ip->isLastCopy())) ||
       (status == 1 && abs(ip->pdgId()) == 13 && (ip->isPromptFinalState() || ip->isLastCopy())) ||
       (status == 1 && (abs(ip->pdgId()) == 12 || abs(ip->pdgId()) == 14 || abs(ip->pdgId()) == 16)) ||
       (status == 1 && ( abs(ip->pdgId()) >= 11 && abs(ip->pdgId()) <= 16 ) && ip->pt() > 3.0)  ||
       (status < 10 && abs(ip->pdgId()) == 15 && ip->pt() > 3.0);

     bool isEle = ( abs(ip->pdgId()) == 11 && (ip->isPromptFinalState() || ip->isLastCopy()) );
     
     if(photonOrLepton){
     genLep.push_back(ip);
     }
     
     }
   */
   
   ///grab electron collection
   Handle<reco::GsfElectronCollection> theRecoEleCollection;
   iEvent.getByToken(recoelectronCollection_, theRecoEleCollection);
   const reco::GsfElectronCollection theRecoEl = *(theRecoEleCollection.product());

   if (theRecoEleCollection.isValid()) {
     for (uint j = 0; j < theRecoEl.size(); j++) {
       
       DetId seedDetId = (theRecoEl[j].superCluster()->seed()->hitsAndFractions())[0].first;
       bool isBarrel = (seedDetId.subdetId() == EcalBarrel);
       eleE_.push_back(theRecoEl[j].energy());
       elePt_.push_back(theRecoEl[j].pt());
       eleEta_.push_back(theRecoEl[j].eta());
       elePhi_.push_back(theRecoEl[j].phi());
       eleSCEta_.push_back(theRecoEl[j].superCluster()->eta());
       eleSCPhi_.push_back(theRecoEl[j].superCluster()->phi());

       eleSigmaIEtaIEtaFull5x5_    .push_back(theRecoEl[j].full5x5_sigmaIetaIeta());
       eleSigmaIPhiIPhiFull5x5_    .push_back(theRecoEl[j].full5x5_sigmaIphiIphi());
       eleE1x5Full5x5_             .push_back(theRecoEl[j].full5x5_e1x5());
       eleE2x5Full5x5_             .push_back(theRecoEl[j].full5x5_e2x5Max());
       eleE5x5Full5x5_             .push_back(theRecoEl[j].full5x5_e5x5());
       eleR9Full5x5_               .push_back(theRecoEl[j].full5x5_r9());

       eleHoverE_          .push_back(theRecoEl[j].hcalOverEcal());
       eledEtaAtVtx_       .push_back(theRecoEl[j].deltaEtaSuperClusterTrackAtVtx());
       eledPhiAtVtx_       .push_back(theRecoEl[j].deltaPhiSuperClusterTrackAtVtx());
       eledEtaAtCalo_      .push_back(theRecoEl[j].deltaEtaSeedClusterTrackAtCalo());

       eleEoverP_          .push_back(theRecoEl[j].eSuperClusterOverP());
       eleEoverPout_       .push_back(theRecoEl[j].eEleClusterOverPout());
       eleBrem_            .push_back(theRecoEl[j].fbrem());
    
       if (theRecoEl[j].ecalEnergy() == 0)   eleEoverPInv_.push_back(1e30);
       else if (!std::isfinite(theRecoEl[j].ecalEnergy()))  eleEoverPInv_.push_back(1e30);
       else  eleEoverPInv_.push_back((1.0 - theRecoEl[j].eSuperClusterOverP())/theRecoEl[j].ecalEnergy());
       
       reco::GsfElectron::PflowIsolationVariables pfIso = theRecoEl[j].pfIsolationVariables();
       elePFChIso_         .push_back(pfIso.sumChargedHadronPt);
       elePFPhoIso_        .push_back(pfIso.sumPhotonEt);
       elePFNeuIso_        .push_back(pfIso.sumNeutralHadronEt);
	elePFPUIso_         .push_back(pfIso.sumPUPt);
	
	
       ///find the matrix of crystals in 5x5 array around the crystal
       //https://cmssdt.cern.ch/lxr/source/RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h#0869
       std::vector<DetId> v_id = noZS::EcalClusterTools::matrixDetId(caloTopology_.product(), seedDetId, 2);
       nCrys_.push_back((int) v_id.size());

       std::vector<double> hitsEnergy;
       std::vector<double> hitsEnergyWeight;
       std::vector<double> hitsAmp;
       std::vector<double> hitsAmpError;
       std::vector<std::vector<double>> hitsOOAmp;
       std::vector<std::vector<double>> hitsPUSubSamples;
       //std::vector<double> hitsThr;

       //std::vector<std::vector<double>> hitsAmplitudes = getTimeSamplesAroundEle(v_id, pEBDigi, pEEDigi, EBRecHits, EERecHits, thresholds, hitsEnergy, hitsThr);
       std::vector<std::vector<double>> hitsAmplitudes = getTimeSamplesAroundEle(v_id, pEBDigi, pEEDigi, EBRecHits, EERecHits, EBRecHitsWeight, EERecHitsWeight, EBUncalibRecHit, EEUncalibRecHit,
										 hitsEnergy, hitsEnergyWeight,
										 hitsAmp, hitsAmpError, hitsOOAmp,
										 hitsPUSubSamples);
       /*
       std::cout<<"Electron pT "<<theRecoEl[j].pt()<<std::endl;

       for(int ic=0; ic<(int)hitsAmplitudes.size(); ic++){

	 const auto& outerVec = hitsAmplitudes[ic];
	 std::cout<<"hit E : hit WE : "<<hitsEnergy[ic]<<" "<<hitsEnergyWeight[ic]<<std::endl;

	 for (size_t j = 0; j < outerVec.size(); ++j) {

	   std::cout<<outerVec[j]<<std::endl;
	 }
       }
       */
       
       hitsAmplitudes_.push_back(hitsAmplitudes);
       
       hitsEnergy_.push_back(hitsEnergy);
       hitsEnergyWeight_.push_back(hitsEnergyWeight);
       //hitsThr_ = hitsThr;
       hitsMFAmp_.push_back(hitsAmp);
       hitsMFAmpErr_.push_back(hitsAmpError);
       hitsMFOOTAmp_.push_back(hitsOOAmp);
       hitsAmplitudesPUSub_.push_back(hitsPUSubSamples);
       e5x5_.push_back(theRecoEl[j].e5x5());

       /*
       ///SJ
       for(int ioo=0; ioo<(int)hitsOOAmp.size(); ioo++){
	 
	 std::cout<<"SJ!!! Printing OOA in the calling function "<<std::endl;
	 for(int ibx=0; ibx<(int)hitsOOAmp[ioo].size(); ibx++){
	   std::cout<<"SJ!!! ibx : amp "<<ibx<<" "<<hitsOOAmp[ioo][ibx]<<std::endl;
	 }
       }
       */
       
       /*gendR_ = 999;
       std::vector<reco::GenParticle>::const_iterator genMatch = getGenMatch(genLep, theRecoEl[j], gendR_);
       
       if(gendR_ < 999){
	 genPt_ = genMatch->pt();
	 genEta_ = genMatch->eta();
	 genPhi_ = genMatch->phi();
	 genE_ = genMatch->energy();
	 genStatus_ = genMatch->status();
       }
       else{
	 genPt_ = -99;
	 genEta_ = -99;
	 genPhi_ = -99;
	 genE_ = -99;
	 genStatus_ = -99;
       }
       */
       
       /* int icrys = 0;
       int isample = 0;
       for (const auto& array : hitsAmplitudes) {
	 isample = 0;
	 for (const auto& element : array) {

	   //std::cout<<"Conent of hitsAmplitudes_["<<icrys<<"]["<<isample<<"] is "<<hitsAmplitudes_[icrys][isample]<<std::endl;
	   
	   isample++;
	 }
	 icrys++;
       }
       */

     }///end of electron loop
   
     //// check
     // Function to print a vector of vectors of vectors of doubles
     /*
     // Iterate over each element of the outer vector
     for (const auto& outerVec : hitsAmplitudes_) {
       // Iterate over each element of the middle vector
       for (const auto& middleVec : outerVec) {
	 // Iterate over each element of the inner vector
	 for (const auto& innerVec : middleVec) {
	   // Print each double value
	   std::cout << innerVec << " "<<std::endl;
	   std::cout << std::endl; // Newline after printing inner vector
	 }
	 std::cout << std::endl; // Extra newline after printing middle vector
       }
       std::cout << std::endl; // Extra newline after printing outer vector
     }
*/

     /*
     ///SJ
     std::cout<<"In that event, info abt the OOT amp "<<std::endl;
     for(int iel=0; iel<(int)hitsMFOOTAmp_.size(); iel++){
       for(int ioo=0; ioo<(int)hitsMFOOTAmp_[iel].size(); ioo++){
	 
	 std::cout<<"SJ!!!after full loop; Printing OOA in the calling function "<<std::endl;
	 for(int ibx=0; ibx<(int)hitsMFOOTAmp_[iel][ioo].size(); ibx++){
	   std::cout<<"SJ!!!after full loop ibx : amp "<<ibx<<" "<<hitsMFOOTAmp_[iel][ioo][ibx]<<std::endl;
	 }
       }
     }
     */

     
     treeEle->Fill();
   }//if (electronHandle.isValid())



   ////cant run on photons from ALCARECO - will have to run those modules

   /*Handle<reco::PhotonCollection> theRecoPhotonCollection;
   iEvent.getByToken(recophotonCollection_, theRecoPhotonCollection);
   const reco::PhotonCollection theRecoPh = *(theRecoPhotonCollection.product());
   
   if (theRecoPhotonCollection.isValid()) {
     for (uint j = 0; j < theRecoPh.size(); j++){
       
       DetId seedDetId = (theRecoPh[j].superCluster()->seed()->hitsAndFractions())[0].first;
       bool isBarrel = (seedDetId.subdetId() == EcalBarrel);
       
       if(isBarrel){
	 //https://cmssdt.cern.ch/lxr/source/CaloOnlineTools/EcalTools/plugins/EcalCosmicsHists.cc#0109
	 EcalDigiCollection::const_iterator thisdigi = pEBDigi->find(seedDetId);
	 if (thisdigi == pEBDigi->end()){
	   std::cout<<"ECALTimeSampleAnalyser!!!  WARNING: seedDetId not found in the pEBDigi collection!"<<std::endl;
	 }
	 else{
	   EBDataFrame df(*thisdigi);
	   
	   for (unsigned int i = 0; i < (*thisdigi).size(); ++i) {
	     EcalMGPASample samp_crystal(df.sample(i));
	   }//loop over time samples
	 }//else when the id is found
       }//if(isBarrel)
       
       if(!isBarrel){
	 //https://cmssdt.cern.ch/lxr/source/CaloOnlineTools/EcalTools/plugins/EcalCosmicsHists.cc#0109
	 EcalDigiCollection::const_iterator thisdigi = pEEDigi->find(seedDetId);
	 if (thisdigi == pEEDigi->end()){
	   std::cout<<"ECALTimeSampleAnalyser!!!  WARNING: seedDetId not found in the EEDigi collection!"<<std::endl;
	 }
	 else{
	   EEDataFrame df(*thisdigi);
	   
	   for (unsigned int i = 0; i < (*thisdigi).size(); ++i) {
	     EcalMGPASample samp_crystal(df.sample(i));
	   }//loop over time samples
	 }//else when the id is found
       }//if(isBarrel)
       

     }///end of loop over photons
   }//if (photonHandle.isValid()) 
   /////////////////
   */



   /*
   if (pEBDigi) {
     for (EBDigiCollection::const_iterator digiItr = pEBDigi->begin(); digiItr != pEBDigi->end();
	  ++digiItr) {  // Loop on EB crystals
       EBDetId id_crystal(digiItr->id());
       EBDataFrame df(*digiItr);

       
       int etaG = id_crystal.ieta();  // global
       int phiG = id_crystal.iphi();  // global

       int etaL;  // local
       int phiL;  // local
       std::pair<int, int> LocalCoord = MEEBGeom::localCoord(etaG, phiG);
       etaL = LocalCoord.first;
       phiL = LocalCoord.second;
       eta = etaG;
       phi = phiG;
       side = MEEBGeom::side(etaG, phiG);
       EcalElectronicsId elecid_crystal = TheMapping.getElectronicsId(id_crystal);
       towerID = elecid_crystal.towerId();
       int strip = elecid_crystal.stripId();
       int xtal = elecid_crystal.xtalId();
       channelID = 5 * (strip - 1) + xtal - 1;  // FIXME
       int module = MEEBGeom::lmmod(etaG, phiG);
       int iMod = module - 1;
       assert(module >= *min_element(modules.begin(), modules.end()) &&
	      module <= *max_element(modules.begin(), modules.end()));
     

       std::cout<<"Size of hte digi sample "<<(*digiItr).size()<<std::endl;
       for (unsigned int i = 0; i < (*digiItr).size(); ++i) {
	 EcalMGPASample samp_crystal(df.sample(i));
	 std::cout<<"ADC of "<<i<<"th sample is "<<samp_crystal.adc()<<" and gain is "<<samp_crystal.gainId()<<std::endl;
	 
	 
	 adc[i] = samp_crystal.adc();
	 adcG[i] = samp_crystal.gainId();
	 if (i == 0)
	 adcgain = adcG[i];
	 if (i > 0)
	 adcgain = TMath::Max(adcG[i], adcgain);
	
       }
       
       // Remove pedestal
       //====================
       for (dsum = 0., dsum1 = 0., k = 0; k < _presample; k++) {
	 dsum += adc[k];
	 if (k < _presample - 1)
	   dsum1 += adc[k];
       }
       bl = dsum / ((double)_presample);
       for (val_max = 0., k = 0; k < _nsamples; k++) {
	 yrange[k] = adc[k] - bl;
	 if (yrange[k] > val_max) {
	   val_max = yrange[k];
	 }
       }
       
     }
   }
   */
}

   
// ------------ method called once each job just before starting event loop  ------------
void
ECALTimeSampleAnalyser::beginJob()
{
  edm::Service<TFileService> fs;
  treeEle    = fs->make<TTree>("EventTreeEle", "Event data");

  //
  treeEle->Branch("run",                  &run_);
  treeEle->Branch("event",                &event_);
  treeEle->Branch("lumis",                &lumis_);
  treeEle->Branch("eleE",                    &eleE_);
  treeEle->Branch("elePt",                   &elePt_);
  treeEle->Branch("eleEta",                  &eleEta_);
  treeEle->Branch("eleSCEta",                  &eleSCEta_);
  treeEle->Branch("elePhi",                  &elePhi_);
  treeEle->Branch("eleSCPhi",                 &eleSCPhi_);
  treeEle->Branch("hitsAmplitudes",         &hitsAmplitudes_);
  treeEle->Branch("hitsEnergy",         &hitsEnergy_);
  treeEle->Branch("hitsEnergyWeight",         &hitsEnergyWeight_);

  treeEle->Branch("hitsAmplitudesPUSub",         &hitsAmplitudesPUSub_);
  treeEle->Branch("hitsMFAmp",         &hitsMFAmp_);
  treeEle->Branch("hitsMFAmpErr",         &hitsMFAmpErr_);
  treeEle->Branch("hitsMFOOTAmp",         &hitsMFOOTAmp_);
  

  treeEle->Branch("eleHoverE",               &eleHoverE_);
  treeEle->Branch("eleEoverP",               &eleEoverP_);
  treeEle->Branch("eleEoverPout",            &eleEoverPout_);
  treeEle->Branch("eleEoverPInv",            &eleEoverPInv_);
  treeEle->Branch("eleBrem",                 &eleBrem_);
  treeEle->Branch("eledEtaAtVtx",            &eledEtaAtVtx_);
  treeEle->Branch("eledPhiAtVtx",            &eledPhiAtVtx_);
  treeEle->Branch("eledEtaAtCalo",           &eledEtaAtCalo_);
  treeEle->Branch("eleSigmaIEtaIEtaFull5x5", &eleSigmaIEtaIEtaFull5x5_);
  treeEle->Branch("eleSigmaIPhiIPhiFull5x5", &eleSigmaIPhiIPhiFull5x5_);

  treeEle->Branch("elePFChIso",              &elePFChIso_);
  treeEle->Branch("elePFPhoIso",             &elePFPhoIso_);
  treeEle->Branch("elePFNeuIso",             &elePFNeuIso_);
  treeEle->Branch("elePFPUIso",              &elePFPUIso_);
  treeEle->Branch("eleE1x5Full5x5",          &eleE1x5Full5x5_);
  treeEle->Branch("eleE2x5Full5x5",          &eleE2x5Full5x5_);
  treeEle->Branch("eleE5x5Full5x5",          &eleE5x5Full5x5_);
  treeEle->Branch("eleR9Full5x5",                &eleR9Full5x5_);
  
  
  //treeEle->Branch("hitsThr",         &hitsThr_);
  treeEle->Branch("nsamples",         &nsamples_);
  treeEle->Branch("nCrys",         &nCrys_);
  treeEle->Branch("e5x5",         &e5x5_);
  treeEle->Branch("nVtx",                 &nVtx_);
  treeEle->Branch("rho",                  &rho_);
  
  /*treeEle->Branch("genPt",         &genPt_);
  treeEle->Branch("genEta",         &genEta_);
  treeEle->Branch("genPhi",         &genPhi_);
  treeEle->Branch("genE",         &genE_);
  treeEle->Branch("genStatus",         &genStatus_);
  treeEle->Branch("gendR",         &gendR_);
  */
  //treeEle->Branch("hitsAmplitudes",         hitsAmplitudes_, "hitsAmplitudes_[100][nsamples]");
}

// ------------ method called once each job just after ending the event loop  ------------
void
ECALTimeSampleAnalyser::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ECALTimeSampleAnalyser::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}


std::vector<std::vector<double>> ECALTimeSampleAnalyser::getTimeSamplesAroundEle(std::vector<DetId> v_id, 
										 edm::Handle<EBDigiCollection> pEBDigi, 
										 edm::Handle<EEDigiCollection> pEEDigi, 
										 const EcalRecHitCollection* EBRecHits, 
										 const EcalRecHitCollection* EERecHits,
										 const EcalRecHitCollection* EBRecHitsWeight, 
										 const EcalRecHitCollection* EERecHitsWeight,
										 const EBUncalibratedRecHitCollection* EBUncalibRecHit,
										 const EEUncalibratedRecHitCollection* EEUncalibRecHit,
										 //const EcalPFRecHitThresholds* thresholds,
										 std::vector<double> &hitsEnergy,
										 std::vector<double> &hitsEnergyWeight,
										 std::vector<double> &hitsAmp,
										 std::vector<double> &hitsAmpError,
										 std::vector<std::vector<double>> &hitsOOAmp,
										 std::vector<std::vector<double>> &hitsPUSubSamples
										 ){
                                                                                 //std::vector<double> &hitsThr){
  
  //// follow ecalmultifit algo as linked in teh analyse function above and do the pedestal subtraction. link also given below
  ///https://cmssdt.cern.ch/lxr/source/RecoLocalCalo/EcalRecAlgos/src/EcalUncalibRecHitMultiFitAlgo.cc#0064
  
  ////https://cmssdt.cern.ch/lxr/source/RecoLocalCalo/EcalRecProducers/plugins/EcalUncalibRecHitWorkerMultiFit.cc#0512

  //using sampleType = std::array<double, NMAXSAMPLES>;
  std::vector<std::vector<double>> hitTimeSamples;
  
   FullSampleVector fullpulse(FullSampleVector::Zero());
   
    for (const auto& id : v_id) {

      /*
      float rhThres = -99.0;
      if (thresholds != nullptr) {
        rhThres = (*thresholds)[id];  // access PFRechit thresholds for noise cleaning
      }
      
      hitsThr.push_back(rhThres);
      */
      
      bool isBarrel = (id.subdetId() == EcalBarrel); 
      const EcalPedestals::Item* aped = nullptr;
      const EcalMGPAGainRatio* aGain = nullptr;
      
      std::vector<double> amplitudes(NMAXSAMPLES);
      std::vector<double> ootAmp_(activeBX.size());
      std::vector<double> PUSubSamples(NMAXSAMPLES);
	    
      for(int isample=0; isample<NMAXSAMPLES; isample++){
	amplitudes[isample] = -99.;
	PUSubSamples[isample] = -99.;
      }

      for(int isample=0; isample<(int)activeBX.size(); isample++){ 
	ootAmp_[isample] = -99.;
      }
      
      double rechitEn = -99;
      double rechitEnWeight = -99;
      double amp_, ampError_;
      
      const EcalPulseShapes::Item* aPulse = nullptr;
      //std::cout<<" "<<std::endl;
      //std::cout<<"SJ!!! this ID is in EB "<<isBarrel<<std::endl;
      if(isBarrel){
	
	//https://cmssdt.cern.ch/lxr/source/CaloOnlineTools/EcalTools/plugins/EcalCosmicsHists.cc#

	EBUncalibratedRecHitCollection::const_iterator itunc =  EBUncalibRecHit->find(id);
	if(itunc != EBUncalibRecHit->end()){
	  if ((itunc->checkFlag(EcalUncalibratedRecHit::kPoorReco))){
	    amp_ = -99;
	    ampError_ = -99;
	    for(int ibx=0; ibx<(int)activeBX.size(); ibx++){
	      ootAmp_[ibx] = -99;
	    }
	    
	  }else{
	    amp_ = itunc->amplitude();
	    ampError_ = itunc->amplitudeError();
	    for(int ibx=0; ibx<(int)activeBX.size(); ibx++){
	      ootAmp_[ibx] = itunc->outOfTimeAmplitude(ibx);
	      //std::cout<<"SJ!!! ibx : oot "<<ibx<<" "<<ootAmp_[ibx]<<std::endl;
	    }
	  }
	} else {
	  amp_ = -99;
	  ampError_ = -99;
	  for(int ibx=0; ibx<(int)activeBX.size(); ibx++){
	    ootAmp_[ibx] = -99;
	  }
	}
	
	  for (int i = 0; i < EcalPulseShape::TEMPLATESAMPLES; ++i){
	    fullpulse(i + 7) = -99;
	  }
	  
	  //std::cout<<"SJ!!! size of ootAmp in EB "<<ootAmp_.size()<<std::endl;
	
      
	/// Nominal method
	EcalRecHitCollection::const_iterator it = EBRecHits->find(id);
	
	if (it != EBRecHits->end()) {
	  if ((it->checkFlag(EcalRecHit::kTowerRecovered) || it->checkFlag(EcalRecHit::kWeird) ||
	       (it->detid().subdetId() == EcalBarrel && it->checkFlag(EcalRecHit::kDiWeird)))){
	    rechitEn = 0.0;
	  }
	  else{
	    rechitEn = it->energy();
	  }
	} else {
	  rechitEn =  0;
	}
      

	/////Weight method
	EcalRecHitCollection::const_iterator itw = EBRecHitsWeight->find(id);
	
	if (itw != EBRecHitsWeight->end()) {
	  if ((itw->checkFlag(EcalRecHit::kTowerRecovered) || itw->checkFlag(EcalRecHit::kWeird) ||
	       (itw->detid().subdetId() == EcalBarrel && itw->checkFlag(EcalRecHit::kDiWeird)))){
	    rechitEnWeight = 0.0;
	  
	  }
	  else{
	    rechitEnWeight = itw->energy();
	  }
	} else {
	  rechitEnWeight =  0;
	}
	
	//std::cout<<"ECAL rechit energy-   multiFit : weight "<<rechitEn<<" "<<rechitEnWeight<<std::endl;

	
	unsigned int hashedIndex = EBDetId(id).hashedIndex();
	aped = &peds->barrel(hashedIndex);
	aGain = &gains->barrel(hashedIndex);

	aPulse = &pulseshapes->barrel(hashedIndex);
	
	///https://cmssdt.cern.ch/lxr/source/RecoLocalCalo/EcalRecProducers/plugins/EcalUncalibRecHitWorkerMultiFit.cc#0465
	for (int i = 0; i < EcalPulseShape::TEMPLATESAMPLES; ++i){
	  fullpulse(i + 7) = aPulse->pdfval[i];
	}
	

	EcalDigiCollection::const_iterator thisdigi = pEBDigi->find(id);
	if (thisdigi == pEBDigi->end()){
	  //std::cout<<"ECALTimeSampleAnalyser!!!  WARNING: seedDetId not found in the pEBDigi collection!"<<std::endl;
	
	  for(int isample=0; isample<NMAXSAMPLES; isample++){
	    amplitudes[isample] = -99.;
	    PUSubSamples[isample] = -99;    
	  }
	  hitTimeSamples.push_back(amplitudes);
	  hitsPUSubSamples.push_back(PUSubSamples);
	}
	
	else{
	  EBDataFrame df(*thisdigi);

	  float pulsenorm = 0.;
	  for (unsigned int i = 0; i < (*thisdigi).size(); i++) {
	    EcalMGPASample samp_crystal(df.sample(i));
	    //std::cout<<"ADC of "<<i<<"th sample is "<<samp_crystal.adc()<<" and gain is "<<samp_crystal.gainId()<<std::endl;
	    
	    double amplitude = 0.;
	    int gainId = samp_crystal.gainId();
	    
	    double pedestal = 0.;
	    double gainratio = 1.;
	    
	    
	    if (gainId == 0 || gainId == 3) {
	      pedestal = aped->mean_x1;
	      gainratio = aGain->gain6Over1() * aGain->gain12Over6();
	    } else if (gainId == 1) {
	      pedestal = aped->mean_x12;
	      gainratio = 1.;
	    } else if (gainId == 2) {
	      pedestal = aped->mean_x6;
	      gainratio = aGain->gain12Over6();
	    }
	  /// static pedestal here - no multifit at the moment
	    amplitude = ((double)(samp_crystal.adc()) - pedestal) * gainratio;
	    if (gainId == 0) {
	      //saturation
	      amplitude = (4095. - pedestal) * gainratio;
	    }
	    
	    amplitudes[i] = amplitude;
	    PUSubSamples[i] = amplitudes[i];
	    pulsenorm += fullpulse(i);
	    //std::cout<<"amplitude after sub "<<amplitude<<std::endl;
	  }///loop over time samples
	  hitTimeSamples.push_back(amplitudes);
	  
	
	  //std::cout<<"======"<<std::endl;
	  ///PU subtracted (as estimated from multifit subtracted sample)
	  int ipulse = -1;

	  for (auto const& amplit : ootAmp_) {
	    ipulse++;
	    int bxp3 = ipulse - 2;
	    int firstsamplet = std::max(0, bxp3);
	    int offset = 7 - bxp3;

	    //std::cout<<"bxp3 : firstsamplet : offset : "<<bxp3<<" "<<firstsamplet<<" "<<offset<<std::endl;
	    for (unsigned int isample = firstsamplet; isample < NMAXSAMPLES; ++isample) {
	      auto const pulse = fullpulse(isample + offset);
	      //std::cout<<"i : amplitudes : amplit : pulse : pulsenorm "<<isample<<" "<<amplitudes[isample]<<" "<<amplit<<" "<<pulse<<" "<<pulsenorm<<std::endl;
	      PUSubSamples[isample] = std::max(0., amplitudes[isample] - amplit * pulse / pulsenorm);
	      //std::cout<<"PU subtracted amp till here "<<PUSubSamples[isample]<<std::endl;
	    }
	  }//for (auto const& amplit : amplitudes)
	  
	}//else when the id is found
	/*
	for(int i=0; i<NMAXSAMPLES; i++){
	  std::cout<<"before PU : PU subt final isample "<<i<<" is "<<amplitudes[i]<<" "<<PUSubSamples[i]<<std::endl;
	}
	*/
	hitsPUSubSamples.push_back(PUSubSamples);
		    
      }//if(isBarrel)
      
      if(!isBarrel){
	//https://cmssdt.cern.ch/lxr/source/CaloOnlineTools/EcalTools/plugins/EcalCosmicsHists.cc#
	
	EEUncalibratedRecHitCollection::const_iterator itunc =  EEUncalibRecHit->find(id);
	if(itunc != EEUncalibRecHit->end()){
	  if ((itunc->checkFlag(EcalUncalibratedRecHit::kPoorReco))){
	    amp_ = -99;
	    ampError_ = -99;
	    for(int ibx=0; ibx<(int)activeBX.size(); ibx++){
	      ootAmp_[ibx] = -99;
	    }
	    
	  }else{
	    amp_ = itunc->amplitude();
	    ampError_ = itunc->amplitudeError();
	    for(int ibx=0; ibx<(int)activeBX.size(); ibx++){
	      ootAmp_[ibx] = itunc->outOfTimeAmplitude(ibx);
	    }
	  }
	} else {
	    amp_ = -99;
	    ampError_ = -99;
	    for(int ibx=0; ibx<(int)activeBX.size(); ibx++){
	      ootAmp_[ibx] = -99;
	    }
	}

	for (int i = 0; i < EcalPulseShape::TEMPLATESAMPLES; ++i){
	  fullpulse(i + 7) = -99;
	}
	
	//std::cout<<"SJ!!! size of ootAmp in EE "<<ootAmp_.size()<<std::endl;
	
	/// nominal method
	EcalRecHitCollection::const_iterator it = EERecHits->find(id);
	if (it != EERecHits->end()) {
	  if ((it->checkFlag(EcalRecHit::kTowerRecovered) || it->checkFlag(EcalRecHit::kWeird) ||
	       (it->detid().subdetId() == EcalEndcap && it->checkFlag(EcalRecHit::kDiWeird)))){
	    rechitEn = 0.0;
	  }
	  else{
	    rechitEn = it->energy();
	  }
	} else {
	  rechitEn =  0;
	}
	/////Weight method
	EcalRecHitCollection::const_iterator itw = EERecHitsWeight->find(id);
	
	if (itw != EERecHitsWeight->end()) {
	  if ((itw->checkFlag(EcalRecHit::kTowerRecovered) || itw->checkFlag(EcalRecHit::kWeird) ||
	       (itw->detid().subdetId() == EcalBarrel && itw->checkFlag(EcalRecHit::kDiWeird)))){
	    rechitEnWeight = 0.0;
	  }
	  else{
	    rechitEnWeight = itw->energy();
	  }
	} else {
	  rechitEnWeight =  0;
	}
	
	unsigned int hashedIndex = EEDetId(id).hashedIndex();
	aped = &peds->endcap(hashedIndex);
	aGain = &gains->endcap(hashedIndex);

	aPulse = &pulseshapes->endcap(hashedIndex);
	
	///https://cmssdt.cern.ch/lxr/source/RecoLocalCalo/EcalRecProducers/plugins/EcalUncalibRecHitWorkerMultiFit.cc#0465
	for (int i = 0; i < EcalPulseShape::TEMPLATESAMPLES; ++i){
	  fullpulse(i + 7) = aPulse->pdfval[i];
	}
      
	EcalDigiCollection::const_iterator thisdigi = pEEDigi->find(id);
	if (thisdigi == pEEDigi->end()){
	  //std::cout<<"ECALTimeSampleAnalyser!!!  WARNING: seedDetId not found in the pEEDigi collection!"<<std::endl;
	  for(int isample=0; isample<NMAXSAMPLES; isample++){
	    amplitudes[isample] = -99.;
	    PUSubSamples[isample] = -99;
	  }
	  hitTimeSamples.push_back(amplitudes);
	  hitsPUSubSamples.push_back(PUSubSamples);
	}
	
	else{
	  EEDataFrame df(*thisdigi);

	  float pulsenorm = 0.;
	  for (unsigned int i = 0; i < (*thisdigi).size(); i++) {
	    EcalMGPASample samp_crystal(df.sample(i));
	    //std::cout<<"ADC of "<<i<<"th sample is "<<samp_crystal.adc()<<" and gain is "<<samp_crystal.gainId()<<std::endl;
	    
	    double amplitude = 0.;
	    int gainId = samp_crystal.gainId();
	    
	    double pedestal = 0.;
	    double gainratio = 1.;
	    
	    
	    if (gainId == 0 || gainId == 3) {
	      pedestal = aped->mean_x1;
	      gainratio = aGain->gain6Over1() * aGain->gain12Over6();
	    } else if (gainId == 1) {
	      pedestal = aped->mean_x12;
	      gainratio = 1.;
	    } else if (gainId == 2) {
	      pedestal = aped->mean_x6;
	      gainratio = aGain->gain12Over6();
	    }
	    /// static pedestal here - no multifit at the moment
	    amplitude = ((double)(samp_crystal.adc()) - pedestal) * gainratio;
	    if (gainId == 0) {
	      //saturation
	    amplitude = (4095. - pedestal) * gainratio;
	    }
	    
	    amplitudes[i] = amplitude;
	    PUSubSamples[i] = amplitudes[i];
	    //std::cout<<"amplitude after sub "<<amplitude<<std::endl;
	    pulsenorm += fullpulse(i);
	  }//loop over time samples
	  hitTimeSamples.push_back(amplitudes);

	  ///PU subtracted (as estimated from multifit subtracted sample)
	  int ipulse = -1;
	  for (auto const& amplit : ootAmp_) {
	    ipulse++;
	    int bxp3 = ipulse - 2;
	    int firstsamplet = std::max(0, bxp3);
	    int offset = 7 - bxp3;
	    
	    for (unsigned int isample = firstsamplet; isample < NMAXSAMPLES; ++isample) {
	      auto const pulse = fullpulse(isample + offset);
	      PUSubSamples[isample] = std::max(0., amplitudes[isample] - amplit * pulse / pulsenorm);
	    }
	  }//for (auto const& amplit : amplitudes)

	}//else when the id is found
	hitsPUSubSamples.push_back(PUSubSamples);
	
      }//if(!isBarrel)
      hitsEnergy.push_back(rechitEn);

      hitsEnergyWeight.push_back(rechitEnWeight);
    
      hitsAmp.push_back(amp_);
      hitsAmpError.push_back(ampError_);
      hitsOOAmp.push_back(ootAmp_);

      //std::cout<<"# time samples in PUSubSamples : amp : OOTAmp "<<PUSubSamples.size()<<" "<<amplitudes.size()<<" "<<ootAmp_.size()<<std::endl;
    }//for (const auto& id : v_id)
    
      
    
    ///SJ
    /*
    for(int ioo=0; ioo<(int)hitsOOAmp.size(); ioo++){
      
      std::cout<<"SJ!!! Printing OOA in the called function "<<std::endl;
      for(int ibx=0; ibx<(int)hitsOOAmp[ioo].size(); ibx++){
	std::cout<<"SJ!!! ibx : amp "<<ibx<<" "<<hitsOOAmp[ioo][ibx]<<std::endl;
      }
    }
    */

	
    return hitTimeSamples;
    
}


//genmatch 
std::vector<reco::GenParticle>::const_iterator  ECALTimeSampleAnalyser::getGenMatch(std::vector<std::vector<reco::GenParticle>::const_iterator> genLep, reco::GsfElectron gsfele, double &dRmin){

  std::vector<reco::GenParticle>::const_iterator genMatch;
    TLorentzVector *ele = new TLorentzVector();
    ele->SetPtEtaPhiM(gsfele.pt(), gsfele.eta(), gsfele.phi(), 0.511/1000.);

    dRmin = 999;
    
  for(auto& lep : genLep){
    TLorentzVector *gen = new TLorentzVector(); 
    gen->SetPtEtaPhiM(lep->pt(), lep->eta(), lep->phi(), lep->mass());
    double deltaR = ele->DeltaR(*gen);
    
    if(deltaR < dRmin){
      //std::cout<<"dR at the moment "<<deltaR<<std::endl;
      dRmin = deltaR;
      genMatch = lep;
    }
  }
  
  //std::cout<<"===final dR "<<dRmin<<std::endl;
  return genMatch;
}


//define this as a plug-in
DEFINE_FWK_MODULE(ECALTimeSampleAnalyser);
