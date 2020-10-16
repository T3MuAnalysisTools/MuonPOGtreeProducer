//////////////////////////////////////
// Ntuplizer that fills muon_pog trees
//////////////////////////////////////

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h" 
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/Conversion.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"

#include "DataFormats/Common/interface/View.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"

#include "DataFormats/Candidate/interface/CompositePtrCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"

#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/L1Trigger/interface/Muon.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenStatusFlags.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 

#include "DataFormats/Scalers/interface/LumiScalers.h"
#include "DataFormats/Luminosity/interface/LumiDetails.h"

#include "DataFormats/GeometryVector/interface/VectorUtil.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"


#include "MuonPOGtreeProducer/Tools/src/MuonPogTree.h"
#include "MuonPOGtreeProducer/Tools/src/PDG_Vars.h"
#include "Math/LorentzVector.h"
#include "Math/MatrixFunctions.h"
#include <Math/SVector.h>
#include <Math/SMatrix.h>
#include "TMatrixT.h"
#include "TTree.h"
#include "TLorentzVector.h"

#include <algorithm>
#include <iostream>

typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3>> SMatrixSym3D;
typedef ROOT::Math::SVector<double, 3> SVector3;

bool DEBUG = false;

#define PHI_MASS_WINDOW_LOW 0.96
#define PHI_MASS_WINDOW_HIGH 1.04

using namespace edm;
using namespace reco;
using namespace std;

class MuonPogTreeProducer : public edm::EDAnalyzer 
{
   public:

      MuonPogTreeProducer(const edm::ParameterSet &);

      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void beginJob();
      virtual void endJob();

   private:

      void fillGenInfo(const edm::Handle<std::vector<PileupSummaryInfo> > &,
            const  edm::Handle<GenEventInfoProduct> &);

      void fillGenParticles(const edm::Handle<reco::GenParticleCollection> &);

      void fillHlt(const edm::Handle<edm::TriggerResults> &, 
            const edm::Handle<trigger::TriggerEvent> &,
            const edm::TriggerNames &);

      void fillHlt(const edm::Handle<edm::TriggerResults> &,
            const edm::TriggerNames&,
            const edm::Handle<std::vector<pat::TriggerObjectStandAlone> > &);

      void fillPV(const edm::Handle<std::vector<reco::Vertex> > &);


      Int_t fillMuons(const edm::Handle<edm::View<reco::Muon> > &,
            const edm::Handle<std::vector<reco::Vertex> > &,
            const edm::Handle<reco::BeamSpot> &);

      Bool_t MuonMatchedKs(const reco::Muon& ,
            const edm::Handle<std::vector<reco::VertexCompositePtrCandidate>> &);

      void fillSimMuonInfo(const edm::Handle<edm::ValueMap<reco::MuonSimInfo>>& );

      void fillL1(const edm::Handle<l1t::MuonBxCollection> &);

      void fillKsVertices(const edm::Handle<std::vector<reco::VertexCompositePtrCandidate> > &,
            const edm::Handle<std::vector<reco::Vertex> >&,
            const edm::Handle<reco::BeamSpot>& );

      void fillKsVertices(const edm::Handle<std::vector<reco::VertexCompositeCandidate> >&,
            const edm::Handle<std::vector<reco::Vertex> >&, 
            const edm::Handle<reco::BeamSpot>& );

      void fillPhiVertices(const edm::EventSetup&,const edm::Handle<std::vector<reco::Track>>& ,
            const edm::Handle<std::vector<reco::Vertex> >& ,
            const edm::Handle<reco::BeamSpot>& );

      void fillPhiVertices(const edm::Handle<std::vector<reco::VertexCompositePtrCandidate> >&,
            const edm::Handle<std::vector<reco::Vertex> >&,
            const edm::Handle<reco::BeamSpot>& );

      float getVhitsComb(const reco::Muon& );

      // returns false in case the match is for a RPC chamber
      bool getMuonChamberId(DetId & id, muon_pog::MuonDetType & det, Int_t & r, Int_t & phi, Int_t & eta) const ;

      edm::EDGetTokenT<edm::TriggerResults> trigResultsToken_;
      edm::EDGetTokenT<trigger::TriggerEvent> trigSummaryToken_;
      edm::EDGetTokenT<pat::TriggerObjectStandAlone> triggerObjectToken_;

      std::string trigFilterCut_;
      std::string trigPathCut_;

      edm::EDGetTokenT<edm::View<reco::Muon> > muonToken_;
      edm::EDGetTokenT<std::vector<reco::Track> > trackToken_;
      edm::EDGetTokenT<edm::ValueMap<reco::MuonSimInfo> > muonSimToken_;
      edm::EDGetTokenT<std::vector<reco::Vertex> > primaryVertexToken_;
      edm::EDGetTokenT<std::vector<reco::VertexCompositePtrCandidate> > secondaryKsVertexToken_;
      edm::EDGetTokenT<std::vector<reco::VertexCompositeCandidate> > secondaryKsVertexTokenAOD_;
      edm::EDGetTokenT<std::vector<reco::VertexCompositePtrCandidate> >secondaryVertexToken_;
      edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;

      edm::EDGetTokenT<reco::PFMETCollection> pfMetToken_;
      edm::EDGetTokenT<reco::PFMETCollection> pfChMetToken_;
      edm::EDGetTokenT<reco::CaloMETCollection> caloMetToken_;

      edm::EDGetTokenT<reco::GenParticleCollection> genToken_;
      edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileUpInfoToken_;
      edm::EDGetTokenT<GenEventInfoProduct> genInfoToken_;
      edm::EDGetTokenT<edm::ValueMap<reco::MuonSimInfo>> simMuonToken_;

      edm::EDGetTokenT<LumiScalersCollection> scalersToken_;
      edm::EDGetTokenT<l1t::MuonBxCollection> l1Token_;

      Float_t m_minMuPtCut;
      Int_t m_minNMuCut;
      Bool_t miniAODRun;
      Bool_t doKsVertices;
      Bool_t doPhiVertices;

      muon_pog::Event event_;
      muon_pog::EventId eventId_;
      std::map<std::string,TTree*> tree_;

};


MuonPogTreeProducer::MuonPogTreeProducer( const edm::ParameterSet & cfg )
{

   // Input collections
   edm::InputTag tag = cfg.getUntrackedParameter<edm::InputTag>("TrigResultsTag", edm::InputTag("TriggerResults::HLT"));
   if (tag.label() != "none") trigResultsToken_ = consumes<edm::TriggerResults>(tag);

   tag = cfg.getUntrackedParameter<edm::InputTag>("TrigSummaryTag", edm::InputTag("hltTriggerSummaryAOD::HLT")); 
   if (tag.label() != "none") trigSummaryToken_ =consumes<trigger::TriggerEvent>(tag);

   tag = cfg.getUntrackedParameter<edm::InputTag>("TriggerObjectStandAlone", edm::InputTag("slimmedPatTrigger::RECO"));
   if (tag.label() != "none") triggerObjectToken_ =consumes<pat::TriggerObjectStandAlone>(tag);

   trigFilterCut_ = cfg.getUntrackedParameter<std::string>("TrigFilterCut", std::string("all"));
   trigPathCut_ = cfg.getUntrackedParameter<std::string>("TrigPathCut", std::string("all"));

   tag = cfg.getUntrackedParameter<edm::InputTag>("MuonTag", edm::InputTag("muons"));
   if (tag.label() != "none") muonToken_ = consumes<edm::View<reco::Muon> >(tag);

   tag = cfg.getUntrackedParameter<edm::InputTag>("TrackTag", edm::InputTag("pfTracks"));
   if (tag.label() != "none") trackToken_ = consumes<std::vector<reco::Track> >(tag); 

   tag = cfg.getUntrackedParameter<edm::InputTag>("PrimaryVertexTag", edm::InputTag("offlinePrimaryVertices"));
   if (tag.label() != "none") primaryVertexToken_ = consumes<std::vector<reco::Vertex> >(tag);

   tag = cfg.getUntrackedParameter<edm::InputTag>("SecondaryKsVertexTag", edm::InputTag("slimmedKshortVertices"));
   if (tag.label()!= "none") secondaryKsVertexToken_ = consumes<std::vector<reco::VertexCompositePtrCandidate> >(tag);

   tag = cfg.getUntrackedParameter<edm::InputTag>("SecondaryKsVertexTag", edm::InputTag("generalV0Candidates:Kshort:RECO"));
   if (tag.label()!= "none") secondaryKsVertexTokenAOD_ = consumes<std::vector<reco::VertexCompositeCandidate> >(tag);

   tag = cfg.getUntrackedParameter<edm::InputTag>("InclusiveSVTag", edm::InputTag("slimmedSecondaryVertices"));
   if (tag.label()!="none") secondaryVertexToken_ = consumes<std::vector<reco::VertexCompositePtrCandidate> >(tag);

   tag = cfg.getUntrackedParameter<edm::InputTag>("BeamSpotTag", edm::InputTag("offlineBeamSpot"));
   if (tag.label() != "none") beamSpotToken_ = consumes<reco::BeamSpot>(tag);

   tag = cfg.getUntrackedParameter<edm::InputTag>("PFMetTag", edm::InputTag("pfMet"));
   if (tag.label() != "none") pfMetToken_ = consumes<reco::PFMETCollection>(tag);

   tag = cfg.getUntrackedParameter<edm::InputTag>("PFChMetTag", edm::InputTag("pfChMet"));
   if (tag.label() != "none") pfChMetToken_ = consumes<reco::PFMETCollection>(tag);

   tag = cfg.getUntrackedParameter<edm::InputTag>("CaloMetTag", edm::InputTag("caloMet"));
   if (tag.label() != "none") caloMetToken_ = consumes<reco::CaloMETCollection>(tag); 

   tag = cfg.getUntrackedParameter<edm::InputTag>("GenTag", edm::InputTag("prunedGenParticles"));
   if (tag.label() != "none") genToken_ = consumes<reco::GenParticleCollection>(tag);

   tag = cfg.getUntrackedParameter<edm::InputTag>("PileUpInfoTag", edm::InputTag("pileupInfo"));
   if (tag.label() != "none") pileUpInfoToken_ = consumes<std::vector<PileupSummaryInfo> >(tag);

   tag = cfg.getUntrackedParameter<edm::InputTag>("GenInfoTag", edm::InputTag("generator"));
   if (tag.label() != "none") genInfoToken_ = consumes<GenEventInfoProduct>(tag);  

   tag = cfg.getUntrackedParameter<edm::InputTag>("ScalersTag", edm::InputTag("scalersRawToDigi"));
   if (tag.label() != "none") scalersToken_ = consumes<LumiScalersCollection>(tag);

   tag = cfg.getUntrackedParameter<edm::InputTag>("l1MuonsTag", edm::InputTag("gmtStage2Digis:Muon:"));
   if (tag.label() != "none") l1Token_ = consumes<l1t::MuonBxCollection>(tag);

   tag = cfg.getUntrackedParameter<edm::InputTag>("simMuonTag", edm::InputTag("muonSimClassifier"));
   if (tag.label() != "none") simMuonToken_ = consumes<edm::ValueMap<reco::MuonSimInfo>>(tag);

   m_minMuPtCut = cfg.getUntrackedParameter<double>("MinMuPtCut", 0.);
   m_minNMuCut  = cfg.getUntrackedParameter<int>("MinNMuCut",  0);
   miniAODRun = cfg.getUntrackedParameter<bool>("miniAODRun", false);
   doKsVertices = cfg.getUntrackedParameter<bool>("doKsVertices", false);
   doPhiVertices = cfg.getUntrackedParameter<bool>("doPhiVertices", false);
}


void MuonPogTreeProducer::beginJob() 
{

   edm::Service<TFileService> fs;
   tree_["muPogTree"] = fs->make<TTree>("MUONPOGTREE","Muon POG Tree");

   tree_["muPogTree"]->Branch("event",&event_);
   tree_["muPogTree"]->Branch("eventId",&eventId_);
}



void MuonPogTreeProducer::endJob() 
{

}


void MuonPogTreeProducer::analyze (const edm::Event & ev, const edm::EventSetup & iSetup)
{

   // Clearing branch variables
   // and setting default values
   event_.hlt.triggers.clear();
   event_.hlt.objects.clear();
   event_.l1muons.clear();

   event_.genParticles.clear();
   event_.genInfos.clear();
   event_.muons.clear();
   event_.sim_muons.clear();
   event_.kshorts.clear();
   event_.phis.clear();

   event_.mets.pfMet   = -999; 
   event_.mets.pfChMet = -999; 
   event_.mets.caloMet = -999; 


   for (unsigned int ix=0; ix<3; ++ix) {
      event_.primaryVertex[ix] = 0.;
      for (unsigned int iy=0; iy<3; ++iy) {
         event_.cov_primaryVertex[ix][iy] = 0.;
      }
   }
   event_.nVtx = -1;

   // Fill general information
   // run, luminosity block, event
   event_.runNumber = ev.id().run();
   event_.luminosityBlockNumber = ev.id().luminosityBlock();
   event_.eventNumber = ev.id().event();

   // Fill GEN pile up information
   if (!ev.isRealData()) 
   {
      if (!pileUpInfoToken_.isUninitialized() &&
            !genInfoToken_.isUninitialized()) 
      {
         edm::Handle<std::vector<PileupSummaryInfo> > puInfo;
         edm::Handle<GenEventInfoProduct> genInfo;

         if (ev.getByToken(pileUpInfoToken_, puInfo) &&
               ev.getByToken(genInfoToken_, genInfo) ) 
            fillGenInfo(puInfo,genInfo);
         else 
            edm::LogError("") << "[MuonPogTreeProducer]: Pile-Up Info collection does not exist !!!";
      }      
   }


   // Fill GEN particles information
   if (!ev.isRealData()) 
   {
      if (!genToken_.isUninitialized() ) 
      { 
         edm::Handle<reco::GenParticleCollection> genParticles;
         if (ev.getByToken(genToken_, genParticles)){ 
            fillGenParticles(genParticles);
         }
         else 
            edm::LogError("") << ">>> GEN collection does not exist !!!";
      }
   }

   if (ev.isRealData()) 
   {

      event_.bxId  = ev.bunchCrossing();
      event_.orbit = ev.orbitNumber();

      if (!scalersToken_.isUninitialized()) 
      { 
         edm::Handle<LumiScalersCollection> lumiScalers;
         if (ev.getByToken(scalersToken_, lumiScalers) && 
               lumiScalers->size() > 0 ) 
            event_.instLumi  = lumiScalers->begin()->instantLumi();
         else 
            edm::LogError("") << ">>> Scaler collection does not exist !!!";
      }
   }

   // Fill trigger information
   if (!trigResultsToken_.isUninitialized() &&
         !trigSummaryToken_.isUninitialized()) 
   {

      edm::Handle<edm::TriggerResults> triggerResults;
      edm::Handle<trigger::TriggerEvent> triggerEvent;
      edm::Handle<std::vector<pat::TriggerObjectStandAlone>> triggerObjects;

      if (miniAODRun){
         if (ev.getByToken(trigResultsToken_, triggerResults) &&
               ev.getByToken(triggerObjectToken_, triggerObjects)) 
            fillHlt(triggerResults, ev.triggerNames(*triggerResults), triggerObjects);
         else 
            edm::LogError("") << "[MuonPogTreeProducer]: Trigger collections do not exist !!!";
      }
      else {
         if (ev.getByToken(trigResultsToken_, triggerResults) &&
               ev.getByToken(trigSummaryToken_, triggerEvent)) 
            fillHlt(triggerResults, triggerEvent,ev.triggerNames(*triggerResults));
         else 
            edm::LogError("") << "[MuonPogTreeProducer]: Trigger collections do not exist !!!";
      }
   }


   // Fill vertex information
   edm::Handle<std::vector<reco::Vertex> > vertexes;
   edm::Handle<std::vector<reco::VertexCompositePtrCandidate> > secondVertexes;

   if(!primaryVertexToken_.isUninitialized()) 
   {
      if (ev.getByToken(primaryVertexToken_, vertexes)){
         fillPV(vertexes);
      }
      else 
         edm::LogError("") << "[MuonPogTreeProducer]: Vertex collection does not exist !!!";
   }
   /*
      if(!secondaryVertexToken_.isUninitialized()){
      if(!ev.getByToken(secondaryVertexToken_, secondVertexes))
      edm::LogError("") << "[MuonPogTreeProducer]: Secondary vertex collection does not exist !!!";
      }
      */
   // Get beam spot for muons
   edm::Handle<reco::BeamSpot> beamSpot;
   if (!beamSpotToken_.isUninitialized() ) 
   { 
      if (!ev.getByToken(beamSpotToken_, beamSpot)) 
         edm::LogError("") << "[MuonPogTreeProducer]: Beam spot collection not found !!!";
   }

   // Fill (raw) MET information: PF, PF charged, Calo    
   edm::Handle<reco::PFMETCollection> pfMet; 
   if(!pfMetToken_.isUninitialized()) 
   { 
      if (!ev.getByToken(pfMetToken_, pfMet)) 
         edm::LogError("") << "[MuonPogTreeProducer] PFMet collection does not exist !!!"; 
      else { 
         const reco::PFMET &iPfMet = (*pfMet)[0]; 
         event_.mets.pfMet = iPfMet.et(); 
      } 
   } 

   edm::Handle<reco::PFMETCollection> pfChMet; 
   if(!pfChMetToken_.isUninitialized()) 
   { 
      if (!ev.getByToken(pfChMetToken_, pfChMet)) 
         edm::LogError("") << "[MuonPogTreeProducer] PFChMet collection does not exist !!!"; 
      else { 
         const reco::PFMET &iPfChMet = (*pfChMet)[0]; 
         event_.mets.pfChMet = iPfChMet.et(); 
      } 
   } 

   edm::Handle<reco::CaloMETCollection> caloMet; 
   if(!caloMetToken_.isUninitialized()) 
   { 
      if (!ev.getByToken(caloMetToken_, caloMet)) 
         edm::LogError("") << "[MuonPogTreeProducer] CaloMet collection does not exist !!!"; 
      else { 
         const reco::CaloMET &iCaloMet = (*caloMet)[0]; 
         event_.mets.caloMet = iCaloMet.et(); 
      } 
   } 

   // Get muons  
   edm::Handle<edm::View<reco::Muon> > muons;
   edm::Handle<reco::VertexCompositePtrCandidateCollection> kshorts;
   edm::Handle<reco::VertexCompositeCandidateCollection> kshorts_aod;

   if (!muonToken_.isUninitialized() ) 
   { 
      if (!ev.getByToken(muonToken_, muons)) 
         edm::LogError("") << "[MuonPogTreeProducer] Muon collection does not exist !!!";
   }

   if(!secondaryKsVertexToken_.isUninitialized() ){
      if (miniAODRun){
         if (!ev.getByToken(secondaryKsVertexToken_, kshorts)) 
            edm::LogError("") << "[MuonPogTreeProducer] Kshort Collection does not exist !!!";
      }
      else if (!ev.getByToken(secondaryKsVertexTokenAOD_, kshorts_aod)) 
         edm::LogError("") << "[MuonPogTreeProducer] Kshort Collection does not exist !!!";
   }

   Int_t nGoodMuons = 0;
   eventId_.maxPTs.clear();

   // Fill muon information
   if (muons.isValid() && vertexes.isValid() && beamSpot.isValid()) 
   {
      nGoodMuons = fillMuons(muons,vertexes,beamSpot);
   }

   eventId_.nMuons = nGoodMuons;

   //Fill Sim muon info
   if (!ev.isRealData()) {
      edm::Handle<edm::ValueMap<reco::MuonSimInfo>> muonSimCollection;
      if (muons.isValid()) 
      { 
         if (!ev.getByToken(simMuonToken_, muonSimCollection)) edm::LogError("") << "[MuonPogTreeProducer] Muon Sim collection does not exist!";
         else fillSimMuonInfo(muonSimCollection);
      }
   }

   //Fill L1 informations
   edm::Handle<l1t::MuonBxCollection> l1s;
   if (!l1Token_.isUninitialized() )
   {
      if (!ev.getByToken(l1Token_, l1s))
         edm::LogError("") << "[MuonPogTreeProducer] L1 muon bx collection does not exist !!!";
      else {
         fillL1(l1s);
      }
   }

   // fill Kshort and phi vertex information
   edm::Handle<std::vector<reco::Track> > tracks;
   if (!trackToken_.isUninitialized() ){
      if (!ev.getByToken(trackToken_, tracks))
         edm::LogError("") << "[MuonPogTreeProducer] Track collection does not exist !!!";
      else {
         if (doKsVertices && kshorts.isValid() && vertexes.isValid() && beamSpot.isValid() && miniAODRun) fillKsVertices(kshorts, vertexes, beamSpot);
         else if (doKsVertices && kshorts_aod.isValid() && vertexes.isValid() && beamSpot.isValid() && !miniAODRun) fillKsVertices(kshorts_aod, vertexes, beamSpot);
         // Fill phi using user defined algorithm
         if (doPhiVertices && tracks.isValid() && vertexes.isValid() && beamSpot.isValid()) fillPhiVertices(iSetup, tracks, vertexes, beamSpot);
         // Fill phi from the secondary vertex collection
         //if (doPhiVertices && tracks.isValid() && vertexes.isValid() && beamSpot.isValid() && secondVertexes.isValid()) fillPhiVertices(secondVertexes, vertexes, beamSpot);
      }
   }

   if (nGoodMuons >= m_minNMuCut)
      tree_["muPogTree"]->Fill();

}


void MuonPogTreeProducer::fillGenInfo(const edm::Handle<std::vector<PileupSummaryInfo> > & puInfo,
      const edm::Handle<GenEventInfoProduct> & gen)
{

   muon_pog::GenInfo genInfo;

   genInfo.trueNumberOfInteractions     = -1.;
   genInfo.actualNumberOfInteractions   = -1.;
   genInfo.genWeight = gen->weight() ;

   std::vector<PileupSummaryInfo>::const_iterator puInfoIt  = puInfo->begin();
   std::vector<PileupSummaryInfo>::const_iterator puInfoEnd = puInfo->end();

   for(; puInfoIt != puInfoEnd; ++puInfoIt) 
   {
      int bx = puInfoIt->getBunchCrossing();

      if(bx == 0) 
      { 
         genInfo.trueNumberOfInteractions   = puInfoIt->getTrueNumInteractions();
         genInfo.actualNumberOfInteractions = puInfoIt->getPU_NumInteractions();
         continue;
      }
   }

   event_.genInfos.push_back(genInfo);

}


void MuonPogTreeProducer::fillGenParticles(const edm::Handle<reco::GenParticleCollection> & genParticles)
{

   unsigned int gensize = genParticles->size();

   // Do not record the initial protons
   for (unsigned int i=0; i<gensize; ++i) 
   {

      const reco::GenParticle& part = genParticles->at(i);

      muon_pog::GenParticle gensel;
      gensel.pdgId = part.pdgId();
      gensel.status = part.status();
      gensel.energy = part.energy();
      gensel.pt = part.pt();
      gensel.eta = part.eta();
      gensel.phi = part.phi();
      gensel.vx = part.vx();
      gensel.vy = part.vy();
      gensel.vz = part.vz();

      // Full set of GenFlags
      gensel.flags.clear();
      reco::GenStatusFlags statusflags = part.statusFlags();
      if (statusflags.flags_.size() == 15)
         for (unsigned int flag = 0; flag < statusflags.flags_.size(); ++flag)
            gensel.flags.push_back(statusflags.flags_[flag]);      

      gensel.mothers.clear();
      unsigned int nMothers = part.numberOfMothers();

      for (unsigned int iMother=0; iMother<nMothers; ++iMother) 
      {
         gensel.mothers.push_back(part.motherRef(iMother)->pdgId());
      }

      // Protect agains bug in genParticles (missing mother => first proton)
      if (i>=2 && nMothers==0) gensel.mothers.push_back(0);

      event_.genParticles.push_back(gensel);
   }

}

void MuonPogTreeProducer::fillHlt(const edm::Handle<edm::TriggerResults> & triggerResults, 
      const edm::Handle<trigger::TriggerEvent> & triggerEvent,
      const edm::TriggerNames & triggerNames)
{    

   for (unsigned int iTrig=0; iTrig<triggerNames.size(); ++iTrig) 
   {

      if (triggerResults->accept(iTrig)) 
      {
         std::string pathName = triggerNames.triggerName(iTrig);
         if (trigPathCut_ == "all" || pathName.find(trigPathCut_) != std::string::npos)
            event_.hlt.triggers.push_back(pathName);
      }
   }

   const trigger::size_type nFilters(triggerEvent->sizeFilters());

   for (trigger::size_type iFilter=0; iFilter!=nFilters; ++iFilter) 
   {

      std::string filterTag = triggerEvent->filterTag(iFilter).encode();

      if (trigFilterCut_ == "all" || filterTag.find(trigFilterCut_) != std::string::npos)
      {

         trigger::Keys objectKeys = triggerEvent->filterKeys(iFilter);
         const trigger::TriggerObjectCollection& triggerObjects(triggerEvent->getObjects());

         for (trigger::size_type iKey=0; iKey<objectKeys.size(); ++iKey) 
         {  
            trigger::size_type objKey = objectKeys.at(iKey);
            const trigger::TriggerObject& triggerObj(triggerObjects[objKey]);

            muon_pog::HLTObject hltObj;

            float trigObjPt = triggerObj.pt();
            float trigObjEta = triggerObj.eta();
            float trigObjPhi = triggerObj.phi();

            hltObj.filterTag = filterTag;

            hltObj.pt  = trigObjPt;
            hltObj.eta = trigObjEta;
            hltObj.phi = trigObjPhi;

            event_.hlt.objects.push_back(hltObj);

         }
      }
   }

}

void MuonPogTreeProducer::fillHlt(const edm::Handle<edm::TriggerResults> & triggerResults, 
      const edm::TriggerNames & triggerNames,
      const edm::Handle<std::vector<pat::TriggerObjectStandAlone> > & triggerObjects)
{  

   for (unsigned int iTrig=0; iTrig<triggerNames.size(); ++iTrig) 
   {

      if (triggerResults->accept(iTrig)) 
      {
         std::string pathName = triggerNames.triggerName(iTrig);
         if (trigPathCut_ == "all" || pathName.find(trigPathCut_) != std::string::npos)
            event_.hlt.triggers.push_back(pathName);
      }
   }


   for (pat::TriggerObjectStandAlone to: *triggerObjects){

      to.unpackPathNames(triggerNames);

      float to_pt = to.pt();
      float to_eta = to.eta();
      float to_phi = to.phi();

      trigger::size_type nFilters = to.filterLabels().size();
      for (trigger::size_type iFilter=0; iFilter<nFilters; iFilter++){

         muon_pog::HLTObject hltObj;

         hltObj.pt = to_pt;
         hltObj.eta = to_eta;
         hltObj.phi = to_phi;

         hltObj.filterTag = to.filterLabels()[iFilter];
         event_.hlt.objects.push_back(hltObj);
      }
   }


}

void MuonPogTreeProducer::fillL1(const edm::Handle<l1t::MuonBxCollection> & l1MuonBxColl)
{

   for (int ibx = l1MuonBxColl->getFirstBX(); ibx <= l1MuonBxColl->getLastBX(); ++ibx) 
   {
      for (auto l1MuIt = l1MuonBxColl->begin(ibx); l1MuIt != l1MuonBxColl->end(ibx); ++l1MuIt)
      {

         muon_pog::L1Muon l1part;
         l1part.pt = l1MuIt->pt();
         l1part.eta = l1MuIt->eta();
         l1part.phi = l1MuIt->phi();
         l1part.charge = l1MuIt->hwChargeValid() ? l1MuIt->charge() : 0;

         l1part.quality = l1MuIt->hwQual();
         l1part.bx = ibx;

         l1part.tfIndex = l1MuIt->tfMuonIndex();

         event_.l1muons.push_back(l1part);

      }
   }
}



void MuonPogTreeProducer::fillPV(const edm::Handle<std::vector<reco::Vertex> > & vertexes)
{

   int nVtx = 0;

   std::vector<reco::Vertex>::const_iterator vertexIt  = vertexes->begin();
   std::vector<reco::Vertex>::const_iterator vertexEnd = vertexes->end();

   for (; vertexIt != vertexEnd; ++vertexIt) 
   {

      const reco::Vertex& vertex = *vertexIt;

      if (!vertex.isValid()) continue;
      ++nVtx;

      if (vertexIt == vertexes->begin()) 
      {
         event_.primaryVertex[0] = vertex.x();
         event_.primaryVertex[1] = vertex.y();
         event_.primaryVertex[2] = vertex.z();

         for (unsigned int ix=0; ix<3; ++ix) 
         {
            for (unsigned int iy=0; iy<3; ++iy) 
            {
               event_.cov_primaryVertex[ix][iy] = vertex.covariance(ix,iy);
            }
         }
      }
   }

   event_.nVtx = nVtx;

}

//void MuonPogTreeProducer::fillSimMuonInfo(const edm::Handle<edm::ValueMap<reco::MuonSimInfo>>& muonSimCollection){

void MuonPogTreeProducer::fillSimMuonInfo(const edm::Handle<edm::ValueMap<reco::MuonSimInfo>>& muonSimCollection){


   muon_pog::simMuon ntupleMu; 

   for (size_t iMu=0; iMu<event_.muons.size(); iMu++){
      unsigned int Muon_index = event_.muons.at(iMu).index;
      reco::MuonSimInfo SimMuon = muonSimCollection->get(Muon_index);
      ntupleMu.index = Muon_index;
      ntupleMu.simPdgId = SimMuon.pdgId;
      ntupleMu.simMotherPdgId = SimMuon.motherPdgId;
      ntupleMu.simFlavour = SimMuon.flavour;
      ntupleMu.simType = SimMuon.primaryClass;
      ntupleMu.simBX = SimMuon.tpBX;
      ntupleMu.pt = SimMuon.p4.pt();
      ntupleMu.eta = SimMuon.p4.eta();
      ntupleMu.phi = SimMuon.p4.phi(); 
   }

   event_.sim_muons.push_back(ntupleMu);
}

Int_t MuonPogTreeProducer::fillMuons(const edm::Handle<edm::View<reco::Muon> > & muons,
      const edm::Handle<std::vector<reco::Vertex> > & vertexes,
      const edm::Handle<reco::BeamSpot> & beamSpot){

   edm::View<reco::Muon>::const_iterator muonIt  = muons->begin();
   edm::View<reco::Muon>::const_iterator muonEnd = muons->end();

   unsigned int Muon_index = 0;

   for (; muonIt != muonEnd; ++muonIt, ++Muon_index) {

      const reco::Muon& mu = (*muonIt);

      bool isGlobal      = mu.isGlobalMuon();
      bool isTracker     = mu.isTrackerMuon();
      bool isTrackerArb  = muon::isGoodMuon(mu, muon::TrackerMuonArbitrated); 
      bool isRPC         = mu.isRPCMuon();
      bool isStandAlone  = mu.isStandAloneMuon();
      bool isPF          = mu.isPFMuon();

      bool hasInnerTrack = !mu.innerTrack().isNull();
      bool hasTunePTrack = !mu.tunePMuonBestTrack().isNull();
      bool hasPickyTrack = !mu.pickyTrack().isNull();
      bool hasDytTrack = !mu.dytTrack().isNull();
      bool hasTpfmsTrack = !mu.tpfmsTrack().isNull();

      if (miniAODRun) {
         hasDytTrack = false; // gaurd against CMSSW error
         hasPickyTrack = false;
         hasDytTrack = false;
         hasTpfmsTrack = false;
      }

      muon_pog::Muon ntupleMu;

      ntupleMu.index = Muon_index;

      ntupleMu.pt     = mu.pt();
      ntupleMu.eta    = mu.eta();
      ntupleMu.phi    = mu.phi();
      ntupleMu.charge = mu.charge();

      ntupleMu.matchedKs = 0;
      ntupleMu.matchedPhi = 0;


      ntupleMu.fits.push_back(muon_pog::MuonFit(mu.pt(),mu.eta(),mu.phi(),
               mu.charge(),mu.muonBestTrack()->ptError()));



      ntupleMu.fits.push_back(muon_pog::MuonFit(hasInnerTrack ? mu.innerTrack()->pt()  : -1000.,
               hasInnerTrack ? mu.innerTrack()->eta() : -1000.,
               hasInnerTrack ? mu.innerTrack()->phi() : -1000.,
               hasInnerTrack ? mu.innerTrack()->charge()  : -1000.,
               hasInnerTrack ? mu.innerTrack()->ptError() : -1000.));


      ntupleMu.fits.push_back(muon_pog::MuonFit(isStandAlone ? mu.outerTrack()->pt()  : -1000.,
               isStandAlone ? mu.outerTrack()->eta() : -1000.,
               isStandAlone ? mu.outerTrack()->phi() : -1000.,
               isStandAlone ? mu.outerTrack()->charge()  : -1000.,
               isStandAlone ? mu.outerTrack()->ptError() : -1000.));


      ntupleMu.fits.push_back(muon_pog::MuonFit(isGlobal ? mu.globalTrack()->pt()  : -1000.,
               isGlobal ? mu.globalTrack()->eta() : -1000.,
               isGlobal ? mu.globalTrack()->phi() : -1000.,
               isGlobal ? mu.globalTrack()->charge()  : -1000.,
               isGlobal ? mu.globalTrack()->ptError() : -1000.));

      ntupleMu.fits.push_back(muon_pog::MuonFit(hasTunePTrack ? mu.tunePMuonBestTrack()->pt()  : -1000.,
               hasTunePTrack ? mu.tunePMuonBestTrack()->eta() : -1000.,
               hasTunePTrack ? mu.tunePMuonBestTrack()->phi() : -1000.,
               hasTunePTrack ? mu.tunePMuonBestTrack()->charge()  : -1000.,
               hasTunePTrack ? mu.tunePMuonBestTrack()->ptError() : -1000.));

      ntupleMu.fits.push_back(muon_pog::MuonFit(hasPickyTrack ? mu.pickyTrack()->pt()  : -1000.,
               hasPickyTrack ? mu.pickyTrack()->eta() : -1000.,
               hasPickyTrack ? mu.pickyTrack()->phi() : -1000.,
               hasPickyTrack ? mu.pickyTrack()->charge()  : -1000.,
               hasPickyTrack ? mu.pickyTrack()->ptError() : -1000.));

      ntupleMu.fits.push_back(muon_pog::MuonFit(hasDytTrack ? mu.dytTrack()->pt()  : -1000.,
               hasDytTrack ? mu.dytTrack()->eta() : -1000.,
               hasDytTrack ? mu.dytTrack()->phi() : -1000.,
               hasDytTrack ? mu.dytTrack()->charge()  : -1000.,
               hasDytTrack ? mu.dytTrack()->ptError() : -1000.));

      ntupleMu.fits.push_back(muon_pog::MuonFit(hasTpfmsTrack ? mu.tpfmsTrack()->pt()  : -1000.,
               hasTpfmsTrack ? mu.tpfmsTrack()->eta() : -1000.,
               hasTpfmsTrack ? mu.tpfmsTrack()->phi() : -1000.,
               hasTpfmsTrack ? mu.tpfmsTrack()->charge()  : -1000.,
               hasTpfmsTrack ? mu.tpfmsTrack()->ptError() : -1000.));


      // Detector Based Isolation
      reco::MuonIsolation detIso03 = mu.isolationR03();

      ntupleMu.trackerIso = detIso03.sumPt;
      ntupleMu.EMCalIso   = detIso03.emEt;
      ntupleMu.HCalIso    = detIso03.hadEt;

      // PF Isolation
      reco::MuonPFIsolation pfIso04 = mu.pfIsolationR04();
      reco::MuonPFIsolation pfIso03 = mu.pfIsolationR03();

      ntupleMu.chargedHadronIso   = pfIso04.sumChargedHadronPt;
      ntupleMu.chargedHadronIsoPU = pfIso04.sumPUPt; 
      ntupleMu.neutralHadronIso   = pfIso04.sumNeutralHadronEt;
      ntupleMu.photonIso          = pfIso04.sumPhotonEt;

      ntupleMu.isGlobal     = isGlobal ? 1 : 0;	
      ntupleMu.isTracker    = isTracker ? 1 : 0;	
      ntupleMu.isTrackerArb = isTrackerArb ? 1 : 0;	
      ntupleMu.isRPC        = isRPC ? 1 : 0;
      ntupleMu.isStandAlone = isStandAlone ? 1 : 0;
      ntupleMu.isPF         = isPF ? 1 : 0;

      ntupleMu.nHitsGlobal     = isGlobal     ? mu.globalTrack()->numberOfValidHits() : -999;	
      ntupleMu.nHitsTracker    = isTracker    ? mu.innerTrack()->numberOfValidHits()  : -999;	
      ntupleMu.nHitsStandAlone = isStandAlone ? mu.outerTrack()->numberOfValidHits()  : -999;

      ntupleMu.glbNormChi2              = isGlobal      ? mu.globalTrack()->normalizedChi2() : -999; 
      ntupleMu.trkNormChi2	             = hasInnerTrack ? mu.innerTrack()->normalizedChi2()  : -999; 
      ntupleMu.trkMuonMatchedStations   = isTracker     ? mu.numberOfMatchedStations()       : -999; 
      ntupleMu.glbMuonValidHits	       = isGlobal      ? mu.globalTrack()->hitPattern().numberOfValidMuonHits()       : -999; 
      ntupleMu.trkPixelValidHits	       = hasInnerTrack ? mu.innerTrack()->hitPattern().numberOfValidPixelHits()       : -999; 
      ntupleMu.trkPixelLayersWithMeas   = hasInnerTrack ? mu.innerTrack()->hitPattern().pixelLayersWithMeasurement()   : -999; 
      ntupleMu.trkTrackerLayersWithMeas = hasInnerTrack ? mu.innerTrack()->hitPattern().trackerLayersWithMeasurement() : -999; 

      ntupleMu.bestMuPtErr              = mu.muonBestTrack()->ptError(); 

      ntupleMu.trkValidHitFrac          = hasInnerTrack           ? mu.innerTrack()->validFraction()       : -999; 
      ntupleMu.trkStaChi2               = isGlobal                ? mu.combinedQuality().chi2LocalPosition : -999; 
      ntupleMu.trkKink                  = isGlobal                ? mu.combinedQuality().trkKink           : -999;
      ntupleMu.isTrkMuOST               = muon::isGoodMuon(mu, muon::TMOneStationTight) ? 1 : 0; 
      ntupleMu.isTrkHP                  = hasInnerTrack && mu.innerTrack()->quality(reco::TrackBase::highPurity) ? 1 : 0; 

      // TrackerMuonId quantities
      ntupleMu.muSegmComp                 = hasInnerTrack ? muon::segmentCompatibility(mu)         : -999; 
      ntupleMu.muCaloComp                 = hasInnerTrack ? muon::caloCompatibility(mu)            : -999;
      ntupleMu.innerNormChi2              = hasInnerTrack ? mu.innerTrack()->normalizedChi2()       : -999;
      ntupleMu.innerNValidHits            = hasInnerTrack ? mu.innerTrack()->hitPattern().numberOfValidTrackerHits() : -999;
      ntupleMu.innerNLostTrackerHits      = hasInnerTrack ? mu.innerTrack()->hitPattern().numberOfLostTrackerHits(HitPattern::TRACK_HITS): -999;
      ntupleMu.innerNInnerLostTrackerHits = hasInnerTrack ? mu.innerTrack()->hitPattern().numberOfLostTrackerHits(HitPattern::MISSING_INNER_HITS): -999;
      ntupleMu.innerNOuterLostTrackerHits = hasInnerTrack ? mu.innerTrack()->hitPattern().numberOfLostTrackerHits(HitPattern::MISSING_OUTER_HITS): -999;
      ntupleMu.innerPtErrPt               = hasInnerTrack ? mu.innerTrack()->ptError()/mu.innerTrack()->pt(): -999;
      ntupleMu.innerNMatches              = hasInnerTrack ? mu.numberOfMatches() : -999;
      ntupleMu.innerNMatchedStations      = hasInnerTrack ? mu.numberOfMatchedStations(): -999;
      ntupleMu.innerExpectedMatchedStations = hasInnerTrack ? mu.expectedNnumberOfMatchedStations(): -999;

      // TrackerMuonId (calo quantities)
      ntupleMu.em = mu.calEnergy().em;
      ntupleMu.emS9 = mu.calEnergy().emS9;
      ntupleMu.emS25 = mu.calEnergy().emS25;
      ntupleMu.had = mu.calEnergy().had;
      ntupleMu.hadS9 = mu.calEnergy().hadS9;

      // Global muon Id quantities
      ntupleMu.updatedSta = isGlobal ? mu.combinedQuality().updatedSta: -999;
      ntupleMu.glbKink = isGlobal ? mu.combinedQuality().glbKink: -999;
      ntupleMu.trkRelChi2 = isGlobal ? mu.combinedQuality().trkRelChi2: -999;
      ntupleMu.staRelChi2 = isGlobal ? mu.combinedQuality().staRelChi2: -999;
      ntupleMu.chi2LocalMomentum = isGlobal ? mu.combinedQuality().chi2LocalMomentum: -999;
      ntupleMu.localDistance = isGlobal ? mu.combinedQuality().localDistance: -999;
      ntupleMu.globalDeltaEtaPhi = isGlobal ? mu.combinedQuality().globalDeltaEtaPhi: -999;
      ntupleMu.tightMatch = isGlobal ? mu.combinedQuality().tightMatch: -999;
      ntupleMu.glbTrackProbability = isGlobal ? mu.combinedQuality().glbTrackProbability: -999;
      ntupleMu.vhitcomb = isGlobal ? getVhitsComb(mu): -999;

      if ( mu.isMatchesValid() && ntupleMu.isTrackerArb )
      {
         for ( reco::MuonChamberMatch match : mu.matches() )
         {
            muon_pog::ChambMatch ntupleMatch;

            if ( getMuonChamberId(match.id,
                     ntupleMatch.type,ntupleMatch.r,
                     ntupleMatch.phi,ntupleMatch.eta)
               )
            {

               ntupleMatch.x = mu.trackX(match.station(),match.detector());
               ntupleMatch.y = mu.trackY(match.station(),match.detector());
               ntupleMatch.dXdZ = mu.trackDxDz(match.station(),match.detector());
               ntupleMatch.dYdZ = mu.trackDyDz(match.station(),match.detector());

               ntupleMatch.errxTk = mu.trackXErr(match.station(),match.detector());
               ntupleMatch.erryTk = mu.trackYErr(match.station(),match.detector());

               ntupleMatch.errDxDzTk = mu.trackDxDzErr(match.station(),match.detector());
               ntupleMatch.errDyDzTk = mu.trackDyDzErr(match.station(),match.detector());

               ntupleMatch.dx = mu.dX(match.station(),match.detector());
               ntupleMatch.dy = mu.dY(match.station(),match.detector());
               ntupleMatch.dDxDz = mu.dDxDz(match.station(),match.detector());
               ntupleMatch.dDyDz = mu.dDxDz(match.station(),match.detector());

               ntupleMatch.errxSeg = mu.segmentXErr(match.station(),match.detector());
               ntupleMatch.errySeg = mu.segmentYErr(match.station(),match.detector());
               ntupleMatch.errDxDzSeg = mu.segmentDxDzErr(match.station(),match.detector());
               ntupleMatch.errDyDzSeg = mu.segmentDyDzErr(match.station(),match.detector());

               ntupleMu.matches.push_back(ntupleMatch);
            }
         }
      }

      ntupleMu.dxyBest  = -999; 
      ntupleMu.dzBest   = -999; 
      ntupleMu.dxyInner = -999; 
      ntupleMu.dzInner  = -999; 

      ntupleMu.isoPflow04 = (pfIso04.sumChargedHadronPt+ 
            std::max(0.,pfIso04.sumPhotonEt+pfIso04.sumNeutralHadronEt - 0.5*pfIso04.sumPUPt)) / mu.pt();

      ntupleMu.isoPflow03 = (pfIso03.sumChargedHadronPt+ 
            std::max(0.,pfIso03.sumPhotonEt+pfIso03.sumNeutralHadronEt - 0.5*pfIso03.sumPUPt)) / mu.pt();

      double dxybs = isGlobal ? mu.globalTrack()->dxy(beamSpot->position()) :
         hasInnerTrack ? mu.innerTrack()->dxy(beamSpot->position()) : -1000;
      double dzbs  = isGlobal ? mu.globalTrack()->dz(beamSpot->position()) :
         hasInnerTrack ? mu.innerTrack()->dz(beamSpot->position()) : -1000;

      double dxy = -1000.;
      double dz  = -1000.;

      ntupleMu.isSoft    = 0;	  
      ntupleMu.isTight   = 0;	  
      ntupleMu.isHighPt  = 0;
      ntupleMu.isLoose   = muon::isLooseMuon(mu)  ? 1 : 0;	  
      ntupleMu.isMedium  = muon::isMediumMuon(mu) ? 1 : 0;	  

      if (vertexes->size() > 0)
      {
         const reco::Vertex & vertex = vertexes->at(0);

         dxy = isGlobal ? mu.globalTrack()->dxy(vertex.position()) :
            hasInnerTrack ? mu.innerTrack()->dxy(vertex.position()) : -1000;
         dz = isGlobal ? mu.globalTrack()->dz(vertex.position()) :
            hasInnerTrack ? mu.innerTrack()->dz(vertex.position()) : -1000;

         ntupleMu.dxyBest  = mu.muonBestTrack()->dxy(vertex.position()); 
         ntupleMu.dzBest   = mu.muonBestTrack()->dz(vertex.position()); 
         if(hasInnerTrack) { 
            ntupleMu.dxyInner = mu.innerTrack()->dxy(vertex.position()); 
            ntupleMu.dzInner  = mu.innerTrack()->dz(vertex.position()); 
         } 

         ntupleMu.isSoft    = muon::isSoftMuon(mu,vertex)   ? 1 : 0;	  
         ntupleMu.isTight   = muon::isTightMuon(mu,vertex)  ? 1 : 0;	  
         ntupleMu.isHighPt  = muon::isHighPtMuon(mu,vertex) ? 1 : 0;

      }

      ntupleMu.dxy    = dxy;
      ntupleMu.dz     = dz;
      ntupleMu.edxy   = isGlobal ? mu.globalTrack()->dxyError() : hasInnerTrack ? mu.innerTrack()->dxyError() : -1000;
      ntupleMu.edz    = isGlobal ? mu.globalTrack()->dzError()  : hasInnerTrack ? mu.innerTrack()->dzError() : -1000;

      ntupleMu.dxybs  = dxybs;
      ntupleMu.dzbs   = dzbs;

      if(mu.isTimeValid()) { 
         ntupleMu.muonTimeDof = mu.time().nDof; 
         ntupleMu.muonTime    = mu.time().timeAtIpInOut; 
         ntupleMu.muonTimeErr = mu.time().timeAtIpInOutErr; 
      } 
      else { 
         ntupleMu.muonTimeDof = -999; 
         ntupleMu.muonTime    = -999; 
         ntupleMu.muonTimeErr = -999; 
      } 

      if(mu.rpcTime().nDof > 0) { 
         ntupleMu.muonRpcTimeDof = mu.rpcTime().nDof; 
         ntupleMu.muonRpcTime    = mu.rpcTime().timeAtIpInOut; 
         ntupleMu.muonRpcTimeErr = mu.rpcTime().timeAtIpInOutErr; 
      } 
      else { 
         ntupleMu.muonRpcTimeDof = -999; 
         ntupleMu.muonRpcTime    = -999; 
         ntupleMu.muonRpcTimeErr = -999; 
      } 

      // asking for a TRK or GLB muon with minimal pT cut
      // ignoring STA muons in this logic
      if ( m_minMuPtCut < 0 ||
            (
             (isTracker || isGlobal || isStandAlone) &&
             (ntupleMu.fitPt(muon_pog::MuonFitType::DEFAULT) > m_minMuPtCut ||
              ntupleMu.fitPt(muon_pog::MuonFitType::GLB)     > m_minMuPtCut ||
              ntupleMu.fitPt(muon_pog::MuonFitType::TUNEP)   > m_minMuPtCut ||
              ntupleMu.fitPt(muon_pog::MuonFitType::INNER)   > m_minMuPtCut ||
              ntupleMu.fitPt(muon_pog::MuonFitType::STA)     > m_minMuPtCut ||
              ntupleMu.fitPt(muon_pog::MuonFitType::PICKY)   > m_minMuPtCut ||
              ntupleMu.fitPt(muon_pog::MuonFitType::DYT)     > m_minMuPtCut ||
              ntupleMu.fitPt(muon_pog::MuonFitType::TPFMS)   > m_minMuPtCut)
            )
         )
      {
         event_.muons.push_back(ntupleMu);

         std::vector<Float_t> PTs = {ntupleMu.fitPt(muon_pog::MuonFitType::DEFAULT),
            ntupleMu.fitPt(muon_pog::MuonFitType::GLB),
            ntupleMu.fitPt(muon_pog::MuonFitType::TUNEP),
            ntupleMu.fitPt(muon_pog::MuonFitType::INNER),
            ntupleMu.fitPt(muon_pog::MuonFitType::PICKY),
            ntupleMu.fitPt(muon_pog::MuonFitType::DYT),
            ntupleMu.fitPt(muon_pog::MuonFitType::TPFMS)};
         eventId_.maxPTs.push_back(*std::max_element(PTs.begin(), PTs.end()));
      }

   }

   return event_.muons.size();

}

float MuonPogTreeProducer::getVhitsComb(const reco::Muon& mu){

   unsigned int dt1(0),dt2(0),dt3(0),dt4(0);
   unsigned int rpc1(0),rpc2(0),rpc3(0),rpc4(0);
   unsigned int csc1(0),csc2(0),csc3(0),csc4(0);
   float comb(0);
   const reco::HitPattern &pattern = mu.globalTrack()->hitPattern();
   for (int i=0;i<pattern.numberOfAllHits(reco::HitPattern::TRACK_HITS);i++)
   { 
      uint32_t hit = pattern.getHitPattern(reco::HitPattern::TRACK_HITS,i);
      if (pattern.validHitFilter(hit) != 1) {continue;}
      if (pattern.getMuonStation(hit) == 1)
      { 
         if (pattern.muonDTHitFilter(hit))  dt1++;
         if (pattern.muonRPCHitFilter(hit)) rpc1++;
         if (pattern.muonCSCHitFilter(hit)) csc1++;
      }
      else if (pattern.getMuonStation(hit) == 2)
      { 
         if (pattern.muonDTHitFilter(hit))  dt2++;
         if (pattern.muonRPCHitFilter(hit)) rpc2++;
         if (pattern.muonCSCHitFilter(hit)) csc2++;
      }
      else if (pattern.getMuonStation(hit) == 3)
      { 
         if (pattern.muonDTHitFilter(hit))  dt3++;
         if (pattern.muonRPCHitFilter(hit)) rpc3++;
         if (pattern.muonCSCHitFilter(hit)) csc3++;
      }
      else if (pattern.getMuonStation(hit) == 4)
      { 
         if (pattern.muonDTHitFilter(hit))  dt4++;
         if (pattern.muonRPCHitFilter(hit)) rpc4++;
         if (pattern.muonCSCHitFilter(hit)) csc4++;
      }    
   }
   comb = (dt1+dt2+dt3+dt4)/2. + (rpc1+rpc2+rpc3+rpc4);
   csc1>6 ? comb+=6 : comb+=csc1;
   csc2>6 ? comb+=6 : comb+=csc2;
   csc3>6 ? comb+=6 : comb+=csc3;
   csc4>6 ? comb+=6 : comb+=csc4; 

   return comb;
}

bool MuonPogTreeProducer::getMuonChamberId(DetId & id, muon_pog::MuonDetType & det,
      Int_t & r, Int_t & phi, Int_t & eta) const
{
   if (id.det() == DetId::Muon && id.subdetId() == MuonSubdetId::DT)
   {
      DTChamberId dtId(id.rawId());  

      det = muon_pog::MuonDetType::DT;
      r   = dtId.station();
      phi = dtId.sector();
      eta = dtId.wheel();

      return true;
   }

   if (id.det() == DetId::Muon && id.subdetId() == MuonSubdetId::CSC)
   {
      CSCDetId cscId(id.rawId());

      det = muon_pog::MuonDetType::CSC;
      r   = cscId.station() * cscId.zendcap();
      phi = cscId.chamber();
      eta = cscId.ring();

      return true;
   }

   return false;
}

void MuonPogTreeProducer::fillKsVertices(const edm::Handle<std::vector<reco::VertexCompositeCandidate>>& kshorts,
      const edm::Handle<std::vector<reco::Vertex>>& vertexes,
      const edm::Handle<reco::BeamSpot>& beamSpot){

   // Get beam spot 
   const reco::BeamSpot& bs = *beamSpot;
   VertexState BSstate(bs);

   // get primary vertex
   const reco::Vertex& pv = vertexes->at(0);
   TVector3 pvVec(pv.position().x(),
         pv.position().y(),
         pv.position().z());

   std::vector<reco::VertexCompositeCandidate>::const_iterator ksIt = kshorts->begin();
   std::vector<reco::VertexCompositeCandidate>::const_iterator ksEnd = kshorts->end();

   unsigned int ksIndex = 0;

   for (; ksIt != ksEnd; ++ksIt, ++ksIndex){

      const reco::VertexCompositeCandidate& ks = (*ksIt);
      if (ks.numberOfDaughters()!=2) continue;

      const reco::Candidate& pi1 = *ks.daughter(0);
      const reco::Candidate& pi2 = *ks.daughter(1);

      Float_t pi1_eta = pi1.eta();
      Float_t pi1_phi = pi1.phi();
      Float_t pi2_eta = pi2.eta();
      Float_t pi2_phi = pi2.phi();

      TVector3 p_ks;
      p_ks.SetXYZ(ks.px(), ks.py(), ks.pz());

      //Bool_t phiMatchedToMuon = false;
      int muon1_index = -1;
      int muon2_index = -1;

      for (size_t iMu=0; iMu<event_.muons.size(); iMu++){

         Float_t mu_eta = event_.muons.at(iMu).eta;
         Float_t mu_phi = event_.muons.at(iMu).phi;
         Float_t mu_pt = event_.muons.at(iMu).pt;

         float dr_pi1_mu = reco::deltaR(pi1_eta, pi1_phi, mu_eta, mu_phi);
         float dr_pi2_mu = reco::deltaR(pi2_eta, pi2_phi, mu_eta, mu_phi);

         if ( dr_pi1_mu < 0.01 && fabs(1-pi1.pt()/mu_pt) < 0.1) {
            event_.muons.at(iMu).matchedKs = 1; // muon matched to track 1
            muon1_index = iMu;
         }
         if ( dr_pi2_mu < 0.01 && fabs(1-pi2.pt()/mu_pt) < 0.1) {
            event_.muons.at(iMu).matchedKs = 1; // muon matched to track 2
            muon2_index = iMu;
         }
      }

      math::Error<3>::type covarianceMatrix;

      for (int i = 0; i < 3; i++){
         for (int j = 0; j < 3; j++) {
            covarianceMatrix(i, j) = ks.vertexCovariance(i, j);
            covarianceMatrix(j, i) = ks.vertexCovariance(i, j);
         }
      }

      reco::Vertex ksVertex(ks.vertex(), covarianceMatrix, ks.vertexChi2(), ks.vertexNdof(), ks.numberOfDaughters());

      TVector3 svVec(ks.vertex().x(),
            ks.vertex().y(),
            ks.vertex().z());

      VertexDistanceXY vertexTool2D;
      VertexDistance3D vertexTool3D;

      Measurement1D BSSV_distanceXY, PVSV_distanceXY, PVSV_distanceXYZ;
      BSSV_distanceXY = vertexTool2D.distance(BSstate, ksVertex);
      PVSV_distanceXY = vertexTool2D.distance(pv, ksVertex);
      PVSV_distanceXYZ = vertexTool3D.distance(pv, ksVertex);


      double Lxy = PVSV_distanceXY.value();
      double vertexSig = PVSV_distanceXY.significance();
      double totalChi2 = ks.vertexChi2();
      double normChi2 = ks.vertexNormalizedChi2();
      double cosTheta = (pvVec-svVec).Dot(p_ks)/((pvVec-svVec).Mag()*p_ks.Mag());

      // Constraints on the vertex fit
      if (totalChi2<0.) continue;
      if (normChi2>3) continue;
      if (vertexSig<5) continue;
      if (Lxy>4) continue;

      muon_pog::CompositeVertex ks_;

      ks_.mass = ks.p4().M();
      ks_.reco_pdgId = 312;
      ks_.bssv_dxy = BSSV_distanceXY.value();
      ks_.bssv_dxyErr = BSSV_distanceXY.error();
      ks_.bssv_dxySig = BSSV_distanceXY.significance();
      ks_.pvsv_dxy = Lxy;
      ks_.pvsv_dxyError = PVSV_distanceXY.error();
      ks_.pvsv_dxySig = vertexSig;
      ks_.pvsv_3d = PVSV_distanceXYZ.value();
      ks_.pvsv_3dError = PVSV_distanceXYZ.error();
      ks_.pvsv_3dSig = PVSV_distanceXYZ.significance();
      ks_.pvsv_normChi2 = normChi2;
      ks_.pvsv_totalChi2 = totalChi2;
      ks_.pvsv_cosTheta = cosTheta;
      ks_.muon1_index = muon1_index;
      ks_.muon2_index = muon2_index;

      event_.kshorts.push_back(ks_);
   }
}

void MuonPogTreeProducer::fillKsVertices(const edm::Handle<std::vector<reco::VertexCompositePtrCandidate>>& kshorts,
      const edm::Handle<std::vector<reco::Vertex>>& vertexes,
      const edm::Handle<reco::BeamSpot>& beamSpot){

   // Get beam spot 
   const reco::BeamSpot& bs = *beamSpot;
   VertexState BSstate(bs);

   // get primary vertex
   const reco::Vertex& pv = vertexes->at(0);
   TVector3 pvVec(pv.position().x(),
         pv.position().y(),
         pv.position().z());


   std::vector<reco::VertexCompositePtrCandidate>::const_iterator ksIt = kshorts->begin();
   std::vector<reco::VertexCompositePtrCandidate>::const_iterator ksEnd = kshorts->end();

   unsigned int ksIndex = 0;

   for (; ksIt != ksEnd; ++ksIt, ++ksIndex){

      const reco::VertexCompositePtrCandidate& ks = (*ksIt);
      if (ks.numberOfDaughters()!=2) continue;

      const reco::Candidate& pi1 = *ks.daughter(0);
      const reco::Candidate& pi2 = *ks.daughter(1);

      Float_t pi1_eta = pi1.eta();
      Float_t pi1_phi = pi1.phi();
      Float_t pi2_eta = pi2.eta();
      Float_t pi2_phi = pi2.phi();

      TVector3 p_ks;
      p_ks.SetXYZ(ks.px(), ks.py(), ks.pz());

      //Bool_t phiMatchedToMuon = false;
      int muon1_index = -1;
      int muon2_index = -1;

      for (size_t iMu=0; iMu<event_.muons.size(); iMu++){

         Float_t mu_eta = event_.muons.at(iMu).eta;
         Float_t mu_phi = event_.muons.at(iMu).phi;
         Float_t mu_pt = event_.muons.at(iMu).pt;

         float dr_pi1_mu = reco::deltaR(pi1_eta, pi1_phi, mu_eta, mu_phi);
         float dr_pi2_mu = reco::deltaR(pi2_eta, pi2_phi, mu_eta, mu_phi);

         if ( dr_pi1_mu < 0.01 && fabs(1-pi1.pt()/mu_pt) < 0.1) {
            event_.muons.at(iMu).matchedKs = 1; // muon matched to track 1
            muon1_index = iMu;
         }
         if ( dr_pi2_mu < 0.01 && fabs(1-pi2.pt()/mu_pt) < 0.1) {
            event_.muons.at(iMu).matchedKs = 1; // muon matched to track 2
            muon2_index = iMu;
         }
      }

      math::Error<3>::type covarianceMatrix;

      for (int i = 0; i < 3; i++){
         for (int j = 0; j < 3; j++) {
            covarianceMatrix(i, j) = ks.vertexCovariance(i, j);
            covarianceMatrix(j, i) = ks.vertexCovariance(i, j);
         }
      }

      reco::Vertex ksVertex(ks.vertex(), covarianceMatrix, ks.vertexChi2(), ks.vertexNdof(), ks.numberOfDaughters());


      TVector3 svVec(ks.vertex().x(),
            ks.vertex().y(),
            ks.vertex().z());

      // find distances with respect to the beam spot


      VertexDistanceXY vertexTool2D;
      VertexDistance3D vertexTool3D;

      Measurement1D BSSV_distanceXY, PVSV_distanceXY, PVSV_distanceXYZ;
      BSSV_distanceXY = vertexTool2D.distance(BSstate, ksVertex);
      PVSV_distanceXY = vertexTool2D.distance(pv, ksVertex);
      PVSV_distanceXYZ = vertexTool3D.distance(pv, ksVertex);

      double Lxy = PVSV_distanceXY.value();
      double vertexSig = PVSV_distanceXY.significance();
      double totalChi2 = ks.vertexChi2();
      double normChi2 = ks.vertexNormalizedChi2();
      double cosTheta = (pvVec-svVec).Dot(p_ks)/((pvVec-svVec).Mag()*p_ks.Mag());


      // Constraints on the vertex fit
      if (totalChi2<0.) continue;
      if (normChi2>3) continue;
      if (vertexSig<5) continue;
      if (Lxy>4) continue;

      muon_pog::CompositeVertex ks_;

      ks_.mass = ks.p4().M();
      ks_.reco_pdgId = 312;
      ks_.bssv_dxy = BSSV_distanceXY.value();
      ks_.bssv_dxyErr = BSSV_distanceXY.error();
      ks_.bssv_dxySig = BSSV_distanceXY.significance();
      ks_.pvsv_dxy = Lxy;
      ks_.pvsv_dxyError = PVSV_distanceXY.error();
      ks_.pvsv_dxySig = vertexSig;
      ks_.pvsv_3d = PVSV_distanceXYZ.value();
      ks_.pvsv_3dError = PVSV_distanceXYZ.error();
      ks_.pvsv_3dSig = PVSV_distanceXYZ.significance();
      ks_.pvsv_normChi2 = normChi2;
      ks_.pvsv_totalChi2 = totalChi2;
      ks_.pvsv_cosTheta = cosTheta;
      ks_.muon1_index = muon1_index;
      ks_.muon2_index = muon2_index;

      event_.kshorts.push_back(ks_);
   }
}

void MuonPogTreeProducer::fillPhiVertices(const edm::Handle<std::vector<reco::VertexCompositePtrCandidate>>& secVertexes,
      const edm::Handle<std::vector<reco::Vertex>>& vertexes,
      const edm::Handle<reco::BeamSpot>& beamSpot){

   // Get beam spot position
   const reco::BeamSpot& bs = *beamSpot;

   VertexState BSstate(bs);

   if(!BSstate.isValid()) return;

   // get primary vertex
   const reco::Vertex& pv = vertexes->at(0);
   TVector3 pvVec(pv.position().x(),
         pv.position().y(),
         pv.position().z());


   std::vector<reco::VertexCompositePtrCandidate>::const_iterator phiIt = secVertexes->begin();
   std::vector<reco::VertexCompositePtrCandidate>::const_iterator phiEnd = secVertexes->end();

   unsigned int phiIndex = 0;

   for (; phiIt != phiEnd; ++phiIt, ++phiIndex){

      const reco::VertexCompositePtrCandidate& phi = (*phiIt);
      if (phi.numberOfDaughters()!=2) continue;

      const reco::Candidate& pi1 = *phi.daughter(0);
      const reco::Candidate& pi2 = *phi.daughter(1);

      Float_t pi1_eta = pi1.eta();
      Float_t pi1_phi = pi1.phi();
      Float_t pi2_eta = pi2.eta();
      Float_t pi2_phi = pi2.phi();

      TVector3 p_phi;
      p_phi.SetXYZ(phi.px(), phi.py(), phi.pz());

      //Bool_t phiMatchedToMuon = false;
      int muon1_index = -1;
      int muon2_index = -1;

      for (size_t iMu=0; iMu<event_.muons.size(); iMu++){

         Float_t mu_eta = event_.muons.at(iMu).eta;
         Float_t mu_phi = event_.muons.at(iMu).phi;
         Float_t mu_pt = event_.muons.at(iMu).pt;

         float dr_pi1_mu = reco::deltaR(pi1_eta, pi1_phi, mu_eta, mu_phi);
         float dr_pi2_mu = reco::deltaR(pi2_eta, pi2_phi, mu_eta, mu_phi);

         if ( dr_pi1_mu < 0.01 && fabs(1-pi1.pt()/mu_pt) < 0.1) {
            event_.muons.at(iMu).matchedKs = 1; // muon matched to track 1
            muon1_index = iMu;
         }
         if ( dr_pi2_mu < 0.01 && fabs(1-pi2.pt()/mu_pt) < 0.1) {
            event_.muons.at(iMu).matchedKs = 1; // muon matched to track 2
            muon2_index = iMu;
         }
      }

      math::Error<3>::type covarianceMatrix;

      for (int i = 0; i < 3; i++){
         for (int j = 0; j < 3; j++) {
            covarianceMatrix(i, j) = phi.vertexCovariance(i, j);
            covarianceMatrix(j, i) = phi.vertexCovariance(i, j);
         }
      }

      reco::Vertex phiVertex(phi.vertex(), covarianceMatrix, phi.vertexChi2(), phi.vertexNdof(), phi.numberOfDaughters());

      TVector3 svVec(phi.vertex().x(),
            phi.vertex().y(),
            phi.vertex().z());

      // find distances with respect to the beam spot


      // Compute vertex distances
      VertexDistanceXY vertexTool2D;
      VertexDistance3D vertexTool3D;

      Measurement1D BSSV_distanceXY, PVSV_distanceXY, PVSV_distanceXYZ;
      BSSV_distanceXY = vertexTool2D.distance(BSstate, phiVertex);
      PVSV_distanceXY = vertexTool2D.distance(pv, phiVertex);
      PVSV_distanceXYZ = vertexTool3D.distance(pv, phiVertex);


      double Lxy = PVSV_distanceXY.value();
      double vertexSig = PVSV_distanceXY.significance();
      double totalChi2 = phi.vertexChi2();
      double normChi2 = phi.vertexNormalizedChi2();
      double cosTheta = (pvVec-svVec).Dot(p_phi)/((pvVec-svVec).Mag()*p_phi.Mag());


      // Constraints on the vertex fit
      if (totalChi2<0.) continue;
      if (normChi2>3) continue;
      if (vertexSig<5) continue;
      if (Lxy>4) continue;
      if ( phi.p4().M() < PHI_MASS_WINDOW_LOW || phi.p4().M() > PHI_MASS_WINDOW_HIGH ) continue;

      muon_pog::CompositeVertex phi_;

      phi_.mass = phi.p4().M();
      phi_.reco_pdgId = 312;
      phi_.bssv_dxy = BSSV_distanceXY.value();
      phi_.bssv_dxyErr = BSSV_distanceXY.error();
      phi_.bssv_dxySig = BSSV_distanceXY.significance();
      phi_.pvsv_dxy = Lxy;
      phi_.pvsv_dxyError = PVSV_distanceXY.error();
      phi_.pvsv_dxySig = vertexSig;
      phi_.pvsv_3d = PVSV_distanceXYZ.value();
      phi_.pvsv_3dError = PVSV_distanceXYZ.error();
      phi_.pvsv_3dSig = PVSV_distanceXYZ.significance();
      phi_.pvsv_normChi2 = normChi2;
      phi_.pvsv_totalChi2 = totalChi2;
      phi_.pvsv_cosTheta = cosTheta;
      phi_.muon1_index = muon1_index;
      phi_.muon2_index = muon2_index;

      event_.phis.push_back(phi_);
   }
}

void MuonPogTreeProducer::fillPhiVertices(const edm::EventSetup& iSetup, const edm::Handle<std::vector<reco::Track>>& tracks,
      const edm::Handle<std::vector<reco::Vertex>>& vertexes,
      const edm::Handle<reco::BeamSpot>& beamSpot){

   ESHandle<TransientTrackBuilder> bField;
   iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", bField);

   std::vector<reco::Track>::const_iterator trackIt = tracks->begin();
   std::vector<reco::Track>::const_iterator trackEnd = tracks->end();

   const reco::Vertex& pv = vertexes->at(0);
   TVector3 pvVec(pv.position().x(), pv.position().y(), pv.position().z());

   const reco::BeamSpot& bs = *beamSpot;

   VertexState BSstate(bs);

   math::XYZPoint bsPoint = bs.position();

   TVector3 bsVec(bsPoint.x(), bsPoint.y(), bsPoint.z());
   VertexDistanceXY vertexTool2D;
   VertexDistance3D vertexTool3D;

   unsigned int iTrack(0), jTrack(0);

   for (; trackIt != trackEnd; ++trackIt, iTrack++){

      const reco::Track& track1 = (*trackIt);
      double track1_ipError = track1.dxyError();
      double track1_ipValue = std::abs(track1.dxy(bsPoint));

      //if (abs(track1.eta())>2.4) continue;
      //if (track1.hitPattern().trackerLayersWithMeasurement()<4 ||
      //      track1.hitPattern().pixelLayersWithMeasurement()<2) continue;
      if ( !track1.quality(reco::TrackBase::loose) ) continue;
      if (track1.pt()<1.0) continue;
      if ( (track1_ipValue/track1_ipError) < -999.0 ) continue;
      if (track1.numberOfValidHits()<6) continue;
      if (track1.normalizedChi2()>5) continue;

      for (jTrack=iTrack+1; jTrack<tracks->size(); ++jTrack) {

         const reco::Track& track2 = tracks->at(jTrack);
         double track2_ipError = track2.dxyError();
         double track2_ipValue = std::abs(track2.dxy(bsPoint));

         // both tracks must be within the fiducial volume,
         // have at least one hit in the pixel detector,
         // have at least 6 hits in the tracker and must have opposite charge
         // In addition, the IP significance is required to be more than 0.5
         // and the distance of closest approach between two tracks must be less than 1 cm
         // calculate distance of closest approach

         //if (abs(track2.eta())>2.4) continue;
         //if (abs(track1.dz(bsPoint)-track2.dz(bsPoint))>0.5) continue;
         //if (track2.hitPattern().trackerLayersWithMeasurement()<4 ||
         //      track2.hitPattern().pixelLayersWithMeasurement()<2) continue;
         if ( !track2.quality(reco::TrackBase::loose) ) continue;
         if ( (track2_ipValue/track2_ipError) < -999.0 ) continue;
         if (track2.pt()<1.0) continue;
         if (track2.numberOfValidHits()<6) continue;
         if (track2.normalizedChi2()>5) continue;
         if(track1.charge()==track2.charge()) continue;

         // fit the common vertex with two tracks
         std::vector<TransientTrack> vTTracks;
         vTTracks.push_back(bField->build(track1));
         vTTracks.push_back(bField->build(track2));

         TransientVertex transVtx;
         KalmanVertexFitter kvf(true);

         bool FitOk(true);

         try{
            transVtx = kvf.vertex(vTTracks);
         } catch (...){
            FitOk = false;
         }
         if (!transVtx.isValid()) FitOk = false;
         if (!transVtx.hasRefittedTracks()) FitOk = false;
         if (transVtx.refittedTracks().size()!=vTTracks.size()) FitOk = false;

         if (!FitOk) continue;

         reco::Vertex theSecondaryVertex = transVtx;
         math::XYZPoint phiVertexPoint = math::XYZPoint(transVtx.position().x(),
               transVtx.position().y(),
               transVtx.position().z());

         GlobalPoint phiVertexGlobalPoint(transVtx.position().x(),
               transVtx.position().y(),
               transVtx.position().z());

         ClosestApproachInRPhi dca_tracks;

         dca_tracks.calculate(vTTracks[0].initialFreeState(), vTTracks[1].initialFreeState());

         if (!dca_tracks.status()) continue;
         if (abs(dca_tracks.distance()) > 1.0) continue;

         GlobalPoint cxPt = dca_tracks.crossingPoint();
         if (sqrt(cxPt.x() * cxPt.x() + cxPt.y() * cxPt.y()) > 120. || std::abs(cxPt.z()) > 300.) continue;

         // the tracks should at least point in the same quadrant
         TrajectoryStateClosestToPoint const& track1_TSCP = transVtx.refittedTracks().at(0).trajectoryStateClosestToPoint(cxPt);
         TrajectoryStateClosestToPoint const& track2_TSCP = transVtx.refittedTracks().at(1).trajectoryStateClosestToPoint(cxPt);

         if (!track1_TSCP.isValid() || !track2_TSCP.isValid()) continue;
         if (track1_TSCP.momentum().dot(track2_TSCP.momentum()) < 0) continue;

         // Find the trajectory states closest to the vertex point using refitted tracks 
         //TrajectoryStateClosestToPoint const& track1_TSVP = transVtx.refittedTracks().at(0).trajectoryStateClosestToPoint(phiVertexGlobalPoint);
         //TrajectoryStateClosestToPoint const& track2_TSVP = transVtx.refittedTracks().at(1).trajectoryStateClosestToPoint(phiVertexGlobalPoint);

         if (!track1_TSCP.isValid() || !track2_TSCP.isValid()) continue;

         TLorentzVector track1_p4(track1_TSCP.momentum().x(), 
               track1_TSCP.momentum().y(), 
               track1_TSCP.momentum().z(), 
               sqrt(track1_TSCP.momentum().mag2()+pow(PDG_Vars::Kp_mass(),2)));

         TLorentzVector track2_p4(track2_TSCP.momentum().x(), 
               track2_TSCP.momentum().y(), 
               track2_TSCP.momentum().z(), 
               sqrt(track2_TSCP.momentum().mag2()+pow(PDG_Vars::Kp_mass(),2)));

         TLorentzVector phi_p4 = track1_p4+track2_p4;

         TVector3 p_phi(phi_p4.Px(), phi_p4.Py(), phi_p4.Pz());

         math::Error<3>::type covarianceMatrix;
         for (int i = 0; i <3; i++){
            for (int j = 0; j < 3; j++) {
               covarianceMatrix(i, j) = theSecondaryVertex.covariance(i, j);
               covarianceMatrix(j, i) = theSecondaryVertex.covariance(i, j);
            }
         }

         math::XYZPoint svPoint(transVtx.position().x(), transVtx.position().y(), transVtx.position().z());
         TVector3 svVec(svPoint.x(), svPoint.y(), svPoint.z());

         SVector3 distVecXY(svPoint.x()-bsPoint.x(), svPoint.y()-bsPoint.y(), 0.);
         reco::Vertex phiVertex = Vertex(phiVertexPoint, covarianceMatrix, transVtx.totalChiSquared(), transVtx.degreesOfFreedom(), 2);

         // calculate vertex distances
         Measurement1D BSSV_distanceXY, PVSV_distanceXY, PVSV_distanceXYZ;

         BSSV_distanceXY = vertexTool2D.distance(BSstate, phiVertex);
         PVSV_distanceXY = vertexTool2D.distance(pv, phiVertex);
         PVSV_distanceXYZ = vertexTool3D.distance(pv, phiVertex);

         double distXY = BSSV_distanceXY.value();
         double sigmaXY = BSSV_distanceXY.error();

         double Lxy = distXY;
         double vertexSig = distXY/sigmaXY;
         double totalChi2 = phiVertex.chi2();
         double normChi2 = phiVertex.normalizedChi2();
         double cosTheta = (svVec-bsVec).Dot(p_phi)/((svVec-bsVec).Mag()*p_phi.Mag());

         // Constraints on the vertex fit
         if (totalChi2 < 0.0) continue;
         if (normChi2 > 3.0) continue;
         if (vertexSig < -999.0) continue;
         if (Lxy > 4.0) continue;

         if ( phi_p4.M() < PHI_MASS_WINDOW_LOW || phi_p4.M() > PHI_MASS_WINDOW_HIGH ) continue;

         //Bool_t kMatchedToMuon = false;
         int muon1_index = -1;
         int muon2_index = -1;

         TVector3 transverseP(phi_p4.Px(), phi_p4.Py(), 0.);
         double angleXY = transverseP.Dot(svVec-bsVec)/(distXY*transverseP.Mag());
         if (angleXY < 0.998) continue;

         for (size_t iMu=0; iMu<event_.muons.size(); iMu++){
            TLorentzVector mu_p4;
            mu_p4.SetPtEtaPhiM((double)event_.muons.at(iMu).pt, (double)event_.muons.at(iMu).eta, (double)event_.muons.at(iMu).phi, PDG_Vars::Muon_mass());
            if (mu_p4.DeltaR(track1_p4)<0.01 && abs(track1.pt()-event_.muons.at(iMu).pt)<0.1*track1.pt()) {
               //kMatchedToMuon = true;
               event_.muons.at(iMu).matchedPhi = 1;
               muon1_index = iMu;
            }
            if (mu_p4.DeltaR(track2_p4)<0.01 && abs(track2.pt()-event_.muons.at(iMu).pt)<0.1*track2.pt()) {
               //kMatchedToMuon = true;
               event_.muons.at(iMu).matchedPhi = 1;
               muon2_index = iMu;
            }
         }

         // find angle with respect to the beam spot
         muon_pog::CompositeVertex phi_;

         phi_.mass = phi_p4.M();
         phi_.reco_pdgId = 333;
         phi_.bssv_dxy = distXY;
         phi_.bssv_dxyErr = sigmaXY;
         phi_.bssv_dxySig = vertexSig;
         phi_.pvsv_dxy = PVSV_distanceXY.value();
         phi_.pvsv_dxySig = PVSV_distanceXY.significance();
         phi_.pvsv_dxyError = PVSV_distanceXY.error();
         phi_.pvsv_3d = PVSV_distanceXYZ.value();
         phi_.pvsv_3dError = PVSV_distanceXYZ.error();
         phi_.pvsv_3dSig = PVSV_distanceXYZ.significance();
         phi_.pvsv_normChi2 = normChi2;
         phi_.pvsv_totalChi2 = totalChi2;
         phi_.pvsv_cosTheta = cosTheta;
         phi_.muon1_index = muon1_index;
         phi_.muon2_index = muon2_index;

         event_.phis.push_back(phi_);
         /*
         // printout for debugging
         cout<<"=========================="<<endl;
         cout<<"dxy (vertexTool2D) = "<<Lxy<<endl;
         cout<<"dxy V0Producer = "<<distXY<<endl;
         cout<<"dxyError (vertexTool2D) = "<<vertexTool2D.distance(BSstate, phiVertex).error()<<endl;
         cout<<"dxyError V0Producer = "<<sigmaXY<<endl;
         cout<<"dxySig (vertexTool2D) = "<<vertexSig<<endl;
         cout<<"dxySig (V0Producer) = "<<distXY/sigmaXY<<endl;
         cout<<"=========================="<<endl;
         */
      }
   }
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(MuonPogTreeProducer);
