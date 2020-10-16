// -*- C++ -*-
//
// Package:    TrackCollectionProducer
// Documentation:
// For MiniAOD analysis produce tracks from PF candidates and lost tracks


// system include files
#include <memory>
#include <algorithm>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "MuonPOGtreeProducer/Tools/interface/TrackCollectionProducer.h"
#include "TFile.h"
#include "TH1.h"
#include <vector>
#include <iostream>

using namespace edm;
using namespace std;
using namespace reco;

TrackCollectionProducer::TrackCollectionProducer(const edm::ParameterSet& iConfig)
{    
    edm::InputTag tag;

    if (DEBUG) cout<<"Initializing PFCandidate Tag"<<endl;
    tag = iConfig.getUntrackedParameter<edm::InputTag>("PFCandidateTag", edm::InputTag("packedPFCandidates::RECO") ); 
    pfCandidateToken_ = consumes<std::vector<pat::PackedCandidate> >(tag);
    
    if (DEBUG) cout<<"Initializing LostTrack Tag"<<endl;
    tag = iConfig.getUntrackedParameter<edm::InputTag>("LostTrackTag", edm::InputTag("lostTracks::RECO") );
    lostTrackToken_ = consumes<std::vector<pat::PackedCandidate> >(tag);

    produces<std::vector<reco::Track>>("pfTracks");
}


TrackCollectionProducer::~TrackCollectionProducer()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called on each new Event  ------------
void
TrackCollectionProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    
   
   edm::Handle<std::vector<pat::PackedCandidate> > PFCands;
   edm::Handle<std::vector<pat::PackedCandidate> > lostTracks;

   tracks = std::unique_ptr<std::vector<reco::Track> >(new std::vector<reco::Track>);

   unsigned int nTracks(0);

   if (!pfCandidateToken_.isUninitialized()){
      if (iEvent.getByToken(pfCandidateToken_, PFCands)){
        if (PFCands.isValid()) nTracks += fillPFTracks(PFCands); 
      }
      else edm::LogError("") << "[TrackCollectionProducer]: ParticleFlow Candidates not found!!!";
      
   }

   if (!lostTrackToken_.isUninitialized()){
      if (iEvent.getByToken(lostTrackToken_, lostTracks)){
         //if (lostTracks.isValid()) nTracks += fillLostTracks(lostTracks);
      }
      else edm::LogError("") << "[TrackCollectionProducer]: Lost tracks not found!!!";
   }


   if (DEBUG) cout<<" ================== Number of tracks : "<<nTracks<<" ================"<<endl;
   iEvent.put(std::move(tracks), "pfTracks");
}

unsigned int
TrackCollectionProducer::fillPFTracks(const edm::Handle<std::vector<pat::PackedCandidate> >& PFCands){

   unsigned int iTrack(0);

   vector<pat::PackedCandidate>::const_iterator candIt = PFCands->begin();
   vector<pat::PackedCandidate>::const_iterator candEnd = PFCands->end();

   for (; candIt != candEnd; ++candIt, ++iTrack) {
      
      const pat::PackedCandidate& cand = *candIt; 
      if (cand.bestTrack()!=nullptr){   
         auto track = *cand.bestTrack();
         tracks.get()->push_back(track);
      }
   }

   return iTrack;
}

unsigned int
TrackCollectionProducer::fillLostTracks(const edm::Handle<std::vector<pat::PackedCandidate> >& lostTracks){

   unsigned int iTrack(0);

   vector<pat::PackedCandidate>::const_iterator candIt = lostTracks->begin();
   vector<pat::PackedCandidate>::const_iterator candEnd = lostTracks->end();

   for (; candIt != candEnd; ++candIt, ++iTrack) {
      
      const pat::PackedCandidate& cand = *candIt; 
      if (cand.bestTrack()!=nullptr){   
         auto track = *cand.bestTrack();
         tracks.get()->push_back(track);
      }
   }
   
   return iTrack;
}

//define this as a plug-in
DEFINE_FWK_MODULE(TrackCollectionProducer);
