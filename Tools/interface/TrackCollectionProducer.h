#ifndef MuonPogTreeProducer_Tools_TrackCollectionProducer_H
#define MuonPogTreeProducer_Tools_TrackCollectionProducer_H
/* \class TrackCollectionProducer
 *
 * \author Bhargav Joshi, UF
 *
 */
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "DataFormats/Candidate/interface/Candidate.h"

class TrackCollectionProducer : public edm::EDProducer {
   
   public:
      explicit TrackCollectionProducer(const edm::ParameterSet &);
       ~TrackCollectionProducer();

   private:
      virtual void produce(edm::Event&, const edm::EventSetup&) override;

      unsigned int fillPFTracks(const edm::Handle<std::vector<pat::PackedCandidate> > &);
      unsigned int fillLostTracks(const edm::Handle<std::vector<pat::PackedCandidate> > &);
  
      // ----------member data ---------------------------
      edm::EDGetTokenT<std::vector<pat::PackedCandidate> > pfCandidateToken_;
      edm::EDGetTokenT<std::vector<pat::PackedCandidate> > lostTrackToken_;
      edm::EDPutTokenT<std::vector<reco::Track> > trackToken_;
      
      std::unique_ptr<std::vector<reco::Track>> tracks;
      
      bool DEBUG = false;
      
};
#endif
