// -*- C++ -*-
//
// Package:    GGAMuMuAnalysis/GenAnalysis
// Class:      GenAnalyzer
// 
/**\class GenAnalyzer GenAnalyzer.cc GGAMuMuAnalysis/GenAnalysis/plugins/GenAnalyzer.cc

 Description: gen-level studies

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rachel Yohay
//         Created:  Wed, 04 Nov 2015 16:50:51 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "Tools/Common/interface/Common.h"

#include "TH2F.h"
#include "TFile.h"

//
// class declaration
//

class GenAnalyzer : public edm::EDAnalyzer {
   public:
      explicit GenAnalyzer(const edm::ParameterSet&);
      ~GenAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

  //delete memory
  void reset(const bool);

  //recursively plot the PDG ID of all daughters with a given status
  void plotPDGID(const reco::GenParticleRef&, const int, TH1F*);

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------

  std::string outFileName_;
  TFile* out_;

  edm::InputTag genMatchedMuonTag_;

  edm::InputTag genParticleTag_;

  TH2F* mu2PTVsMu1PT_;

  TH1F* status1TauDaughterPDGIDs_;
  TH1F* status2TauDaughterPDGIDs_;
  TH1F* statusNon1Or2TauDaughterPDGIDs_;
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
GenAnalyzer::GenAnalyzer(const edm::ParameterSet& iConfig) :
  outFileName_(iConfig.getParameter<std::string>("outFileName")),
  // genMatchedMuonTag_(iConfig.getParameter<edm::InputTag>("genMatchedMuonTag")),
  genParticleTag_(iConfig.getParameter<edm::InputTag>("genParticleTag"))
{
   //now do what ever initialization is needed

  reset(false);
}


GenAnalyzer::~GenAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

  reset(true);
}


//
// member functions
//

// ------------ method called for each event  ------------
void
GenAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // edm::Handle<reco::MuonRefVector> pGenMatchedMuons;
  // iEvent.getByLabel(genMatchedMuonTag_, pGenMatchedMuons);

  edm::Handle<reco::GenParticleCollection> pGenParticles;
  iEvent.getByLabel(genParticleTag_, pGenParticles);

  // //create a collection of gen-matched muon pointers sorted in ascending order of pT
  // std::vector<reco::Muon*> genMatchedMuonPtrs;
  // if (pGenMatchedMuons.isValid()) {
  //   for (reco::MuonRefVector::const_iterator iGenMatchedMuon = pGenMatchedMuons->begin(); 
  // 	 iGenMatchedMuon != pGenMatchedMuons->end(); ++iGenMatchedMuon) {
  //     genMatchedMuonPtrs.push_back(const_cast<reco::Muon*>((*iGenMatchedMuon).get()));
  //   }
  // }
  // Common::sortByPT(genMatchedMuonPtrs);

  // //define the highest and second-highest pT gen-matched muons
  // const unsigned int nGenMatchedMuonPtrs = genMatchedMuonPtrs.size();
  // const reco::Muon* genMatchedMu1 = 
  //   const_cast<const reco::Muon*>(genMatchedMuonPtrs[nGenMatchedMuonPtrs - 2]);
  // const reco::Muon* genMatchedMu2 = 
  //   const_cast<const reco::Muon*>(genMatchedMuonPtrs[nGenMatchedMuonPtrs - 1]);

  // mu2PTVsMu1PT_->Fill(genMatchedMu2->pt(), genMatchedMu1->pt());

  //plot PDG IDs of all status 1, 2, and non-1-or-2 tau daughters
  for (reco::GenParticleCollection::const_iterator iGenParticle = pGenParticles->begin(); 
       iGenParticle != pGenParticles->end(); ++iGenParticle) {
    if ((fabs(iGenParticle->pdgId()) == GenTauDecayID::TAUPDGID) && (iGenParticle->status() == 2)) {
      plotPDGID(reco::GenParticleRef(pGenParticles, iGenParticle - pGenParticles->begin()), 1, 
		status1TauDaughterPDGIDs_);
      plotPDGID(reco::GenParticleRef(pGenParticles, iGenParticle - pGenParticles->begin()), 2, 
		status2TauDaughterPDGIDs_);
      plotPDGID(reco::GenParticleRef(pGenParticles, iGenParticle - pGenParticles->begin()), 0, 
		statusNon1Or2TauDaughterPDGIDs_);
    }
  }
}


// ------------ method called once each job just before starting event loop  ------------
void 
GenAnalyzer::beginJob()
{
  out_ = new TFile(outFileName_.c_str(), "RECREATE");

  mu2PTVsMu1PT_ = new TH2F("mu2PTVsMu1PT", ";p_{T}^{1} (GeV);p_{T}^{2} (GeV)", 
			   25, 0.0, 50.0, 25, 0.0, 50.0);

  status1TauDaughterPDGIDs_ = 
    new TH1F("status1TauDaughterPDGIDs", ";PDG ID;", 101000, 1.0, 101001.0);
  status2TauDaughterPDGIDs_ = 
    new TH1F("status2TauDaughterPDGIDs", ";PDG ID;", 101000, 1.0, 101001.0);
  statusNon1Or2TauDaughterPDGIDs_ = 
    new TH1F("statusNon1Or2TauDaughterPDGIDs", ";PDG ID;", 101000, 1.0, 101001.0);

  mu2PTVsMu1PT_->Sumw2();
  status1TauDaughterPDGIDs_->Sumw2();
  status2TauDaughterPDGIDs_->Sumw2();
  statusNon1Or2TauDaughterPDGIDs_->Sumw2();
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GenAnalyzer::endJob() 
{
  out_->cd();

  TCanvas mu2PTVsMu1PTCanvas("mu2PTVsMu1PTCanvas", "", 600, 600);

  TCanvas status1TauDaughterPDGIDsCanvas("status1TauDaughterPDGIDsCanvas", "", 600, 600);
  TCanvas status2TauDaughterPDGIDsCanvas("status2TauDaughterPDGIDsCanvas", "", 600, 600);
  TCanvas 
    statusNon1Or2TauDaughterPDGIDsCanvas("statusNon1Or2TauDaughterPDGIDsCanvas", "", 600, 600);

  Common::draw2DHistograms(mu2PTVsMu1PTCanvas, mu2PTVsMu1PT_);

  Common::draw1DHistograms(status1TauDaughterPDGIDsCanvas, status1TauDaughterPDGIDs_);
  Common::draw1DHistograms(status2TauDaughterPDGIDsCanvas, status2TauDaughterPDGIDs_);
  Common::draw1DHistograms(statusNon1Or2TauDaughterPDGIDsCanvas, statusNon1Or2TauDaughterPDGIDs_);

  for (Int_t iBin = 1; iBin <= (status1TauDaughterPDGIDs_->GetNbinsX() + 1); ++iBin) {
    const unsigned int nParticles = status1TauDaughterPDGIDs_->GetBinContent(iBin);
    if (nParticles > 0) {
      std::cout << "Number of status 1 particles with PDG ID ";
      std::cout << status1TauDaughterPDGIDs_->GetBinLowEdge(iBin) << ": " << nParticles;
      std::cout << std::endl;
    }
  }
  for (Int_t iBin = 1; iBin <= (status2TauDaughterPDGIDs_->GetNbinsX() + 1); ++iBin) {
    const unsigned int nParticles = status2TauDaughterPDGIDs_->GetBinContent(iBin);
    if (nParticles > 0) {
      std::cout << "Number of status 2 particles with PDG ID ";
      std::cout << status2TauDaughterPDGIDs_->GetBinLowEdge(iBin) << ": " << nParticles;
      std::cout << std::endl;
    }
  }
  for (Int_t iBin = 1; iBin <= (statusNon1Or2TauDaughterPDGIDs_->GetNbinsX() + 1); ++iBin) {
    const unsigned int nParticles = statusNon1Or2TauDaughterPDGIDs_->GetBinContent(iBin);
    if (nParticles > 0) {
      std::cout << "Number of status non-1-or-2 particles with PDG ID ";
      std::cout << statusNon1Or2TauDaughterPDGIDs_->GetBinLowEdge(iBin) << ": " << nParticles;
      std::cout << std::endl;
    }
  }

  out_->cd();

  mu2PTVsMu1PTCanvas.Write();

  status1TauDaughterPDGIDsCanvas.Write();
  status2TauDaughterPDGIDsCanvas.Write();
  statusNon1Or2TauDaughterPDGIDsCanvas.Write();

  out_->Write();
  out_->Close();
}

void GenAnalyzer::reset(const bool doDelete)
{
  if (doDelete && (out_ != NULL)) delete out_;
  out_ = NULL;

  if (doDelete && (mu2PTVsMu1PT_ != NULL)) delete mu2PTVsMu1PT_;
  mu2PTVsMu1PT_ = NULL;

  if (doDelete && (status1TauDaughterPDGIDs_ != NULL)) delete status1TauDaughterPDGIDs_;
  status1TauDaughterPDGIDs_ = NULL;
  if (doDelete && (status2TauDaughterPDGIDs_ != NULL)) delete status2TauDaughterPDGIDs_;
  status2TauDaughterPDGIDs_ = NULL;
  if (doDelete && (statusNon1Or2TauDaughterPDGIDs_ != NULL)) delete statusNon1Or2TauDaughterPDGIDs_;
  statusNon1Or2TauDaughterPDGIDs_ = NULL;
}

// ------------ method called when starting to processes a run  ------------
/*
void 
GenAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
GenAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
GenAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
GenAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

 void GenAnalyzer::plotPDGID(const reco::GenParticleRef& momRef, const int status, TH1F* hist)
{
  for (unsigned int iDaughter = 0; iDaughter < momRef->numberOfDaughters(); 
       ++iDaughter) {
    reco::GenParticleRef kidRef = momRef->daughterRef(iDaughter);
    const unsigned int absDaughterPDGID = fabs(kidRef->pdgId());
    const int daughterStatus = kidRef->status();
    if (((status != 0) && (daughterStatus == status)) || 
	((status == 0) && (daughterStatus != 1) && (daughterStatus != 2))) {
      hist->Fill(absDaughterPDGID);
    }
    if (kidRef->status() != 1) plotPDGID(kidRef, status, hist);
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenAnalyzer);
