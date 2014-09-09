// -*- C++ -*-
//
// Package:    ntuple_maker
// Class:      ntuple_maker
// 
/**\class ntuple_maker ntuple_maker.cc ntuple_maker/ntuple_maker/plugins/ntuple_maker.cc

 Description: makes smaller ntuples from miniaod

 Implementation:
 writes a flat ttree
*/
//
// Original Author:  Andrew Michael Levin
//         Created:  Fri, 25 Jul 2014 16:58:43 GMT
// $Id$
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

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TFile.h"
#include "TTree.h"
#include "Math/LorentzVector.h"

//
// class declaration
//



class ntuple_maker : public edm::EDAnalyzer {
public:
  typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;
  
      explicit ntuple_maker(const edm::ParameterSet&);
      ~ntuple_maker();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  enum Cuts {
    Lep1FullSelection  = 1UL<<2, 
    Lep2FullSelection  = 1UL<<9 
  };


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------

  // ----------member data ---------------------------
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<pat::MuonCollection> muonToken_;
  edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
  edm::EDGetTokenT<pat::TauCollection> tauToken_;
  edm::EDGetTokenT<pat::PhotonCollection> photonToken_;
  edm::EDGetTokenT<pat::JetCollection> jetToken_;
  edm::EDGetTokenT<pat::JetCollection> fatjetToken_;
  edm::EDGetTokenT<pat::METCollection> metToken_;

  UInt_t cuts;
  UInt_t event;
  UInt_t run;
  UInt_t lumi;
  UInt_t nvtx;
  TTree * tree;
  Float_t jetpt;
  Float_t jet1pujetid;
  Float_t jet1btag;
  Float_t jet1btagincl;
  Float_t jet2pujetid;
  Float_t jet2btag;
  Float_t jet2btagincl;
  Float_t metpt;
  Float_t metphi;
  Float_t metsumet;
  Float_t metgenmetpt;
  Float_t metptshiftup;
  Float_t metptshiftdown;
  LorentzVector jet1;
  LorentzVector jet2;
  LorentzVector lep1;
  LorentzVector lep2;
  Int_t lep1id;
  Int_t lep2id;
  Int_t lep1q;
  Int_t lep2q;

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
ntuple_maker::ntuple_maker(const edm::ParameterSet& iConfig):
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
  tauToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
  photonToken_(consumes<pat::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
  jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
  fatjetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("fatjets"))),
  metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets")))
{
  //now do what ever initialization is needed

}


ntuple_maker::~ntuple_maker()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ntuple_maker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  cuts = 0;

   using namespace edm;

   run=iEvent.eventAuxiliary().run(); 
   lumi=iEvent.eventAuxiliary().luminosityBlock();
   event=iEvent.eventAuxiliary().event();

   edm::Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vtxToken_, vertices);
   if (vertices->empty()) return; // skip the event if no PV found
   const reco::Vertex &PV = vertices->front();

   nvtx=vertices->size();

   edm::Handle<pat::MuonCollection> muons;
   iEvent.getByToken(muonToken_, muons);

   Bool_t found_muon = kFALSE;

   for (const pat::Muon &mu : *muons) {
     if (mu.pt() < 10 || !mu.isLooseMuon() || mu.chargedHadronIso()/mu.pt() > 0.3)
       continue;

     Bool_t lep1passfullselection = kTRUE;

     Float_t relative_isolation = (mu.pfIsolationR04().sumChargedHadronPt+ std::max(0.0,mu.pfIsolationR04().sumNeutralHadronEt + mu.pfIsolationR04().sumPhotonEt - 0.5 * mu.pfIsolationR04().sumChargedHadronPt))/mu.pt();
     if (!mu.isTightMuon(PV) || relative_isolation > 0.12) 
       lep1passfullselection = kFALSE;

     if (lep1passfullselection)
       cuts = cuts | Lep1FullSelection;

     lep1 = mu.p4();
     lep1q = mu.charge();
     lep1id = mu.pdgId();
     found_muon=kTRUE;
     break;

     //printf("muon with pt %4.1f, dz(PV) %+5.3f, POG loose id %d, tight id %d\n",
     //	    mu.pt(), mu.muonBestTrack()->dz(PV.position()), mu.isLooseMuon(), mu.isTightMuon(PV));
   }

   if(!found_muon)
     return;

   Bool_t found_electron = kFALSE;

   edm::Handle<pat::ElectronCollection> electrons;
   iEvent.getByToken(electronToken_, electrons);

   for (const pat::Electron &el : *electrons) {
     if (el.pt() < 10 || el.chargedHadronIso()/el.pt() > 0.3) 
       continue;
     //if (el.chargedHadronIso()/el.pt() > 0.3)
     //  continue;

     reco::GsfElectron::PflowIsolationVariables pfIso = el.pfIsolationVariables();

     Float_t absiso = pfIso.sumChargedHadronPt + std::max(0.0 , pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - 0.5 * pfIso.sumPUPt );
     Float_t relIsoWithDBeta = absiso/el.pt();
     Float_t ooEmooP = 0;

     if( el.ecalEnergy() == 0 ){
       std::cout << "Electron energy is zero!" << std::endl;
       ooEmooP = 1e30;
     }
     else if (!std::isfinite(el.ecalEnergy())){
       std::cout << "Electron energy is not finite!" << std::endl;
       ooEmooP = 1e30;
     }
     else
       ooEmooP = abs(1.0/el.ecalEnergy() - el.eSuperClusterOverP()/el.ecalEnergy() );
       

     Bool_t lep2passfullselection = kTRUE;

     //medium working point from here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
     if(el.superCluster()->eta() < 2.5 && el.superCluster()->eta() > 1.479 ){
       if (abs(el.deltaEtaSuperClusterTrackAtVtx()) > 0.0108)
	 lep2passfullselection = kFALSE;
       if ( abs(el.deltaPhiSuperClusterTrackAtVtx()) > 0.0455)
	 lep2passfullselection = kFALSE;
       if (el.full5x5_sigmaIetaIeta() > 0.0318)
	 lep2passfullselection = kFALSE;
       if ( el.hcalOverEcal() > 0.097)
	 lep2passfullselection = kFALSE;
       if( abs((-1) * el.gsfTrack()->dxy(PV.position())) > 0.0845)
	 lep2passfullselection = kFALSE;
       if(  abs(el.gsfTrack()->dz( PV.position() )) > 0.7523)
	 lep2passfullselection = kFALSE;
       if(abs(ooEmooP) > 0.1201)
	 lep2passfullselection = kFALSE;
       if(relIsoWithDBeta > 0.254)
	 lep2passfullselection = kFALSE;
       if(!el.passConversionVeto())
	 lep2passfullselection = kFALSE;
       if(el.gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() > 1)
	 lep2passfullselection = kFALSE;
     } else if (el.superCluster()->eta() < 1.479) {
       if (abs(el.deltaEtaSuperClusterTrackAtVtx()) > 0.0106)
	 lep2passfullselection = kFALSE;
       if ( abs(el.deltaPhiSuperClusterTrackAtVtx()) > 0.0323)
	 lep2passfullselection = kFALSE;
       if (el.full5x5_sigmaIetaIeta() > 0.0107)
	 lep2passfullselection = kFALSE;
       if ( el.hcalOverEcal() > 0.067)
	 lep2passfullselection = kFALSE;
       if( abs((-1) * el.gsfTrack()->dxy(PV.position())) > 0.0131)
	 lep2passfullselection = kFALSE;
       if(  abs(el.gsfTrack()->dz( PV.position() )) > 0.22310)
	 lep2passfullselection = kFALSE;
       if(abs(ooEmooP) > 0.1043)
	 lep2passfullselection = kFALSE;
       if(relIsoWithDBeta > 0.2179)
	 lep2passfullselection = kFALSE;
       if(!el.passConversionVeto())
	 lep2passfullselection = kFALSE;
       if(el.gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() > 1)
	 lep2passfullselection = kFALSE;
     } else
       continue;
     
     if (lep2passfullselection)
       cuts = cuts | Lep2FullSelection;

     lep2= el.p4();
     lep2q = el.charge();
     lep2id = el.pdgId();
     found_electron=true;
     break;

     //printf("elec with pt %4.1f, supercluster eta %+5.3f, sigmaIetaIeta %.3f (%.3f with full5x5 shower shapes), lost hits %d, pass conv veto %d\n",
     //	    el.pt(), el.superCluster()->eta(), el.sigmaIetaIeta(), el.full5x5_sigmaIetaIeta(), el.gsfTrack()->trackerExpectedHitsInner().numberOfLostHits(), el.passConversionVeto());
   }

  if(!found_electron)
     return;

   edm::Handle<pat::PhotonCollection> photons;
   iEvent.getByToken(photonToken_, photons);
   for (const pat::Photon &pho : *photons) {
     if (pho.pt() < 20 or pho.chargedHadronIso()/pho.pt() > 0.3) continue;
     printf("phot with pt %4.1f, supercluster eta %+5.3f, sigmaIetaIeta %.3f (%.3f with full5x5 shower shapes)\n",
	    pho.pt(), pho.superCluster()->eta(), pho.sigmaIetaIeta(), pho.full5x5_sigmaIetaIeta());
   }


   edm::Handle<pat::TauCollection> taus;
   iEvent.getByToken(tauToken_, taus);
   for (const pat::Tau &tau : *taus) {
     if (tau.pt() < 20) continue;
     printf("tau  with pt %4.1f, dxy signif %.1f, ID(byMediumCombinedIsolationDeltaBetaCorr3Hits) %.1f, lead candidate pt %.1f, pdgId %d \n",
	    tau.pt(), tau.dxy_Sig(), tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits"), tau.leadCand()->pt(), tau.leadCand()->pdgId());
   }


   edm::Handle<pat::JetCollection> jets;
   iEvent.getByToken(jetToken_, jets);

   if (jets->size() < 2) 
     return;
   else if ((*jets)[1].pt() < 20)
     return;
   jet1=(*jets)[0].p4();
   jet2=(*jets)[1].p4();   
   jet1btag = std::max(0.f,(*jets)[0].bDiscriminator("combinedSecondaryVertexBJetTags"));
   jet1btagincl = std::max(0.f,(*jets)[0].bDiscriminator("combinedInclusiveSecondaryVertexBJetTags"));
   jet2btag = std::max(0.f,(*jets)[1].bDiscriminator("combinedSecondaryVertexBJetTags"));
   jet2btagincl = std::max(0.f,(*jets)[1].bDiscriminator("combinedInclusiveSecondaryVertexBJetTags"));
   jet1pujetid = (*jets)[0].userFloat("pileupJetId:fullDiscriminant");
   jet2pujetid = (*jets)[1].userFloat("pileupJetId:fullDiscriminant");


   for (const pat::Jet &j : *jets) {
     if (j.pt() < 20) continue;

     if( std::max(0.f,j.bDiscriminator("combinedSecondaryVertexBJetTags")) > 0.5)
       return;

   }


   edm::Handle<pat::JetCollection> fatjets;
   iEvent.getByToken(fatjetToken_, fatjets);
   for (const pat::Jet &j : *fatjets) {
     printf("AK8j with pt %5.1f (raw pt %5.1f), eta %+4.2f, mass %5.1f ungroomed, %5.1f pruned, %5.1f trimmed, %5.1f filtered. CMS TopTagger %.1f\n",
            j.pt(), j.pt()*j.jecFactor("Uncorrected"), j.eta(), j.mass(), j.userFloat("ak8PFJetsCHSPrunedLinks"), j.userFloat("ak8PFJetsCHSTrimmedLinks"), j.userFloat("ak8PFJetsCHSFilteredLinks"), j.userFloat("cmsTopTagPFJetsCHSLinksAK8"));
   }

   
 
   edm::Handle<pat::METCollection> mets;
   iEvent.getByToken(metToken_, mets);
   const pat::MET &met = mets->front();
   printf("MET: pt %5.1f, phi %+4.2f, sumEt (%.1f). genMET %.1f. MET with JES up/down: %.1f/%.1f\n",
	  met.pt(), met.phi(), met.sumEt(),
	  met.genMET()->pt(),
	  met.shiftedPt(pat::MET::JetEnUp), met.shiftedPt(pat::MET::JetEnDown));

   printf("\n");

   metpt = met.pt();
   metphi = met.phi();
   metsumet = met.sumEt();
   metgenmetpt = met.genMET()->pt();
   metptshiftup = met.shiftedPt(pat::MET::JetEnUp);
   metptshiftdown = met.shiftedPt(pat::MET::JetEnDown);

   tree->Fill();

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
ntuple_maker::beginJob()
{

  edm::Service<TFileService> fs;

  tree = fs->make<TTree>( "events"  , "events");

  tree->Branch("cuts",&cuts);

  tree->Branch("event",&event);
  tree->Branch("lumi",&lumi);
  tree->Branch("run",&run);

  tree->Branch("metpt",&metpt);
  tree->Branch("metphi",&metphi);
  tree->Branch("metsumet",&metsumet);
  tree->Branch("metgenmetpt",&metgenmetpt);
  tree->Branch("metptshiftup",&metptshiftup);
  tree->Branch("metptshiftdown",&metptshiftdown);
  tree->Branch("jet1",&jet1);
  tree->Branch("jet2",&jet2);
  tree->Branch("jet1pujetid",&jet1pujetid);
  tree->Branch("jet2pujetid",&jet2pujetid);
  tree->Branch("jet1btag",&jet1btag);
  tree->Branch("jet2btag",&jet2btag);
  tree->Branch("jet1btagincl",&jet1btagincl);
  tree->Branch("jet2btagincl",&jet2btagincl);
  tree->Branch("lep1",&lep1);
  tree->Branch("lep2",&lep2);
  tree->Branch("nvtx",&nvtx);
  tree->Branch("lep1q",&lep1q);
  tree->Branch("lep1id",&lep1id);
  tree->Branch("lep2q",&lep2q);
  tree->Branch("lep2id",&lep2id);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
ntuple_maker::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------

/*
void 
ntuple_maker::beginRun(edm::Run const&, edm::EventSetup const&)
{



}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
ntuple_maker::endRun(edm::Run const&, edm::EventSetup const&)
{



}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
ntuple_maker::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
ntuple_maker::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ntuple_maker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ntuple_maker);
