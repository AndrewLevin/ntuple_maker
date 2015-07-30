#ifndef lhe_and_gen_H
#define lhe_and_gen_H

#include <memory>
#include <vector>
using namespace std;

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// TFile Service
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

// MINIAOD INCLUDES
//

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

//VALUE MAP
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

// GEN
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"


#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

//PU
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

// TRIGGER
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

// Try to add something
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"

// TRIGGER
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"

// ROOT
#include "TTree.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

class lhe_and_gen
{
 public:
  lhe_and_gen();
  int analyze(const edm::Event& iEvent, LorentzVector & lep1, LorentzVector & lep2, std::vector<int> & pdf_weight_indices, int qcd_weight_up_index, int qcd_weight_down_index);
  void defineBranches(TTree * tree);

  edm::EDGetTokenT<LHEEventProduct> lheEvtToken_;
  edm::EDGetTokenT<edm::View<reco::GenParticle> > prunedGenToken_;
  edm::EDGetTokenT<edm::View<pat::PackedGenParticle> > packedGenToken_;

  LorentzVector lep1_nearestparton_4mom;
  LorentzVector lep2_nearestparton_4mom;
  Int_t lep1_nearestparton_pdgid;
  Int_t lep2_nearestparton_pdgid;
  Int_t lep1_matching_real_gen_lepton_pdgid;
  Int_t lep2_matching_real_gen_lepton_pdgid;
  Int_t lep1_matching_real_gen_lepton_q;
  Int_t lep2_matching_real_gen_lepton_q;
  Float_t lep1_nearest_gen_electron_dr;
  Float_t lep2_nearest_gen_electron_dr;
  Float_t lep1_nearest_gen_muon_dr;
  Float_t lep2_nearest_gen_muon_dr;

  std::vector<Float_t> * lhe_weights;  
  std::vector<Float_t> * pdf_weights;  
  Float_t lhe_weight_orig;
  Float_t qcd_pdf_weight_orig;
  Float_t qcd_weight_up;
  Float_t qcd_weight_down;

  bool syscalcinfo_;

};

#endif
