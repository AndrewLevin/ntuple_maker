#include "ntuple_maker/ntuple_maker/interface/lhe_and_gen.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

void lhe_and_gen::defineBranches(TTree *tree){

  tree->Branch("lep1_nearestparton_4mom",&lep1_nearestparton_4mom);
  tree->Branch("lep1_nearestparton_pdgid",&lep1_nearestparton_pdgid);
  tree->Branch("lep2_nearestparton_4mom",&lep2_nearestparton_4mom);
  tree->Branch("lep2_nearestparton_pdgid",&lep2_nearestparton_pdgid);
  tree->Branch("lep1_matching_real_gen_lepton_pdgid",&lep1_matching_real_gen_lepton_pdgid);
  tree->Branch("lep2_matching_real_gen_lepton_pdgid",&lep2_matching_real_gen_lepton_pdgid);
  tree->Branch("lep1_matching_real_gen_lepton_q",&lep1_matching_real_gen_lepton_q);
  tree->Branch("lep2_matching_real_gen_lepton_q",&lep2_matching_real_gen_lepton_q);
  tree->Branch("lep1_nearest_gen_electron_dr",&lep1_nearest_gen_electron_dr);
  tree->Branch("lep2_nearest_gen_electron_dr",&lep2_nearest_gen_electron_dr);
  tree->Branch("lep1_nearest_gen_muon_dr",&lep1_nearest_gen_muon_dr);
  tree->Branch("lep2_nearest_gen_muon_dr",&lep2_nearest_gen_muon_dr);

  tree->Branch("lhe_weight_orig",&lhe_weight_orig);
  tree->Branch("lhe_weights",&lhe_weights);



}

reco::GenParticle find_nearest_parton(LorentzVector &lepton_4mom, edm::Handle<edm::View<reco::GenParticle> > pruned)
{

  Float_t minpartondR = std::numeric_limits<Float_t>::max();

  reco::GenParticle nearest_parton;

  for(size_t j=0; j<pruned->size();j++){
    
    if ( abs((*pruned)[j].pdgId()) == 1 || abs((*pruned)[j].pdgId()) == 2 || abs((*pruned)[j].pdgId()) == 3 || abs((*pruned)[j].pdgId()) == 4 || abs((*pruned)[j].pdgId()) == 5 || abs((*pruned)[j].pdgId()) == 21 ){
      
      //status 23 means outgoing: http://home.thep.lu.se/~torbjorn/pythia81html/ParticleProperties.html
      if (((*pruned)[j].status() == 23 || (*pruned)[j].status() == 71) && reco::deltaR((*pruned)[j].p4(),lepton_4mom) < minpartondR){
	minpartondR = reco::deltaR((*pruned)[j].p4(),lepton_4mom);
	nearest_parton = (*pruned)[j];
      }
    }
  }

  return nearest_parton;

};

Float_t find_nearest_gen_lepton(LorentzVector &lepton_4mom, edm::Handle<edm::View<pat::PackedGenParticle> > packed, std::string flavor)
{

  assert(flavor == "muon" || flavor == "electron");

  Int_t abspdgid;

  if (flavor == "muon")
    abspdgid = 13;
  else
    abspdgid = 11;

  Float_t minleptondR = std::numeric_limits<Float_t>::max();

  for(size_t j=0; j<packed->size();j++){
    
    if ( abs((*packed)[j].pdgId()) == abspdgid ){
      
      if ((*packed)[j].status() == 1 && reco::deltaR((*packed)[j].p4(),lepton_4mom) < minleptondR){
	minleptondR = reco::deltaR((*packed)[j].p4(),lepton_4mom);
      }
    }
  }

  return minleptondR;

};

lhe_and_gen::lhe_and_gen()
{

}

int lhe_and_gen::analyze(const edm::Event& iEvent, LorentzVector & lep1, LorentzVector & lep2){


  edm::Handle<LHEEventProduct> hLheEvt;
  iEvent.getByToken(lheEvtToken_,hLheEvt);

  lhe_weights = new std::vector<Float_t>();

  lhe_weight_orig = hLheEvt->originalXWGTUP();

  for(unsigned int i = 0; i < hLheEvt->weights().size(); i++){
    lhe_weights->push_back(hLheEvt->weights()[i].wgt);
  }

  edm::Handle<edm::View<pat::PackedGenParticle> > packed;
  iEvent.getByToken(packedGenToken_,packed);

  lep1_nearest_gen_electron_dr = std::numeric_limits<Float_t>::max();
  lep2_nearest_gen_electron_dr = std::numeric_limits<Float_t>::max();
  lep1_nearest_gen_muon_dr = std::numeric_limits<Float_t>::max();
  lep2_nearest_gen_muon_dr = std::numeric_limits<Float_t>::max();

  lep1_nearest_gen_muon_dr = find_nearest_gen_lepton(lep1, packed, "muon");
  lep1_nearest_gen_electron_dr = find_nearest_gen_lepton(lep1, packed,"electron");
  lep2_nearest_gen_muon_dr = find_nearest_gen_lepton(lep1, packed, "muon");
  lep2_nearest_gen_electron_dr = find_nearest_gen_lepton(lep2, packed, "electron");

  edm::Handle<edm::View<reco::GenParticle> > pruned;
  iEvent.getByToken(prunedGenToken_,pruned);
  
  reco::GenParticle lep1_nearest_parton;
  reco::GenParticle lep2_nearest_parton;

  lep1_nearest_parton=find_nearest_parton(lep1,pruned);
  lep2_nearest_parton=find_nearest_parton(lep2,pruned);

  lep1_matching_real_gen_lepton_pdgid=0;
  lep1_matching_real_gen_lepton_q=0;
  lep2_matching_real_gen_lepton_pdgid=0;
  lep2_matching_real_gen_lepton_q=0;
  
  for(size_t j=0; j<pruned->size();j++){
    if ( (abs((*pruned)[j].pdgId()) == 11 || abs((*pruned)[j].pdgId()) == 13) && reco::deltaR((*pruned)[j].p4(),lep1) < 0.3 ){
      
      if (! (*pruned)[j].mother())
	continue;
      
      const reco::GenParticle * ancestor = dynamic_cast<const reco::GenParticle *>((*pruned)[j].mother());
      
      while(abs(ancestor->pdgId()) == 15){
	ancestor = dynamic_cast<const reco::GenParticle *>(ancestor->mother());
      }
      
      if (abs(ancestor->pdgId()) != 24)
	continue;
      
      lep1_matching_real_gen_lepton_pdgid=(*pruned)[j].pdgId();
      lep1_matching_real_gen_lepton_q=(*pruned)[j].charge();
    }
  }
  
  for(size_t j=0; j<pruned->size();j++){
    if ( (abs((*pruned)[j].pdgId()) == 11 || abs((*pruned)[j].pdgId()) == 13) && reco::deltaR((*pruned)[j].p4(),lep2) < 0.3 ){
      
      if (! (*pruned)[j].mother())
	continue;
      
      const reco::GenParticle * ancestor = dynamic_cast<const reco::GenParticle *>((*pruned)[j].mother());
      
      while(abs(ancestor->pdgId()) == 15){
	ancestor = dynamic_cast<const reco::GenParticle *>(ancestor->mother());
      }
      
      if (abs(ancestor->pdgId()) != 24)
	continue;
      
      lep2_matching_real_gen_lepton_pdgid=(*pruned)[j].pdgId();
      lep2_matching_real_gen_lepton_q=(*pruned)[j].charge();
    }
  }
  
  lep1_nearestparton_4mom = lep1_nearest_parton.p4();
  lep2_nearestparton_4mom = lep2_nearest_parton.p4();
  lep1_nearestparton_pdgid = lep1_nearest_parton.pdgId();
  lep2_nearestparton_pdgid = lep2_nearest_parton.pdgId();

  return 0;

}

