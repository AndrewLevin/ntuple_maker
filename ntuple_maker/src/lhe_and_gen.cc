#include "ntuple_maker/ntuple_maker/interface/lhe_and_gen.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

void lhe_and_gen::defineBranches(TTree *tree){

  initrwgt_header_tree_->Branch("initrwgt_header_line",&initrwgt_header_line_);

  slha_header_tree_->Branch("slha_header_line",&slha_header_line_);

  if (lhe_lepton_info_){
    tree->Branch("lhelep1pdgid",&lhelep1pdgid);
    tree->Branch("lhelep2pdgid",&lhelep2pdgid);
  }

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

  tree->Branch("gen_weight",&gen_weight);

  tree->Branch("lhe_weight_orig",&lhe_weight_orig);

  tree->Branch("qcd_pdf_weight_orig",&qcd_pdf_weight_orig);

  tree->Branch("mgreweight_weights",&mgreweight_weights);

  tree->Branch("pdf_weights",&pdf_weights);

  tree->Branch("qcd_weight_mur1muf2",&qcd_weight_mur1muf2);
  tree->Branch("qcd_weight_mur1muf0p5",&qcd_weight_mur1muf0p5);
  tree->Branch("qcd_weight_mur2muf1",&qcd_weight_mur2muf1);
  tree->Branch("qcd_weight_mur2muf2",&qcd_weight_mur2muf2);
  tree->Branch("qcd_weight_mur2muf0p5",&qcd_weight_mur2muf0p5);
  tree->Branch("qcd_weight_mur0p5muf1",&qcd_weight_mur0p5muf1);
  tree->Branch("qcd_weight_mur0p5muf2",&qcd_weight_mur0p5muf2);
  tree->Branch("qcd_weight_mur0p5muf0p5",&qcd_weight_mur0p5muf0p5);


  //  tree->Branch("qcd_weight_down",&qcd_weight_down);

  //  tree->Branch("qcd_weight_up",&qcd_weight_up);

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

lhe_and_gen::lhe_and_gen():
  qcd_weight_mur1muf2_index(std::numeric_limits<int>::max()),
  qcd_weight_mur1muf0p5_index(std::numeric_limits<int>::max()),
  qcd_weight_mur2muf1_index(std::numeric_limits<int>::max()),
  qcd_weight_mur2muf2_index(std::numeric_limits<int>::max()),
  qcd_weight_mur2muf0p5_index(std::numeric_limits<int>::max()),
  qcd_weight_mur0p5muf1_index(std::numeric_limits<int>::max()),
  qcd_weight_mur0p5muf2_index(std::numeric_limits<int>::max()),
  qcd_weight_mur0p5muf0p5_index(std::numeric_limits<int>::max())
{

}

int lhe_and_gen::analyze(const edm::Event& iEvent, LorentzVector & lep1, LorentzVector & lep2){

  mgreweight_weights = new std::vector<Float_t>();
  pdf_weights = new std::vector<Float_t>();

  if (! isMC_)
    return 0;

  if(mgreweightinfo_ || syscalcinfo_ || lhe_lepton_info_){

  edm::Handle<LHEEventProduct> hLheEvt;
  iEvent.getByToken(lheEvtToken_,hLheEvt);

  edm::Handle<GenEventInfoProduct> hGenEvt;
  iEvent.getByToken(genEvtToken_,hGenEvt);

  //assert(hGenEvt->weight() == hLheEvt->originalXWGTUP());

  gen_weight = hGenEvt->weight();

  lhe_weight_orig = hLheEvt->originalXWGTUP();

  std::vector<unsigned int> lhe_lepton_indices;

  if (lhe_lepton_info_){
    for(unsigned int i = 0; i < hLheEvt->hepeup().IDUP.size(); i++){

      if (abs(hLheEvt->hepeup().IDUP[i]) == 11 || abs(hLheEvt->hepeup().IDUP[i]) == 13 || abs(hLheEvt->hepeup().IDUP[i]) == 15){
	
	lhe_lepton_indices.push_back(i);
      }
    }
    
    assert(lhe_lepton_indices.size() == 2);

    lhelep1pdgid = hLheEvt->hepeup().IDUP[lhe_lepton_indices[0]];
    lhelep2pdgid = hLheEvt->hepeup().IDUP[lhe_lepton_indices[1]];



  }

  if (syscalcinfo_){

    assert(qcd_weight_mur1muf2_index != std::numeric_limits<int>::max());
    assert(qcd_weight_mur1muf0p5_index != std::numeric_limits<int>::max());
    assert(qcd_weight_mur2muf1_index != std::numeric_limits<int>::max());
    assert(qcd_weight_mur2muf2_index != std::numeric_limits<int>::max());
    assert(qcd_weight_mur2muf0p5_index != std::numeric_limits<int>::max());
    assert(qcd_weight_mur0p5muf1_index != std::numeric_limits<int>::max());
    assert(qcd_weight_mur0p5muf2_index != std::numeric_limits<int>::max());
    assert(qcd_weight_mur0p5muf0p5_index != std::numeric_limits<int>::max());

    qcd_pdf_weight_orig = hLheEvt->weights()[0].wgt;

    qcd_weight_mur1muf2 = hLheEvt->weights()[qcd_weight_mur1muf2_index].wgt;
    qcd_weight_mur1muf0p5 = hLheEvt->weights()[qcd_weight_mur1muf0p5_index].wgt;
    qcd_weight_mur2muf1 = hLheEvt->weights()[qcd_weight_mur2muf1_index].wgt;
    qcd_weight_mur2muf2 = hLheEvt->weights()[qcd_weight_mur2muf2_index].wgt;
    qcd_weight_mur2muf0p5 = hLheEvt->weights()[qcd_weight_mur2muf0p5_index].wgt;
    qcd_weight_mur0p5muf1 = hLheEvt->weights()[qcd_weight_mur0p5muf1_index].wgt;
    qcd_weight_mur0p5muf2 = hLheEvt->weights()[qcd_weight_mur0p5muf2_index].wgt;
    qcd_weight_mur0p5muf0p5 = hLheEvt->weights()[qcd_weight_mur0p5muf0p5_index].wgt;

    for (unsigned int i = 0; i < pdf_weights_indices.size(); i++){
      pdf_weights->push_back(hLheEvt->weights()[pdf_weights_indices[i]].wgt);
    }

  }

  if (mgreweightinfo_) {

    for(unsigned int i = 0; i < mgreweight_weights_indices.size(); i++){
      mgreweight_weights->push_back(hLheEvt->weights()[mgreweight_weights_indices[i]].wgt);
    }

  }

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

void 
lhe_and_gen::beginRun(edm::Run const& iRun)
{

  if(mgreweightinfo_ || syscalcinfo_){

  edm::Handle<LHERunInfoProduct> hLheRun;
  iRun.getByLabel(lheRunInfoLabel_,hLheRun);


  for ( LHERunInfoProduct::headers_const_iterator lheruniter = hLheRun.product()->headers_begin(); lheruniter != hLheRun.product()->headers_end(); lheruniter++ ) {
    
    if (lheruniter->tag() != "slha")
      continue;



    bool leading= true;

   for ( LHERunInfoProduct::Header::const_iterator iter = lheruniter->begin(); iter != lheruniter->end(); iter++ ) {

      slha_header_line_ = (*iter);

      //remove trailing empty line
      if (slha_header_line_ == "\n" && leading){
	leading = false;
	continue;
      }

      leading= false;

      slha_header_tree_->Fill();

   }

   slha_header_line_ = "END_SLHA_HEADER\n";

   slha_header_tree_->Fill();

   

  }


  for ( LHERunInfoProduct::headers_const_iterator lheruniter = hLheRun.product()->headers_begin(); lheruniter != hLheRun.product()->headers_end(); lheruniter++ ) {

    if (lheruniter->tag() != "initrwgt")
      continue;

    bool leading = true;

    for ( LHERunInfoProduct::Header::const_iterator iter = lheruniter->begin(); iter != lheruniter->end(); iter++ ) {

      std::cout << "(*iter) = " << (*iter) << std::endl;

      initrwgt_header_line_ = (*iter);

      //skip leading empty line
      if (initrwgt_header_line_ == "\n" && leading){
	leading = false;
	continue;
      }

      leading = false;

      initrwgt_header_tree_->Fill();

    }

    initrwgt_header_line_ = "END_INITRWGT_HEADER\n";

    initrwgt_header_tree_->Fill();

  }

  // use the lines below for /WLLJJ_WToLNu_MLL-4To60_EWK_TuneCUETP8M1_13TeV_madgraph-madspin-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/MINIAODSIM which does not have the initrwgt header because  madspin drops it for some reason

  /*

  qcd_weight_mur1muf2_index = 2 - 1;
  qcd_weight_mur1muf0p5_index = 3 - 1;
  qcd_weight_mur2muf1_index = 4 - 1;
  qcd_weight_mur2muf2_index= 5 - 1;
  qcd_weight_mur2muf0p5_index = 6 - 1;
  qcd_weight_mur0p5muf1_index= 7 - 1;
  qcd_weight_mur0p5muf2_index= 8 - 1;
  qcd_weight_mur0p5muf0p5_index= 9 - 1;

  for (int i = 10; i < 111; i++)
    pdf_weights_indices.push_back(i - 1);

  */

  for ( LHERunInfoProduct::headers_const_iterator lheruniter = hLheRun.product()->headers_begin(); lheruniter != hLheRun.product()->headers_end(); lheruniter++ ) {
    
    if (lheruniter->tag() != "initrwgt")
      continue;

    bool in_NNPDF23_lo_as_0130_qed = false;
    bool in_mg_reweighting = false;

    int weight_index = -1;

    for ( LHERunInfoProduct::Header::const_iterator iter = lheruniter->begin(); iter != lheruniter->end(); iter++ ) {

     if (  (*iter).find("weight id") != std::string::npos || (*iter).find("weight MUF=") != std::string::npos )
	weight_index++;

      if ( (*iter).find("mur=1 muf=2") != std::string::npos || (*iter).find("muR=0.10000E+01 muF=0.20000E+01") != std::string::npos){

	assert(weight_index != -1);

	qcd_weight_mur1muf2_index = weight_index;
	continue;
      }

      if ( (*iter).find("mur=2 muf=2") != std::string::npos || (*iter).find("muR=0.20000E+01 muF=0.20000E+01") != std::string::npos){

	assert(weight_index != -1);
	qcd_weight_mur2muf2_index=weight_index;

	continue;
      }

      if ( (*iter).find("mur=2 muf=1") != std::string::npos || (*iter).find("muR=0.20000E+01 muF=0.10000E+01") != std::string::npos ){

	assert(weight_index != -1);
	qcd_weight_mur2muf1_index = weight_index;
	continue;
      }

      if ( (*iter).find("mur=2 muf=0.5") != std::string::npos || (*iter).find("muR=0.20000E+01 muF=0.50000E+00") != std::string::npos  ){

	assert(weight_index != -1);

	qcd_weight_mur2muf0p5_index=weight_index;
	continue;
      }

      if ( (*iter).find("mur=0.5 muf=2") != std::string::npos || (*iter).find("muR=0.50000E+00 muF=0.20000E+01") != std::string::npos){

	assert(weight_index != -1);
	qcd_weight_mur0p5muf2_index=weight_index;
	continue;
      }

      if ( (*iter).find("mur=0.5 muf=0.5") != std::string::npos || (*iter).find("muR=0.50000E+00 muF=0.50000E+00") != std::string::npos){

	assert(weight_index != -1);
	qcd_weight_mur0p5muf0p5_index=weight_index;
	continue;
      }

      if ( (*iter).find("mur=1 muf=0.5") != std::string::npos || (*iter).find("muR=0.10000E+01 muF=0.50000E+00") != std::string::npos ){

	assert(weight_index != -1);
	qcd_weight_mur1muf0p5_index=weight_index;
	continue;
      }

      if ( (*iter).find("mur=0.5 muf=1") != std::string::npos || (*iter).find("muR=0.50000E+00 muF=0.10000E+01") != std::string::npos ){

	assert(weight_index != -1);
	qcd_weight_mur0p5muf1_index=weight_index;
	continue;
      }

      if ( ((*iter).find("\"NNPDF30_lo_as_0130.LHgrid\"") != std::string::npos || (*iter).find("\"NNPDF30_lo_as_0130\"") != std::string::npos) && (*iter).find("weightgroup") != std::string::npos && (*iter).find("hessian") != std::string::npos){

	in_NNPDF23_lo_as_0130_qed = true;
	continue;
      }

      if ( (*iter).find("mg_reweighting") != std::string::npos){
	in_mg_reweighting = true;
	continue;
      }

      if (in_mg_reweighting) {
	if ( (*iter).find("/weightgroup") != std::string::npos){
	  in_mg_reweighting = false;
	  continue;
	}


	//the mg reweight weights can take more than one line

	assert((*iter).find("set param_card") != std::string::npos || (*iter).find("</weight>") != std::string::npos);

	assert(weight_index != -1);

	if ((*iter).find("weight id") != std::string::npos)
	  mgreweight_weights_indices.push_back(weight_index);

      }

      if (in_NNPDF23_lo_as_0130_qed){
	if ( (*iter).find("/weightgroup") != std::string::npos){

	  in_NNPDF23_lo_as_0130_qed = false;
	  continue;
	}


	assert((*iter).find("Member") != std::string::npos);

	//std::cout << (*iter).substr((*iter).find("<weight id=")+12,(*iter).find(">Member") - (*iter).find("<weight id=") - 13) << std::endl;

	assert(weight_index != -1);

	pdf_weights_indices.push_back(weight_index);

      }

      //std::cout << (*iter) << std::endl;
    }

  }

  }

}
