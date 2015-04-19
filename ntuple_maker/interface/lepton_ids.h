#ifndef NTUPLEMAKER_LEPTONIDS_H
#define NTUPLEMAKER_LEPTONIDS_H

#include <memory>

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

inline Bool_t passTightMuonId(const pat::Muon & muon, const reco::Vertex &vtx) {

  if(!muon.isPFMuon() || !muon.isGlobalMuon() ) return false;

  bool muID = muon.isGlobalMuon() && muon.globalTrack()->normalizedChi2()<10. && muon.globalTrack()->hitPattern().numberOfValidMuonHits() >0 && (muon.numberOfMatchedStations() > 1);
    
  
  bool hits = muon.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 &&
    muon.innerTrack()->hitPattern().numberOfValidPixelHits() > 0; 

  
  bool ip = fabs(muon.muonBestTrack()->dxy(vtx.position())) < 0.2 && fabs(muon.muonBestTrack()->dz(vtx.position())) < 0.5;
  
  return muID && hits && ip;



}

inline Bool_t passLooseMuonId(const pat::Muon & muon, const reco::Vertex &vtx) {

  if(!muon.isPFMuon() || !muon.isGlobalMuon() ) return false;

  bool muID = muon.isGlobalMuon() && muon.globalTrack()->normalizedChi2()<10. && muon.globalTrack()->hitPattern().numberOfValidMuonHits() >0 && (muon.numberOfMatchedStations() > 1);
    
  
  bool hits = muon.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 &&
    muon.innerTrack()->hitPattern().numberOfValidPixelHits() > 0; 

  
  bool ip = fabs(muon.muonBestTrack()->dxy(vtx.position())) < 0.2 && fabs(muon.muonBestTrack()->dz(vtx.position())) < 0.5;
  
  return muID && hits && ip;

}

inline Bool_t passTightElectronId(const pat::Electron & el, const reco::Vertex &PV) {

  Bool_t pass = kFALSE;

  if (!el.chargeInfo().isGsfCtfScPixConsistent)
    return kFALSE;

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
       ooEmooP = fabs(1.0/el.ecalEnergy() - el.eSuperClusterOverP()/el.ecalEnergy() );

     //medium working point from here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
     if(el.superCluster()->eta() < 2.5 && el.superCluster()->eta() > 1.479 ){
       if(
	  (fabs(el.deltaEtaSuperClusterTrackAtVtx()) < 0.009285)
	  &&
	  ( fabs(el.deltaPhiSuperClusterTrackAtVtx()) < 0.042447)
	  &&
	  (el.full5x5_sigmaIetaIeta() < 0.029524)
	  &&
	  ( el.hcalOverEcal() < 0.104263)
	  &&
	  ( fabs((-1) * el.gsfTrack()->dxy(PV.position())) < 0.051682)
	  &&
	  (  fabs(el.gsfTrack()->dz( PV.position() )) < 0.180720)
	  &&
	  (fabs(ooEmooP) < 0.137468)
	  &&
	  (relIsoWithDBeta < 0.116708)
	  &&
	  (el.passConversionVeto())
	  &&
	  (el.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) <= 1 )
	  )
	 pass = kTRUE;
     } else if (el.superCluster()->eta() < 1.479) {
       if(
	  (fabs(el.deltaEtaSuperClusterTrackAtVtx()) < 0.007641)
	  &&
	  ( fabs(el.deltaPhiSuperClusterTrackAtVtx()) < 0.032643)
	  &&
	  (el.full5x5_sigmaIetaIeta() < 0.010399)
	  &&
	  ( el.hcalOverEcal() < 0.060662)
	  &&
	  ( fabs((-1) * el.gsfTrack()->dxy(PV.position())) < 0.011811)
	  &&
	  (  fabs(el.gsfTrack()->dz( PV.position() )) < 0.070775)
	  &&
	  (fabs(ooEmooP) < 0.153897)
	  &&
	  (relIsoWithDBeta < 0.097213)
	  &&
	  (el.passConversionVeto())
	  &&
	  (el.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) <= 1  )
	  )
	 pass = kTRUE;
     } 

  return pass;

}    

inline Bool_t passLooseElectronId(const pat::Electron & el, const reco::Vertex &PV) {

  Bool_t pass = kFALSE;

  if (!el.chargeInfo().isGsfCtfScPixConsistent)
    return kFALSE;

     reco::GsfElectron::PflowIsolationVariables pfIso = el.pfIsolationVariables();

     Float_t absiso = pfIso.sumChargedHadronPt + std::max(0.0 , pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - 0.5 * pfIso.sumPUPt );
     Float_t relIsoWithDBeta = absiso/el.pt();
     Float_t ooEmooP = 0;

     if( el.ecalEnergy() == 0 ){
       std::cout << "el.pt() = " << el.pt() << std::endl;
       std::cout << "Electron energy is zero!" << std::endl;
       ooEmooP = 1e30;
     }
     else if (!std::isfinite(el.ecalEnergy())){
       std::cout << "Electron energy is not finite!" << std::endl;
       ooEmooP = 1e30;
     }
     else
       ooEmooP = fabs(1.0/el.ecalEnergy() - el.eSuperClusterOverP()/el.ecalEnergy() );

     //medium working point from here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
     if(el.superCluster()->eta() < 2.5 && el.superCluster()->eta() > 1.479 ){
       if(
	  (fabs(el.deltaEtaSuperClusterTrackAtVtx()) < 0.009285)
	  &&
	  ( fabs(el.deltaPhiSuperClusterTrackAtVtx()) < 0.042447)
	  &&
	  (el.full5x5_sigmaIetaIeta() < 0.029524)
	  &&
	  ( el.hcalOverEcal() < 0.104263)
	  &&
	  //	  ( fabs((-1) * el.gsfTrack()->dxy(PV.position())) < 0.051682)
	  //	  &&
	  (  fabs(el.gsfTrack()->dz( PV.position() )) < 0.180720)
	  &&
	  (fabs(ooEmooP) < 0.137468)
	  &&
	  (relIsoWithDBeta < 0.116708*6)
	  &&
	  (el.passConversionVeto())
	  &&
	  (el.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) <= 1 )
	  )
	 pass = kTRUE;
     } else if (el.superCluster()->eta() < 1.479) {
       if(
	  (fabs(el.deltaEtaSuperClusterTrackAtVtx()) < 0.007641)
	  &&
	  ( fabs(el.deltaPhiSuperClusterTrackAtVtx()) < 0.032643)
	  &&
	  (el.full5x5_sigmaIetaIeta() < 0.010399)
	  &&
	  ( el.hcalOverEcal() < 0.060662)
	  &&
	  //	  ( fabs((-1) * el.gsfTrack()->dxy(PV.position())) < 0.011811)
	  //	  &&
	  (  fabs(el.gsfTrack()->dz( PV.position() )) < 0.070775)
	  &&
	  (fabs(ooEmooP) < 0.153897)
	  &&
	  (relIsoWithDBeta < 0.097213*6)
	  &&
	  (el.passConversionVeto())
	  &&
	  (el.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) <= 1  )
	  )
	 pass = kTRUE;
     } 

  return pass;

}    


#endif
