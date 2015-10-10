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

inline Bool_t passSoftMuonId(const pat::Muon & muon, const reco::Vertex &vtx){

  return muon.isSoftMuon(vtx);

}

inline Bool_t passWLLJJVetoMuonId(const pat::Muon & muon, const reco::Vertex &vtx){

  Float_t relative_isolation = ( muon.pfIsolationR04().sumChargedHadronPt+ std::max(0.0,muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5 * muon.pfIsolationR04().sumPUPt) )/muon.pt();  

  if (relative_isolation < 0.1 && fabs(muon.muonBestTrack()->dxy(vtx.position())) < 0.02 && fabs(muon.muonBestTrack()->dz(vtx.position())) < 0.1)
    return true;
  else
    return false;


}

inline Bool_t passWLLJJVetoElectronId(const pat::Electron & el, const reco::Vertex &PV){

  reco::GsfElectron::PflowIsolationVariables pfIso = el.pfIsolationVariables();
  
  Float_t absiso = pfIso.sumChargedHadronPt + std::max(0.0 , pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - 0.5 * pfIso.sumPUPt );
  Float_t relIsoWithDBeta = absiso/el.pt();
  
  if (relIsoWithDBeta < 0.1 && fabs(el.gsfTrack()->dz( PV.position() )) < 0.02 && fabs((-1) * el.gsfTrack()->dxy(PV.position())) < 0.1 )
    return true;
  else
    return false;

}


inline Bool_t passLooseMuonId(const pat::Muon & muon, const reco::Vertex &vtx) {

  if(!muon.isPFMuon() || !muon.isGlobalMuon() ) return false;

  bool muID = muon.isGlobalMuon() && muon.globalTrack()->normalizedChi2()<10. && muon.globalTrack()->hitPattern().numberOfValidMuonHits() >0 && (muon.numberOfMatchedStations() > 1);
  
  bool hits = muon.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 &&
    muon.innerTrack()->hitPattern().numberOfValidPixelHits() > 0; 
  
  bool ip = fabs(muon.muonBestTrack()->dxy(vtx.position())) < 0.2 && fabs(muon.muonBestTrack()->dz(vtx.position())) < 0.5;

  return muID && hits && ip;

}


inline Bool_t passLooseMuonSelection(const pat::Muon & muon, const reco::Vertex &vtx) {

  Float_t relative_isolation = ( muon.pfIsolationR04().sumChargedHadronPt+ std::max(0.0,muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5 * muon.pfIsolationR04().sumPUPt) )/muon.pt();

  return passLooseMuonId(muon,vtx) && relative_isolation < 1.0;

}


inline Bool_t passVeryLooseMuonSelection(const pat::Muon & muon, const reco::Vertex &vtx ){

  return passLooseMuonSelection(muon,vtx);

}

inline Bool_t passTightMuonIdV1(const pat::Muon & muon, const reco::Vertex &vtx) {

  if(!muon.isPFMuon() || !muon.isGlobalMuon() ) return false;

  bool muID = muon.isGlobalMuon() && muon.globalTrack()->normalizedChi2()<10. && muon.globalTrack()->hitPattern().numberOfValidMuonHits() >0 && (muon.numberOfMatchedStations() > 1);
  
  bool hits = muon.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 &&
    muon.innerTrack()->hitPattern().numberOfValidPixelHits() > 0; 

  //tighten the d0 cut from 0.2 to 0.02 because the muon POG ids consider muons from bs as real
  
  bool ip = fabs(muon.muonBestTrack()->dxy(vtx.position())) < 0.02 && fabs(muon.muonBestTrack()->dz(vtx.position())) < 0.5;
  
  return muID && hits && ip;

}


inline Bool_t passTightMuonSelectionV1(const pat::Muon & muon, const reco::Vertex &vtx) {

  Float_t relative_isolation = ( muon.pfIsolationR04().sumChargedHadronPt+ std::max(0.0,muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5 * muon.pfIsolationR04().sumPUPt) )/muon.pt();

  return passTightMuonIdV1(muon,vtx) && relative_isolation < 0.2;

}

inline Bool_t passTightMuonSelectionV2(const pat::Muon & muon, const reco::Vertex &vtx) {

  Float_t relative_isolation = ( muon.pfIsolationR04().sumChargedHadronPt+ std::max(0.0,muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5 * muon.pfIsolationR04().sumPUPt) )/muon.pt();

  return passTightMuonIdV1(muon,vtx) && relative_isolation < 0.12;

}

inline Bool_t passTightMuonSelectionV3(const pat::Muon & muon, const reco::Vertex &vtx) {

  Float_t relative_isolation = ( muon.pfIsolationR04().sumChargedHadronPt+ std::max(0.0,muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5 * muon.pfIsolationR04().sumPUPt) )/muon.pt();

  return passTightMuonIdV1(muon,vtx) && relative_isolation < 0.4;

}


inline Bool_t passTightElectronSelection(const pat::Electron & el, const reco::Vertex &PV) {

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
	  (fabs(el.deltaEtaSuperClusterTrackAtVtx()) < 0.00733 )
	  &&
	  ( fabs(el.deltaPhiSuperClusterTrackAtVtx()) < 0.114)
	  &&
	  (el.full5x5_sigmaIetaIeta() < 0.0283)
	  &&
	  ( el.hcalOverEcal() < 0.0678)
	  &&
	  ( fabs((-1) * el.gsfTrack()->dxy(PV.position())) < 0.0739 )
	  &&
	  (  fabs(el.gsfTrack()->dz( PV.position() )) < 0.602 )
	  &&
	  (fabs(ooEmooP) < 0.0898)
	  &&
	  (relIsoWithDBeta < 0.0678)
	  &&
	  (el.passConversionVeto())
	  &&
	  (el.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) <= 1 )
	  )
	 pass = kTRUE;
     } else if (el.superCluster()->eta() < 1.479) {
       if(
	  (fabs(el.deltaEtaSuperClusterTrackAtVtx()) < 0.0103)
	  &&
	  ( fabs(el.deltaPhiSuperClusterTrackAtVtx()) < 0.0336)
	  &&
	  (el.full5x5_sigmaIetaIeta() < 0.0101)
	  &&
	  ( el.hcalOverEcal() < 0.0876)
	  &&
	  ( fabs((-1) * el.gsfTrack()->dxy(PV.position())) < 0.0118)
	  &&
	  (  fabs(el.gsfTrack()->dz( PV.position() )) < 0.373)
	  &&
	  (fabs(ooEmooP) < 0.0174)
	  &&
	  (relIsoWithDBeta < 0.0766)
	  &&
	  (el.passConversionVeto())
	  &&
	  (el.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) <= 2  )
	  )
	 pass = kTRUE;
     } 

  return pass;

}    

inline Bool_t passLooseElectronSelectionV1(const pat::Electron & el, const reco::Vertex &PV) {

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
	  (fabs(el.deltaEtaSuperClusterTrackAtVtx()) < 0.00733 )
	  &&
	  ( fabs(el.deltaPhiSuperClusterTrackAtVtx()) < 0.114)
	  &&
	  (el.full5x5_sigmaIetaIeta() < 0.0283)
	  &&
	  ( el.hcalOverEcal() < 0.0678)
	  &&
	  //	  ( fabs((-1) * el.gsfTrack()->dxy(PV.position())) < 0.0739 )
	  //	  &&
	  (  fabs(el.gsfTrack()->dz( PV.position() )) < 0.602 )
	  &&
	  (fabs(ooEmooP) < 0.0898)
	  &&
	  (relIsoWithDBeta < 0.0678*6)
	  &&
	  (el.passConversionVeto())
	  &&
	  (el.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) <= 1 )
	  )
	 pass = kTRUE;
     } else if (el.superCluster()->eta() < 1.479) {
       if(
	  (fabs(el.deltaEtaSuperClusterTrackAtVtx()) < 0.0103)
	  &&
	  ( fabs(el.deltaPhiSuperClusterTrackAtVtx()) < 0.0336)
	  &&
	  (el.full5x5_sigmaIetaIeta() < 0.0101)
	  &&
	  ( el.hcalOverEcal() < 0.0876)
	  &&
	  //	  ( fabs((-1) * el.gsfTrack()->dxy(PV.position())) < 0.0118)
	  //	  &&
	  (  fabs(el.gsfTrack()->dz( PV.position() )) < 0.373)
	  &&
	  (fabs(ooEmooP) < 0.0174)
	  &&
	  (relIsoWithDBeta < 0.0766*6)
	  &&
	  (el.passConversionVeto())
	  &&
	  (el.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) <= 2  )
	  )
	 pass = kTRUE;
     } 

  return pass;

}    

inline Bool_t passLooseElectronSelectionV2(const pat::Electron & el, const reco::Vertex &PV) {

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
	  (fabs(el.deltaEtaSuperClusterTrackAtVtx()) < 0.00733 )
	  &&
	  ( fabs(el.deltaPhiSuperClusterTrackAtVtx()) < 0.114)
	  &&
	  (el.full5x5_sigmaIetaIeta() < 0.0283)
	  &&
	  ( el.hcalOverEcal() < 0.0678)
	  &&
	  //	  ( fabs((-1) * el.gsfTrack()->dxy(PV.position())) < 0.0739 )
	  //	  &&
	  (  fabs(el.gsfTrack()->dz( PV.position() )) < 0.602 )
	  &&
	  (fabs(ooEmooP) < 0.0898)
	  &&
	  (relIsoWithDBeta < 0.0678*10)
	  &&
	  (el.passConversionVeto())
	  &&
	  (el.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) <= 1 )
	  )
	 pass = kTRUE;
     } else if (el.superCluster()->eta() < 1.479) {
       if(
	  (fabs(el.deltaEtaSuperClusterTrackAtVtx()) < 0.0103)
	  &&
	  ( fabs(el.deltaPhiSuperClusterTrackAtVtx()) < 0.0336)
	  &&
	  (el.full5x5_sigmaIetaIeta() < 0.0101)
	  &&
	  ( el.hcalOverEcal() < 0.0876)
	  &&
	  //	  ( fabs((-1) * el.gsfTrack()->dxy(PV.position())) < 0.0118)
	  //	  &&
	  (  fabs(el.gsfTrack()->dz( PV.position() )) < 0.373)
	  &&
	  (fabs(ooEmooP) < 0.0174)
	  &&
	  (relIsoWithDBeta < 0.0766*10)
	  &&
	  (el.passConversionVeto())
	  &&
	  (el.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) <= 2  )
	  )
	 pass = kTRUE;
     } 



  return pass;

}    

inline Bool_t passVeryLooseElectronSelection(const pat::Electron & el, const reco::Vertex &PV) {

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
	  (fabs(el.deltaEtaSuperClusterTrackAtVtx()) < 0.00733 )
	  &&
	  ( fabs(el.deltaPhiSuperClusterTrackAtVtx()) < 0.114)
	  &&
	  (el.full5x5_sigmaIetaIeta() < 0.0283)
	  &&
	  ( el.hcalOverEcal() < 0.0678)
	  &&
	  //	  ( fabs((-1) * el.gsfTrack()->dxy(PV.position())) < 0.0739 )
	  //	  &&
	  (  fabs(el.gsfTrack()->dz( PV.position() )) < 0.602 )
	  &&
	  (fabs(ooEmooP) < 0.0898)
	  &&
	  (relIsoWithDBeta < 0.0678*10)
	  &&
	  (el.passConversionVeto())
	  &&
	  (el.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) <= 1 )
	  )
	 pass = kTRUE;
     } else if (el.superCluster()->eta() < 1.479) {
       if(
	  (fabs(el.deltaEtaSuperClusterTrackAtVtx()) < 0.0103)
	  &&
	  ( fabs(el.deltaPhiSuperClusterTrackAtVtx()) < 0.0336)
	  &&
	  (el.full5x5_sigmaIetaIeta() < 0.0101)
	  &&
	  ( el.hcalOverEcal() < 0.0876)
	  &&
	  //	  ( fabs((-1) * el.gsfTrack()->dxy(PV.position())) < 0.0118)
	  //	  &&
	  (  fabs(el.gsfTrack()->dz( PV.position() )) < 0.373)
	  &&
	  (fabs(ooEmooP) < 0.0174)
	  &&
	  (relIsoWithDBeta < 0.0766*10)
	  &&
	  (el.passConversionVeto())
	  &&
	  (el.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) <= 2  )
	  )
	 pass = kTRUE;
     } 

  return pass;

}    


#endif

