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



inline Float_t muon_isolation(const pat::Muon & muon, const reco::Vertex &vtx) {

  return ( muon.pfIsolationR04().sumChargedHadronPt+ std::max(0.0,muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5 * muon.pfIsolationR04().sumPUPt) )/muon.pt();  

}

inline Float_t electron_isolation(const pat::Electron & el, const reco::Vertex &PV, const double &rho) {

     reco::GsfElectron::PflowIsolationVariables pfIso = el.pfIsolationVariables();

     //see here https://indico.cern.ch/event/369239/contribution/4/attachments/1134761/1623262/talk_effective_areas_25ns.pdf slide 12

     double EffectiveArea = 0;

     float abseta = fabs(fabs(el.superCluster()->eta()));
     if (abseta >= 0.0 && abseta < 1.0 ) EffectiveArea = 0.1752;
     if (abseta >= 1.0 && abseta < 1.479 ) EffectiveArea = 0.1862;
     if (abseta >= 1.479 && abseta < 2.0 ) EffectiveArea = 0.1411;
     if (abseta >= 2.0 && abseta < 2.2 ) EffectiveArea = 0.1534;
     if (abseta >= 2.2 && abseta < 2.3 ) EffectiveArea = 0.1903;
     if (abseta >= 2.3 && abseta < 2.4 ) EffectiveArea = 0.2243;
     if (abseta >= 2.4 && abseta < 5.0 ) EffectiveArea = 0.2687;

     //EffectiveArea = 0;

     Float_t absiso = pfIso.sumChargedHadronPt + std::max(0.0 , pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - rho * EffectiveArea );
     Float_t relIso = absiso/el.pt();

     return relIso;

}

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
  
  bool ip = fabs(muon.muonBestTrack()->dxy(vtx.position())) < 0.02 && fabs(muon.muonBestTrack()->dz(vtx.position())) < 0.5;

  return muID && hits && ip;

}


inline Bool_t passLooseMuonSelectionV1(const pat::Muon & muon, const reco::Vertex &vtx) {

  Float_t relative_isolation = ( muon.pfIsolationR04().sumChargedHadronPt+ std::max(0.0,muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5 * muon.pfIsolationR04().sumPUPt) )/muon.pt();

  return passLooseMuonId(muon,vtx) && relative_isolation < 1.0;

}

inline Bool_t passLooseMuonSelectionV2(const pat::Muon & muon, const reco::Vertex &vtx) {

  if(!muon.isGlobalMuon() ) return false;

  return fabs(muon.muonBestTrack()->dxy(vtx.position())) < 0.5 && fabs(muon.muonBestTrack()->dz(vtx.position())) < 1;

}

inline Bool_t passLooseMuonSelectionV3(const pat::Muon & muon, const reco::Vertex &vtx) {

  if(!muon.isPFMuon() || !muon.isGlobalMuon() ) return false;

  bool muID = muon.isGlobalMuon() && muon.globalTrack()->normalizedChi2()<10. && muon.globalTrack()->hitPattern().numberOfValidMuonHits() >0 && (muon.numberOfMatchedStations() > 1);
  
  bool hits = muon.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 &&
    muon.innerTrack()->hitPattern().numberOfValidPixelHits() > 0; 

  bool ip = fabs(muon.muonBestTrack()->dxy(vtx.position())) < 0.5 && fabs(muon.muonBestTrack()->dz(vtx.position())) < 1;

  return muID && hits && ip;

}

inline Bool_t passLooseMuonSelectionV4(const pat::Muon & muon, const reco::Vertex &vtx) {

  if(!muon.isPFMuon() || !muon.isGlobalMuon() ) return false;

  bool muID = muon.isGlobalMuon() && muon.globalTrack()->normalizedChi2()<10. && muon.globalTrack()->hitPattern().numberOfValidMuonHits() >0 && (muon.numberOfMatchedStations() > 1);
  
  bool hits = muon.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 &&
    muon.innerTrack()->hitPattern().numberOfValidPixelHits() > 0; 

  bool ip =  fabs(muon.muonBestTrack()->dxy(vtx.position())) < 0.5 && fabs(muon.muonBestTrack()->dz(vtx.position())) < 1;

  Float_t relative_isolation = ( muon.pfIsolationR04().sumChargedHadronPt+ std::max(0.0,muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5 * muon.pfIsolationR04().sumPUPt) )/muon.pt();

  bool iso = relative_isolation < 1.0;

  return muID && hits && ip && iso;  

}

inline Bool_t passLooseMuonSelectionV5(const pat::Muon & muon, const reco::Vertex &vtx) {

  return fabs(muon.muonBestTrack()->dxy(vtx.position())) < 0.5 && fabs(muon.muonBestTrack()->dz(vtx.position())) < 1;

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

//should be the same as the member function is istightmuon()
inline Bool_t passTightMuonIdV2(const pat::Muon & muon, const reco::Vertex &vtx) {

  if(!muon.isPFMuon() || !muon.isGlobalMuon() ) return false;

  bool muID = muon.isGlobalMuon() && muon.globalTrack()->normalizedChi2()<10. && muon.globalTrack()->hitPattern().numberOfValidMuonHits() >0 && (muon.numberOfMatchedStations() > 1);
  
  bool hits = muon.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 &&
    muon.innerTrack()->hitPattern().numberOfValidPixelHits() > 0; 

  bool ip = fabs(muon.muonBestTrack()->dxy(vtx.position())) < 0.2 && fabs(muon.muonBestTrack()->dz(vtx.position())) < 0.5;
  
  return muID && hits && ip;

}


inline Bool_t passTightMuonSelectionV1(const pat::Muon & muon, const reco::Vertex &vtx) {

  Float_t relative_isolation = ( muon.pfIsolationR04().sumChargedHadronPt+ std::max(0.0,muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5 * muon.pfIsolationR04().sumPUPt) )/muon.pt();

  return passTightMuonIdV1(muon,vtx) && relative_isolation < 0.15;

}

inline Bool_t passTightMuonSelectionV2(const pat::Muon & muon, const reco::Vertex &vtx) {

  Float_t relative_isolation = ( muon.pfIsolationR04().sumChargedHadronPt+ std::max(0.0,muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5 * muon.pfIsolationR04().sumPUPt) )/muon.pt();

  return passTightMuonIdV1(muon,vtx) && relative_isolation < 0.25;

}

inline Bool_t passTightMuonSelectionV3(const pat::Muon & muon, const reco::Vertex &vtx) {

  Float_t relative_isolation = ( muon.pfIsolationR04().sumChargedHadronPt+ std::max(0.0,muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5 * muon.pfIsolationR04().sumPUPt) )/muon.pt();

  return passTightMuonIdV1(muon,vtx) && relative_isolation < 0.4;

}


inline Bool_t passTightElectronSelectionV1(const pat::Electron & el, const reco::Vertex &PV, const double &rho) {

  Bool_t pass = kFALSE;

  if (!el.chargeInfo().isGsfCtfScPixConsistent)
    return kFALSE;

     reco::GsfElectron::PflowIsolationVariables pfIso = el.pfIsolationVariables();

     //see here https://indico.cern.ch/event/369239/contribution/4/attachments/1134761/1623262/talk_effective_areas_25ns.pdf slide 12

     double EffectiveArea = 0;

     float abseta = fabs(fabs(el.superCluster()->eta()));
     if (abseta >= 0.0 && abseta < 1.0 ) EffectiveArea = 0.1752;
     if (abseta >= 1.0 && abseta < 1.479 ) EffectiveArea = 0.1862;
     if (abseta >= 1.479 && abseta < 2.0 ) EffectiveArea = 0.1411;
     if (abseta >= 2.0 && abseta < 2.2 ) EffectiveArea = 0.1534;
     if (abseta >= 2.2 && abseta < 2.3 ) EffectiveArea = 0.1903;
     if (abseta >= 2.3 && abseta < 2.4 ) EffectiveArea = 0.2243;
     if (abseta >= 2.4 && abseta < 5.0 ) EffectiveArea = 0.2687;

     //EffectiveArea = 0;

     Float_t absiso = pfIso.sumChargedHadronPt + std::max(0.0 , pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - rho * EffectiveArea );
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
     if(fabs(el.superCluster()->eta()) < 2.5 && fabs(el.superCluster()->eta()) > 1.479 ){
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
     } else if (fabs(el.superCluster()->eta()) < 1.479) {
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

inline Bool_t passTightElectronSelectionV2(const pat::Electron & el, const reco::Vertex &PV, const double &rho) {

  Bool_t pass = kFALSE;

  if (!el.chargeInfo().isGsfCtfScPixConsistent)
    return kFALSE;


     reco::GsfElectron::PflowIsolationVariables pfIso = el.pfIsolationVariables();

     //see here https://indico.cern.ch/event/369239/contribution/4/attachments/1134761/1623262/talk_effective_areas_25ns.pdf slide 12

     double EffectiveArea = 0;

     //float abseta = fabs(el.eta());
     float abseta = fabs(fabs(el.superCluster()->eta()));

     if (abseta >= 0.0 && abseta < 1.0 ) EffectiveArea = 0.1752;
     if (abseta >= 1.0 && abseta < 1.479 ) EffectiveArea = 0.1862;
     if (abseta >= 1.479 && abseta < 2.0 ) EffectiveArea = 0.1411;
     if (abseta >= 2.0 && abseta < 2.2 ) EffectiveArea = 0.1534;
     if (abseta >= 2.2 && abseta < 2.3 ) EffectiveArea = 0.1903;
     if (abseta >= 2.3 && abseta < 2.4 ) EffectiveArea = 0.2243;
     if (abseta >= 2.4 && abseta < 5.0 ) EffectiveArea = 0.2687;

     Float_t absiso = pfIso.sumChargedHadronPt + std::max(0.0 , pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - rho * EffectiveArea );
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

     /*

     std::cout << "rho = " << rho << std::endl;
     std::cout << el.superCluster()->eta() << std::endl;
     std::cout << fabs(el.deltaEtaSuperClusterTrackAtVtx()) << std::endl;
     std::cout << fabs(el.deltaPhiSuperClusterTrackAtVtx()) << std::endl;
     std::cout << el.full5x5_sigmaIetaIeta() << std::endl;
     std::cout << el.hcalOverEcal() << std::endl;
     std::cout << fabs((-1) * el.gsfTrack()->dxy(PV.position())) << std::endl;
     std::cout << fabs(el.gsfTrack()->dz( PV.position() )) << std::endl;
     std::cout << fabs(ooEmooP) << std::endl;
     std::cout << relIsoWithDBeta  << std::endl;
     std::cout << el.passConversionVeto() << std::endl;
     std::cout << el.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) << std::endl;

     */

     //tight working point from here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
     if(fabs(el.superCluster()->eta()) < 2.5 && fabs(el.superCluster()->eta()) > 1.479 ){
       if(
	  (fabs(el.deltaEtaSuperClusterTrackAtVtx()) < 0.00724 )
	  &&
	  ( fabs(el.deltaPhiSuperClusterTrackAtVtx()) < 0.0918)
	  &&
	  (el.full5x5_sigmaIetaIeta() < 0.0279)
	  &&
	  ( el.hcalOverEcal() < 0.0646)
	  &&
	  ( fabs((-1) * el.gsfTrack()->dxy(PV.position())) <  0.0351)
	  &&
	  (  fabs(el.gsfTrack()->dz( PV.position() )) <  0.471)
	  &&
	  (fabs(ooEmooP) < 0.00999)
	  &&
	  (relIsoWithDBeta < 0.0646)
	  &&
	  (el.passConversionVeto())
	  &&
	  (el.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) <= 1 )
	  )
	 pass = kTRUE;
     } else if (fabs(el.superCluster()->eta()) < 1.479) {
       if(
	  (fabs(el.deltaEtaSuperClusterTrackAtVtx()) < 0.00926)
	  &&
	  ( fabs(el.deltaPhiSuperClusterTrackAtVtx()) < 0.0336)
	  &&
	  (el.full5x5_sigmaIetaIeta() < 0.0101)
	  &&
	  ( el.hcalOverEcal() < 0.0597)
	  &&
	  ( fabs((-1) * el.gsfTrack()->dxy(PV.position())) < 0.0111)
	  &&
	  (  fabs(el.gsfTrack()->dz( PV.position() )) < 0.0466)
	  &&
	  (fabs(ooEmooP) < 0.012)
	  &&
	  (relIsoWithDBeta < 0.0354)
	  &&
	  (el.passConversionVeto())
	  &&
	  (el.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) <= 2 )
	  )
	 pass = kTRUE;
     } 

  return pass;

}    

inline Bool_t passTightElectronSelectionV3(const pat::Electron & el, const reco::Vertex &PV, const double &rho) {

  return el.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-tight");

}    

inline Bool_t passLooseElectronSelectionV1(const pat::Electron & el, const reco::Vertex &PV, const double &rho) {

  Bool_t pass = kFALSE;

  if (!el.chargeInfo().isGsfCtfScPixConsistent)
    return kFALSE;

     reco::GsfElectron::PflowIsolationVariables pfIso = el.pfIsolationVariables();

     double EffectiveArea = 0;

     //float abseta = fabs(el.eta());

     float abseta = fabs(fabs(el.superCluster()->eta()));


     if (abseta >= 0.0 && abseta < 1.0 ) EffectiveArea = 0.1752;
     if (abseta >= 1.0 && abseta < 1.479 ) EffectiveArea = 0.1862;
     if (abseta >= 1.479 && abseta < 2.0 ) EffectiveArea = 0.1411;
     if (abseta >= 2.0 && abseta < 2.2 ) EffectiveArea = 0.1534;
     if (abseta >= 2.2 && abseta < 2.3 ) EffectiveArea = 0.1903;
     if (abseta >= 2.3 && abseta < 2.4 ) EffectiveArea = 0.2243;
     if (abseta >= 2.4 && abseta < 5.0 ) EffectiveArea = 0.2687;

     Float_t absiso = pfIso.sumChargedHadronPt + std::max(0.0 , pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - rho * EffectiveArea );
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
     if(fabs(el.superCluster()->eta()) < 2.5 && fabs(el.superCluster()->eta()) > 1.479 ){
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
     } else if (fabs(el.superCluster()->eta()) < 1.479) {
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

inline Bool_t passLooseElectronSelectionV2(const pat::Electron & el, const reco::Vertex &PV, const double &rho) {

  Bool_t pass = kFALSE;

  if (!el.chargeInfo().isGsfCtfScPixConsistent)
    return kFALSE;

     reco::GsfElectron::PflowIsolationVariables pfIso = el.pfIsolationVariables();

     double EffectiveArea = 0;

     //float abseta = fabs(el.eta());

     float abseta = fabs(fabs(el.superCluster()->eta()));


     if (abseta >= 0.0 && abseta < 1.0 ) EffectiveArea = 0.1752;
     if (abseta >= 1.0 && abseta < 1.479 ) EffectiveArea = 0.1862;
     if (abseta >= 1.479 && abseta < 2.0 ) EffectiveArea = 0.1411;
     if (abseta >= 2.0 && abseta < 2.2 ) EffectiveArea = 0.1534;
     if (abseta >= 2.2 && abseta < 2.3 ) EffectiveArea = 0.1903;
     if (abseta >= 2.3 && abseta < 2.4 ) EffectiveArea = 0.2243;
     if (abseta >= 2.4 && abseta < 5.0 ) EffectiveArea = 0.2687;

     Float_t absiso = pfIso.sumChargedHadronPt + std::max(0.0 , pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - rho * EffectiveArea );
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
     if(fabs(el.superCluster()->eta()) < 2.5 && fabs(el.superCluster()->eta()) > 1.479 ){
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
     } else if (fabs(el.superCluster()->eta()) < 1.479) {
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

inline Bool_t passLooseElectronSelectionV3(const pat::Electron & el, const reco::Vertex &PV, const double &rho) {

  Bool_t pass = kFALSE;

  if (!el.chargeInfo().isGsfCtfScPixConsistent)
    return kFALSE;

     reco::GsfElectron::PflowIsolationVariables pfIso = el.pfIsolationVariables();

     //see here https://indico.cern.ch/event/369239/contribution/4/attachments/1134761/1623262/talk_effective_areas_25ns.pdf slide 12

     double EffectiveArea = 0;

     //float abseta = fabs(el.eta());
     float abseta = fabs(fabs(el.superCluster()->eta()));

     if (abseta >= 0.0 && abseta < 1.0 ) EffectiveArea = 0.1752;
     if (abseta >= 1.0 && abseta < 1.479 ) EffectiveArea = 0.1862;
     if (abseta >= 1.479 && abseta < 2.0 ) EffectiveArea = 0.1411;
     if (abseta >= 2.0 && abseta < 2.2 ) EffectiveArea = 0.1534;
     if (abseta >= 2.2 && abseta < 2.3 ) EffectiveArea = 0.1903;
     if (abseta >= 2.3 && abseta < 2.4 ) EffectiveArea = 0.2243;
     if (abseta >= 2.4 && abseta < 5.0 ) EffectiveArea = 0.2687;

     Float_t absiso = pfIso.sumChargedHadronPt + std::max(0.0 , pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - rho * EffectiveArea );
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

     //tight working point from here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
     if(fabs(el.superCluster()->eta()) < 2.5 && fabs(el.superCluster()->eta()) > 1.479 ){
       if(
	  (fabs(el.deltaEtaSuperClusterTrackAtVtx()) < 0.00724 )
	  &&
	  ( fabs(el.deltaPhiSuperClusterTrackAtVtx()) < 0.0918)
	  &&
	  (el.full5x5_sigmaIetaIeta() < 0.0279)
	  &&
	  ( el.hcalOverEcal() < 0.0646)
	  &&
	  //	  ( fabs((-1) * el.gsfTrack()->dxy(PV.position())) <  0.0351)
	  //	  &&
	  (  fabs(el.gsfTrack()->dz( PV.position() )) <  0.471)
	  &&
	  (fabs(ooEmooP) < 0.00999)
	  &&
	  (relIsoWithDBeta < 10*0.0646)
	  &&
	  (el.passConversionVeto())
	  &&
	  (el.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) <= 1 )
	  )
	 pass = kTRUE;
     } else if (fabs(el.superCluster()->eta()) < 1.479) {
       if(
	  (fabs(el.deltaEtaSuperClusterTrackAtVtx()) < 0.00926)
	  &&
	  ( fabs(el.deltaPhiSuperClusterTrackAtVtx()) < 0.0336)
	  &&
	  (el.full5x5_sigmaIetaIeta() < 0.0101)
	  &&
	  ( el.hcalOverEcal() < 0.0597)
	  &&
	  //	  ( fabs((-1) * el.gsfTrack()->dxy(PV.position())) < 0.0111)
	  //	  &&
	  (  fabs(el.gsfTrack()->dz( PV.position() )) < 0.0466)
	  &&
	  (fabs(ooEmooP) < 0.012)
	  &&
	  (relIsoWithDBeta < 10*0.0354)
	  &&
	  (el.passConversionVeto())
	  &&
	  (el.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) <= 2 )
	  )
	 pass = kTRUE;
     } 

  return pass;

}    

inline Bool_t passLooseElectronSelectionV4(const pat::Electron & el, const reco::Vertex &PV, const double &rho) {

  Bool_t pass = kFALSE;

  if (!el.chargeInfo().isGsfCtfScPixConsistent)
    return kFALSE;

     reco::GsfElectron::PflowIsolationVariables pfIso = el.pfIsolationVariables();

     //see here https://indico.cern.ch/event/369239/contribution/4/attachments/1134761/1623262/talk_effective_areas_25ns.pdf slide 12

     double EffectiveArea = 0;

     //float abseta = fabs(el.eta());
     float abseta = fabs(fabs(el.superCluster()->eta()));

     if (abseta >= 0.0 && abseta < 1.0 ) EffectiveArea = 0.1752;
     if (abseta >= 1.0 && abseta < 1.479 ) EffectiveArea = 0.1862;
     if (abseta >= 1.479 && abseta < 2.0 ) EffectiveArea = 0.1411;
     if (abseta >= 2.0 && abseta < 2.2 ) EffectiveArea = 0.1534;
     if (abseta >= 2.2 && abseta < 2.3 ) EffectiveArea = 0.1903;
     if (abseta >= 2.3 && abseta < 2.4 ) EffectiveArea = 0.2243;
     if (abseta >= 2.4 && abseta < 5.0 ) EffectiveArea = 0.2687;

     Float_t absiso = pfIso.sumChargedHadronPt + std::max(0.0 , pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - rho * EffectiveArea );
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

     //tight working point from here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
     if(fabs(el.superCluster()->eta()) < 2.5 && fabs(el.superCluster()->eta()) > 1.479 ){
       if(
	  (fabs(el.deltaEtaSuperClusterTrackAtVtx()) < 0.00724 )
	  &&
	  ( fabs(el.deltaPhiSuperClusterTrackAtVtx()) < 0.0918)
	  &&
	  (el.full5x5_sigmaIetaIeta() < 0.0279)
	  &&
	  ( el.hcalOverEcal() < 0.0646)
	  &&
	  //	  ( fabs((-1) * el.gsfTrack()->dxy(PV.position())) <  0.0351)
	  //	  &&
	  (  fabs(el.gsfTrack()->dz( PV.position() )) <  0.471)
	  &&
	  (fabs(ooEmooP) < 0.00999)
	  &&
	  (relIsoWithDBeta < 20*0.0646)
	  &&
	  (el.passConversionVeto())
	  &&
	  (el.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) <= 1 )
	  )
	 pass = kTRUE;
     } else if (fabs(el.superCluster()->eta()) < 1.479) {
       if(
	  (fabs(el.deltaEtaSuperClusterTrackAtVtx()) < 0.00926)
	  &&
	  ( fabs(el.deltaPhiSuperClusterTrackAtVtx()) < 0.0336)
	  &&
	  (el.full5x5_sigmaIetaIeta() < 0.0101)
	  &&
	  ( el.hcalOverEcal() < 0.0597)
	  &&
	  //	  ( fabs((-1) * el.gsfTrack()->dxy(PV.position())) < 0.0111)
	  //	  &&
	  (  fabs(el.gsfTrack()->dz( PV.position() )) < 0.0466)
	  &&
	  (fabs(ooEmooP) < 0.012)
	  &&
	  (relIsoWithDBeta < 20*0.0354)
	  &&
	  (el.passConversionVeto())
	  &&
	  (el.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) <= 2 )
	  )
	 pass = kTRUE;
     } 

  return pass;

}    

inline Bool_t passLooseElectronSelectionV5(const pat::Electron & el, const reco::Vertex &PV, const double &rho) {

  Bool_t pass = kFALSE;

  //  std::cout << "el.chargeInfo().isGsfCtfScPixConsistent = " << el.chargeInfo().isGsfCtfScPixConsistent << std::endl;

  if (!el.chargeInfo().isGsfCtfScPixConsistent)
    return kFALSE;

  //  std::cout << "fabs((-1) * el.gsfTrack()->dxy(PV.position())) = " << fabs((-1) * el.gsfTrack()->dxy(PV.position())) << std::endl;
  //  std::cout << "fabs(el.gsfTrack()->dz( PV.position() )) = " << fabs(el.gsfTrack()->dz( PV.position() )) << std::endl;
  //  std::cout << "el.passConversionVeto() = " << el.passConversionVeto() << std::endl;

  if(fabs(el.superCluster()->eta()) < 2.5 && fabs(el.superCluster()->eta()) > 1.479 ){
    if(
       ( fabs((-1) * el.gsfTrack()->dxy(PV.position())) <  0.5)
       &&
       (  fabs(el.gsfTrack()->dz( PV.position() )) <  1)
       &&
       (el.passConversionVeto())
       )
      pass = kTRUE;
  } else if (fabs(el.superCluster()->eta()) < 1.479) {
    if(
       ( fabs((-1) * el.gsfTrack()->dxy(PV.position())) <  0.5)
       &&
       (  fabs(el.gsfTrack()->dz( PV.position() )) <  1)
       &&
       (el.passConversionVeto())
       )
      pass = kTRUE;
  } 

  return pass;

}    


inline Bool_t passVeryLooseElectronSelection(const pat::Electron & el, const reco::Vertex &PV, const double &rho) {

  return passLooseElectronSelectionV5(el,PV,rho) || passLooseElectronSelectionV4(el,PV,rho) || passLooseElectronSelectionV2(el,PV,rho);;
  //return passTightElectronSelectionV3(el,PV,rho);

}    

inline Bool_t passVeryLooseMuonSelection(const pat::Muon & muon, const reco::Vertex &vtx ){


  //return passTightMuonSelectionV1(muon,vtx);

  return passLooseMuonSelectionV1(muon,vtx) || passLooseMuonSelectionV2(muon,vtx) || passLooseMuonSelectionV3(muon,vtx) || passLooseMuonSelectionV4(muon,vtx) || passLooseMuonSelectionV5(muon,vtx);

}


#endif

