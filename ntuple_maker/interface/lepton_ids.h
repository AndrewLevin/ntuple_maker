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



inline Bool_t passTightMuonIdV1(const pat::Muon & muon, const reco::Vertex &vtx) {

  if(!muon.isPFMuon() || !muon.isGlobalMuon() ) return false;

  bool muID = muon.isGlobalMuon() && muon.globalTrack()->normalizedChi2()<10. && muon.globalTrack()->hitPattern().numberOfValidMuonHits() >0 && (muon.numberOfMatchedStations() > 1);
  
  bool hits = muon.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 &&
    muon.innerTrack()->hitPattern().numberOfValidPixelHits() > 0; 

  //tighten the d0 cut from 0.2 to 0.02 because the muon POG ids consider muons from bs as real

  bool ip = fabs(muon.muonBestTrack()->dxy(vtx.position())) < 0.02 && fabs(muon.muonBestTrack()->dz(vtx.position())) < 0.1;

  return muID && hits && ip;

}

//should be the same as the member function is istightmuon()
inline Bool_t passTightMuonIdV2(const pat::Muon & muon, const reco::Vertex &vtx) {

  if(!muon.isPFMuon() || !muon.isGlobalMuon() ) return false;

  bool muID = muon.isGlobalMuon() && muon.globalTrack()->normalizedChi2()<10. && muon.globalTrack()->hitPattern().numberOfValidMuonHits() >0 && (muon.numberOfMatchedStations() > 1);
  
  bool hits = muon.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 &&
    muon.innerTrack()->hitPattern().numberOfValidPixelHits() > 0; 

  bool ip = fabs(muon.muonBestTrack()->dxy(vtx.position())) < 0.2 && fabs(muon.muonBestTrack()->dz(vtx.position())) < 0.1;
  
  return muID && hits && ip;

}


inline Bool_t passTightMuonSelectionV1(const pat::Muon & muon, const reco::Vertex &vtx) {

  Float_t relative_isolation = ( muon.pfIsolationR04().sumChargedHadronPt+ std::max(0.0,muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5 * muon.pfIsolationR04().sumPUPt) )/muon.pt();

  //  std::cout << "relative_isolation = " << relative_isolation << std::endl;

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


inline Bool_t passLooseMuonId(const pat::Muon & muon, const reco::Vertex &vtx) {

  if(!muon.isPFMuon()) return false;

  if(!muon.isTrackerMuon() && !muon.isGlobalMuon() ) return false;

  //  bool muID = muon.isGlobalMuon() && muon.globalTrack()->normalizedChi2()<10. && muon.globalTrack()->hitPattern().numberOfValidMuonHits() >0 && (muon.numberOfMatchedStations() > 1);
  
  //  bool hits = muon.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 &&
  //    muon.innerTrack()->hitPattern().numberOfValidPixelHits() > 0; 
  
    //  bool ip = fabs(muon.muonBestTrack()->dxy(vtx.position())) < 0.02 && fabs(muon.muonBestTrack()->dz(vtx.position())) < 0.5;

  return true;

}


inline Bool_t passLooseMuonSelectionV1(const pat::Muon & muon, const reco::Vertex &vtx) {

  //  Float_t relative_isolation = ( muon.pfIsolationR04().sumChargedHadronPt+ std::max(0.0,muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5 * muon.pfIsolationR04().sumPUPt) )/muon.pt();

  return passLooseMuonId(muon,vtx); // && relative_isolation < 1.0;

}

inline Bool_t passLooseMuonSelectionV2(const pat::Muon & muon, const reco::Vertex &vtx) {

  if(!muon.isGlobalMuon() ) return false;

  return fabs(muon.muonBestTrack()->dxy(vtx.position())) < 0.5 && fabs(muon.muonBestTrack()->dz(vtx.position())) < 1;

}

//Guillelmo's loose muon selection
inline Bool_t passLooseMuonSelectionV3(const pat::Muon & muon, const reco::Vertex &vtx) {

  Float_t relative_total_isolation = ( muon.pfIsolationR04().sumChargedHadronPt+ std::max(0.0,muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5 * muon.pfIsolationR04().sumPUPt) )/muon.pt();

  //from https://github.com/MiT-HEP/NeroProducer/blob/master/Nero/src/NeroLeptons.cpp#L94
  bool iso = ((relative_total_isolation < 0.4) && (muon.isolationR03().sumPt/muon.pt() < 0.4));

  //bool iso = (relative_total_isolation < 0.4);

  //bool iso = (muon.isolationR03().sumPt/muon.pt() < 0.4);

  return iso && passTightMuonIdV1(muon,vtx);

}

//Guillelmo's loose muon selection
inline Bool_t passLooseMuonSelectionV4(const pat::Muon & muon, const reco::Vertex &vtx) {

  if(!muon.isPFMuon() || !(muon.isGlobalMuon() || muon.isTrackerMuon()) ) return false;

  Float_t relative_total_isolation = ( muon.pfIsolationR04().sumChargedHadronPt+ std::max(0.0,muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5 * muon.pfIsolationR04().sumPUPt) )/muon.pt();

  //from https://github.com/MiT-HEP/NeroProducer/blob/master/Nero/src/NeroLeptons.cpp#L94
  //bool iso = ((relative_total_isolation < 0.4) && (muon.isolationR03().sumPt/muon.pt() < 0.4));

  //bool iso = (relative_total_isolation < 0.4);

  bool iso = (muon.isolationR03().sumPt/muon.pt() < 0.4);

  return iso;  

}

inline Bool_t passLooseMuonSelectionV5(const pat::Muon & muon, const reco::Vertex &vtx) {

  return fabs(muon.muonBestTrack()->dxy(vtx.position())) < 0.5 && fabs(muon.muonBestTrack()->dz(vtx.position())) < 1;

}




inline Bool_t passTightElectronSelectionV1(const pat::Electron & el, const reco::Vertex &PV, const double &rho) {

  Bool_t pass = kFALSE;

  if (!el.chargeInfo().isGsfCtfScPixConsistent)
    return kFALSE;

     reco::GsfElectron::PflowIsolationVariables pfIso = el.pfIsolationVariables();

     //see here https://indico.cern.ch/event/369239/contribution/4/attachments/1134761/1623262/talk_effective_areas_25ns.pdf slide 12

     double EffectiveArea = 0;

     float abseta = fabs(fabs(el.superCluster()->eta()));
     if (abseta >= 0.0 && abseta < 1.0 ) EffectiveArea = 0.1703;
     if (abseta >= 1.0 && abseta < 1.479 ) EffectiveArea = 0.1715;
     if (abseta >= 1.479 && abseta < 2.0 ) EffectiveArea = 0.1213;
     if (abseta >= 2.0 && abseta < 2.2 ) EffectiveArea = 0.1230;
     if (abseta >= 2.2 && abseta < 2.3 ) EffectiveArea = 0.1635;
     if (abseta >= 2.3 && abseta < 2.4 ) EffectiveArea = 0.1937;
     if (abseta >= 2.4 && abseta < 5.0 ) EffectiveArea = 0.2393;

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
	  (fabs(el.deltaEtaSuperClusterTrackAtVtx() - el.superCluster()->eta() + el.superCluster()->seed()->eta()) < 0.00609 )
	  &&
	  ( fabs(el.deltaPhiSuperClusterTrackAtVtx()) < 0.045)
	  &&
	  (el.full5x5_sigmaIetaIeta() < 0.0298)
	  &&
	  ( el.hcalOverEcal() < 0.0878)
	  &&
	  ( fabs((-1) * el.gsfTrack()->dxy(PV.position())) <  0.1)
	  &&
	  (  fabs(el.gsfTrack()->dz( PV.position() )) <  0.2)
	  &&
	  (fabs(ooEmooP) < 0.0298)
	  &&
	  (relIsoWithDBeta < 0.0821)
	  &&
	  (el.passConversionVeto())
	  &&
	  (el.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) <=  1)
	  )
	 pass = kTRUE;
     } else if (fabs(el.superCluster()->eta()) < 1.479) {
       if(
	  (fabs(el.deltaEtaSuperClusterTrackAtVtx() - el.superCluster()->eta() + el.superCluster()->seed()->eta()) < 0.00311)
	  &&
	  ( fabs(el.deltaPhiSuperClusterTrackAtVtx()) < 0.103)
	  &&
	  (el.full5x5_sigmaIetaIeta() < 0.00998)
	  &&
	  ( el.hcalOverEcal() < 0.253)
	  &&
	  ( fabs((-1) * el.gsfTrack()->dxy(PV.position())) < 0.05 )
	  &&
	  (  fabs(el.gsfTrack()->dz( PV.position() )) < 0.1 )
	  &&
	  (fabs(ooEmooP) < 0.134)
	  &&
	  (relIsoWithDBeta < 0.0695)
	  &&
	  (el.passConversionVeto())
	  &&
	  (el.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) <= 1  )
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
     if (abseta >= 0.0 && abseta < 1.0 ) EffectiveArea = 0.1703;
     if (abseta >= 1.0 && abseta < 1.479 ) EffectiveArea = 0.1715;
     if (abseta >= 1.479 && abseta < 2.0 ) EffectiveArea = 0.1213;
     if (abseta >= 2.0 && abseta < 2.2 ) EffectiveArea = 0.1230;
     if (abseta >= 2.2 && abseta < 2.3 ) EffectiveArea = 0.1635;
     if (abseta >= 2.3 && abseta < 2.4 ) EffectiveArea = 0.1937;
     if (abseta >= 2.4 && abseta < 5.0 ) EffectiveArea = 0.2393;

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
	  (fabs(el.deltaEtaSuperClusterTrackAtVtx() - el.superCluster()->eta() + el.superCluster()->seed()->eta()) < 0.00605)
	  &&
	  ( fabs(el.deltaPhiSuperClusterTrackAtVtx()) < 0.0394)
	  &&
	  (el.full5x5_sigmaIetaIeta() < 0.0292)
	  &&
	  ( el.hcalOverEcal() < 0.0641)
	  &&
	  ( fabs((-1) * el.gsfTrack()->dxy(PV.position())) < 0.1)
	  &&
	  (  fabs(el.gsfTrack()->dz( PV.position() )) < 0.2)
	  &&
	  (fabs(ooEmooP) < 0.0129)
	  &&
	  (relIsoWithDBeta < 0.0571)
	  &&
	  (el.passConversionVeto())
	  &&
	  (el.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) <= 1 )
	  )
	 pass = kTRUE;
     } else if (fabs(el.superCluster()->eta()) < 1.479) {
       if(
	  (fabs(el.deltaEtaSuperClusterTrackAtVtx() - el.superCluster()->eta() + el.superCluster()->seed()->eta()) <  0.00308)
	  &&
	  ( fabs(el.deltaPhiSuperClusterTrackAtVtx()) < 0.0816)
	  &&
	  (el.full5x5_sigmaIetaIeta() < 0.00998)
	  &&
	  ( el.hcalOverEcal() < 0.0414)
	  &&
	  ( fabs((-1) * el.gsfTrack()->dxy(PV.position())) < 0.05 )
	  &&
	  (  fabs(el.gsfTrack()->dz( PV.position() )) < 0.1 )
	  &&
	  (fabs(ooEmooP) < 0.0129)
	  &&
	  (relIsoWithDBeta < 0.0588)
	  &&
	  (el.passConversionVeto())
	  &&
	  (el.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) <= 1)
	  )
	 pass = kTRUE;
     } 


  return pass;

}    

inline Bool_t passTightElectronSelectionV3(const pat::Electron & el, const reco::Vertex &PV, const double &rho) {
  
  return el.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-tight");

}    

inline Bool_t passTightElectronSelectionV4(const pat::Electron & el, const reco::Vertex &PV, const double &rho) {

  Bool_t pass = kFALSE;

  if (!el.chargeInfo().isGsfCtfScPixConsistent)
    return kFALSE;


     reco::GsfElectron::PflowIsolationVariables pfIso = el.pfIsolationVariables();

     //see here https://indico.cern.ch/event/369239/contribution/4/attachments/1134761/1623262/talk_effective_areas_25ns.pdf slide 12

     double EffectiveArea = 0;

     //float abseta = fabs(el.eta());
     float abseta = fabs(fabs(el.superCluster()->eta()));
     if (abseta >= 0.0 && abseta < 1.0 ) EffectiveArea = 0.1703;
     if (abseta >= 1.0 && abseta < 1.479 ) EffectiveArea = 0.1715;
     if (abseta >= 1.479 && abseta < 2.0 ) EffectiveArea = 0.1213;
     if (abseta >= 2.0 && abseta < 2.2 ) EffectiveArea = 0.1230;
     if (abseta >= 2.2 && abseta < 2.3 ) EffectiveArea = 0.1635;
     if (abseta >= 2.3 && abseta < 2.4 ) EffectiveArea = 0.1937;
     if (abseta >= 2.4 && abseta < 5.0 ) EffectiveArea = 0.2393;

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
	  (fabs(el.deltaEtaSuperClusterTrackAtVtx() - el.superCluster()->eta() + el.superCluster()->seed()->eta()) < 0.00605)
	  &&
	  ( fabs(el.deltaPhiSuperClusterTrackAtVtx()) < 0.0394)
	  &&
	  (el.full5x5_sigmaIetaIeta() < 0.0292)
	  &&
	  ( el.hadronicOverEm() < 0.0641)
	  &&
	  //	  ( fabs((-1) * el.gsfTrack()->dxy(PV.position())) < 0.1)
	  //	  &&
	  //	  (  fabs(el.gsfTrack()->dz( PV.position() )) < 0.2)
	  //	  &&
	  (fabs(ooEmooP) < 0.0129)
	  &&
	  (relIsoWithDBeta < 0.0571)
	  &&
	  (el.passConversionVeto())
	  &&
	  (el.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) <= 0)
	  )
	 pass = kTRUE;
     } else if (fabs(el.superCluster()->eta()) < 1.479) {
       if(
	  (fabs(el.deltaEtaSuperClusterTrackAtVtx() - el.superCluster()->eta() + el.superCluster()->seed()->eta()) <  0.00308)
	  &&
	  ( fabs(el.deltaPhiSuperClusterTrackAtVtx()) < 0.0816)
	  &&
	  (el.full5x5_sigmaIetaIeta() < 0.00998)
	  &&
	  ( el.hadronicOverEm() < 0.0414)
	  &&
	  //	  ( fabs((-1) * el.gsfTrack()->dxy(PV.position())) < 0.05 )
	  //	  &&
	  //	  (  fabs(el.gsfTrack()->dz( PV.position() )) < 0.1 )
	  //	  &&
	  (fabs(ooEmooP) < 0.0129)
	  &&
	  (relIsoWithDBeta < 0.0588)
	  &&
	  (el.passConversionVeto())
	  &&
	  (el.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) <= 0) 
	  )
	 pass = kTRUE;
     } 

     //     std::cout << "el.gsfTrack()->numberOfLostHits () = " << el.gsfTrack()->numberOfLostHits () << std::endl;

     //     if (el.gsfTrack()->numberOfLostHits () != 0)
     //       pass = false;


     //     std::cout << "el.gsfTrack()->dz( PV.position() ) = " << el.gsfTrack()->dz( PV.position() ) << std::endl;

     //     std::cout << "el.superCluster()->eta() = " << el.superCluster()->eta() << std::endl;
     //     std::cout << "el.eta() = " << el.eta() << std::endl;
     //     std::cout << "el.pt() = " << el.pt() << std::endl;
     //     std::cout << "pass = " << pass << std::endl;

  return pass;

}    


//Guillelmo's electron selection
inline Bool_t passTightElectronSelectionV5(const pat::Electron & el, const reco::Vertex &PV, const double &rho) {

  Bool_t pass = kFALSE;

  if (!el.chargeInfo().isGsfCtfScPixConsistent)
    return kFALSE;


     reco::GsfElectron::PflowIsolationVariables pfIso = el.pfIsolationVariables();

     //see here https://indico.cern.ch/event/369239/contribution/4/attachments/1134761/1623262/talk_effective_areas_25ns.pdf slide 12

     double EffectiveArea = 0;

     //float abseta = fabs(el.eta());
     float abseta = fabs(fabs(el.superCluster()->eta()));
     if (abseta >= 0.0 && abseta < 1.0 ) EffectiveArea = 0.1703;
     if (abseta >= 1.0 && abseta < 1.479 ) EffectiveArea = 0.1715;
     if (abseta >= 1.479 && abseta < 2.0 ) EffectiveArea = 0.1213;
     if (abseta >= 2.0 && abseta < 2.2 ) EffectiveArea = 0.1230;
     if (abseta >= 2.2 && abseta < 2.3 ) EffectiveArea = 0.1635;
     if (abseta >= 2.3 && abseta < 2.4 ) EffectiveArea = 0.1937;
     if (abseta >= 2.4 && abseta < 5.0 ) EffectiveArea = 0.2393;

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
	  (fabs(el.deltaEtaSuperClusterTrackAtVtx() - el.superCluster()->eta() + el.superCluster()->seed()->eta()) < 0.00605)
	  &&
	  ( fabs(el.deltaPhiSuperClusterTrackAtVtx()) < 0.0394)
	  &&
	  (el.full5x5_sigmaIetaIeta() < 0.0292)
	  &&
	  ( el.hadronicOverEm() < 0.0641)
	  &&
	  (fabs(ooEmooP) < 0.0129)
	  &&
	  (relIsoWithDBeta < 0.0571)
	  &&
	  (el.passConversionVeto())
	  &&
	  (el.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) <= 1 )
	  )
	 pass = kTRUE;
     } else if (fabs(el.superCluster()->eta()) < 1.479) {
       if(
	  (fabs(el.deltaEtaSuperClusterTrackAtVtx() - el.superCluster()->eta() + el.superCluster()->seed()->eta()) <  0.00308)
	  &&
	  ( fabs(el.deltaPhiSuperClusterTrackAtVtx()) < 0.0816)
	  &&
	  (el.full5x5_sigmaIetaIeta() < 0.00998)
	  &&
	  ( el.hadronicOverEm() < 0.0414)
	  &&
	  (fabs(ooEmooP) < 0.0129)
	  &&
	  (relIsoWithDBeta < 0.0588)
	  &&
	  (el.passConversionVeto())
	  &&
	  (el.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) <= 1) 
	  )
	 pass = kTRUE;
     } 

  return pass;

}    

inline Bool_t passLooseElectronSelectionV1(const pat::Electron & el, const reco::Vertex &PV, const double &rho, const double &rhoHLTElectronSelection) {

  Bool_t pass = kFALSE;

  //  if (!el.chargeInfo().isGsfCtfScPixConsistent)
  //    return kFALSE;


     reco::GsfElectron::PflowIsolationVariables pfIso = el.pfIsolationVariables();

     //see here https://indico.cern.ch/event/369239/contribution/4/attachments/1134761/1623262/talk_effective_areas_25ns.pdf slide 12

     double EffectiveArea = 0;

     //float abseta = fabs(el.eta());
     float abseta = fabs(fabs(el.superCluster()->eta()));
     if (abseta >= 0.0 && abseta < 1.0 ) EffectiveArea = 0.1703;
     if (abseta >= 1.0 && abseta < 1.479 ) EffectiveArea = 0.1715;
     if (abseta >= 1.479 && abseta < 2.0 ) EffectiveArea = 0.1213;
     if (abseta >= 2.0 && abseta < 2.2 ) EffectiveArea = 0.1230;
     if (abseta >= 2.2 && abseta < 2.3 ) EffectiveArea = 0.1635;
     if (abseta >= 2.3 && abseta < 2.4 ) EffectiveArea = 0.1937;
     if (abseta >= 2.4 && abseta < 5.0 ) EffectiveArea = 0.2393;

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

     //loose working point from here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
     if(fabs(el.superCluster()->eta()) < 2.5 && fabs(el.superCluster()->eta()) > 1.479 ){
       if(
	  (fabs(el.deltaEtaSuperClusterTrackAtVtx() - el.superCluster()->eta() + el.superCluster()->seed()->eta()) < 0.00868 )
	  &&
	  ( fabs(el.deltaPhiSuperClusterTrackAtVtx()) < 0.213)
	  &&
	  (el.full5x5_sigmaIetaIeta() < 0.0314)
	  &&
	  ( el.hadronicOverEm() < 0.101)
	  &&
	  //	  ( fabs((-1) * el.gsfTrack()->dxy(PV.position())) < 0.05 )
	  //	  &&
	  //	  (  fabs(el.gsfTrack()->dz( PV.position() )) < 0.1 )
	  //	  &&
	  (fabs(ooEmooP) < 0.14)
	  &&
	  (relIsoWithDBeta < 0.107)
	  &&
	  (el.passConversionVeto())
	  &&
	  (el.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) <= 1) 
	  )
	 pass = kTRUE;


     } else if (fabs(el.superCluster()->eta()) < 1.479) {
       if(
	  (fabs(el.deltaEtaSuperClusterTrackAtVtx() - el.superCluster()->eta() + el.superCluster()->seed()->eta()) < 0.00477)
	  &&
	  ( fabs(el.deltaPhiSuperClusterTrackAtVtx()) < 0.222)
	  &&
	  (el.full5x5_sigmaIetaIeta() < 0.011)
	  &&
	  ( el.hadronicOverEm() < 0.298)
	  &&
	  //	  ( fabs((-1) * el.gsfTrack()->dxy(PV.position())) < 0.1)
	  //	  &&
	  //	  (  fabs(el.gsfTrack()->dz( PV.position() )) < 0.2)
	  //	  &&
	  (fabs(ooEmooP) < 0.241)
	  &&
	  (relIsoWithDBeta < 0.0994)
	  &&
	  (el.passConversionVeto())
	  &&
	  (el.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) <= 1)
	  )
	 pass = kTRUE;
     } 

  return pass;

}    

//Guillelmo's version of the HLT safe electron ID
inline Bool_t passLooseElectronSelectionV2(const pat::Electron & el, const reco::Vertex &PV, const double &rho, const double &rhoHLTElectronSelection) {

  Bool_t pass = kFALSE;

     //see here https://indico.cern.ch/event/369239/contribution/4/attachments/1134761/1623262/talk_effective_areas_25ns.pdf slide 12

     double EAecal = 0;
     double EAhcal = 0;

     //float abseta = fabs(el.eta());
     float abseta = fabs(fabs(el.superCluster()->eta()));
     if (abseta >= 0.0 && abseta < 1.479 ) EAecal = 0.165;
     if (abseta >= 1.479 && abseta < 5.0 ) EAecal = 0.132;
     if (abseta >= 0.0 && abseta < 1.479 ) EAhcal = 0.060;
     if (abseta >= 1.479 && abseta < 5.0 ) EAhcal = 0.131;

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
     
     std::cout << "el.full5x5_sigmaIetaIeta() = " << el.full5x5_sigmaIetaIeta() << std::endl;
     std::cout << "el.hadronicOverEm() = " << el.hadronicOverEm() << std::endl;
     std::cout << "fabs(ooEmooP) = " << fabs(ooEmooP) << std::endl;
     std::cout << "std::max(0.0, el.ecalPFClusterIso() - rhoHLTElectronSelection*EAecal)/el.pt() = " << std::max(0.0, el.ecalPFClusterIso() - rhoHLTElectronSelection*EAecal)/el.pt() << std::endl;
     std::cout << "std::max(0.0, el.hcalPFClusterIso() - rhoHLTElectronSelection*EAhcal)/el.pt() = " << std::max(0.0, el.hcalPFClusterIso() - rhoHLTElectronSelection*EAhcal)/el.pt() << std::endl;
     std::cout << "el.dr03TkSumPt()/el.pt() = " << el.dr03TkSumPt()/el.pt() << std::endl;
     std::cout << "std::abs(el.gsfTrack().isNonnull() ? el.gsfTrack()->normalizedChi2() : std::numeric_limits<float>::max()) = " << std::abs(el.gsfTrack().isNonnull() ? el.gsfTrack()->normalizedChi2() : std::numeric_limits<float>::max()) << std::endl;

     */

     //     std::cout << "el.electronID(\"cutBasedElectronID-Spring15-25ns-V1-standalone-loose\") = " << el.electronID("cutBasedElectronID-Spring15-25ns-V1-standalone-loose") << std::endl;
     
     if(fabs(el.superCluster()->eta()) < 2.5 && fabs(el.superCluster()->eta()) > 1.479 ){
       if(
	  (el.full5x5_sigmaIetaIeta() < 0.031)
	  &&
	  ( el.hadronicOverEm() < 0.065)
	  &&
	  (fabs(ooEmooP) < 0.013)
	  &&
	  (std::max(0.0, el.ecalPFClusterIso() - rhoHLTElectronSelection*EAecal)/el.pt() < 0.120)
	  &&
	  (std::max(0.0, el.hcalPFClusterIso() - rhoHLTElectronSelection*EAhcal)/el.pt() < 0.120)
	  &&
	  (el.dr03TkSumPt()/el.pt() < 0.08)
          &&
	  ((el.gsfTrack().isNonnull() ? el.gsfTrack()->normalizedChi2() : std::numeric_limits<float>::max()) < 3.0)
	  )
	 pass = kTRUE;
     } else if (fabs(el.superCluster()->eta()) < 1.479) {
       if(
	  (el.full5x5_sigmaIetaIeta() < 0.011)
	  &&
	  (fabs(el.deltaEtaSuperClusterTrackAtVtx() - el.superCluster()->eta() + el.superCluster()->seed()->eta()) <  0.004)
	  &&
	  ( fabs(el.deltaPhiSuperClusterTrackAtVtx()) < 0.020)
	  &&
	  ( el.hadronicOverEm() < 0.060)
	  &&
	  (fabs(ooEmooP) < 0.013)
	  &&
	  (std::max(0.0, el.ecalPFClusterIso() - rhoHLTElectronSelection*EAecal)/el.pt() < 0.160)
	  &&
	  (std::max(0.0, el.hcalPFClusterIso() - rhoHLTElectronSelection*EAhcal)/el.pt() < 0.120)
	  &&
	  (el.dr03TkSumPt()/el.pt() < 0.08)
	  )
	 pass = kTRUE;
     } 

  pass = pass && el.passConversionVeto();

  return pass;

}    

//double electron HLT safe electron selection + triple charge + impact parameter cuts
//see here https://twiki.cern.ch/twiki/bin/view/CMS/ChangesEGMHLTAlgo2014#Double_Electron_cuts and https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF#ID_IP_ISO_AN1 and https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
inline Bool_t passLooseElectronSelectionV3(const pat::Electron & el, const reco::Vertex &PV, const double &rho, const double &rhoHLTElectronSelection) {

  Bool_t pass = kFALSE;

  if (!el.chargeInfo().isGsfCtfScPixConsistent)
    return kFALSE;

     //see here https://indico.cern.ch/event/369239/contribution/4/attachments/1134761/1623262/talk_effective_areas_25ns.pdf slide 12

     double EAecal = 0;
     double EAhcal = 0;

     //float abseta = fabs(el.eta());
     float abseta = fabs(fabs(el.superCluster()->eta()));
     if (abseta >= 0.0 && abseta < 1.479 ) EAecal = 0.165;
     if (abseta >= 1.479 && abseta < 5.0 ) EAecal = 0.132;
     if (abseta >= 0.0 && abseta < 1.479 ) EAhcal = 0.060;
     if (abseta >= 1.479 && abseta < 5.0 ) EAhcal = 0.131;

     if(fabs(el.superCluster()->eta()) < 2.5 && fabs(el.superCluster()->eta()) > 1.479 ){
       if(
	  (el.full5x5_sigmaIetaIeta() < 0.035)
	  &&
	  ( el.hadronicOverEm() < 0.13)
	  &&
	  (el.ecalPFClusterIso()/el.pt() < 0.5)
	  &&
	  (el.hcalPFClusterIso()/el.pt() < 0.3)
	  &&
	  (el.dr03TkSumPt()/el.pt() < 0.2)
	  &&
	  (std::max(0.0, el.ecalPFClusterIso() - rhoHLTElectronSelection*EAecal)/el.pt() < 0.132)
	  &&
	  (std::max(0.0, el.hcalPFClusterIso() - rhoHLTElectronSelection*EAhcal)/el.pt() < 0.131)
	  )
	 pass = kTRUE;
     } else if (fabs(el.superCluster()->eta()) < 1.479) {
       if(
	  (el.full5x5_sigmaIetaIeta() < 0.013)
	  &&
	  (fabs(el.deltaEtaSuperClusterTrackAtVtx() - el.superCluster()->eta() + el.superCluster()->seed()->eta()) <  0.01)
	  &&
	  ( fabs(el.deltaPhiSuperClusterTrackAtVtx()) < 0.07)
	  &&
	  ( el.hadronicOverEm() < 0.13)
	  &&
	  (el.ecalPFClusterIso()/el.pt() < 0.5)
	  &&
	  (el.hcalPFClusterIso()/el.pt() < 0.3)
	  &&
	  (el.dr03TkSumPt()/el.pt() < 0.2)
	  &&
	  (std::max(0.0, el.ecalPFClusterIso() - rhoHLTElectronSelection*EAecal)/el.pt() < 0.165)
	  &&
	  (std::max(0.0, el.hcalPFClusterIso() - rhoHLTElectronSelection*EAhcal)/el.pt() < 0.06)
	  )
	 pass = kTRUE;
     } 

  pass = pass && el.passConversionVeto() && ( fabs((-1) * el.gsfTrack()->dxy(PV.position())) <  0.5) && (  fabs(el.gsfTrack()->dz( PV.position() )) <  1);

  return pass;

}    


//single electron HLT safe electron selection + triple charge + impact parameter cuts
//see here https://twiki.cern.ch/twiki/bin/view/CMS/ChangesEGMHLTAlgo2014#Double_Electron_cuts and https://twiki.cern.ch/twiki/bin/view/CMS/SUSLeptonSF#ID_IP_ISO_AN1 and https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
inline Bool_t passLooseElectronSelectionV4(const pat::Electron & el, const reco::Vertex &PV, const double &rho, const double &rhoHLTElectronSelection) {

  Bool_t pass = kFALSE;

  if (!el.chargeInfo().isGsfCtfScPixConsistent)
    return kFALSE;

     //see here https://indico.cern.ch/event/369239/contribution/4/attachments/1134761/1623262/talk_effective_areas_25ns.pdf slide 12

     double EAecal = 0;
     double EAhcal = 0;

     //float abseta = fabs(el.eta());
     float abseta = fabs(fabs(el.superCluster()->eta()));
     if (abseta >= 0.0 && abseta < 1.479 ) EAecal = 0.165;
     if (abseta >= 1.479 && abseta < 5.0 ) EAecal = 0.132;
     if (abseta >= 0.0 && abseta < 1.479 ) EAhcal = 0.060;
     if (abseta >= 1.479 && abseta < 5.0 ) EAhcal = 0.131;

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

     std::cout << "el.full5x5_sigmaIetaIeta() = " << el.full5x5_sigmaIetaIeta() << std::endl;
     std::cout << "el.hadronicOverEm() = " << el.hadronicOverEm() << std::endl;
     std::cout << "fabs(ooEmooP) = " << fabs(ooEmooP) << std::endl;
     std::cout << "std::max(0.0, el.ecalPFClusterIso() - rhoHLTElectronSelection*EAecal)/el.pt() = " << std::max(0.0, el.ecalPFClusterIso() - rhoHLTElectronSelection*EAecal)/el.pt() << std::endl;
     std::cout << "std::max(0.0, el.hcalPFClusterIso() - rhoHLTElectronSelection*EAhcal)/el.pt() = " << std::max(0.0, el.hcalPFClusterIso() - rhoHLTElectronSelection*EAhcal)/el.pt() << std::endl;
     std::cout << "el.dr03TkSumPt()/el.pt() = " << el.dr03TkSumPt()/el.pt() << std::endl;
     std::cout << "(el.gsfTrack().isNonnull() ? el.gsfTrack()->normalizedChi2() : std::numeric_limits<float>::max()) = " << (el.gsfTrack().isNonnull() ? el.gsfTrack()->normalizedChi2() : std::numeric_limits<float>::max()) << std::endl;

     */

     if(fabs(el.superCluster()->eta()) < 2.5 && fabs(el.superCluster()->eta()) > 1.479 ){
       if(
	  (el.full5x5_sigmaIetaIeta() < 0.031)
	  &&
	  ( el.hadronicOverEm() < 0.065)
	  &&
	  (fabs(ooEmooP) < 0.013)
	  &&
	  (std::max(0.0, el.ecalPFClusterIso() - rhoHLTElectronSelection*EAecal)/el.pt() < 0.120)
	  &&
	  (std::max(0.0, el.hcalPFClusterIso() - rhoHLTElectronSelection*EAhcal)/el.pt() < 0.120)
	  &&
	  (el.dr03TkSumPt()/el.pt() < 0.08)
          &&
	  ((el.gsfTrack().isNonnull() ? el.gsfTrack()->normalizedChi2() : std::numeric_limits<float>::max()) < 3.0)
	  )
	 pass = kTRUE;
     } else if (fabs(el.superCluster()->eta()) < 1.479) {
       if(
	  (el.full5x5_sigmaIetaIeta() < 0.011)
	  &&
	  (fabs(el.deltaEtaSuperClusterTrackAtVtx() - el.superCluster()->eta() + el.superCluster()->seed()->eta()) <  0.004)
	  &&
	  ( fabs(el.deltaPhiSuperClusterTrackAtVtx()) < 0.020)
	  &&
	  ( el.hadronicOverEm() < 0.060)
	  &&
	  (fabs(ooEmooP) < 0.013)
	  &&
	  (std::max(0.0, el.ecalPFClusterIso() - rhoHLTElectronSelection*EAecal)/el.pt() < 0.160)
	  &&
	  (std::max(0.0, el.hcalPFClusterIso() - rhoHLTElectronSelection*EAhcal)/el.pt() < 0.120)
	  &&
	  (el.dr03TkSumPt()/el.pt() < 0.08)
	  )
	 pass = kTRUE;
     } 

  pass = pass && el.passConversionVeto() && ( fabs((-1) * el.gsfTrack()->dxy(PV.position())) <  0.5) && (  fabs(el.gsfTrack()->dz( PV.position() )) <  1);

  return pass;

}    

inline Bool_t passLooseElectronSelectionV5(const pat::Electron & el, const reco::Vertex &PV, const double &rho, const double &rhoHLTElectronSelection) {

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


inline Bool_t passVeryLooseElectronSelection(const pat::Electron & el, const reco::Vertex &PV, const double &rho, const double &rhoHLTElectronSelection) {

  //  return passLooseElectronSelectionV5(el,PV,rho,rhoHLTElectronSelection) || passLooseElectronSelectionV4(el,PV,rho,rhoHLTElectronSelection) || passLooseElectronSelectionV3(el,PV,rho,rhoHLTElectronSelection) || passLooseElectronSelectionV2(el,PV,rho,rhoHLTElectronSelection) || passLooseElectronSelectionV1(el,PV,rho,rhoHLTElectronSelection);
  return passLooseElectronSelectionV2(el,PV,rho,rhoHLTElectronSelection) || passTightElectronSelectionV5(el,PV,rho);
  //return passTightElectronSelectionV3(el,PV,rho);

}    

inline Bool_t passVeryLooseMuonSelection(const pat::Muon & muon, const reco::Vertex &vtx ){


  //return passTightMuonSelectionV1(muon,vtx);

  //  return passLooseMuonSelectionV1(muon,vtx) || passLooseMuonSelectionV2(muon,vtx) || passLooseMuonSelectionV3(muon,vtx) || passLooseMuonSelectionV4(muon,vtx) || passLooseMuonSelectionV5(muon,vtx);
  return passLooseMuonSelectionV3(muon,vtx);

}

inline Bool_t passWLLJJVetoMuonId(const pat::Muon & muon, const reco::Vertex &vtx){

  //  Float_t relative_isolation = ( muon.pfIsolationR04().sumChargedHadronPt+ std::max(0.0,muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5 * muon.pfIsolationR04().sumPUPt) )/muon.pt();  
  //  Float_t relative_isolation = ( muon.pfIsolationR04().sumChargedHadronPt+ std::max(0.0,muon.pfIsolationR04().sumNeutralHadronEt  - 0.5 * muon.pfIsolationR04().sumPUPt) )/muon.pt();  

  if (muon.isSoftMuon(vtx) && fabs(muon.muonBestTrack()->dxy(vtx.position())) < 0.02 && fabs(muon.muonBestTrack()->dz(vtx.position())) < 0.1)
    return true;
  else
    return false;


}

inline Bool_t passWLLJJVetoElectronId(const pat::Electron & el, const reco::Vertex &PV, const float &rho, const double &rhoHLTElectronSelection){

     return passLooseElectronSelectionV2(el, PV, rho, rhoHLTElectronSelection) || passLooseElectronSelectionV1(el, PV, rho, rhoHLTElectronSelection);
}

#endif

