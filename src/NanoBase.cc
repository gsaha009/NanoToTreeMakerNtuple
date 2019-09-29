#include <iostream>
#include <iomanip>
#include <algorithm>
#include <iterator>
#include <functional>
#include <numeric>
#include <string>
#include <climits>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <sstream>

#include "TROOT.h"
#include "TSystem.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TClonesArray.h"
#include "TFile.h"

#include "NanoUtil.h"
#include "NanoBase.h"

using std::cout;
using std::cerr;
using std::endl;
using std::ios;
using std::string;
using std::vector;
using std::map;
using std::pair;
using std::setprecision;
using std::setw;
using std::setiosflags;
using std::resetiosflags;


//***C O N S T R U C T O R***//
#ifdef 0
NanoBase::NanoBase(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("NanoBaseAODProdMc_NANO.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("NanoBaseAODProdMc_NANO.root");
      }
      f->GetObject("Events",tree);
   }
   Init(tree);
}
#endif

NanoBase::NanoBase()
{
  fileList_.clear();
  hmap_.clear();
  if (!fChain) Init(fChain);
}

//***D E S T R C T O R***//
NanoBase::~NanoBase()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}



Int_t NanoBase::GetEntry(Long64_t entry)
{
// Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

Long64_t NanoBase::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void NanoBase::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("luminosityBlock", &luminosityBlock, &b_luminosityBlock);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("btagWeight_CSVV2", &btagWeight_CSVV2, &b_btagWeight_CSVV2);
   fChain->SetBranchAddress("btagWeight_CMVA", &btagWeight_CMVA, &b_btagWeight_CMVA);
   fChain->SetBranchAddress("CaloMET_phi", &CaloMET_phi, &b_CaloMET_phi);
   fChain->SetBranchAddress("CaloMET_pt", &CaloMET_pt, &b_CaloMET_pt);
   fChain->SetBranchAddress("CaloMET_sumEt", &CaloMET_sumEt, &b_CaloMET_sumEt);
   fChain->SetBranchAddress("nElectron", &nElectron, &b_nElectron);
   fChain->SetBranchAddress("Electron_deltaEtaSC", Electron_deltaEtaSC, &b_Electron_deltaEtaSC);
   fChain->SetBranchAddress("Electron_dr03EcalRecHitSumEt", Electron_dr03EcalRecHitSumEt, &b_Electron_dr03EcalRecHitSumEt);
   fChain->SetBranchAddress("Electron_dr03HcalDepth1TowerSumEt", Electron_dr03HcalDepth1TowerSumEt, &b_Electron_dr03HcalDepth1TowerSumEt);
   fChain->SetBranchAddress("Electron_dr03TkSumPt", Electron_dr03TkSumPt, &b_Electron_dr03TkSumPt);
   fChain->SetBranchAddress("Electron_dxy", Electron_dxy, &b_Electron_dxy);
   fChain->SetBranchAddress("Electron_dxyErr", Electron_dxyErr, &b_Electron_dxyErr);
   fChain->SetBranchAddress("Electron_dz", Electron_dz, &b_Electron_dz);
   fChain->SetBranchAddress("Electron_dzErr", Electron_dzErr, &b_Electron_dzErr);
   fChain->SetBranchAddress("Electron_eCorr", Electron_eCorr, &b_Electron_eCorr);
   fChain->SetBranchAddress("Electron_eInvMinusPInv", Electron_eInvMinusPInv, &b_Electron_eInvMinusPInv);
   fChain->SetBranchAddress("Electron_energyErr", Electron_energyErr, &b_Electron_energyErr);
   fChain->SetBranchAddress("Electron_eta", Electron_eta, &b_Electron_eta);
   fChain->SetBranchAddress("Electron_hoe", Electron_hoe, &b_Electron_hoe);
   fChain->SetBranchAddress("Electron_ip3d", Electron_ip3d, &b_Electron_ip3d);
   fChain->SetBranchAddress("Electron_mass", Electron_mass, &b_Electron_mass);
   fChain->SetBranchAddress("Electron_miniPFRelIso_all", Electron_miniPFRelIso_all, &b_Electron_miniPFRelIso_all);
   fChain->SetBranchAddress("Electron_miniPFRelIso_chg", Electron_miniPFRelIso_chg, &b_Electron_miniPFRelIso_chg);
   fChain->SetBranchAddress("Electron_mvaSpring16GP", Electron_mvaSpring16GP, &b_Electron_mvaSpring16GP);
   fChain->SetBranchAddress("Electron_mvaSpring16HZZ", Electron_mvaSpring16HZZ, &b_Electron_mvaSpring16HZZ);
   fChain->SetBranchAddress("Electron_pfRelIso03_all", Electron_pfRelIso03_all, &b_Electron_pfRelIso03_all);
   fChain->SetBranchAddress("Electron_pfRelIso03_chg", Electron_pfRelIso03_chg, &b_Electron_pfRelIso03_chg);
   fChain->SetBranchAddress("Electron_phi", Electron_phi, &b_Electron_phi);
   fChain->SetBranchAddress("Electron_pt", Electron_pt, &b_Electron_pt);
   fChain->SetBranchAddress("Electron_r9", Electron_r9, &b_Electron_r9);
   fChain->SetBranchAddress("Electron_sieie", Electron_sieie, &b_Electron_sieie);
   fChain->SetBranchAddress("Electron_sip3d", Electron_sip3d, &b_Electron_sip3d);
   fChain->SetBranchAddress("Electron_mvaTTH", Electron_mvaTTH, &b_Electron_mvaTTH);
   fChain->SetBranchAddress("Electron_charge", Electron_charge, &b_Electron_charge);
   fChain->SetBranchAddress("Electron_cutBased", Electron_cutBased, &b_Electron_cutBased);
   fChain->SetBranchAddress("Electron_cutBased_HLTPreSel", Electron_cutBased_HLTPreSel, &b_Electron_cutBased_HLTPreSel);
   fChain->SetBranchAddress("Electron_jetIdx", Electron_jetIdx, &b_Electron_jetIdx);
   fChain->SetBranchAddress("Electron_pdgId", Electron_pdgId, &b_Electron_pdgId);
   fChain->SetBranchAddress("Electron_photonIdx", Electron_photonIdx, &b_Electron_photonIdx);
   fChain->SetBranchAddress("Electron_tightCharge", Electron_tightCharge, &b_Electron_tightCharge);
   fChain->SetBranchAddress("Electron_vidNestedWPBitmap", Electron_vidNestedWPBitmap, &b_Electron_vidNestedWPBitmap);
   fChain->SetBranchAddress("Electron_convVeto", Electron_convVeto, &b_Electron_convVeto);
   fChain->SetBranchAddress("Electron_cutBased_HEEP", Electron_cutBased_HEEP, &b_Electron_cutBased_HEEP);
   fChain->SetBranchAddress("Electron_isPFcand", Electron_isPFcand, &b_Electron_isPFcand);
   fChain->SetBranchAddress("Electron_lostHits", Electron_lostHits, &b_Electron_lostHits);
   fChain->SetBranchAddress("Electron_mvaSpring16GP_WP80", Electron_mvaSpring16GP_WP80, &b_Electron_mvaSpring16GP_WP80);
   fChain->SetBranchAddress("Electron_mvaSpring16GP_WP90", Electron_mvaSpring16GP_WP90, &b_Electron_mvaSpring16GP_WP90);
   fChain->SetBranchAddress("Electron_mvaSpring16HZZ_WPL", Electron_mvaSpring16HZZ_WPL, &b_Electron_mvaSpring16HZZ_WPL);
   fChain->SetBranchAddress("Flag_BadChargedCandidateFilter", &Flag_BadChargedCandidateFilter, &b_Flag_BadChargedCandidateFilter);
   fChain->SetBranchAddress("Flag_BadGlobalMuon", &Flag_BadGlobalMuon, &b_Flag_BadGlobalMuon);
   fChain->SetBranchAddress("Flag_BadPFMuonFilter", &Flag_BadPFMuonFilter, &b_Flag_BadPFMuonFilter);
   fChain->SetBranchAddress("Flag_CloneGlobalMuon", &Flag_CloneGlobalMuon, &b_Flag_CloneGlobalMuon);
   fChain->SetBranchAddress("nFatJet", &nFatJet, &b_nFatJet);
   fChain->SetBranchAddress("FatJet_area", FatJet_area, &b_FatJet_area);
   fChain->SetBranchAddress("FatJet_btagCMVA", FatJet_btagCMVA, &b_FatJet_btagCMVA);
   fChain->SetBranchAddress("FatJet_btagCSVV2", FatJet_btagCSVV2, &b_FatJet_btagCSVV2);
   fChain->SetBranchAddress("FatJet_btagDeepB", FatJet_btagDeepB, &b_FatJet_btagDeepB);
   fChain->SetBranchAddress("FatJet_btagHbb", FatJet_btagHbb, &b_FatJet_btagHbb);
   fChain->SetBranchAddress("FatJet_eta", FatJet_eta, &b_FatJet_eta);
   fChain->SetBranchAddress("FatJet_mass", FatJet_mass, &b_FatJet_mass);
   fChain->SetBranchAddress("FatJet_msoftdrop", FatJet_msoftdrop, &b_FatJet_msoftdrop);
   fChain->SetBranchAddress("FatJet_msoftdrop_chs", FatJet_msoftdrop_chs, &b_FatJet_msoftdrop_chs);
   fChain->SetBranchAddress("FatJet_n2b1", FatJet_n2b1, &b_FatJet_n2b1);
   fChain->SetBranchAddress("FatJet_n3b1", FatJet_n3b1, &b_FatJet_n3b1);
   fChain->SetBranchAddress("FatJet_phi", FatJet_phi, &b_FatJet_phi);
   fChain->SetBranchAddress("FatJet_pt", FatJet_pt, &b_FatJet_pt);
   fChain->SetBranchAddress("FatJet_rawFactor", FatJet_rawFactor, &b_FatJet_rawFactor);
   fChain->SetBranchAddress("FatJet_tau1", FatJet_tau1, &b_FatJet_tau1);
   fChain->SetBranchAddress("FatJet_tau2", FatJet_tau2, &b_FatJet_tau2);
   fChain->SetBranchAddress("FatJet_tau3", FatJet_tau3, &b_FatJet_tau3);
   fChain->SetBranchAddress("FatJet_tau4", FatJet_tau4, &b_FatJet_tau4);
   fChain->SetBranchAddress("FatJet_jetId", FatJet_jetId, &b_FatJet_jetId);
   fChain->SetBranchAddress("FatJet_subJetIdx1", FatJet_subJetIdx1, &b_FatJet_subJetIdx1);
   fChain->SetBranchAddress("FatJet_subJetIdx2", FatJet_subJetIdx2, &b_FatJet_subJetIdx2);
   fChain->SetBranchAddress("nGenJetAK8", &nGenJetAK8, &b_nGenJetAK8);
   fChain->SetBranchAddress("GenJetAK8_eta", GenJetAK8_eta, &b_GenJetAK8_eta);
   fChain->SetBranchAddress("GenJetAK8_mass", GenJetAK8_mass, &b_GenJetAK8_mass);
   fChain->SetBranchAddress("GenJetAK8_phi", GenJetAK8_phi, &b_GenJetAK8_phi);
   fChain->SetBranchAddress("GenJetAK8_pt", GenJetAK8_pt, &b_GenJetAK8_pt);
   fChain->SetBranchAddress("nGenJet", &nGenJet, &b_nGenJet);
   fChain->SetBranchAddress("GenJet_eta", GenJet_eta, &b_GenJet_eta);
   fChain->SetBranchAddress("GenJet_mass", GenJet_mass, &b_GenJet_mass);
   fChain->SetBranchAddress("GenJet_phi", GenJet_phi, &b_GenJet_phi);
   fChain->SetBranchAddress("GenJet_pt", GenJet_pt, &b_GenJet_pt);
   fChain->SetBranchAddress("nGenPart", &nGenPart, &b_nGenPart);
   fChain->SetBranchAddress("GenPart_eta", GenPart_eta, &b_GenPart_eta);
   fChain->SetBranchAddress("GenPart_mass", GenPart_mass, &b_GenPart_mass);
   fChain->SetBranchAddress("GenPart_phi", GenPart_phi, &b_GenPart_phi);
   fChain->SetBranchAddress("GenPart_pt", GenPart_pt, &b_GenPart_pt);
   fChain->SetBranchAddress("GenPart_genPartIdxMother", GenPart_genPartIdxMother, &b_GenPart_genPartIdxMother);
   fChain->SetBranchAddress("GenPart_pdgId", GenPart_pdgId, &b_GenPart_pdgId);
   fChain->SetBranchAddress("GenPart_status", GenPart_status, &b_GenPart_status);
   fChain->SetBranchAddress("GenPart_statusFlags", GenPart_statusFlags, &b_GenPart_statusFlags);
   fChain->SetBranchAddress("Generator_binvar", &Generator_binvar, &b_Generator_binvar);
   fChain->SetBranchAddress("Generator_scalePDF", &Generator_scalePDF, &b_Generator_scalePDF);
   fChain->SetBranchAddress("Generator_weight", &Generator_weight, &b_Generator_weight);
   fChain->SetBranchAddress("Generator_x1", &Generator_x1, &b_Generator_x1);
   fChain->SetBranchAddress("Generator_x2", &Generator_x2, &b_Generator_x2);
   fChain->SetBranchAddress("Generator_xpdf1", &Generator_xpdf1, &b_Generator_xpdf1);
   fChain->SetBranchAddress("Generator_xpdf2", &Generator_xpdf2, &b_Generator_xpdf2);
   fChain->SetBranchAddress("Generator_id1", &Generator_id1, &b_Generator_id1);
   fChain->SetBranchAddress("Generator_id2", &Generator_id2, &b_Generator_id2);
   fChain->SetBranchAddress("nGenVisTau", &nGenVisTau, &b_nGenVisTau);
   fChain->SetBranchAddress("GenVisTau_eta", &GenVisTau_eta, &b_GenVisTau_eta);
   fChain->SetBranchAddress("GenVisTau_mass", &GenVisTau_mass, &b_GenVisTau_mass);
   fChain->SetBranchAddress("GenVisTau_phi", &GenVisTau_phi, &b_GenVisTau_phi);
   fChain->SetBranchAddress("GenVisTau_pt", &GenVisTau_pt, &b_GenVisTau_pt);
   fChain->SetBranchAddress("GenVisTau_charge", &GenVisTau_charge, &b_GenVisTau_charge);
   fChain->SetBranchAddress("GenVisTau_genPartIdxMother", &GenVisTau_genPartIdxMother, &b_GenVisTau_genPartIdxMother);
   fChain->SetBranchAddress("GenVisTau_status", &GenVisTau_status, &b_GenVisTau_status);
   fChain->SetBranchAddress("genWeight", &genWeight, &b_genWeight);
   fChain->SetBranchAddress("nPSWeight", &nPSWeight, &b_nPSWeight);
   fChain->SetBranchAddress("PSWeight", PSWeight, &b_PSWeight);
   fChain->SetBranchAddress("nJet", &nJet, &b_nJet);
   fChain->SetBranchAddress("Jet_area", Jet_area, &b_Jet_area);
   fChain->SetBranchAddress("Jet_btagCMVA", Jet_btagCMVA, &b_Jet_btagCMVA);
   fChain->SetBranchAddress("Jet_btagCSVV2", Jet_btagCSVV2, &b_Jet_btagCSVV2);
   fChain->SetBranchAddress("Jet_btagDeepB", Jet_btagDeepB, &b_Jet_btagDeepB);
   fChain->SetBranchAddress("Jet_btagDeepC", Jet_btagDeepC, &b_Jet_btagDeepC);
   fChain->SetBranchAddress("Jet_btagDeepFlavB", Jet_btagDeepFlavB, &b_Jet_btagDeepFlavB);
   fChain->SetBranchAddress("Jet_chEmEF", Jet_chEmEF, &b_Jet_chEmEF);
   fChain->SetBranchAddress("Jet_chHEF", Jet_chHEF, &b_Jet_chHEF);
   fChain->SetBranchAddress("Jet_eta", Jet_eta, &b_Jet_eta);
   fChain->SetBranchAddress("Jet_mass", Jet_mass, &b_Jet_mass);
   fChain->SetBranchAddress("Jet_muEF", Jet_muEF, &b_Jet_muEF);
   fChain->SetBranchAddress("Jet_neEmEF", Jet_neEmEF, &b_Jet_neEmEF);
   fChain->SetBranchAddress("Jet_neHEF", Jet_neHEF, &b_Jet_neHEF);
   fChain->SetBranchAddress("Jet_phi", Jet_phi, &b_Jet_phi);
   fChain->SetBranchAddress("Jet_pt", Jet_pt, &b_Jet_pt);
   fChain->SetBranchAddress("Jet_qgl", Jet_qgl, &b_Jet_qgl);
   fChain->SetBranchAddress("Jet_rawFactor", Jet_rawFactor, &b_Jet_rawFactor);
   fChain->SetBranchAddress("Jet_bRegCorr", Jet_bRegCorr, &b_Jet_bRegCorr);
   fChain->SetBranchAddress("Jet_bRegRes", Jet_bRegRes, &b_Jet_bRegRes);
   fChain->SetBranchAddress("Jet_electronIdx1", Jet_electronIdx1, &b_Jet_electronIdx1);
   fChain->SetBranchAddress("Jet_electronIdx2", Jet_electronIdx2, &b_Jet_electronIdx2);
   fChain->SetBranchAddress("Jet_jetId", Jet_jetId, &b_Jet_jetId);
   fChain->SetBranchAddress("Jet_muonIdx1", Jet_muonIdx1, &b_Jet_muonIdx1);
   fChain->SetBranchAddress("Jet_muonIdx2", Jet_muonIdx2, &b_Jet_muonIdx2);
   fChain->SetBranchAddress("Jet_nConstituents", Jet_nConstituents, &b_Jet_nConstituents);
   fChain->SetBranchAddress("Jet_nElectrons", Jet_nElectrons, &b_Jet_nElectrons);
   fChain->SetBranchAddress("Jet_nMuons", Jet_nMuons, &b_Jet_nMuons);
   fChain->SetBranchAddress("Jet_puId", Jet_puId, &b_Jet_puId);
   fChain->SetBranchAddress("GenMET_phi", &GenMET_phi, &b_GenMET_phi);
   fChain->SetBranchAddress("GenMET_pt", &GenMET_pt, &b_GenMET_pt);
   fChain->SetBranchAddress("MET_MetUnclustEnUpDeltaX", &MET_MetUnclustEnUpDeltaX, &b_MET_MetUnclustEnUpDeltaX);
   fChain->SetBranchAddress("MET_MetUnclustEnUpDeltaY", &MET_MetUnclustEnUpDeltaY, &b_MET_MetUnclustEnUpDeltaY);
   fChain->SetBranchAddress("MET_covXX", &MET_covXX, &b_MET_covXX);
   fChain->SetBranchAddress("MET_covXY", &MET_covXY, &b_MET_covXY);
   fChain->SetBranchAddress("MET_covYY", &MET_covYY, &b_MET_covYY);
   fChain->SetBranchAddress("MET_phi", &MET_phi, &b_MET_phi);
   fChain->SetBranchAddress("MET_pt", &MET_pt, &b_MET_pt);
   fChain->SetBranchAddress("MET_significance", &MET_significance, &b_MET_significance);
   fChain->SetBranchAddress("MET_sumEt", &MET_sumEt, &b_MET_sumEt);
   fChain->SetBranchAddress("nMuon", &nMuon, &b_nMuon);
   fChain->SetBranchAddress("Muon_dxy", Muon_dxy, &b_Muon_dxy);
   fChain->SetBranchAddress("Muon_dxyErr", Muon_dxyErr, &b_Muon_dxyErr);
   fChain->SetBranchAddress("Muon_dz", Muon_dz, &b_Muon_dz);
   fChain->SetBranchAddress("Muon_dzErr", Muon_dzErr, &b_Muon_dzErr);
   fChain->SetBranchAddress("Muon_eta", Muon_eta, &b_Muon_eta);
   fChain->SetBranchAddress("Muon_ip3d", Muon_ip3d, &b_Muon_ip3d);
   fChain->SetBranchAddress("Muon_mass", Muon_mass, &b_Muon_mass);
   fChain->SetBranchAddress("Muon_miniPFRelIso_all", Muon_miniPFRelIso_all, &b_Muon_miniPFRelIso_all);
   fChain->SetBranchAddress("Muon_miniPFRelIso_chg", Muon_miniPFRelIso_chg, &b_Muon_miniPFRelIso_chg);
   fChain->SetBranchAddress("Muon_pfRelIso03_all", Muon_pfRelIso03_all, &b_Muon_pfRelIso03_all);
   fChain->SetBranchAddress("Muon_pfRelIso03_chg", Muon_pfRelIso03_chg, &b_Muon_pfRelIso03_chg);
   fChain->SetBranchAddress("Muon_pfRelIso04_all", Muon_pfRelIso04_all, &b_Muon_pfRelIso04_all);
   fChain->SetBranchAddress("Muon_phi", Muon_phi, &b_Muon_phi);
   fChain->SetBranchAddress("Muon_pt", Muon_pt, &b_Muon_pt);
   fChain->SetBranchAddress("Muon_ptErr", Muon_ptErr, &b_Muon_ptErr);
   fChain->SetBranchAddress("Muon_segmentComp", Muon_segmentComp, &b_Muon_segmentComp);
   fChain->SetBranchAddress("Muon_sip3d", Muon_sip3d, &b_Muon_sip3d);
   fChain->SetBranchAddress("Muon_mvaTTH", Muon_mvaTTH, &b_Muon_mvaTTH);
   fChain->SetBranchAddress("Muon_charge", Muon_charge, &b_Muon_charge);
   fChain->SetBranchAddress("Muon_jetIdx", Muon_jetIdx, &b_Muon_jetIdx);
   fChain->SetBranchAddress("Muon_nStations", Muon_nStations, &b_Muon_nStations);
   fChain->SetBranchAddress("Muon_nTrackerLayers", Muon_nTrackerLayers, &b_Muon_nTrackerLayers);
   fChain->SetBranchAddress("Muon_pdgId", Muon_pdgId, &b_Muon_pdgId);
   fChain->SetBranchAddress("Muon_tightCharge", Muon_tightCharge, &b_Muon_tightCharge);
   fChain->SetBranchAddress("Muon_highPtId", Muon_highPtId, &b_Muon_highPtId);
   fChain->SetBranchAddress("Muon_isPFcand", Muon_isPFcand, &b_Muon_isPFcand);
   fChain->SetBranchAddress("Muon_mediumId", Muon_mediumId, &b_Muon_mediumId);
   fChain->SetBranchAddress("Muon_softId", Muon_softId, &b_Muon_softId);
   fChain->SetBranchAddress("Muon_tightId", Muon_tightId, &b_Muon_tightId);
   fChain->SetBranchAddress("nPhoton", &nPhoton, &b_nPhoton);
   fChain->SetBranchAddress("Photon_eCorr", Photon_eCorr, &b_Photon_eCorr);
   fChain->SetBranchAddress("Photon_energyErr", Photon_energyErr, &b_Photon_energyErr);
   fChain->SetBranchAddress("Photon_eta", Photon_eta, &b_Photon_eta);
   fChain->SetBranchAddress("Photon_hoe", Photon_hoe, &b_Photon_hoe);
   fChain->SetBranchAddress("Photon_mass", Photon_mass, &b_Photon_mass);
   fChain->SetBranchAddress("Photon_mvaID", Photon_mvaID, &b_Photon_mvaID);
   fChain->SetBranchAddress("Photon_pfRelIso03_all", Photon_pfRelIso03_all, &b_Photon_pfRelIso03_all);
   fChain->SetBranchAddress("Photon_pfRelIso03_chg", Photon_pfRelIso03_chg, &b_Photon_pfRelIso03_chg);
   fChain->SetBranchAddress("Photon_phi", Photon_phi, &b_Photon_phi);
   fChain->SetBranchAddress("Photon_pt", Photon_pt, &b_Photon_pt);
   fChain->SetBranchAddress("Photon_r9", Photon_r9, &b_Photon_r9);
   fChain->SetBranchAddress("Photon_sieie", Photon_sieie, &b_Photon_sieie);
   fChain->SetBranchAddress("Photon_charge", Photon_charge, &b_Photon_charge);
   fChain->SetBranchAddress("Photon_cutBased", Photon_cutBased, &b_Photon_cutBased);
   fChain->SetBranchAddress("Photon_electronIdx", Photon_electronIdx, &b_Photon_electronIdx);
   fChain->SetBranchAddress("Photon_jetIdx", Photon_jetIdx, &b_Photon_jetIdx);
   fChain->SetBranchAddress("Photon_pdgId", Photon_pdgId, &b_Photon_pdgId);
   fChain->SetBranchAddress("Photon_vidNestedWPBitmap", Photon_vidNestedWPBitmap, &b_Photon_vidNestedWPBitmap);
   fChain->SetBranchAddress("Photon_electronVeto", Photon_electronVeto, &b_Photon_electronVeto);
   fChain->SetBranchAddress("Photon_isScEtaEB", Photon_isScEtaEB, &b_Photon_isScEtaEB);
   fChain->SetBranchAddress("Photon_isScEtaEE", Photon_isScEtaEE, &b_Photon_isScEtaEE);
   fChain->SetBranchAddress("Photon_mvaID_WP80", Photon_mvaID_WP80, &b_Photon_mvaID_WP80);
   fChain->SetBranchAddress("Photon_mvaID_WP90", Photon_mvaID_WP90, &b_Photon_mvaID_WP90);
   fChain->SetBranchAddress("Photon_pixelSeed", Photon_pixelSeed, &b_Photon_pixelSeed);
   fChain->SetBranchAddress("Pileup_nTrueInt", &Pileup_nTrueInt, &b_Pileup_nTrueInt);
   fChain->SetBranchAddress("Pileup_nPU", &Pileup_nPU, &b_Pileup_nPU);
   fChain->SetBranchAddress("Pileup_sumEOOT", &Pileup_sumEOOT, &b_Pileup_sumEOOT);
   fChain->SetBranchAddress("Pileup_sumLOOT", &Pileup_sumLOOT, &b_Pileup_sumLOOT);
   fChain->SetBranchAddress("PuppiMET_phi", &PuppiMET_phi, &b_PuppiMET_phi);
   fChain->SetBranchAddress("PuppiMET_pt", &PuppiMET_pt, &b_PuppiMET_pt);
   fChain->SetBranchAddress("PuppiMET_sumEt", &PuppiMET_sumEt, &b_PuppiMET_sumEt);
   fChain->SetBranchAddress("RawMET_phi", &RawMET_phi, &b_RawMET_phi);
   fChain->SetBranchAddress("RawMET_pt", &RawMET_pt, &b_RawMET_pt);
   fChain->SetBranchAddress("RawMET_sumEt", &RawMET_sumEt, &b_RawMET_sumEt);
   fChain->SetBranchAddress("fixedGridRhoFastjetAll", &fixedGridRhoFastjetAll, &b_fixedGridRhoFastjetAll);
   fChain->SetBranchAddress("fixedGridRhoFastjetCentralCalo", &fixedGridRhoFastjetCentralCalo, &b_fixedGridRhoFastjetCentralCalo);
   fChain->SetBranchAddress("fixedGridRhoFastjetCentralNeutral", &fixedGridRhoFastjetCentralNeutral, &b_fixedGridRhoFastjetCentralNeutral);
   fChain->SetBranchAddress("nGenDressedLepton", &nGenDressedLepton, &b_nGenDressedLepton);
   fChain->SetBranchAddress("GenDressedLepton_eta", GenDressedLepton_eta, &b_GenDressedLepton_eta);
   fChain->SetBranchAddress("GenDressedLepton_mass", GenDressedLepton_mass, &b_GenDressedLepton_mass);
   fChain->SetBranchAddress("GenDressedLepton_phi", GenDressedLepton_phi, &b_GenDressedLepton_phi);
   fChain->SetBranchAddress("GenDressedLepton_pt", GenDressedLepton_pt, &b_GenDressedLepton_pt);
   fChain->SetBranchAddress("GenDressedLepton_pdgId", GenDressedLepton_pdgId, &b_GenDressedLepton_pdgId);
   fChain->SetBranchAddress("nSoftActivityJet", &nSoftActivityJet, &b_nSoftActivityJet);
   fChain->SetBranchAddress("SoftActivityJet_eta", SoftActivityJet_eta, &b_SoftActivityJet_eta);
   fChain->SetBranchAddress("SoftActivityJet_phi", SoftActivityJet_phi, &b_SoftActivityJet_phi);
   fChain->SetBranchAddress("SoftActivityJet_pt", SoftActivityJet_pt, &b_SoftActivityJet_pt);
   fChain->SetBranchAddress("SoftActivityJetHT", &SoftActivityJetHT, &b_SoftActivityJetHT);
   fChain->SetBranchAddress("SoftActivityJetHT10", &SoftActivityJetHT10, &b_SoftActivityJetHT10);
   fChain->SetBranchAddress("SoftActivityJetHT2", &SoftActivityJetHT2, &b_SoftActivityJetHT2);
   fChain->SetBranchAddress("SoftActivityJetHT5", &SoftActivityJetHT5, &b_SoftActivityJetHT5);
   fChain->SetBranchAddress("SoftActivityJetNjets10", &SoftActivityJetNjets10, &b_SoftActivityJetNjets10);
   fChain->SetBranchAddress("SoftActivityJetNjets2", &SoftActivityJetNjets2, &b_SoftActivityJetNjets2);
   fChain->SetBranchAddress("SoftActivityJetNjets5", &SoftActivityJetNjets5, &b_SoftActivityJetNjets5);
   fChain->SetBranchAddress("nSubJet", &nSubJet, &b_nSubJet);
   fChain->SetBranchAddress("SubJet_btagCMVA", SubJet_btagCMVA, &b_SubJet_btagCMVA);
   fChain->SetBranchAddress("SubJet_btagCSVV2", SubJet_btagCSVV2, &b_SubJet_btagCSVV2);
   fChain->SetBranchAddress("SubJet_btagDeepB", SubJet_btagDeepB, &b_SubJet_btagDeepB);
   fChain->SetBranchAddress("SubJet_eta", SubJet_eta, &b_SubJet_eta);
   fChain->SetBranchAddress("SubJet_mass", SubJet_mass, &b_SubJet_mass);
   fChain->SetBranchAddress("SubJet_n2b1", SubJet_n2b1, &b_SubJet_n2b1);
   fChain->SetBranchAddress("SubJet_n3b1", SubJet_n3b1, &b_SubJet_n3b1);
   fChain->SetBranchAddress("SubJet_phi", SubJet_phi, &b_SubJet_phi);
   fChain->SetBranchAddress("SubJet_pt", SubJet_pt, &b_SubJet_pt);
   fChain->SetBranchAddress("SubJet_tau1", SubJet_tau1, &b_SubJet_tau1);
   fChain->SetBranchAddress("SubJet_tau2", SubJet_tau2, &b_SubJet_tau2);
   fChain->SetBranchAddress("SubJet_tau3", SubJet_tau3, &b_SubJet_tau3);
   fChain->SetBranchAddress("SubJet_tau4", SubJet_tau4, &b_SubJet_tau4);
   fChain->SetBranchAddress("nTau", &nTau, &b_nTau);
   fChain->SetBranchAddress("Tau_chargedIso", Tau_chargedIso, &b_Tau_chargedIso);
   fChain->SetBranchAddress("Tau_dxy", Tau_dxy, &b_Tau_dxy);
   fChain->SetBranchAddress("Tau_dz", Tau_dz, &b_Tau_dz);
   fChain->SetBranchAddress("Tau_eta", Tau_eta, &b_Tau_eta);
   fChain->SetBranchAddress("Tau_leadTkDeltaEta", Tau_leadTkDeltaEta, &b_Tau_leadTkDeltaEta);
   fChain->SetBranchAddress("Tau_leadTkDeltaPhi", Tau_leadTkDeltaPhi, &b_Tau_leadTkDeltaPhi);
   fChain->SetBranchAddress("Tau_leadTkPtOverTauPt", Tau_leadTkPtOverTauPt, &b_Tau_leadTkPtOverTauPt);
   fChain->SetBranchAddress("Tau_mass", Tau_mass, &b_Tau_mass);
   fChain->SetBranchAddress("Tau_neutralIso", Tau_neutralIso, &b_Tau_neutralIso);
   fChain->SetBranchAddress("Tau_phi", Tau_phi, &b_Tau_phi);
   fChain->SetBranchAddress("Tau_photonsOutsideSignalCone", Tau_photonsOutsideSignalCone, &b_Tau_photonsOutsideSignalCone);
   fChain->SetBranchAddress("Tau_pt", Tau_pt, &b_Tau_pt);
   fChain->SetBranchAddress("Tau_puCorr", Tau_puCorr, &b_Tau_puCorr);
   fChain->SetBranchAddress("Tau_rawAntiEle", Tau_rawAntiEle, &b_Tau_rawAntiEle);
   fChain->SetBranchAddress("Tau_rawIso", Tau_rawIso, &b_Tau_rawIso);
   fChain->SetBranchAddress("Tau_rawIsodR03", Tau_rawIsodR03, &b_Tau_rawIsodR03);
   fChain->SetBranchAddress("Tau_rawMVAnewDM", Tau_rawMVAnewDM, &b_Tau_rawMVAnewDM);
   fChain->SetBranchAddress("Tau_rawMVAoldDM", Tau_rawMVAoldDM, &b_Tau_rawMVAoldDM);
   fChain->SetBranchAddress("Tau_rawMVAoldDMdR03", Tau_rawMVAoldDMdR03, &b_Tau_rawMVAoldDMdR03);
   fChain->SetBranchAddress("Tau_charge", Tau_charge, &b_Tau_charge);
   fChain->SetBranchAddress("Tau_decayMode", Tau_decayMode, &b_Tau_decayMode);
   fChain->SetBranchAddress("Tau_jetIdx", Tau_jetIdx, &b_Tau_jetIdx);
   fChain->SetBranchAddress("Tau_rawAntiEleCat", Tau_rawAntiEleCat, &b_Tau_rawAntiEleCat);
   fChain->SetBranchAddress("Tau_idAntiEle", Tau_idAntiEle, &b_Tau_idAntiEle);
   fChain->SetBranchAddress("Tau_idAntiMu", Tau_idAntiMu, &b_Tau_idAntiMu);
   fChain->SetBranchAddress("Tau_idDecayMode", Tau_idDecayMode, &b_Tau_idDecayMode);
   fChain->SetBranchAddress("Tau_idDecayModeNewDMs", Tau_idDecayModeNewDMs, &b_Tau_idDecayModeNewDMs);
   fChain->SetBranchAddress("Tau_idMVAnew", Tau_idMVAnew, &b_Tau_idMVAnew);
   fChain->SetBranchAddress("Tau_idMVAoldDM", Tau_idMVAoldDM, &b_Tau_idMVAoldDM);
   fChain->SetBranchAddress("Tau_idMVAoldDMdR03", Tau_idMVAoldDMdR03, &b_Tau_idMVAoldDMdR03);
   fChain->SetBranchAddress("TkMET_phi", &TkMET_phi, &b_TkMET_phi);
   fChain->SetBranchAddress("TkMET_pt", &TkMET_pt, &b_TkMET_pt);
   fChain->SetBranchAddress("TkMET_sumEt", &TkMET_sumEt, &b_TkMET_sumEt);
   fChain->SetBranchAddress("nTrigObj", &nTrigObj, &b_nTrigObj);
   fChain->SetBranchAddress("TrigObj_pt", TrigObj_pt, &b_TrigObj_pt);
   fChain->SetBranchAddress("TrigObj_eta", TrigObj_eta, &b_TrigObj_eta);
   fChain->SetBranchAddress("TrigObj_phi", TrigObj_phi, &b_TrigObj_phi);
   fChain->SetBranchAddress("TrigObj_l1pt", TrigObj_l1pt, &b_TrigObj_l1pt);
   fChain->SetBranchAddress("TrigObj_l1pt_2", TrigObj_l1pt_2, &b_TrigObj_l1pt_2);
   fChain->SetBranchAddress("TrigObj_l2pt", TrigObj_l2pt, &b_TrigObj_l2pt);
   fChain->SetBranchAddress("TrigObj_id", TrigObj_id, &b_TrigObj_id);
   fChain->SetBranchAddress("TrigObj_l1iso", TrigObj_l1iso, &b_TrigObj_l1iso);
   fChain->SetBranchAddress("TrigObj_l1charge", TrigObj_l1charge, &b_TrigObj_l1charge);
   fChain->SetBranchAddress("TrigObj_filterBits", TrigObj_filterBits, &b_TrigObj_filterBits);
   fChain->SetBranchAddress("genTtbarId", &genTtbarId, &b_genTtbarId);
   fChain->SetBranchAddress("nOtherPV", &nOtherPV, &b_nOtherPV);
   fChain->SetBranchAddress("OtherPV_z", OtherPV_z, &b_OtherPV_z);
   fChain->SetBranchAddress("PV_ndof", &PV_ndof, &b_PV_ndof);
   fChain->SetBranchAddress("PV_x", &PV_x, &b_PV_x);
   fChain->SetBranchAddress("PV_y", &PV_y, &b_PV_y);
   fChain->SetBranchAddress("PV_z", &PV_z, &b_PV_z);
   fChain->SetBranchAddress("PV_chi2", &PV_chi2, &b_PV_chi2);
   fChain->SetBranchAddress("PV_score", &PV_score, &b_PV_score);
   fChain->SetBranchAddress("PV_npvs", &PV_npvs, &b_PV_npvs);
   fChain->SetBranchAddress("PV_npvsGood", &PV_npvsGood, &b_PV_npvsGood);
   fChain->SetBranchAddress("nSV", &nSV, &b_nSV);
   fChain->SetBranchAddress("SV_dlen", SV_dlen, &b_SV_dlen);
   fChain->SetBranchAddress("SV_dlenSig", SV_dlenSig, &b_SV_dlenSig);
   fChain->SetBranchAddress("SV_pAngle", SV_pAngle, &b_SV_pAngle);
   fChain->SetBranchAddress("Electron_genPartIdx", Electron_genPartIdx, &b_Electron_genPartIdx);
   fChain->SetBranchAddress("Electron_genPartFlav", Electron_genPartFlav, &b_Electron_genPartFlav);
   fChain->SetBranchAddress("GenJetAK8_partonFlavour", GenJetAK8_partonFlavour, &b_GenJetAK8_partonFlavour);
   fChain->SetBranchAddress("GenJetAK8_hadronFlavour", GenJetAK8_hadronFlavour, &b_GenJetAK8_hadronFlavour);
   fChain->SetBranchAddress("GenJet_partonFlavour", GenJet_partonFlavour, &b_GenJet_partonFlavour);
   fChain->SetBranchAddress("GenJet_hadronFlavour", GenJet_hadronFlavour, &b_GenJet_hadronFlavour);
   fChain->SetBranchAddress("Jet_genJetIdx", Jet_genJetIdx, &b_Jet_genJetIdx);
   fChain->SetBranchAddress("Jet_hadronFlavour", Jet_hadronFlavour, &b_Jet_hadronFlavour);
   fChain->SetBranchAddress("Jet_partonFlavour", Jet_partonFlavour, &b_Jet_partonFlavour);
   fChain->SetBranchAddress("Muon_genPartIdx", Muon_genPartIdx, &b_Muon_genPartIdx);
   fChain->SetBranchAddress("Muon_genPartFlav", Muon_genPartFlav, &b_Muon_genPartFlav);
   fChain->SetBranchAddress("Photon_genPartIdx", Photon_genPartIdx, &b_Photon_genPartIdx);
   fChain->SetBranchAddress("Photon_genPartFlav", Photon_genPartFlav, &b_Photon_genPartFlav);
   fChain->SetBranchAddress("MET_fiducialGenPhi", &MET_fiducialGenPhi, &b_MET_fiducialGenPhi);
   fChain->SetBranchAddress("MET_fiducialGenPt", &MET_fiducialGenPt, &b_MET_fiducialGenPt);
   fChain->SetBranchAddress("Electron_cleanmask", Electron_cleanmask, &b_Electron_cleanmask);
   fChain->SetBranchAddress("Jet_cleanmask", Jet_cleanmask, &b_Jet_cleanmask);
   fChain->SetBranchAddress("Muon_cleanmask", Muon_cleanmask, &b_Muon_cleanmask);
   fChain->SetBranchAddress("Photon_cleanmask", Photon_cleanmask, &b_Photon_cleanmask);
   fChain->SetBranchAddress("Tau_cleanmask", Tau_cleanmask, &b_Tau_cleanmask);
   fChain->SetBranchAddress("SV_chi2", SV_chi2, &b_SV_chi2);
   fChain->SetBranchAddress("SV_eta", SV_eta, &b_SV_eta);
   fChain->SetBranchAddress("SV_mass", SV_mass, &b_SV_mass);
   fChain->SetBranchAddress("SV_ndof", SV_ndof, &b_SV_ndof);
   fChain->SetBranchAddress("SV_phi", SV_phi, &b_SV_phi);
   fChain->SetBranchAddress("SV_pt", SV_pt, &b_SV_pt);
   fChain->SetBranchAddress("SV_x", SV_x, &b_SV_x);
   fChain->SetBranchAddress("SV_y", SV_y, &b_SV_y);
   fChain->SetBranchAddress("SV_z", SV_z, &b_SV_z);
   fChain->SetBranchAddress("Tau_genPartIdx", Tau_genPartIdx, &b_Tau_genPartIdx);
   fChain->SetBranchAddress("Tau_genPartFlav", Tau_genPartFlav, &b_Tau_genPartFlav);
   fChain->SetBranchAddress("L1simulation_step", &L1simulation_step, &b_L1simulation_step);
   fChain->SetBranchAddress("HLTriggerFirstPath", &HLTriggerFirstPath, &b_HLTriggerFirstPath);
   fChain->SetBranchAddress("HLT_AK8PFJet360_TrimMass30", &HLT_AK8PFJet360_TrimMass30, &b_HLT_AK8PFJet360_TrimMass30);
   fChain->SetBranchAddress("HLT_AK8PFHT700_TrimR0p1PT0p03Mass50", &HLT_AK8PFHT700_TrimR0p1PT0p03Mass50, &b_HLT_AK8PFHT700_TrimR0p1PT0p03Mass50);
   fChain->SetBranchAddress("HLT_AK8PFHT650_TrimR0p1PT0p03Mass50", &HLT_AK8PFHT650_TrimR0p1PT0p03Mass50, &b_HLT_AK8PFHT650_TrimR0p1PT0p03Mass50);
   fChain->SetBranchAddress("HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20", &HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20, &b_HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20);
   fChain->SetBranchAddress("HLT_CaloJet500_NoJetID", &HLT_CaloJet500_NoJetID, &b_HLT_CaloJet500_NoJetID);
   fChain->SetBranchAddress("HLT_Dimuon13_PsiPrime", &HLT_Dimuon13_PsiPrime, &b_HLT_Dimuon13_PsiPrime);
   fChain->SetBranchAddress("HLT_Dimuon13_Upsilon", &HLT_Dimuon13_Upsilon, &b_HLT_Dimuon13_Upsilon);
   fChain->SetBranchAddress("HLT_Dimuon20_Jpsi", &HLT_Dimuon20_Jpsi, &b_HLT_Dimuon20_Jpsi);
   fChain->SetBranchAddress("HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf", &HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf, &b_HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf);
   fChain->SetBranchAddress("HLT_DoubleEle33_CaloIdL", &HLT_DoubleEle33_CaloIdL, &b_HLT_DoubleEle33_CaloIdL);
   fChain->SetBranchAddress("HLT_DoubleEle33_CaloIdL_MW", &HLT_DoubleEle33_CaloIdL_MW, &b_HLT_DoubleEle33_CaloIdL_MW);
   fChain->SetBranchAddress("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW", &HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW, &b_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW);
   fChain->SetBranchAddress("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL", &HLT_DoubleEle33_CaloIdL_GsfTrkIdVL, &b_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL);
   fChain->SetBranchAddress("HLT_DoubleEle37_Ele27_CaloIdL_GsfTrkIdVL", &HLT_DoubleEle37_Ele27_CaloIdL_GsfTrkIdVL, &b_HLT_DoubleEle37_Ele27_CaloIdL_GsfTrkIdVL);
   fChain->SetBranchAddress("HLT_DoubleMediumIsoPFTau32_Trk1_eta2p1_Reg", &HLT_DoubleMediumIsoPFTau32_Trk1_eta2p1_Reg, &b_HLT_DoubleMediumIsoPFTau32_Trk1_eta2p1_Reg);
   fChain->SetBranchAddress("HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg", &HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg, &b_HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg);
   fChain->SetBranchAddress("HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg", &HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg, &b_HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg);
   fChain->SetBranchAddress("HLT_DoubleMu33NoFiltersNoVtx", &HLT_DoubleMu33NoFiltersNoVtx, &b_HLT_DoubleMu33NoFiltersNoVtx);
   fChain->SetBranchAddress("HLT_DoubleMu38NoFiltersNoVtx", &HLT_DoubleMu38NoFiltersNoVtx, &b_HLT_DoubleMu38NoFiltersNoVtx);
   fChain->SetBranchAddress("HLT_DoubleMu23NoFiltersNoVtxDisplaced", &HLT_DoubleMu23NoFiltersNoVtxDisplaced, &b_HLT_DoubleMu23NoFiltersNoVtxDisplaced);
   fChain->SetBranchAddress("HLT_DoubleMu28NoFiltersNoVtxDisplaced", &HLT_DoubleMu28NoFiltersNoVtxDisplaced, &b_HLT_DoubleMu28NoFiltersNoVtxDisplaced);
   fChain->SetBranchAddress("HLT_DoubleMu4_3_Bs", &HLT_DoubleMu4_3_Bs, &b_HLT_DoubleMu4_3_Bs);
   fChain->SetBranchAddress("HLT_DoubleMu4_3_Jpsi_Displaced", &HLT_DoubleMu4_3_Jpsi_Displaced, &b_HLT_DoubleMu4_3_Jpsi_Displaced);
   fChain->SetBranchAddress("HLT_DoubleMu4_JpsiTrk_Displaced", &HLT_DoubleMu4_JpsiTrk_Displaced, &b_HLT_DoubleMu4_JpsiTrk_Displaced);
   fChain->SetBranchAddress("HLT_DoubleMu4_LowMassNonResonantTrk_Displaced", &HLT_DoubleMu4_LowMassNonResonantTrk_Displaced, &b_HLT_DoubleMu4_LowMassNonResonantTrk_Displaced);
   fChain->SetBranchAddress("HLT_DoubleMu3_Trk_Tau3mu", &HLT_DoubleMu3_Trk_Tau3mu, &b_HLT_DoubleMu3_Trk_Tau3mu);
   fChain->SetBranchAddress("HLT_DoubleMu4_PsiPrimeTrk_Displaced", &HLT_DoubleMu4_PsiPrimeTrk_Displaced, &b_HLT_DoubleMu4_PsiPrimeTrk_Displaced);
   fChain->SetBranchAddress("HLT_Mu7p5_L2Mu2_Jpsi", &HLT_Mu7p5_L2Mu2_Jpsi, &b_HLT_Mu7p5_L2Mu2_Jpsi);
   fChain->SetBranchAddress("HLT_Mu7p5_L2Mu2_Upsilon", &HLT_Mu7p5_L2Mu2_Upsilon, &b_HLT_Mu7p5_L2Mu2_Upsilon);
   fChain->SetBranchAddress("HLT_Mu7p5_Track2_Jpsi", &HLT_Mu7p5_Track2_Jpsi, &b_HLT_Mu7p5_Track2_Jpsi);
   fChain->SetBranchAddress("HLT_Mu7p5_Track3p5_Jpsi", &HLT_Mu7p5_Track3p5_Jpsi, &b_HLT_Mu7p5_Track3p5_Jpsi);
   fChain->SetBranchAddress("HLT_Mu7p5_Track7_Jpsi", &HLT_Mu7p5_Track7_Jpsi, &b_HLT_Mu7p5_Track7_Jpsi);
   fChain->SetBranchAddress("HLT_Mu7p5_Track2_Upsilon", &HLT_Mu7p5_Track2_Upsilon, &b_HLT_Mu7p5_Track2_Upsilon);
   fChain->SetBranchAddress("HLT_Mu7p5_Track3p5_Upsilon", &HLT_Mu7p5_Track3p5_Upsilon, &b_HLT_Mu7p5_Track3p5_Upsilon);
   fChain->SetBranchAddress("HLT_Mu7p5_Track7_Upsilon", &HLT_Mu7p5_Track7_Upsilon, &b_HLT_Mu7p5_Track7_Upsilon);
   fChain->SetBranchAddress("HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing", &HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing, &b_HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing);
   fChain->SetBranchAddress("HLT_Dimuon0er16_Jpsi_NoVertexing", &HLT_Dimuon0er16_Jpsi_NoVertexing, &b_HLT_Dimuon0er16_Jpsi_NoVertexing);
   fChain->SetBranchAddress("HLT_Dimuon6_Jpsi_NoVertexing", &HLT_Dimuon6_Jpsi_NoVertexing, &b_HLT_Dimuon6_Jpsi_NoVertexing);
   fChain->SetBranchAddress("HLT_DoublePhoton60", &HLT_DoublePhoton60, &b_HLT_DoublePhoton60);
   fChain->SetBranchAddress("HLT_DoublePhoton85", &HLT_DoublePhoton85, &b_HLT_DoublePhoton85);
   fChain->SetBranchAddress("HLT_Ele22_eta2p1_WPLoose_Gsf", &HLT_Ele22_eta2p1_WPLoose_Gsf, &b_HLT_Ele22_eta2p1_WPLoose_Gsf);
   fChain->SetBranchAddress("HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1", &HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1, &b_HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1);
   fChain->SetBranchAddress("HLT_Ele23_WPLoose_Gsf", &HLT_Ele23_WPLoose_Gsf, &b_HLT_Ele23_WPLoose_Gsf);
   fChain->SetBranchAddress("HLT_Ele23_WPLoose_Gsf_WHbbBoost", &HLT_Ele23_WPLoose_Gsf_WHbbBoost, &b_HLT_Ele23_WPLoose_Gsf_WHbbBoost);
   fChain->SetBranchAddress("HLT_Ele24_eta2p1_WPLoose_Gsf", &HLT_Ele24_eta2p1_WPLoose_Gsf, &b_HLT_Ele24_eta2p1_WPLoose_Gsf);
   fChain->SetBranchAddress("HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20", &HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20, &b_HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20);
   fChain->SetBranchAddress("HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1", &HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1, &b_HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1);
   fChain->SetBranchAddress("HLT_Ele25_WPTight_Gsf", &HLT_Ele25_WPTight_Gsf, &b_HLT_Ele25_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele25_eta2p1_WPLoose_Gsf", &HLT_Ele25_eta2p1_WPLoose_Gsf, &b_HLT_Ele25_eta2p1_WPLoose_Gsf);
   fChain->SetBranchAddress("HLT_Ele25_eta2p1_WPTight_Gsf", &HLT_Ele25_eta2p1_WPTight_Gsf, &b_HLT_Ele25_eta2p1_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele27_WPLoose_Gsf", &HLT_Ele27_WPLoose_Gsf, &b_HLT_Ele27_WPLoose_Gsf);
   fChain->SetBranchAddress("HLT_Ele27_WPLoose_Gsf_WHbbBoost", &HLT_Ele27_WPLoose_Gsf_WHbbBoost, &b_HLT_Ele27_WPLoose_Gsf_WHbbBoost);
   fChain->SetBranchAddress("HLT_Ele27_WPTight_Gsf", &HLT_Ele27_WPTight_Gsf, &b_HLT_Ele27_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele27_eta2p1_WPLoose_Gsf", &HLT_Ele27_eta2p1_WPLoose_Gsf, &b_HLT_Ele27_eta2p1_WPLoose_Gsf);
   fChain->SetBranchAddress("HLT_Ele27_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1", &HLT_Ele27_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1, &b_HLT_Ele27_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1);
   fChain->SetBranchAddress("HLT_Ele27_eta2p1_WPLoose_Gsf_DoubleMediumIsoPFTau32_Trk1_eta2p1_Reg", &HLT_Ele27_eta2p1_WPLoose_Gsf_DoubleMediumIsoPFTau32_Trk1_eta2p1_Reg, &b_HLT_Ele27_eta2p1_WPLoose_Gsf_DoubleMediumIsoPFTau32_Trk1_eta2p1_Reg);
   fChain->SetBranchAddress("HLT_Ele27_eta2p1_WPLoose_Gsf_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg", &HLT_Ele27_eta2p1_WPLoose_Gsf_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg, &b_HLT_Ele27_eta2p1_WPLoose_Gsf_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg);
   fChain->SetBranchAddress("HLT_Ele27_eta2p1_WPLoose_Gsf_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg", &HLT_Ele27_eta2p1_WPLoose_Gsf_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg, &b_HLT_Ele27_eta2p1_WPLoose_Gsf_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg);
   fChain->SetBranchAddress("HLT_Ele27_eta2p1_WPTight_Gsf", &HLT_Ele27_eta2p1_WPTight_Gsf, &b_HLT_Ele27_eta2p1_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele32_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1", &HLT_Ele32_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1, &b_HLT_Ele32_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1);
   fChain->SetBranchAddress("HLT_Ele32_eta2p1_WPTight_Gsf", &HLT_Ele32_eta2p1_WPTight_Gsf, &b_HLT_Ele32_eta2p1_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele35_WPLoose_Gsf", &HLT_Ele35_WPLoose_Gsf, &b_HLT_Ele35_WPLoose_Gsf);
   fChain->SetBranchAddress("HLT_Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50", &HLT_Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50, &b_HLT_Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50);
   fChain->SetBranchAddress("HLT_Ele45_WPLoose_Gsf", &HLT_Ele45_WPLoose_Gsf, &b_HLT_Ele45_WPLoose_Gsf);
   fChain->SetBranchAddress("HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50", &HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50, &b_HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50);
   fChain->SetBranchAddress("HLT_Ele105_CaloIdVT_GsfTrkIdT", &HLT_Ele105_CaloIdVT_GsfTrkIdT, &b_HLT_Ele105_CaloIdVT_GsfTrkIdT);
   fChain->SetBranchAddress("HLT_Ele30WP60_SC4_Mass55", &HLT_Ele30WP60_SC4_Mass55, &b_HLT_Ele30WP60_SC4_Mass55);
   fChain->SetBranchAddress("HLT_Ele30WP60_Ele8_Mass55", &HLT_Ele30WP60_Ele8_Mass55, &b_HLT_Ele30WP60_Ele8_Mass55);
   fChain->SetBranchAddress("HLT_HT200", &HLT_HT200, &b_HLT_HT200);
   fChain->SetBranchAddress("HLT_HT275", &HLT_HT275, &b_HLT_HT275);
   fChain->SetBranchAddress("HLT_HT325", &HLT_HT325, &b_HLT_HT325);
   fChain->SetBranchAddress("HLT_HT425", &HLT_HT425, &b_HLT_HT425);
   fChain->SetBranchAddress("HLT_HT575", &HLT_HT575, &b_HLT_HT575);
   fChain->SetBranchAddress("HLT_HT410to430", &HLT_HT410to430, &b_HLT_HT410to430);
   fChain->SetBranchAddress("HLT_HT430to450", &HLT_HT430to450, &b_HLT_HT430to450);
   fChain->SetBranchAddress("HLT_HT450to470", &HLT_HT450to470, &b_HLT_HT450to470);
   fChain->SetBranchAddress("HLT_HT470to500", &HLT_HT470to500, &b_HLT_HT470to500);
   fChain->SetBranchAddress("HLT_HT500to550", &HLT_HT500to550, &b_HLT_HT500to550);
   fChain->SetBranchAddress("HLT_HT550to650", &HLT_HT550to650, &b_HLT_HT550to650);
   fChain->SetBranchAddress("HLT_HT650", &HLT_HT650, &b_HLT_HT650);
   fChain->SetBranchAddress("HLT_Mu16_eta2p1_MET30", &HLT_Mu16_eta2p1_MET30, &b_HLT_Mu16_eta2p1_MET30);
   fChain->SetBranchAddress("HLT_IsoMu16_eta2p1_MET30", &HLT_IsoMu16_eta2p1_MET30, &b_HLT_IsoMu16_eta2p1_MET30);
   fChain->SetBranchAddress("HLT_IsoMu16_eta2p1_MET30_LooseIsoPFTau50_Trk30_eta2p1", &HLT_IsoMu16_eta2p1_MET30_LooseIsoPFTau50_Trk30_eta2p1, &b_HLT_IsoMu16_eta2p1_MET30_LooseIsoPFTau50_Trk30_eta2p1);
   fChain->SetBranchAddress("HLT_IsoMu17_eta2p1_LooseIsoPFTau20", &HLT_IsoMu17_eta2p1_LooseIsoPFTau20, &b_HLT_IsoMu17_eta2p1_LooseIsoPFTau20);
   fChain->SetBranchAddress("HLT_IsoMu17_eta2p1_LooseIsoPFTau20_SingleL1", &HLT_IsoMu17_eta2p1_LooseIsoPFTau20_SingleL1, &b_HLT_IsoMu17_eta2p1_LooseIsoPFTau20_SingleL1);
   fChain->SetBranchAddress("HLT_DoubleIsoMu17_eta2p1", &HLT_DoubleIsoMu17_eta2p1, &b_HLT_DoubleIsoMu17_eta2p1);
   fChain->SetBranchAddress("HLT_DoubleIsoMu17_eta2p1_noDzCut", &HLT_DoubleIsoMu17_eta2p1_noDzCut, &b_HLT_DoubleIsoMu17_eta2p1_noDzCut);
   fChain->SetBranchAddress("HLT_IsoMu18", &HLT_IsoMu18, &b_HLT_IsoMu18);
   fChain->SetBranchAddress("HLT_IsoMu19_eta2p1_LooseIsoPFTau20", &HLT_IsoMu19_eta2p1_LooseIsoPFTau20, &b_HLT_IsoMu19_eta2p1_LooseIsoPFTau20);
   fChain->SetBranchAddress("HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1", &HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1, &b_HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1);
   fChain->SetBranchAddress("HLT_IsoMu19_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg", &HLT_IsoMu19_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg, &b_HLT_IsoMu19_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg);
   fChain->SetBranchAddress("HLT_IsoMu20", &HLT_IsoMu20, &b_HLT_IsoMu20);
   fChain->SetBranchAddress("HLT_IsoMu21_eta2p1_LooseIsoPFTau20_SingleL1", &HLT_IsoMu21_eta2p1_LooseIsoPFTau20_SingleL1, &b_HLT_IsoMu21_eta2p1_LooseIsoPFTau20_SingleL1);
   fChain->SetBranchAddress("HLT_IsoMu21_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg", &HLT_IsoMu21_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg, &b_HLT_IsoMu21_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg);
   fChain->SetBranchAddress("HLT_IsoMu22", &HLT_IsoMu22, &b_HLT_IsoMu22);
   fChain->SetBranchAddress("HLT_IsoMu22_eta2p1", &HLT_IsoMu22_eta2p1, &b_HLT_IsoMu22_eta2p1);
   fChain->SetBranchAddress("HLT_IsoMu24", &HLT_IsoMu24, &b_HLT_IsoMu24);
   fChain->SetBranchAddress("HLT_IsoMu27", &HLT_IsoMu27, &b_HLT_IsoMu27);
   fChain->SetBranchAddress("HLT_IsoTkMu18", &HLT_IsoTkMu18, &b_HLT_IsoTkMu18);
   fChain->SetBranchAddress("HLT_IsoTkMu20", &HLT_IsoTkMu20, &b_HLT_IsoTkMu20);
   fChain->SetBranchAddress("HLT_IsoTkMu22", &HLT_IsoTkMu22, &b_HLT_IsoTkMu22);
   fChain->SetBranchAddress("HLT_IsoTkMu22_eta2p1", &HLT_IsoTkMu22_eta2p1, &b_HLT_IsoTkMu22_eta2p1);
   fChain->SetBranchAddress("HLT_IsoTkMu24", &HLT_IsoTkMu24, &b_HLT_IsoTkMu24);
   fChain->SetBranchAddress("HLT_IsoTkMu27", &HLT_IsoTkMu27, &b_HLT_IsoTkMu27);
   fChain->SetBranchAddress("HLT_JetE30_NoBPTX3BX", &HLT_JetE30_NoBPTX3BX, &b_HLT_JetE30_NoBPTX3BX);
   fChain->SetBranchAddress("HLT_JetE30_NoBPTX", &HLT_JetE30_NoBPTX, &b_HLT_JetE30_NoBPTX);
   fChain->SetBranchAddress("HLT_JetE50_NoBPTX3BX", &HLT_JetE50_NoBPTX3BX, &b_HLT_JetE50_NoBPTX3BX);
   fChain->SetBranchAddress("HLT_JetE70_NoBPTX3BX", &HLT_JetE70_NoBPTX3BX, &b_HLT_JetE70_NoBPTX3BX);
   fChain->SetBranchAddress("HLT_L1SingleMu18", &HLT_L1SingleMu18, &b_HLT_L1SingleMu18);
   fChain->SetBranchAddress("HLT_L2Mu10", &HLT_L2Mu10, &b_HLT_L2Mu10);
   fChain->SetBranchAddress("HLT_L1SingleMuOpen", &HLT_L1SingleMuOpen, &b_HLT_L1SingleMuOpen);
   fChain->SetBranchAddress("HLT_L1SingleMuOpen_DT", &HLT_L1SingleMuOpen_DT, &b_HLT_L1SingleMuOpen_DT);
   fChain->SetBranchAddress("HLT_L2DoubleMu23_NoVertex", &HLT_L2DoubleMu23_NoVertex, &b_HLT_L2DoubleMu23_NoVertex);
   fChain->SetBranchAddress("HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10", &HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10, &b_HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10);
   fChain->SetBranchAddress("HLT_L2DoubleMu38_NoVertex_2Cha_Angle2p5_Mass10", &HLT_L2DoubleMu38_NoVertex_2Cha_Angle2p5_Mass10, &b_HLT_L2DoubleMu38_NoVertex_2Cha_Angle2p5_Mass10);
   fChain->SetBranchAddress("HLT_L2Mu10_NoVertex_NoBPTX3BX", &HLT_L2Mu10_NoVertex_NoBPTX3BX, &b_HLT_L2Mu10_NoVertex_NoBPTX3BX);
   fChain->SetBranchAddress("HLT_L2Mu10_NoVertex_NoBPTX", &HLT_L2Mu10_NoVertex_NoBPTX, &b_HLT_L2Mu10_NoVertex_NoBPTX);
   fChain->SetBranchAddress("HLT_L2Mu35_NoVertex_3Sta_NoBPTX3BX", &HLT_L2Mu35_NoVertex_3Sta_NoBPTX3BX, &b_HLT_L2Mu35_NoVertex_3Sta_NoBPTX3BX);
   fChain->SetBranchAddress("HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX", &HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX, &b_HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX);
   fChain->SetBranchAddress("HLT_LooseIsoPFTau50_Trk30_eta2p1", &HLT_LooseIsoPFTau50_Trk30_eta2p1, &b_HLT_LooseIsoPFTau50_Trk30_eta2p1);
   fChain->SetBranchAddress("HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80", &HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80, &b_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80);
   fChain->SetBranchAddress("HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90", &HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90, &b_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90);
   fChain->SetBranchAddress("HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110", &HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110, &b_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110);
   fChain->SetBranchAddress("HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120", &HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120, &b_HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120);
   fChain->SetBranchAddress("HLT_PFTau120_eta2p1", &HLT_PFTau120_eta2p1, &b_HLT_PFTau120_eta2p1);
   fChain->SetBranchAddress("HLT_VLooseIsoPFTau120_Trk50_eta2p1", &HLT_VLooseIsoPFTau120_Trk50_eta2p1, &b_HLT_VLooseIsoPFTau120_Trk50_eta2p1);
   fChain->SetBranchAddress("HLT_VLooseIsoPFTau140_Trk50_eta2p1", &HLT_VLooseIsoPFTau140_Trk50_eta2p1, &b_HLT_VLooseIsoPFTau140_Trk50_eta2p1);
   fChain->SetBranchAddress("HLT_Mu17_Mu8", &HLT_Mu17_Mu8, &b_HLT_Mu17_Mu8);
   fChain->SetBranchAddress("HLT_Mu17_Mu8_DZ", &HLT_Mu17_Mu8_DZ, &b_HLT_Mu17_Mu8_DZ);
   fChain->SetBranchAddress("HLT_Mu17_Mu8_SameSign", &HLT_Mu17_Mu8_SameSign, &b_HLT_Mu17_Mu8_SameSign);
   fChain->SetBranchAddress("HLT_Mu17_Mu8_SameSign_DZ", &HLT_Mu17_Mu8_SameSign_DZ, &b_HLT_Mu17_Mu8_SameSign_DZ);
   fChain->SetBranchAddress("HLT_Mu20_Mu10", &HLT_Mu20_Mu10, &b_HLT_Mu20_Mu10);
   fChain->SetBranchAddress("HLT_Mu20_Mu10_DZ", &HLT_Mu20_Mu10_DZ, &b_HLT_Mu20_Mu10_DZ);
   fChain->SetBranchAddress("HLT_Mu20_Mu10_SameSign", &HLT_Mu20_Mu10_SameSign, &b_HLT_Mu20_Mu10_SameSign);
   fChain->SetBranchAddress("HLT_Mu20_Mu10_SameSign_DZ", &HLT_Mu20_Mu10_SameSign_DZ, &b_HLT_Mu20_Mu10_SameSign_DZ);
   fChain->SetBranchAddress("HLT_Mu17_TkMu8_DZ", &HLT_Mu17_TkMu8_DZ, &b_HLT_Mu17_TkMu8_DZ);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ", &HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ, &b_HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL", &HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL, &b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ", &HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ, &b_HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ);
   fChain->SetBranchAddress("HLT_Mu25_TkMu0_dEta18_Onia", &HLT_Mu25_TkMu0_dEta18_Onia, &b_HLT_Mu25_TkMu0_dEta18_Onia);
   fChain->SetBranchAddress("HLT_Mu27_TkMu8", &HLT_Mu27_TkMu8, &b_HLT_Mu27_TkMu8);
   fChain->SetBranchAddress("HLT_Mu30_TkMu11", &HLT_Mu30_TkMu11, &b_HLT_Mu30_TkMu11);
   fChain->SetBranchAddress("HLT_Mu30_eta2p1_PFJet150_PFJet50", &HLT_Mu30_eta2p1_PFJet150_PFJet50, &b_HLT_Mu30_eta2p1_PFJet150_PFJet50);
   fChain->SetBranchAddress("HLT_Mu40_TkMu11", &HLT_Mu40_TkMu11, &b_HLT_Mu40_TkMu11);
   fChain->SetBranchAddress("HLT_Mu40_eta2p1_PFJet200_PFJet50", &HLT_Mu40_eta2p1_PFJet200_PFJet50, &b_HLT_Mu40_eta2p1_PFJet200_PFJet50);
   fChain->SetBranchAddress("HLT_Mu20", &HLT_Mu20, &b_HLT_Mu20);
   fChain->SetBranchAddress("HLT_TkMu20", &HLT_TkMu20, &b_HLT_TkMu20);
   fChain->SetBranchAddress("HLT_Mu24_eta2p1", &HLT_Mu24_eta2p1, &b_HLT_Mu24_eta2p1);
   fChain->SetBranchAddress("HLT_TkMu24_eta2p1", &HLT_TkMu24_eta2p1, &b_HLT_TkMu24_eta2p1);
   fChain->SetBranchAddress("HLT_Mu27", &HLT_Mu27, &b_HLT_Mu27);
   fChain->SetBranchAddress("HLT_TkMu27", &HLT_TkMu27, &b_HLT_TkMu27);
   fChain->SetBranchAddress("HLT_Mu45_eta2p1", &HLT_Mu45_eta2p1, &b_HLT_Mu45_eta2p1);
   fChain->SetBranchAddress("HLT_Mu50", &HLT_Mu50, &b_HLT_Mu50);
   fChain->SetBranchAddress("HLT_TkMu50", &HLT_TkMu50, &b_HLT_TkMu50);
   fChain->SetBranchAddress("HLT_Mu38NoFiltersNoVtx_Photon38_CaloIdL", &HLT_Mu38NoFiltersNoVtx_Photon38_CaloIdL, &b_HLT_Mu38NoFiltersNoVtx_Photon38_CaloIdL);
   fChain->SetBranchAddress("HLT_Mu42NoFiltersNoVtx_Photon42_CaloIdL", &HLT_Mu42NoFiltersNoVtx_Photon42_CaloIdL, &b_HLT_Mu42NoFiltersNoVtx_Photon42_CaloIdL);
   fChain->SetBranchAddress("HLT_Mu28NoFiltersNoVtxDisplaced_Photon28_CaloIdL", &HLT_Mu28NoFiltersNoVtxDisplaced_Photon28_CaloIdL, &b_HLT_Mu28NoFiltersNoVtxDisplaced_Photon28_CaloIdL);
   fChain->SetBranchAddress("HLT_Mu33NoFiltersNoVtxDisplaced_Photon33_CaloIdL", &HLT_Mu33NoFiltersNoVtxDisplaced_Photon33_CaloIdL, &b_HLT_Mu33NoFiltersNoVtxDisplaced_Photon33_CaloIdL);
   fChain->SetBranchAddress("HLT_Mu23NoFiltersNoVtx_Photon23_CaloIdL", &HLT_Mu23NoFiltersNoVtx_Photon23_CaloIdL, &b_HLT_Mu23NoFiltersNoVtx_Photon23_CaloIdL);
   fChain->SetBranchAddress("HLT_DoubleMu18NoFiltersNoVtx", &HLT_DoubleMu18NoFiltersNoVtx, &b_HLT_DoubleMu18NoFiltersNoVtx);
   fChain->SetBranchAddress("HLT_Mu33NoFiltersNoVtxDisplaced_DisplacedJet50_Tight", &HLT_Mu33NoFiltersNoVtxDisplaced_DisplacedJet50_Tight, &b_HLT_Mu33NoFiltersNoVtxDisplaced_DisplacedJet50_Tight);
   fChain->SetBranchAddress("HLT_Mu33NoFiltersNoVtxDisplaced_DisplacedJet50_Loose", &HLT_Mu33NoFiltersNoVtxDisplaced_DisplacedJet50_Loose, &b_HLT_Mu33NoFiltersNoVtxDisplaced_DisplacedJet50_Loose);
   fChain->SetBranchAddress("HLT_Mu28NoFiltersNoVtx_DisplacedJet40_Loose", &HLT_Mu28NoFiltersNoVtx_DisplacedJet40_Loose, &b_HLT_Mu28NoFiltersNoVtx_DisplacedJet40_Loose);
   fChain->SetBranchAddress("HLT_Mu38NoFiltersNoVtxDisplaced_DisplacedJet60_Tight", &HLT_Mu38NoFiltersNoVtxDisplaced_DisplacedJet60_Tight, &b_HLT_Mu38NoFiltersNoVtxDisplaced_DisplacedJet60_Tight);
   fChain->SetBranchAddress("HLT_Mu38NoFiltersNoVtxDisplaced_DisplacedJet60_Loose", &HLT_Mu38NoFiltersNoVtxDisplaced_DisplacedJet60_Loose, &b_HLT_Mu38NoFiltersNoVtxDisplaced_DisplacedJet60_Loose);
   fChain->SetBranchAddress("HLT_Mu38NoFiltersNoVtx_DisplacedJet60_Loose", &HLT_Mu38NoFiltersNoVtx_DisplacedJet60_Loose, &b_HLT_Mu38NoFiltersNoVtx_DisplacedJet60_Loose);
   fChain->SetBranchAddress("HLT_Mu28NoFiltersNoVtx_CentralCaloJet40", &HLT_Mu28NoFiltersNoVtx_CentralCaloJet40, &b_HLT_Mu28NoFiltersNoVtx_CentralCaloJet40);
   fChain->SetBranchAddress("HLT_PFHT300_PFMET100", &HLT_PFHT300_PFMET100, &b_HLT_PFHT300_PFMET100);
   fChain->SetBranchAddress("HLT_PFHT300_PFMET110", &HLT_PFHT300_PFMET110, &b_HLT_PFHT300_PFMET110);
   fChain->SetBranchAddress("HLT_PFHT550_4JetPt50", &HLT_PFHT550_4JetPt50, &b_HLT_PFHT550_4JetPt50);
   fChain->SetBranchAddress("HLT_PFHT650_4JetPt50", &HLT_PFHT650_4JetPt50, &b_HLT_PFHT650_4JetPt50);
   fChain->SetBranchAddress("HLT_PFHT750_4JetPt50", &HLT_PFHT750_4JetPt50, &b_HLT_PFHT750_4JetPt50);
   fChain->SetBranchAddress("HLT_PFJet15_NoCaloMatched", &HLT_PFJet15_NoCaloMatched, &b_HLT_PFJet15_NoCaloMatched);
   fChain->SetBranchAddress("HLT_PFJet25_NoCaloMatched", &HLT_PFJet25_NoCaloMatched, &b_HLT_PFJet25_NoCaloMatched);
   fChain->SetBranchAddress("HLT_DiPFJet15_NoCaloMatched", &HLT_DiPFJet15_NoCaloMatched, &b_HLT_DiPFJet15_NoCaloMatched);
   fChain->SetBranchAddress("HLT_DiPFJet25_NoCaloMatched", &HLT_DiPFJet25_NoCaloMatched, &b_HLT_DiPFJet25_NoCaloMatched);
   fChain->SetBranchAddress("HLT_DiPFJet15_FBEta3_NoCaloMatched", &HLT_DiPFJet15_FBEta3_NoCaloMatched, &b_HLT_DiPFJet15_FBEta3_NoCaloMatched);
   fChain->SetBranchAddress("HLT_DiPFJet25_FBEta3_NoCaloMatched", &HLT_DiPFJet25_FBEta3_NoCaloMatched, &b_HLT_DiPFJet25_FBEta3_NoCaloMatched);
   fChain->SetBranchAddress("HLT_DiPFJetAve15_HFJEC", &HLT_DiPFJetAve15_HFJEC, &b_HLT_DiPFJetAve15_HFJEC);
   fChain->SetBranchAddress("HLT_DiPFJetAve25_HFJEC", &HLT_DiPFJetAve25_HFJEC, &b_HLT_DiPFJetAve25_HFJEC);
   fChain->SetBranchAddress("HLT_DiPFJetAve35_HFJEC", &HLT_DiPFJetAve35_HFJEC, &b_HLT_DiPFJetAve35_HFJEC);
   fChain->SetBranchAddress("HLT_AK8PFJet40", &HLT_AK8PFJet40, &b_HLT_AK8PFJet40);
   fChain->SetBranchAddress("HLT_AK8PFJet60", &HLT_AK8PFJet60, &b_HLT_AK8PFJet60);
   fChain->SetBranchAddress("HLT_AK8PFJet80", &HLT_AK8PFJet80, &b_HLT_AK8PFJet80);
   fChain->SetBranchAddress("HLT_AK8PFJet140", &HLT_AK8PFJet140, &b_HLT_AK8PFJet140);
   fChain->SetBranchAddress("HLT_AK8PFJet200", &HLT_AK8PFJet200, &b_HLT_AK8PFJet200);
   fChain->SetBranchAddress("HLT_AK8PFJet260", &HLT_AK8PFJet260, &b_HLT_AK8PFJet260);
   fChain->SetBranchAddress("HLT_AK8PFJet320", &HLT_AK8PFJet320, &b_HLT_AK8PFJet320);
   fChain->SetBranchAddress("HLT_AK8PFJet400", &HLT_AK8PFJet400, &b_HLT_AK8PFJet400);
   fChain->SetBranchAddress("HLT_AK8PFJet450", &HLT_AK8PFJet450, &b_HLT_AK8PFJet450);
   fChain->SetBranchAddress("HLT_AK8PFJet500", &HLT_AK8PFJet500, &b_HLT_AK8PFJet500);
   fChain->SetBranchAddress("HLT_PFJet40", &HLT_PFJet40, &b_HLT_PFJet40);
   fChain->SetBranchAddress("HLT_PFJet60", &HLT_PFJet60, &b_HLT_PFJet60);
   fChain->SetBranchAddress("HLT_PFJet80", &HLT_PFJet80, &b_HLT_PFJet80);
   fChain->SetBranchAddress("HLT_PFJet140", &HLT_PFJet140, &b_HLT_PFJet140);
   fChain->SetBranchAddress("HLT_PFJet200", &HLT_PFJet200, &b_HLT_PFJet200);
   fChain->SetBranchAddress("HLT_PFJet260", &HLT_PFJet260, &b_HLT_PFJet260);
   fChain->SetBranchAddress("HLT_PFJet320", &HLT_PFJet320, &b_HLT_PFJet320);
   fChain->SetBranchAddress("HLT_PFJet400", &HLT_PFJet400, &b_HLT_PFJet400);
   fChain->SetBranchAddress("HLT_PFJet450", &HLT_PFJet450, &b_HLT_PFJet450);
   fChain->SetBranchAddress("HLT_PFJet500", &HLT_PFJet500, &b_HLT_PFJet500);
   fChain->SetBranchAddress("HLT_DiPFJetAve40", &HLT_DiPFJetAve40, &b_HLT_DiPFJetAve40);
   fChain->SetBranchAddress("HLT_DiPFJetAve60", &HLT_DiPFJetAve60, &b_HLT_DiPFJetAve60);
   fChain->SetBranchAddress("HLT_DiPFJetAve80", &HLT_DiPFJetAve80, &b_HLT_DiPFJetAve80);
   fChain->SetBranchAddress("HLT_DiPFJetAve140", &HLT_DiPFJetAve140, &b_HLT_DiPFJetAve140);
   fChain->SetBranchAddress("HLT_DiPFJetAve200", &HLT_DiPFJetAve200, &b_HLT_DiPFJetAve200);
   fChain->SetBranchAddress("HLT_DiPFJetAve260", &HLT_DiPFJetAve260, &b_HLT_DiPFJetAve260);
   fChain->SetBranchAddress("HLT_DiPFJetAve320", &HLT_DiPFJetAve320, &b_HLT_DiPFJetAve320);
   fChain->SetBranchAddress("HLT_DiPFJetAve400", &HLT_DiPFJetAve400, &b_HLT_DiPFJetAve400);
   fChain->SetBranchAddress("HLT_DiPFJetAve500", &HLT_DiPFJetAve500, &b_HLT_DiPFJetAve500);
   fChain->SetBranchAddress("HLT_DiPFJetAve60_HFJEC", &HLT_DiPFJetAve60_HFJEC, &b_HLT_DiPFJetAve60_HFJEC);
   fChain->SetBranchAddress("HLT_DiPFJetAve80_HFJEC", &HLT_DiPFJetAve80_HFJEC, &b_HLT_DiPFJetAve80_HFJEC);
   fChain->SetBranchAddress("HLT_DiPFJetAve100_HFJEC", &HLT_DiPFJetAve100_HFJEC, &b_HLT_DiPFJetAve100_HFJEC);
   fChain->SetBranchAddress("HLT_DiPFJetAve160_HFJEC", &HLT_DiPFJetAve160_HFJEC, &b_HLT_DiPFJetAve160_HFJEC);
   fChain->SetBranchAddress("HLT_DiPFJetAve220_HFJEC", &HLT_DiPFJetAve220_HFJEC, &b_HLT_DiPFJetAve220_HFJEC);
   fChain->SetBranchAddress("HLT_DiPFJetAve300_HFJEC", &HLT_DiPFJetAve300_HFJEC, &b_HLT_DiPFJetAve300_HFJEC);
   fChain->SetBranchAddress("HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu140", &HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu140, &b_HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu140);
   fChain->SetBranchAddress("HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu80", &HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu80, &b_HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu80);
   fChain->SetBranchAddress("HLT_DiCentralPFJet55_PFMET110", &HLT_DiCentralPFJet55_PFMET110, &b_HLT_DiCentralPFJet55_PFMET110);
   fChain->SetBranchAddress("HLT_DiCentralPFJet170", &HLT_DiCentralPFJet170, &b_HLT_DiCentralPFJet170);
   fChain->SetBranchAddress("HLT_SingleCentralPFJet170_CFMax0p1", &HLT_SingleCentralPFJet170_CFMax0p1, &b_HLT_SingleCentralPFJet170_CFMax0p1);
   fChain->SetBranchAddress("HLT_DiCentralPFJet170_CFMax0p1", &HLT_DiCentralPFJet170_CFMax0p1, &b_HLT_DiCentralPFJet170_CFMax0p1);
   fChain->SetBranchAddress("HLT_DiCentralPFJet220_CFMax0p3", &HLT_DiCentralPFJet220_CFMax0p3, &b_HLT_DiCentralPFJet220_CFMax0p3);
   fChain->SetBranchAddress("HLT_DiCentralPFJet330_CFMax0p5", &HLT_DiCentralPFJet330_CFMax0p5, &b_HLT_DiCentralPFJet330_CFMax0p5);
   fChain->SetBranchAddress("HLT_DiCentralPFJet430", &HLT_DiCentralPFJet430, &b_HLT_DiCentralPFJet430);
   fChain->SetBranchAddress("HLT_PFHT125", &HLT_PFHT125, &b_HLT_PFHT125);
   fChain->SetBranchAddress("HLT_PFHT200", &HLT_PFHT200, &b_HLT_PFHT200);
   fChain->SetBranchAddress("HLT_PFHT250", &HLT_PFHT250, &b_HLT_PFHT250);
   fChain->SetBranchAddress("HLT_PFHT300", &HLT_PFHT300, &b_HLT_PFHT300);
   fChain->SetBranchAddress("HLT_PFHT350", &HLT_PFHT350, &b_HLT_PFHT350);
   fChain->SetBranchAddress("HLT_PFHT400", &HLT_PFHT400, &b_HLT_PFHT400);
   fChain->SetBranchAddress("HLT_PFHT475", &HLT_PFHT475, &b_HLT_PFHT475);
   fChain->SetBranchAddress("HLT_PFHT600", &HLT_PFHT600, &b_HLT_PFHT600);
   fChain->SetBranchAddress("HLT_PFHT650", &HLT_PFHT650, &b_HLT_PFHT650);
   fChain->SetBranchAddress("HLT_PFHT800", &HLT_PFHT800, &b_HLT_PFHT800);
   fChain->SetBranchAddress("HLT_PFHT900", &HLT_PFHT900, &b_HLT_PFHT900);
   fChain->SetBranchAddress("HLT_PFHT200_PFAlphaT0p51", &HLT_PFHT200_PFAlphaT0p51, &b_HLT_PFHT200_PFAlphaT0p51);
   fChain->SetBranchAddress("HLT_PFHT200_DiPFJetAve90_PFAlphaT0p57", &HLT_PFHT200_DiPFJetAve90_PFAlphaT0p57, &b_HLT_PFHT200_DiPFJetAve90_PFAlphaT0p57);
   fChain->SetBranchAddress("HLT_PFHT200_DiPFJetAve90_PFAlphaT0p63", &HLT_PFHT200_DiPFJetAve90_PFAlphaT0p63, &b_HLT_PFHT200_DiPFJetAve90_PFAlphaT0p63);
   fChain->SetBranchAddress("HLT_PFHT250_DiPFJetAve90_PFAlphaT0p55", &HLT_PFHT250_DiPFJetAve90_PFAlphaT0p55, &b_HLT_PFHT250_DiPFJetAve90_PFAlphaT0p55);
   fChain->SetBranchAddress("HLT_PFHT250_DiPFJetAve90_PFAlphaT0p58", &HLT_PFHT250_DiPFJetAve90_PFAlphaT0p58, &b_HLT_PFHT250_DiPFJetAve90_PFAlphaT0p58);
   fChain->SetBranchAddress("HLT_PFHT300_DiPFJetAve90_PFAlphaT0p53", &HLT_PFHT300_DiPFJetAve90_PFAlphaT0p53, &b_HLT_PFHT300_DiPFJetAve90_PFAlphaT0p53);
   fChain->SetBranchAddress("HLT_PFHT300_DiPFJetAve90_PFAlphaT0p54", &HLT_PFHT300_DiPFJetAve90_PFAlphaT0p54, &b_HLT_PFHT300_DiPFJetAve90_PFAlphaT0p54);
   fChain->SetBranchAddress("HLT_PFHT350_DiPFJetAve90_PFAlphaT0p52", &HLT_PFHT350_DiPFJetAve90_PFAlphaT0p52, &b_HLT_PFHT350_DiPFJetAve90_PFAlphaT0p52);
   fChain->SetBranchAddress("HLT_PFHT350_DiPFJetAve90_PFAlphaT0p53", &HLT_PFHT350_DiPFJetAve90_PFAlphaT0p53, &b_HLT_PFHT350_DiPFJetAve90_PFAlphaT0p53);
   fChain->SetBranchAddress("HLT_PFHT400_DiPFJetAve90_PFAlphaT0p51", &HLT_PFHT400_DiPFJetAve90_PFAlphaT0p51, &b_HLT_PFHT400_DiPFJetAve90_PFAlphaT0p51);
   fChain->SetBranchAddress("HLT_PFHT400_DiPFJetAve90_PFAlphaT0p52", &HLT_PFHT400_DiPFJetAve90_PFAlphaT0p52, &b_HLT_PFHT400_DiPFJetAve90_PFAlphaT0p52);
   fChain->SetBranchAddress("HLT_MET60_IsoTrk35_Loose", &HLT_MET60_IsoTrk35_Loose, &b_HLT_MET60_IsoTrk35_Loose);
   fChain->SetBranchAddress("HLT_MET75_IsoTrk50", &HLT_MET75_IsoTrk50, &b_HLT_MET75_IsoTrk50);
   fChain->SetBranchAddress("HLT_MET90_IsoTrk50", &HLT_MET90_IsoTrk50, &b_HLT_MET90_IsoTrk50);
   fChain->SetBranchAddress("HLT_PFMET120_BTagCSV_p067", &HLT_PFMET120_BTagCSV_p067, &b_HLT_PFMET120_BTagCSV_p067);
   fChain->SetBranchAddress("HLT_PFMET120_Mu5", &HLT_PFMET120_Mu5, &b_HLT_PFMET120_Mu5);
   fChain->SetBranchAddress("HLT_PFMET170_NoiseCleaned", &HLT_PFMET170_NoiseCleaned, &b_HLT_PFMET170_NoiseCleaned);
   fChain->SetBranchAddress("HLT_PFMET170_HBHECleaned", &HLT_PFMET170_HBHECleaned, &b_HLT_PFMET170_HBHECleaned);
   fChain->SetBranchAddress("HLT_PFMET170_JetIdCleaned", &HLT_PFMET170_JetIdCleaned, &b_HLT_PFMET170_JetIdCleaned);
   fChain->SetBranchAddress("HLT_PFMET170_NotCleaned", &HLT_PFMET170_NotCleaned, &b_HLT_PFMET170_NotCleaned);
   fChain->SetBranchAddress("HLT_PFMET170_BeamHaloCleaned", &HLT_PFMET170_BeamHaloCleaned, &b_HLT_PFMET170_BeamHaloCleaned);
   fChain->SetBranchAddress("HLT_PFMET90_PFMHT90_IDTight", &HLT_PFMET90_PFMHT90_IDTight, &b_HLT_PFMET90_PFMHT90_IDTight);
   fChain->SetBranchAddress("HLT_PFMET100_PFMHT100_IDTight", &HLT_PFMET100_PFMHT100_IDTight, &b_HLT_PFMET100_PFMHT100_IDTight);
   fChain->SetBranchAddress("HLT_PFMET110_PFMHT110_IDTight", &HLT_PFMET110_PFMHT110_IDTight, &b_HLT_PFMET110_PFMHT110_IDTight);
   fChain->SetBranchAddress("HLT_PFMET120_PFMHT120_IDTight", &HLT_PFMET120_PFMHT120_IDTight, &b_HLT_PFMET120_PFMHT120_IDTight);
   fChain->SetBranchAddress("HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight_BTagCSV_p067", &HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight_BTagCSV_p067, &b_HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight_BTagCSV_p067);
   fChain->SetBranchAddress("HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight", &HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight, &b_HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight);
   fChain->SetBranchAddress("HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq200", &HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq200, &b_HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq200);
   fChain->SetBranchAddress("HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq460", &HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq460, &b_HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq460);
   fChain->SetBranchAddress("HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq240", &HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq240, &b_HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq240);
   fChain->SetBranchAddress("HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq500", &HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq500, &b_HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq500);
   fChain->SetBranchAddress("HLT_QuadPFJet_VBF", &HLT_QuadPFJet_VBF, &b_HLT_QuadPFJet_VBF);
   fChain->SetBranchAddress("HLT_L1_TripleJet_VBF", &HLT_L1_TripleJet_VBF, &b_HLT_L1_TripleJet_VBF);
   fChain->SetBranchAddress("HLT_QuadJet45_TripleBTagCSV_p087", &HLT_QuadJet45_TripleBTagCSV_p087, &b_HLT_QuadJet45_TripleBTagCSV_p087);
   fChain->SetBranchAddress("HLT_QuadJet45_DoubleBTagCSV_p087", &HLT_QuadJet45_DoubleBTagCSV_p087, &b_HLT_QuadJet45_DoubleBTagCSV_p087);
   fChain->SetBranchAddress("HLT_DoubleJet90_Double30_TripleBTagCSV_p087", &HLT_DoubleJet90_Double30_TripleBTagCSV_p087, &b_HLT_DoubleJet90_Double30_TripleBTagCSV_p087);
   fChain->SetBranchAddress("HLT_DoubleJet90_Double30_DoubleBTagCSV_p087", &HLT_DoubleJet90_Double30_DoubleBTagCSV_p087, &b_HLT_DoubleJet90_Double30_DoubleBTagCSV_p087);
   fChain->SetBranchAddress("HLT_DoubleJetsC100_DoubleBTagCSV_p026_DoublePFJetsC160", &HLT_DoubleJetsC100_DoubleBTagCSV_p026_DoublePFJetsC160, &b_HLT_DoubleJetsC100_DoubleBTagCSV_p026_DoublePFJetsC160);
   fChain->SetBranchAddress("HLT_DoubleJetsC100_DoubleBTagCSV_p014_DoublePFJetsC100MaxDeta1p6", &HLT_DoubleJetsC100_DoubleBTagCSV_p014_DoublePFJetsC100MaxDeta1p6, &b_HLT_DoubleJetsC100_DoubleBTagCSV_p014_DoublePFJetsC100MaxDeta1p6);
   fChain->SetBranchAddress("HLT_DoubleJetsC112_DoubleBTagCSV_p026_DoublePFJetsC172", &HLT_DoubleJetsC112_DoubleBTagCSV_p026_DoublePFJetsC172, &b_HLT_DoubleJetsC112_DoubleBTagCSV_p026_DoublePFJetsC172);
   fChain->SetBranchAddress("HLT_DoubleJetsC112_DoubleBTagCSV_p014_DoublePFJetsC112MaxDeta1p6", &HLT_DoubleJetsC112_DoubleBTagCSV_p014_DoublePFJetsC112MaxDeta1p6, &b_HLT_DoubleJetsC112_DoubleBTagCSV_p014_DoublePFJetsC112MaxDeta1p6);
   fChain->SetBranchAddress("HLT_DoubleJetsC100_SingleBTagCSV_p026", &HLT_DoubleJetsC100_SingleBTagCSV_p026, &b_HLT_DoubleJetsC100_SingleBTagCSV_p026);
   fChain->SetBranchAddress("HLT_DoubleJetsC100_SingleBTagCSV_p014", &HLT_DoubleJetsC100_SingleBTagCSV_p014, &b_HLT_DoubleJetsC100_SingleBTagCSV_p014);
   fChain->SetBranchAddress("HLT_DoubleJetsC100_SingleBTagCSV_p026_SinglePFJetC350", &HLT_DoubleJetsC100_SingleBTagCSV_p026_SinglePFJetC350, &b_HLT_DoubleJetsC100_SingleBTagCSV_p026_SinglePFJetC350);
   fChain->SetBranchAddress("HLT_DoubleJetsC100_SingleBTagCSV_p014_SinglePFJetC350", &HLT_DoubleJetsC100_SingleBTagCSV_p014_SinglePFJetC350, &b_HLT_DoubleJetsC100_SingleBTagCSV_p014_SinglePFJetC350);
   fChain->SetBranchAddress("HLT_Photon135_PFMET100", &HLT_Photon135_PFMET100, &b_HLT_Photon135_PFMET100);
   fChain->SetBranchAddress("HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_PFMET40", &HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_PFMET40, &b_HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_PFMET40);
   fChain->SetBranchAddress("HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_VBF", &HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_VBF, &b_HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_VBF);
   fChain->SetBranchAddress("HLT_Photon250_NoHE", &HLT_Photon250_NoHE, &b_HLT_Photon250_NoHE);
   fChain->SetBranchAddress("HLT_Photon300_NoHE", &HLT_Photon300_NoHE, &b_HLT_Photon300_NoHE);
   fChain->SetBranchAddress("HLT_Photon26_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon16_AND_HE10_R9Id65_Eta2_Mass60", &HLT_Photon26_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon16_AND_HE10_R9Id65_Eta2_Mass60, &b_HLT_Photon26_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon16_AND_HE10_R9Id65_Eta2_Mass60);
   fChain->SetBranchAddress("HLT_Photon36_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon22_AND_HE10_R9Id65_Eta2_Mass15", &HLT_Photon36_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon22_AND_HE10_R9Id65_Eta2_Mass15, &b_HLT_Photon36_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon22_AND_HE10_R9Id65_Eta2_Mass15);
   fChain->SetBranchAddress("HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_PFMET40", &HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_PFMET40, &b_HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_PFMET40);
   fChain->SetBranchAddress("HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_VBF", &HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_VBF, &b_HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_VBF);
   fChain->SetBranchAddress("HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_PFMET40", &HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_PFMET40, &b_HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_PFMET40);
   fChain->SetBranchAddress("HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_VBF", &HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_VBF, &b_HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_VBF);
   fChain->SetBranchAddress("HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_PFMET40", &HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_PFMET40, &b_HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_PFMET40);
   fChain->SetBranchAddress("HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_VBF", &HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_VBF, &b_HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_VBF);
   fChain->SetBranchAddress("HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_PFMET40", &HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_PFMET40, &b_HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_PFMET40);
   fChain->SetBranchAddress("HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_VBF", &HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_VBF, &b_HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_VBF);
   fChain->SetBranchAddress("HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_PFMET40", &HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_PFMET40, &b_HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_PFMET40);
   fChain->SetBranchAddress("HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_VBF", &HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_VBF, &b_HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_VBF);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL", &HLT_Mu8_TrkIsoVVL, &b_HLT_Mu8_TrkIsoVVL);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL", &HLT_Mu17_TrkIsoVVL, &b_HLT_Mu17_TrkIsoVVL);
   fChain->SetBranchAddress("HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30", &HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30, &b_HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30);
   fChain->SetBranchAddress("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30", &HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30, &b_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30);
   fChain->SetBranchAddress("HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30", &HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30, &b_HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30);
   fChain->SetBranchAddress("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30", &HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30, &b_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30);
   fChain->SetBranchAddress("HLT_BTagMu_DiJet20_Mu5", &HLT_BTagMu_DiJet20_Mu5, &b_HLT_BTagMu_DiJet20_Mu5);
   fChain->SetBranchAddress("HLT_BTagMu_DiJet40_Mu5", &HLT_BTagMu_DiJet40_Mu5, &b_HLT_BTagMu_DiJet40_Mu5);
   fChain->SetBranchAddress("HLT_BTagMu_DiJet70_Mu5", &HLT_BTagMu_DiJet70_Mu5, &b_HLT_BTagMu_DiJet70_Mu5);
   fChain->SetBranchAddress("HLT_BTagMu_DiJet110_Mu5", &HLT_BTagMu_DiJet110_Mu5, &b_HLT_BTagMu_DiJet110_Mu5);
   fChain->SetBranchAddress("HLT_BTagMu_DiJet170_Mu5", &HLT_BTagMu_DiJet170_Mu5, &b_HLT_BTagMu_DiJet170_Mu5);
   fChain->SetBranchAddress("HLT_BTagMu_Jet300_Mu5", &HLT_BTagMu_Jet300_Mu5, &b_HLT_BTagMu_Jet300_Mu5);
   fChain->SetBranchAddress("HLT_BTagMu_AK8Jet300_Mu5", &HLT_BTagMu_AK8Jet300_Mu5, &b_HLT_BTagMu_AK8Jet300_Mu5);
   fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ);
   fChain->SetBranchAddress("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ", &HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ, &b_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ);
   fChain->SetBranchAddress("HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL", &HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL, &b_HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL", &HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL, &b_HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL", &HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL, &b_HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL);
   fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL", &HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL, &b_HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL);
   fChain->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL", &HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL, &b_HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL", &HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL, &b_HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL);
   fChain->SetBranchAddress("HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL", &HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL, &b_HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL);
   fChain->SetBranchAddress("HLT_Mu37_Ele27_CaloIdL_GsfTrkIdVL", &HLT_Mu37_Ele27_CaloIdL_GsfTrkIdVL, &b_HLT_Mu37_Ele27_CaloIdL_GsfTrkIdVL);
   fChain->SetBranchAddress("HLT_Mu27_Ele37_CaloIdL_GsfTrkIdVL", &HLT_Mu27_Ele37_CaloIdL_GsfTrkIdVL, &b_HLT_Mu27_Ele37_CaloIdL_GsfTrkIdVL);
   fChain->SetBranchAddress("HLT_Mu8_DiEle12_CaloIdL_TrackIdL", &HLT_Mu8_DiEle12_CaloIdL_TrackIdL, &b_HLT_Mu8_DiEle12_CaloIdL_TrackIdL);
   fChain->SetBranchAddress("HLT_Mu12_Photon25_CaloIdL", &HLT_Mu12_Photon25_CaloIdL, &b_HLT_Mu12_Photon25_CaloIdL);
   fChain->SetBranchAddress("HLT_Mu12_Photon25_CaloIdL_L1ISO", &HLT_Mu12_Photon25_CaloIdL_L1ISO, &b_HLT_Mu12_Photon25_CaloIdL_L1ISO);
   fChain->SetBranchAddress("HLT_Mu12_Photon25_CaloIdL_L1OR", &HLT_Mu12_Photon25_CaloIdL_L1OR, &b_HLT_Mu12_Photon25_CaloIdL_L1OR);
   fChain->SetBranchAddress("HLT_Mu17_Photon22_CaloIdL_L1ISO", &HLT_Mu17_Photon22_CaloIdL_L1ISO, &b_HLT_Mu17_Photon22_CaloIdL_L1ISO);
   fChain->SetBranchAddress("HLT_Mu17_Photon30_CaloIdL_L1ISO", &HLT_Mu17_Photon30_CaloIdL_L1ISO, &b_HLT_Mu17_Photon30_CaloIdL_L1ISO);
   fChain->SetBranchAddress("HLT_Mu17_Photon35_CaloIdL_L1ISO", &HLT_Mu17_Photon35_CaloIdL_L1ISO, &b_HLT_Mu17_Photon35_CaloIdL_L1ISO);
   fChain->SetBranchAddress("HLT_DiMu9_Ele9_CaloIdL_TrackIdL", &HLT_DiMu9_Ele9_CaloIdL_TrackIdL, &b_HLT_DiMu9_Ele9_CaloIdL_TrackIdL);
   fChain->SetBranchAddress("HLT_TripleMu_5_3_3", &HLT_TripleMu_5_3_3, &b_HLT_TripleMu_5_3_3);
   fChain->SetBranchAddress("HLT_TripleMu_12_10_5", &HLT_TripleMu_12_10_5, &b_HLT_TripleMu_12_10_5);
   fChain->SetBranchAddress("HLT_Mu3er_PFHT140_PFMET125", &HLT_Mu3er_PFHT140_PFMET125, &b_HLT_Mu3er_PFHT140_PFMET125);
   fChain->SetBranchAddress("HLT_Mu6_PFHT200_PFMET80_BTagCSV_p067", &HLT_Mu6_PFHT200_PFMET80_BTagCSV_p067, &b_HLT_Mu6_PFHT200_PFMET80_BTagCSV_p067);
   fChain->SetBranchAddress("HLT_Mu6_PFHT200_PFMET100", &HLT_Mu6_PFHT200_PFMET100, &b_HLT_Mu6_PFHT200_PFMET100);
   fChain->SetBranchAddress("HLT_Mu14er_PFMET100", &HLT_Mu14er_PFMET100, &b_HLT_Mu14er_PFMET100);
   fChain->SetBranchAddress("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL", &HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL, &b_HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL);
   fChain->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL", &HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL, &b_HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL);
   fChain->SetBranchAddress("HLT_Ele12_CaloIdL_TrackIdL_IsoVL", &HLT_Ele12_CaloIdL_TrackIdL_IsoVL, &b_HLT_Ele12_CaloIdL_TrackIdL_IsoVL);
   fChain->SetBranchAddress("HLT_Ele17_CaloIdL_GsfTrkIdVL", &HLT_Ele17_CaloIdL_GsfTrkIdVL, &b_HLT_Ele17_CaloIdL_GsfTrkIdVL);
   fChain->SetBranchAddress("HLT_Ele17_CaloIdL_TrackIdL_IsoVL", &HLT_Ele17_CaloIdL_TrackIdL_IsoVL, &b_HLT_Ele17_CaloIdL_TrackIdL_IsoVL);
   fChain->SetBranchAddress("HLT_Ele23_CaloIdL_TrackIdL_IsoVL", &HLT_Ele23_CaloIdL_TrackIdL_IsoVL, &b_HLT_Ele23_CaloIdL_TrackIdL_IsoVL);
   fChain->SetBranchAddress("HLT_AK8DiPFJet280_200_TrimMass30", &HLT_AK8DiPFJet280_200_TrimMass30, &b_HLT_AK8DiPFJet280_200_TrimMass30);
   fChain->SetBranchAddress("HLT_AK8DiPFJet250_200_TrimMass30", &HLT_AK8DiPFJet250_200_TrimMass30, &b_HLT_AK8DiPFJet250_200_TrimMass30);
   fChain->SetBranchAddress("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20", &HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20, &b_HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20);
   fChain->SetBranchAddress("HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20", &HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20, &b_HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20);
   fChain->SetBranchAddress("HLT_PFHT650_WideJetMJJ900DEtaJJ1p5", &HLT_PFHT650_WideJetMJJ900DEtaJJ1p5, &b_HLT_PFHT650_WideJetMJJ900DEtaJJ1p5);
   fChain->SetBranchAddress("HLT_PFHT650_WideJetMJJ950DEtaJJ1p5", &HLT_PFHT650_WideJetMJJ950DEtaJJ1p5, &b_HLT_PFHT650_WideJetMJJ950DEtaJJ1p5);
   fChain->SetBranchAddress("HLT_Photon22", &HLT_Photon22, &b_HLT_Photon22);
   fChain->SetBranchAddress("HLT_Photon30", &HLT_Photon30, &b_HLT_Photon30);
   fChain->SetBranchAddress("HLT_Photon36", &HLT_Photon36, &b_HLT_Photon36);
   fChain->SetBranchAddress("HLT_Photon50", &HLT_Photon50, &b_HLT_Photon50);
   fChain->SetBranchAddress("HLT_Photon75", &HLT_Photon75, &b_HLT_Photon75);
   fChain->SetBranchAddress("HLT_Photon90", &HLT_Photon90, &b_HLT_Photon90);
   fChain->SetBranchAddress("HLT_Photon120", &HLT_Photon120, &b_HLT_Photon120);
   fChain->SetBranchAddress("HLT_Photon175", &HLT_Photon175, &b_HLT_Photon175);
   fChain->SetBranchAddress("HLT_Photon165_HE10", &HLT_Photon165_HE10, &b_HLT_Photon165_HE10);
   fChain->SetBranchAddress("HLT_Photon22_R9Id90_HE10_IsoM", &HLT_Photon22_R9Id90_HE10_IsoM, &b_HLT_Photon22_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("HLT_Photon30_R9Id90_HE10_IsoM", &HLT_Photon30_R9Id90_HE10_IsoM, &b_HLT_Photon30_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("HLT_Photon36_R9Id90_HE10_IsoM", &HLT_Photon36_R9Id90_HE10_IsoM, &b_HLT_Photon36_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("HLT_Photon50_R9Id90_HE10_IsoM", &HLT_Photon50_R9Id90_HE10_IsoM, &b_HLT_Photon50_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("HLT_Photon75_R9Id90_HE10_IsoM", &HLT_Photon75_R9Id90_HE10_IsoM, &b_HLT_Photon75_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("HLT_Photon90_R9Id90_HE10_IsoM", &HLT_Photon90_R9Id90_HE10_IsoM, &b_HLT_Photon90_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("HLT_Photon120_R9Id90_HE10_IsoM", &HLT_Photon120_R9Id90_HE10_IsoM, &b_HLT_Photon120_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("HLT_Photon165_R9Id90_HE10_IsoM", &HLT_Photon165_R9Id90_HE10_IsoM, &b_HLT_Photon165_R9Id90_HE10_IsoM);
   fChain->SetBranchAddress("HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90", &HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90, &b_HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90);
   fChain->SetBranchAddress("HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelSeedMatch_Mass70", &HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelSeedMatch_Mass70, &b_HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelSeedMatch_Mass70);
   fChain->SetBranchAddress("HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55", &HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55, &b_HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55);
   fChain->SetBranchAddress("HLT_Diphoton30_18_Solid_R9Id_AND_IsoCaloId_AND_HE_R9Id_Mass55", &HLT_Diphoton30_18_Solid_R9Id_AND_IsoCaloId_AND_HE_R9Id_Mass55, &b_HLT_Diphoton30_18_Solid_R9Id_AND_IsoCaloId_AND_HE_R9Id_Mass55);
   fChain->SetBranchAddress("HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55", &HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55, &b_HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55);
   fChain->SetBranchAddress("HLT_Dimuon0_Jpsi_Muon", &HLT_Dimuon0_Jpsi_Muon, &b_HLT_Dimuon0_Jpsi_Muon);
   fChain->SetBranchAddress("HLT_Dimuon0_Upsilon_Muon", &HLT_Dimuon0_Upsilon_Muon, &b_HLT_Dimuon0_Upsilon_Muon);
   fChain->SetBranchAddress("HLT_QuadMuon0_Dimuon0_Jpsi", &HLT_QuadMuon0_Dimuon0_Jpsi, &b_HLT_QuadMuon0_Dimuon0_Jpsi);
   fChain->SetBranchAddress("HLT_QuadMuon0_Dimuon0_Upsilon", &HLT_QuadMuon0_Dimuon0_Upsilon, &b_HLT_QuadMuon0_Dimuon0_Upsilon);
   fChain->SetBranchAddress("HLT_Rsq0p25", &HLT_Rsq0p25, &b_HLT_Rsq0p25);
   fChain->SetBranchAddress("HLT_Rsq0p30", &HLT_Rsq0p30, &b_HLT_Rsq0p30);
   fChain->SetBranchAddress("HLT_RsqMR240_Rsq0p09_MR200", &HLT_RsqMR240_Rsq0p09_MR200, &b_HLT_RsqMR240_Rsq0p09_MR200);
   fChain->SetBranchAddress("HLT_RsqMR240_Rsq0p09_MR200_4jet", &HLT_RsqMR240_Rsq0p09_MR200_4jet, &b_HLT_RsqMR240_Rsq0p09_MR200_4jet);
   fChain->SetBranchAddress("HLT_RsqMR270_Rsq0p09_MR200", &HLT_RsqMR270_Rsq0p09_MR200, &b_HLT_RsqMR270_Rsq0p09_MR200);
   fChain->SetBranchAddress("HLT_RsqMR270_Rsq0p09_MR200_4jet", &HLT_RsqMR270_Rsq0p09_MR200_4jet, &b_HLT_RsqMR270_Rsq0p09_MR200_4jet);
   fChain->SetBranchAddress("HLT_Rsq0p02_MR300_TriPFJet80_60_40_BTagCSV_p063_p20_Mbb60_200", &HLT_Rsq0p02_MR300_TriPFJet80_60_40_BTagCSV_p063_p20_Mbb60_200, &b_HLT_Rsq0p02_MR300_TriPFJet80_60_40_BTagCSV_p063_p20_Mbb60_200);
   fChain->SetBranchAddress("HLT_Rsq0p02_MR300_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200", &HLT_Rsq0p02_MR300_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200, &b_HLT_Rsq0p02_MR300_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200);
   fChain->SetBranchAddress("HLT_HT200_DisplacedDijet40_DisplacedTrack", &HLT_HT200_DisplacedDijet40_DisplacedTrack, &b_HLT_HT200_DisplacedDijet40_DisplacedTrack);
   fChain->SetBranchAddress("HLT_HT250_DisplacedDijet40_DisplacedTrack", &HLT_HT250_DisplacedDijet40_DisplacedTrack, &b_HLT_HT250_DisplacedDijet40_DisplacedTrack);
   fChain->SetBranchAddress("HLT_HT350_DisplacedDijet40_DisplacedTrack", &HLT_HT350_DisplacedDijet40_DisplacedTrack, &b_HLT_HT350_DisplacedDijet40_DisplacedTrack);
   fChain->SetBranchAddress("HLT_HT350_DisplacedDijet80_DisplacedTrack", &HLT_HT350_DisplacedDijet80_DisplacedTrack, &b_HLT_HT350_DisplacedDijet80_DisplacedTrack);
   fChain->SetBranchAddress("HLT_HT350_DisplacedDijet80_Tight_DisplacedTrack", &HLT_HT350_DisplacedDijet80_Tight_DisplacedTrack, &b_HLT_HT350_DisplacedDijet80_Tight_DisplacedTrack);
   fChain->SetBranchAddress("HLT_HT350_DisplacedDijet40_Inclusive", &HLT_HT350_DisplacedDijet40_Inclusive, &b_HLT_HT350_DisplacedDijet40_Inclusive);
   fChain->SetBranchAddress("HLT_HT400_DisplacedDijet40_Inclusive", &HLT_HT400_DisplacedDijet40_Inclusive, &b_HLT_HT400_DisplacedDijet40_Inclusive);
   fChain->SetBranchAddress("HLT_HT500_DisplacedDijet40_Inclusive", &HLT_HT500_DisplacedDijet40_Inclusive, &b_HLT_HT500_DisplacedDijet40_Inclusive);
   fChain->SetBranchAddress("HLT_HT550_DisplacedDijet40_Inclusive", &HLT_HT550_DisplacedDijet40_Inclusive, &b_HLT_HT550_DisplacedDijet40_Inclusive);
   fChain->SetBranchAddress("HLT_HT650_DisplacedDijet80_Inclusive", &HLT_HT650_DisplacedDijet80_Inclusive, &b_HLT_HT650_DisplacedDijet80_Inclusive);
   fChain->SetBranchAddress("HLT_HT750_DisplacedDijet80_Inclusive", &HLT_HT750_DisplacedDijet80_Inclusive, &b_HLT_HT750_DisplacedDijet80_Inclusive);
   fChain->SetBranchAddress("HLT_VBF_DisplacedJet40_DisplacedTrack", &HLT_VBF_DisplacedJet40_DisplacedTrack, &b_HLT_VBF_DisplacedJet40_DisplacedTrack);
   fChain->SetBranchAddress("HLT_VBF_DisplacedJet40_DisplacedTrack_2TrackIP2DSig5", &HLT_VBF_DisplacedJet40_DisplacedTrack_2TrackIP2DSig5, &b_HLT_VBF_DisplacedJet40_DisplacedTrack_2TrackIP2DSig5);
   fChain->SetBranchAddress("HLT_VBF_DisplacedJet40_TightID_DisplacedTrack", &HLT_VBF_DisplacedJet40_TightID_DisplacedTrack, &b_HLT_VBF_DisplacedJet40_TightID_DisplacedTrack);
   fChain->SetBranchAddress("HLT_VBF_DisplacedJet40_Hadronic", &HLT_VBF_DisplacedJet40_Hadronic, &b_HLT_VBF_DisplacedJet40_Hadronic);
   fChain->SetBranchAddress("HLT_VBF_DisplacedJet40_Hadronic_2PromptTrack", &HLT_VBF_DisplacedJet40_Hadronic_2PromptTrack, &b_HLT_VBF_DisplacedJet40_Hadronic_2PromptTrack);
   fChain->SetBranchAddress("HLT_VBF_DisplacedJet40_TightID_Hadronic", &HLT_VBF_DisplacedJet40_TightID_Hadronic, &b_HLT_VBF_DisplacedJet40_TightID_Hadronic);
   fChain->SetBranchAddress("HLT_VBF_DisplacedJet40_VTightID_Hadronic", &HLT_VBF_DisplacedJet40_VTightID_Hadronic, &b_HLT_VBF_DisplacedJet40_VTightID_Hadronic);
   fChain->SetBranchAddress("HLT_VBF_DisplacedJet40_VVTightID_Hadronic", &HLT_VBF_DisplacedJet40_VVTightID_Hadronic, &b_HLT_VBF_DisplacedJet40_VVTightID_Hadronic);
   fChain->SetBranchAddress("HLT_VBF_DisplacedJet40_VTightID_DisplacedTrack", &HLT_VBF_DisplacedJet40_VTightID_DisplacedTrack, &b_HLT_VBF_DisplacedJet40_VTightID_DisplacedTrack);
   fChain->SetBranchAddress("HLT_VBF_DisplacedJet40_VVTightID_DisplacedTrack", &HLT_VBF_DisplacedJet40_VVTightID_DisplacedTrack, &b_HLT_VBF_DisplacedJet40_VVTightID_DisplacedTrack);
   fChain->SetBranchAddress("HLT_PFMETNoMu90_PFMHTNoMu90_IDTight", &HLT_PFMETNoMu90_PFMHTNoMu90_IDTight, &b_HLT_PFMETNoMu90_PFMHTNoMu90_IDTight);
   fChain->SetBranchAddress("HLT_PFMETNoMu100_PFMHTNoMu100_IDTight", &HLT_PFMETNoMu100_PFMHTNoMu100_IDTight, &b_HLT_PFMETNoMu100_PFMHTNoMu100_IDTight);
   fChain->SetBranchAddress("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight", &HLT_PFMETNoMu110_PFMHTNoMu110_IDTight, &b_HLT_PFMETNoMu110_PFMHTNoMu110_IDTight);
   fChain->SetBranchAddress("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight, &b_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight);
   fChain->SetBranchAddress("HLT_MonoCentralPFJet80_PFMETNoMu90_PFMHTNoMu90_IDTight", &HLT_MonoCentralPFJet80_PFMETNoMu90_PFMHTNoMu90_IDTight, &b_HLT_MonoCentralPFJet80_PFMETNoMu90_PFMHTNoMu90_IDTight);
   fChain->SetBranchAddress("HLT_MonoCentralPFJet80_PFMETNoMu100_PFMHTNoMu100_IDTight", &HLT_MonoCentralPFJet80_PFMETNoMu100_PFMHTNoMu100_IDTight, &b_HLT_MonoCentralPFJet80_PFMETNoMu100_PFMHTNoMu100_IDTight);
   fChain->SetBranchAddress("HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight", &HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight, &b_HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight);
   fChain->SetBranchAddress("HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight", &HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight, &b_HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight);
   fChain->SetBranchAddress("HLT_Ele27_eta2p1_WPLoose_Gsf_HT200", &HLT_Ele27_eta2p1_WPLoose_Gsf_HT200, &b_HLT_Ele27_eta2p1_WPLoose_Gsf_HT200);
   fChain->SetBranchAddress("HLT_Photon90_CaloIdL_PFHT500", &HLT_Photon90_CaloIdL_PFHT500, &b_HLT_Photon90_CaloIdL_PFHT500);
   fChain->SetBranchAddress("HLT_DoubleMu8_Mass8_PFHT250", &HLT_DoubleMu8_Mass8_PFHT250, &b_HLT_DoubleMu8_Mass8_PFHT250);
   fChain->SetBranchAddress("HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT250", &HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT250, &b_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT250);
   fChain->SetBranchAddress("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT250", &HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT250, &b_HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT250);
   fChain->SetBranchAddress("HLT_DoubleMu8_Mass8_PFHT300", &HLT_DoubleMu8_Mass8_PFHT300, &b_HLT_DoubleMu8_Mass8_PFHT300);
   fChain->SetBranchAddress("HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300", &HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300, &b_HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300);
   fChain->SetBranchAddress("HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300", &HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300, &b_HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300);
   fChain->SetBranchAddress("HLT_Mu10_CentralPFJet30_BTagCSV_p13", &HLT_Mu10_CentralPFJet30_BTagCSV_p13, &b_HLT_Mu10_CentralPFJet30_BTagCSV_p13);
   fChain->SetBranchAddress("HLT_DoubleMu3_PFMET50", &HLT_DoubleMu3_PFMET50, &b_HLT_DoubleMu3_PFMET50);
   fChain->SetBranchAddress("HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV_p13", &HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV_p13, &b_HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV_p13);
   fChain->SetBranchAddress("HLT_Ele15_IsoVVVL_BTagCSV_p067_PFHT400", &HLT_Ele15_IsoVVVL_BTagCSV_p067_PFHT400, &b_HLT_Ele15_IsoVVVL_BTagCSV_p067_PFHT400);
   fChain->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT350_PFMET50", &HLT_Ele15_IsoVVVL_PFHT350_PFMET50, &b_HLT_Ele15_IsoVVVL_PFHT350_PFMET50);
   fChain->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT600", &HLT_Ele15_IsoVVVL_PFHT600, &b_HLT_Ele15_IsoVVVL_PFHT600);
   fChain->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT350", &HLT_Ele15_IsoVVVL_PFHT350, &b_HLT_Ele15_IsoVVVL_PFHT350);
   fChain->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT400_PFMET50", &HLT_Ele15_IsoVVVL_PFHT400_PFMET50, &b_HLT_Ele15_IsoVVVL_PFHT400_PFMET50);
   fChain->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT400", &HLT_Ele15_IsoVVVL_PFHT400, &b_HLT_Ele15_IsoVVVL_PFHT400);
   fChain->SetBranchAddress("HLT_Ele50_IsoVVVL_PFHT400", &HLT_Ele50_IsoVVVL_PFHT400, &b_HLT_Ele50_IsoVVVL_PFHT400);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60", &HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60, &b_HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60);
   fChain->SetBranchAddress("HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60", &HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60, &b_HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60);
   fChain->SetBranchAddress("HLT_Mu15_IsoVVVL_BTagCSV_p067_PFHT400", &HLT_Mu15_IsoVVVL_BTagCSV_p067_PFHT400, &b_HLT_Mu15_IsoVVVL_BTagCSV_p067_PFHT400);
   fChain->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT350_PFMET50", &HLT_Mu15_IsoVVVL_PFHT350_PFMET50, &b_HLT_Mu15_IsoVVVL_PFHT350_PFMET50);
   fChain->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT600", &HLT_Mu15_IsoVVVL_PFHT600, &b_HLT_Mu15_IsoVVVL_PFHT600);
   fChain->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT350", &HLT_Mu15_IsoVVVL_PFHT350, &b_HLT_Mu15_IsoVVVL_PFHT350);
   fChain->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT400_PFMET50", &HLT_Mu15_IsoVVVL_PFHT400_PFMET50, &b_HLT_Mu15_IsoVVVL_PFHT400_PFMET50);
   fChain->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT400", &HLT_Mu15_IsoVVVL_PFHT400, &b_HLT_Mu15_IsoVVVL_PFHT400);
   fChain->SetBranchAddress("HLT_Mu50_IsoVVVL_PFHT400", &HLT_Mu50_IsoVVVL_PFHT400, &b_HLT_Mu50_IsoVVVL_PFHT400);
   fChain->SetBranchAddress("HLT_Dimuon16_Jpsi", &HLT_Dimuon16_Jpsi, &b_HLT_Dimuon16_Jpsi);
   fChain->SetBranchAddress("HLT_Dimuon10_Jpsi_Barrel", &HLT_Dimuon10_Jpsi_Barrel, &b_HLT_Dimuon10_Jpsi_Barrel);
   fChain->SetBranchAddress("HLT_Dimuon8_PsiPrime_Barrel", &HLT_Dimuon8_PsiPrime_Barrel, &b_HLT_Dimuon8_PsiPrime_Barrel);
   fChain->SetBranchAddress("HLT_Dimuon8_Upsilon_Barrel", &HLT_Dimuon8_Upsilon_Barrel, &b_HLT_Dimuon8_Upsilon_Barrel);
   fChain->SetBranchAddress("HLT_Dimuon0_Phi_Barrel", &HLT_Dimuon0_Phi_Barrel, &b_HLT_Dimuon0_Phi_Barrel);
   fChain->SetBranchAddress("HLT_Mu16_TkMu0_dEta18_Onia", &HLT_Mu16_TkMu0_dEta18_Onia, &b_HLT_Mu16_TkMu0_dEta18_Onia);
   fChain->SetBranchAddress("HLT_Mu16_TkMu0_dEta18_Phi", &HLT_Mu16_TkMu0_dEta18_Phi, &b_HLT_Mu16_TkMu0_dEta18_Phi);
   fChain->SetBranchAddress("HLT_TrkMu15_DoubleTrkMu5NoFiltersNoVtx", &HLT_TrkMu15_DoubleTrkMu5NoFiltersNoVtx, &b_HLT_TrkMu15_DoubleTrkMu5NoFiltersNoVtx);
   fChain->SetBranchAddress("HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx", &HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx, &b_HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx);
   fChain->SetBranchAddress("HLT_Mu8", &HLT_Mu8, &b_HLT_Mu8);
   fChain->SetBranchAddress("HLT_Mu17", &HLT_Mu17, &b_HLT_Mu17);
   fChain->SetBranchAddress("HLT_Mu3_PFJet40", &HLT_Mu3_PFJet40, &b_HLT_Mu3_PFJet40);
   fChain->SetBranchAddress("HLT_Ele8_CaloIdM_TrackIdM_PFJet30", &HLT_Ele8_CaloIdM_TrackIdM_PFJet30, &b_HLT_Ele8_CaloIdM_TrackIdM_PFJet30);
   fChain->SetBranchAddress("HLT_Ele12_CaloIdM_TrackIdM_PFJet30", &HLT_Ele12_CaloIdM_TrackIdM_PFJet30, &b_HLT_Ele12_CaloIdM_TrackIdM_PFJet30);
   fChain->SetBranchAddress("HLT_Ele17_CaloIdM_TrackIdM_PFJet30", &HLT_Ele17_CaloIdM_TrackIdM_PFJet30, &b_HLT_Ele17_CaloIdM_TrackIdM_PFJet30);
   fChain->SetBranchAddress("HLT_Ele23_CaloIdM_TrackIdM_PFJet30", &HLT_Ele23_CaloIdM_TrackIdM_PFJet30, &b_HLT_Ele23_CaloIdM_TrackIdM_PFJet30);
   fChain->SetBranchAddress("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet140", &HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet140, &b_HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet140);
   fChain->SetBranchAddress("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165", &HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165, &b_HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165);
   fChain->SetBranchAddress("HLT_PFHT400_SixJet30_DoubleBTagCSV_p056", &HLT_PFHT400_SixJet30_DoubleBTagCSV_p056, &b_HLT_PFHT400_SixJet30_DoubleBTagCSV_p056);
   fChain->SetBranchAddress("HLT_PFHT450_SixJet40_BTagCSV_p056", &HLT_PFHT450_SixJet40_BTagCSV_p056, &b_HLT_PFHT450_SixJet40_BTagCSV_p056);
   fChain->SetBranchAddress("HLT_PFHT400_SixJet30", &HLT_PFHT400_SixJet30, &b_HLT_PFHT400_SixJet30);
   fChain->SetBranchAddress("HLT_PFHT450_SixJet40", &HLT_PFHT450_SixJet40, &b_HLT_PFHT450_SixJet40);
   fChain->SetBranchAddress("HLT_Ele115_CaloIdVT_GsfTrkIdT", &HLT_Ele115_CaloIdVT_GsfTrkIdT, &b_HLT_Ele115_CaloIdVT_GsfTrkIdT);
   fChain->SetBranchAddress("HLT_Mu55", &HLT_Mu55, &b_HLT_Mu55);
   fChain->SetBranchAddress("HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15", &HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15, &b_HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15);
   fChain->SetBranchAddress("HLT_Photon90_CaloIdL_PFHT600", &HLT_Photon90_CaloIdL_PFHT600, &b_HLT_Photon90_CaloIdL_PFHT600);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity60ForEndOfFill", &HLT_PixelTracks_Multiplicity60ForEndOfFill, &b_HLT_PixelTracks_Multiplicity60ForEndOfFill);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity85ForEndOfFill", &HLT_PixelTracks_Multiplicity85ForEndOfFill, &b_HLT_PixelTracks_Multiplicity85ForEndOfFill);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity110ForEndOfFill", &HLT_PixelTracks_Multiplicity110ForEndOfFill, &b_HLT_PixelTracks_Multiplicity110ForEndOfFill);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity135ForEndOfFill", &HLT_PixelTracks_Multiplicity135ForEndOfFill, &b_HLT_PixelTracks_Multiplicity135ForEndOfFill);
   fChain->SetBranchAddress("HLT_PixelTracks_Multiplicity160ForEndOfFill", &HLT_PixelTracks_Multiplicity160ForEndOfFill, &b_HLT_PixelTracks_Multiplicity160ForEndOfFill);
   fChain->SetBranchAddress("HLT_FullTracks_Multiplicity80", &HLT_FullTracks_Multiplicity80, &b_HLT_FullTracks_Multiplicity80);
   fChain->SetBranchAddress("HLT_FullTracks_Multiplicity100", &HLT_FullTracks_Multiplicity100, &b_HLT_FullTracks_Multiplicity100);
   fChain->SetBranchAddress("HLT_FullTracks_Multiplicity130", &HLT_FullTracks_Multiplicity130, &b_HLT_FullTracks_Multiplicity130);
   fChain->SetBranchAddress("HLT_FullTracks_Multiplicity150", &HLT_FullTracks_Multiplicity150, &b_HLT_FullTracks_Multiplicity150);
   fChain->SetBranchAddress("HLT_ECALHT800", &HLT_ECALHT800, &b_HLT_ECALHT800);
   fChain->SetBranchAddress("HLT_DiSC30_18_EIso_AND_HE_Mass70", &HLT_DiSC30_18_EIso_AND_HE_Mass70, &b_HLT_DiSC30_18_EIso_AND_HE_Mass70);
   fChain->SetBranchAddress("HLT_MET200", &HLT_MET200, &b_HLT_MET200);
   fChain->SetBranchAddress("HLT_Ele27_HighEta_Ele20_Mass55", &HLT_Ele27_HighEta_Ele20_Mass55, &b_HLT_Ele27_HighEta_Ele20_Mass55);
   fChain->SetBranchAddress("HLT_L1FatEvents", &HLT_L1FatEvents, &b_HLT_L1FatEvents);
   fChain->SetBranchAddress("HLT_Physics", &HLT_Physics, &b_HLT_Physics);
   fChain->SetBranchAddress("HLT_Physics_part0", &HLT_Physics_part0, &b_HLT_Physics_part0);
   fChain->SetBranchAddress("HLT_Physics_part1", &HLT_Physics_part1, &b_HLT_Physics_part1);
   fChain->SetBranchAddress("HLT_Physics_part2", &HLT_Physics_part2, &b_HLT_Physics_part2);
   fChain->SetBranchAddress("HLT_Physics_part3", &HLT_Physics_part3, &b_HLT_Physics_part3);
   fChain->SetBranchAddress("HLT_Random", &HLT_Random, &b_HLT_Random);
   fChain->SetBranchAddress("HLT_ZeroBias", &HLT_ZeroBias, &b_HLT_ZeroBias);
   fChain->SetBranchAddress("HLT_AK4CaloJet30", &HLT_AK4CaloJet30, &b_HLT_AK4CaloJet30);
   fChain->SetBranchAddress("HLT_AK4CaloJet40", &HLT_AK4CaloJet40, &b_HLT_AK4CaloJet40);
   fChain->SetBranchAddress("HLT_AK4CaloJet50", &HLT_AK4CaloJet50, &b_HLT_AK4CaloJet50);
   fChain->SetBranchAddress("HLT_AK4CaloJet80", &HLT_AK4CaloJet80, &b_HLT_AK4CaloJet80);
   fChain->SetBranchAddress("HLT_AK4CaloJet100", &HLT_AK4CaloJet100, &b_HLT_AK4CaloJet100);
   fChain->SetBranchAddress("HLT_AK4PFJet30", &HLT_AK4PFJet30, &b_HLT_AK4PFJet30);
   fChain->SetBranchAddress("HLT_AK4PFJet50", &HLT_AK4PFJet50, &b_HLT_AK4PFJet50);
   fChain->SetBranchAddress("HLT_AK4PFJet80", &HLT_AK4PFJet80, &b_HLT_AK4PFJet80);
   fChain->SetBranchAddress("HLT_AK4PFJet100", &HLT_AK4PFJet100, &b_HLT_AK4PFJet100);
   fChain->SetBranchAddress("HLT_HISinglePhoton10", &HLT_HISinglePhoton10, &b_HLT_HISinglePhoton10);
   fChain->SetBranchAddress("HLT_HISinglePhoton15", &HLT_HISinglePhoton15, &b_HLT_HISinglePhoton15);
   fChain->SetBranchAddress("HLT_HISinglePhoton20", &HLT_HISinglePhoton20, &b_HLT_HISinglePhoton20);
   fChain->SetBranchAddress("HLT_HISinglePhoton40", &HLT_HISinglePhoton40, &b_HLT_HISinglePhoton40);
   fChain->SetBranchAddress("HLT_HISinglePhoton60", &HLT_HISinglePhoton60, &b_HLT_HISinglePhoton60);
   fChain->SetBranchAddress("HLT_EcalCalibration", &HLT_EcalCalibration, &b_HLT_EcalCalibration);
   fChain->SetBranchAddress("HLT_HcalCalibration", &HLT_HcalCalibration, &b_HLT_HcalCalibration);
   fChain->SetBranchAddress("HLT_GlobalRunHPDNoise", &HLT_GlobalRunHPDNoise, &b_HLT_GlobalRunHPDNoise);
   fChain->SetBranchAddress("HLT_L1BptxMinus", &HLT_L1BptxMinus, &b_HLT_L1BptxMinus);
   fChain->SetBranchAddress("HLT_L1BptxPlus", &HLT_L1BptxPlus, &b_HLT_L1BptxPlus);
   fChain->SetBranchAddress("HLT_L1NotBptxOR", &HLT_L1NotBptxOR, &b_HLT_L1NotBptxOR);
   fChain->SetBranchAddress("HLT_L1BeamGasMinus", &HLT_L1BeamGasMinus, &b_HLT_L1BeamGasMinus);
   fChain->SetBranchAddress("HLT_L1BeamGasPlus", &HLT_L1BeamGasPlus, &b_HLT_L1BeamGasPlus);
   fChain->SetBranchAddress("HLT_L1BptxXOR", &HLT_L1BptxXOR, &b_HLT_L1BptxXOR);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF_OR", &HLT_L1MinimumBiasHF_OR, &b_HLT_L1MinimumBiasHF_OR);
   fChain->SetBranchAddress("HLT_L1MinimumBiasHF_AND", &HLT_L1MinimumBiasHF_AND, &b_HLT_L1MinimumBiasHF_AND);
   fChain->SetBranchAddress("HLT_HcalNZS", &HLT_HcalNZS, &b_HLT_HcalNZS);
   fChain->SetBranchAddress("HLT_HcalPhiSym", &HLT_HcalPhiSym, &b_HLT_HcalPhiSym);
   fChain->SetBranchAddress("HLT_ZeroBias_FirstCollisionAfterAbortGap", &HLT_ZeroBias_FirstCollisionAfterAbortGap, &b_HLT_ZeroBias_FirstCollisionAfterAbortGap);
   fChain->SetBranchAddress("HLT_ZeroBias_FirstCollisionAfterAbortGap_TCDS", &HLT_ZeroBias_FirstCollisionAfterAbortGap_TCDS, &b_HLT_ZeroBias_FirstCollisionAfterAbortGap_TCDS);
   fChain->SetBranchAddress("HLT_ZeroBias_IsolatedBunches", &HLT_ZeroBias_IsolatedBunches, &b_HLT_ZeroBias_IsolatedBunches);
   fChain->SetBranchAddress("HLT_Photon500", &HLT_Photon500, &b_HLT_Photon500);
   fChain->SetBranchAddress("HLT_Photon600", &HLT_Photon600, &b_HLT_Photon600);
   fChain->SetBranchAddress("HLT_Mu300", &HLT_Mu300, &b_HLT_Mu300);
   fChain->SetBranchAddress("HLT_Mu350", &HLT_Mu350, &b_HLT_Mu350);
   fChain->SetBranchAddress("HLT_MET250", &HLT_MET250, &b_HLT_MET250);
   fChain->SetBranchAddress("HLT_MET300", &HLT_MET300, &b_HLT_MET300);
   fChain->SetBranchAddress("HLT_MET600", &HLT_MET600, &b_HLT_MET600);
   fChain->SetBranchAddress("HLT_MET700", &HLT_MET700, &b_HLT_MET700);
   fChain->SetBranchAddress("HLT_PFMET300", &HLT_PFMET300, &b_HLT_PFMET300);
   fChain->SetBranchAddress("HLT_PFMET400", &HLT_PFMET400, &b_HLT_PFMET400);
   fChain->SetBranchAddress("HLT_PFMET500", &HLT_PFMET500, &b_HLT_PFMET500);
   fChain->SetBranchAddress("HLT_PFMET600", &HLT_PFMET600, &b_HLT_PFMET600);
   fChain->SetBranchAddress("HLT_HT2000", &HLT_HT2000, &b_HLT_HT2000);
   fChain->SetBranchAddress("HLT_HT2500", &HLT_HT2500, &b_HLT_HT2500);
   fChain->SetBranchAddress("HLT_IsoTrackHE", &HLT_IsoTrackHE, &b_HLT_IsoTrackHE);
   fChain->SetBranchAddress("HLT_IsoTrackHB", &HLT_IsoTrackHB, &b_HLT_IsoTrackHB);
   fChain->SetBranchAddress("HLTriggerFinalPath", &HLTriggerFinalPath, &b_HLTriggerFinalPath);
   fChain->SetBranchAddress("L1Reco_step", &L1Reco_step, &b_L1Reco_step);
   fChain->SetBranchAddress("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, &b_Flag_HBHENoiseFilter);
   fChain->SetBranchAddress("Flag_HBHENoiseIsoFilter", &Flag_HBHENoiseIsoFilter, &b_Flag_HBHENoiseIsoFilter);
   fChain->SetBranchAddress("Flag_CSCTightHaloFilter", &Flag_CSCTightHaloFilter, &b_Flag_CSCTightHaloFilter);
   fChain->SetBranchAddress("Flag_CSCTightHaloTrkMuUnvetoFilter", &Flag_CSCTightHaloTrkMuUnvetoFilter, &b_Flag_CSCTightHaloTrkMuUnvetoFilter);
   fChain->SetBranchAddress("Flag_CSCTightHalo2015Filter", &Flag_CSCTightHalo2015Filter, &b_Flag_CSCTightHalo2015Filter);
   fChain->SetBranchAddress("Flag_globalTightHalo2016Filter", &Flag_globalTightHalo2016Filter, &b_Flag_globalTightHalo2016Filter);
   fChain->SetBranchAddress("Flag_globalSuperTightHalo2016Filter", &Flag_globalSuperTightHalo2016Filter, &b_Flag_globalSuperTightHalo2016Filter);
   fChain->SetBranchAddress("Flag_HcalStripHaloFilter", &Flag_HcalStripHaloFilter, &b_Flag_HcalStripHaloFilter);
   fChain->SetBranchAddress("Flag_hcalLaserEventFilter", &Flag_hcalLaserEventFilter, &b_Flag_hcalLaserEventFilter);
   fChain->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, &b_Flag_EcalDeadCellTriggerPrimitiveFilter);
   fChain->SetBranchAddress("Flag_EcalDeadCellBoundaryEnergyFilter", &Flag_EcalDeadCellBoundaryEnergyFilter, &b_Flag_EcalDeadCellBoundaryEnergyFilter);
   fChain->SetBranchAddress("Flag_goodVertices", &Flag_goodVertices, &b_Flag_goodVertices);
   fChain->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter, &b_Flag_eeBadScFilter);
   fChain->SetBranchAddress("Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter, &b_Flag_ecalLaserCorrFilter);
   fChain->SetBranchAddress("Flag_trkPOGFilters", &Flag_trkPOGFilters, &b_Flag_trkPOGFilters);
   fChain->SetBranchAddress("Flag_chargedHadronTrackResolutionFilter", &Flag_chargedHadronTrackResolutionFilter, &b_Flag_chargedHadronTrackResolutionFilter);
   fChain->SetBranchAddress("Flag_muonBadTrackFilter", &Flag_muonBadTrackFilter, &b_Flag_muonBadTrackFilter);
//    fChain->SetBranchAddress("Flag_BadChargedCandidateFilter", &Flag_BadChargedCandidateFilter, &b_Flag_BadChargedCandidateFilter);
//    fChain->SetBranchAddress("Flag_BadPFMuonFilter", &Flag_BadPFMuonFilter, &b_Flag_BadPFMuonFilter);
   fChain->SetBranchAddress("Flag_BadChargedCandidateSummer16Filter", &Flag_BadChargedCandidateSummer16Filter, &b_Flag_BadChargedCandidateSummer16Filter);
   fChain->SetBranchAddress("Flag_BadPFMuonSummer16Filter", &Flag_BadPFMuonSummer16Filter, &b_Flag_BadPFMuonSummer16Filter);
   fChain->SetBranchAddress("Flag_trkPOG_manystripclus53X", &Flag_trkPOG_manystripclus53X, &b_Flag_trkPOG_manystripclus53X);
   fChain->SetBranchAddress("Flag_trkPOG_toomanystripclus53X", &Flag_trkPOG_toomanystripclus53X, &b_Flag_trkPOG_toomanystripclus53X);
   fChain->SetBranchAddress("Flag_trkPOG_logErrorTooManyClusters", &Flag_trkPOG_logErrorTooManyClusters, &b_Flag_trkPOG_logErrorTooManyClusters);
   fChain->SetBranchAddress("Flag_METFilters", &Flag_METFilters, &b_Flag_METFilters);
   Notify();
}

Bool_t NanoBase::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void NanoBase::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

Int_t NanoBase::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}


Bool_t NanoBase::ReadJob(const string& jobFile, int& nFiles)
{
  // Open the file containing the datacards            
  ifstream fin(jobFile.c_str(), ios::in);
  if (!fin) {
    cerr << "==> Input File: <<" << jobFile << ">> could not be opened!" << endl;
    return false;
  }

  static constexpr int BUF_SIZE = 256;
  char buf[BUF_SIZE];
  while (fin.getline(buf, BUF_SIZE, '\n')) {  // Pops off the newline character 
    string line(buf);
    if (line.empty() || line == "START") continue;
    // enable '#' and '//' style comments 
    if (line.substr(0,1) == "#" || line.substr(0,2) == "//") continue;
    if (line == "END") break;

    // Split the line into words
    vector<string> tokens;
    NanoUtil::tokenize(line, tokens);
    int vsize = tokens.size();
    assert(vsize > 1);

    const string& key   = tokens.at(0);
    const string& value = tokens.at(1);
    if (key == "dataType") {
      string vtmp(value);
      std::transform(vtmp.begin(), vtmp.end(), vtmp.begin(), ::toupper);
      vector<string> dt;
      NanoUtil::tokenize(vtmp, dt, "#");
      if (dt.size()) {
        isMC_ = (dt.at(0) == "MC") ? true : false;
        if (isMC_ && dt.size() > 1) {
          isSignal_ = (dt.at(1) == "SIGNAL") ? true : false;
        }
      }
    }
    else if (key == "inputFile")
      NanoUtil::buildList(tokens, fileList_);
    else if (key == "maxEvent")
      maxEvt_ = std::stoi(value.c_str());
    else if (key == "outputFile")
      outFile_ = value;
    else if (key == "eventId" && tokens.size() == 4)
      NanoUtil::buildMap(tokens, eventIdMap_);
    else {
      if (0) cout << "==> " << line << endl;
      NanoUtil::storeCuts(tokens, hmap_);
    }
  }
  // Close the file     
  fin.close();
  //  if (!isSignal_) readGenInfo_ = false;
  fChain = new TChain("Events");
  
  // Build the chain of root files
  for (const auto& fname: fileList_) {
    cout << ">>> INFO. Adding input file " << fname << " to TChain " << endl;
    ++nFiles;
    fChain->Add(fname.c_str());
    int nevt = static_cast<int>(fChain->GetEntries());
    if (maxEvt_ > 0 && nevt >= maxEvt_) break;
  }
  if (!nFiles) {
    cerr << ">>> WARN. Input Root file list is empty! exiting ..." << endl;
    return false;
  }
  cout<<nFiles<<endl;
  return true;
}
