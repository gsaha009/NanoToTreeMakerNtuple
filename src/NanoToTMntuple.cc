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
#include <sstream>
#include <utility> 
#include <typeinfo>
#include <memory>

#include "TROOT.h"
#include "TSystem.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TRandom.h"
#include "TStopwatch.h"
#include "TFile.h"
#include "TH1K.h"

#include "NanoToTMntuple.h"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::ostringstream;
using std::vector;
using std::map;
using std::pair;
using std::setprecision;
using std::setw;

using namespace vhtm;

// -----------
// Constructor
// -----------
NanoToTMntuple::NanoToTMntuple()
  : NanoBase()
{
  listEv_      = new std::vector<vhtm::Event>();
  listPVtx_    = new std::vector<vhtm::Vertex>();
  listSVtx_    = new std::vector<vhtm::Vertex>();
  listTau_     = new std::vector<vhtm::Tau>();
  listEle_     = new std::vector<vhtm::Electron>();
  listMu_      = new std::vector<vhtm::Muon>();
  listPhoton_  = new std::vector<vhtm::Photon>();
  listJet_     = new std::vector<vhtm::Jet>();
  listSVtx_    = new std::vector<vhtm::Vertex>();
  listMET_     = new std::vector<vhtm::MET>();

  listTrgObj_  = new std::vector<vhtm::TriggerObject>();

  hltpaths_    = new std::vector<std::string>();
  hltresults_  = new std::vector<bool>();
  hltprescales_= new std::vector<int>();
  flags_       = new std::vector<std::string>();
  flagValues_  = new std::vector<bool>();
  if (isMC()) {
    listGenEv_  = new std::vector<vhtm::GenEvent>();
    listGenP_   = new std::vector<vhtm::GenParticle>();
    listGenJet_ = new std::vector<vhtm::GenJet>();
    listGenMet_ = new std::vector<vhtm::GenMET>();
  }
}

// ----------
// Destructor
// ----------
NanoToTMntuple::~NanoToTMntuple() 
{
}
// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
bool NanoToTMntuple::beginJob() 
{ 
  outf()->cd();
  outf()->mkdir("treeCreator");
  outTree = new TTree("vhtree", "Analysisvhtm");
  assert(outTree);

  if (!outTree) std::cout<<"========No output tree opened======"<<std::endl;
  outTree->Branch("Event", "std::vector<vhtm::Event>", &listEv_, 32000, -1);
  outTree->Branch("PVertex", "std::vector<vhtm::Vertex>", &listPVtx_, 32000, -1);
  outTree->Branch("SVertex", "std::vector<vhtm::Vertex>", &listSVtx_, 32000, -1);

  outTree->Branch("Electron", "std::vector<vhtm::Electron>", &listEle_, 32000, -1);
  outTree->Branch("Muon", "std::vector<vhtm::Muon>", &listMu_, 32000, -1);

  outTree->Branch("Jet", "std::vector<vhtm::Jet>", &listJet_, 32000, -1);
  //outTree->Branch("PackedPFCandidate", "std::vector<vhtm::PackedPFCandidate>", &listPPF_, 32000, -1);
  outTree->Branch("MET", "std::vector<vhtm::MET>", &listMET_, 32000, -1);
  outTree->Branch("RawMET", "std::vector<vhtm::MET>", &listRMET_, 32000, -1);
  outTree->Branch("PuppiMET", "std::vector<vhtm::MET>", &listPMET_, 32000, -1);
  outTree->Branch("Tau", "std::vector<vhtm::Tau>", &listTau_, 32000, -1);
  outTree->Branch("Photon", "std::vector<vhtm::Photon>", &listPhoton_, 32000, -1);
  outTree->Branch("hltpaths", "std::vector<std::string>", &hltpaths_, 32000, -1);
  outTree->Branch("hltresults", "std::vector<bool>", &hltresults_, 32000, -1);
  outTree->Branch("hltprescales", "std::vector<int>", &hltprescales_, 32000, -1);
  outTree->Branch("Flags", "std::vector<std::string>", &flags_, 32000, -1);
  outTree->Branch("FlagValues", "std::vector<bool>", &flagValues_, 32000, -1);
  //new
  if (isMC()) {
    outTree->Branch("GenEvent", "std::vector<vhtm::GenEvent>", &listGenEv_, 32000, -1);
    outTree->Branch("GenParticle", "std::vector<vhtm::GenParticle>", &listGenP_, 32000, -1);
    outTree->Branch("GenJet", "std::vector<vhtm::GenJet>", &listGenJet_, 32000, -1);
    outTree->Branch("GenMET", "std::vector<vhtm::GenMET>", &listGenMet_, 32000, -1);
  }
#ifdef  SKIP_DUPLICATE_ALL
  eventIdStore_.clear();
#endif
  return true;
}

void NanoToTMntuple::clearLists()
{
  // Reset the vector and the nObj variables
  listEv_  -> clear();
  listPVtx_-> clear();
  listSVtx_-> clear();

  listEle_ -> clear();
  listMu_  -> clear();

  listJet_ -> clear();
  //listPPF_->clear();
  listMET_ -> clear();
  listRMET_ -> clear();
  listPMET_ -> clear();
  listTau_ -> clear();
  hltpaths_-> clear();
  hltresults_  -> clear();
  hltprescales_-> clear();
  flags_      -> clear();
  flagValues_ -> clear();
  //new
  if (isMC()) {
    listGenEv_  -> clear();
    listGenP_   -> clear();
    listGenJet_ -> clear();
    listGenMet_ -> clear();
  }
}

void NanoToTMntuple::eventLoop()
{
  size_t nEntries = fChain->GetEntriesFast();
  cout << "** Chain contains " << nEntries << " events" << endl;

  size_t nEvents = (maxEvt_ < 0) ? nEntries : maxEvt_;
  cout << "** Starting Analysis with  " << nEvents << " events" << endl;

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//    Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch


   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      clearLists();

      std::cout<<"Tree No. "<<jentry<<"\t"<<"nEntries ."<<ientry<<std::endl;
      std::cout<<"nMuons :"<<nMuon<<"\t"<<"nElectrons: "<<nElectron<<std::endl;

      //--------------------------
      //Assigning Event Properties
      //--------------------------
      vhtm::Event ev;
      ev.run              = run;
      ev.event            = event;
      ev.lumis            = luminosityBlock;
      ev.btagWeightCMVA   = btagWeight_CMVA; //b-tag event weight for CMVA
      ev.btagWeightCSVV2  = btagWeight_CSVV2; //b-tag event weight for CSVV2
      ev.fixedGridRhoFastjetAll            = fixedGridRhoFastjetAll;
      ev.fixedGridRhoFastjetCentralCalo    = fixedGridRhoFastjetCentralCalo;
      ev.fixedGridRhoFastjetCentralNeutral = fixedGridRhoFastjetCentralNeutral;

      //PileUp
      ev.nPU        = Pileup_nPU; //the number of pileup interactions that have been added to the event in the current bunch crossing
      ev.trueNInt   = Pileup_nTrueInt; //the true mean number of the poisson distribution for this event from which the number of interactions each bunch crossing has been sampled
      ev.sumEOOT    = Pileup_sumEOOT; //number of early out of time pileup
      ev.sumLOOT    = Pileup_sumLOOT; //number of late out of time pileup

      //other PVs
      ev.nvtx       = nOtherPV;
      

      listEv_->push_back(ev);


      //------------------
      //Assigning Flags
      //------------------
      flags_ {
	"Flag_HBHENoiseFilter",
	  "Flag_HBHENoiseIsoFilter",
	  "Flag_CSCTightHaloFilter",
	  "Flag_CSCTightHaloTrkMuUnvetoFilter",
	  "Flag_CSCTightHalo2015Filter",
	  "Flag_globalTightHalo2016Filter",
	  "Flag_globalSuperTightHalo2016Filter",
	  "Flag_HcalStripHaloFilter",
	  "Flag_hcalLaserEventFilter",
	  "Flag_EcalDeadCellTriggerPrimitiveFilter",
	  "Flag_EcalDeadCellBoundaryEnergyFilter",
	  "Flag_goodVertices",
	  "Flag_eeBadScFilter",
	  "Flag_ecalLaserCorrFilter",
	  "Flag_trkPOGFilters",
	  "Flag_chargedHadronTrackResolutionFilter",
	  "Flag_muonBadTrackFilter",
	  "Flag_BadChargedCandidateFilter",
	  "Flag_BadPFMuonFilter",
	  "Flag_BadChargedCandidateSummer16Filter",
	  "Flag_BadPFMuonSummer16Filter",
	  "Flag_trkPOG_manystripclus53X",
	  "Flag_trkPOG_toomanystripclus53X",
	  "Flag_trkPOG_logErrorTooManyClusters",
	  "Flag_METFilters"
	  };
      flagValues_ {
	Flag_HBHENoiseFilter,
	  Flag_HBHENoiseIsoFilter,
	  Flag_CSCTightHaloFilter,
	  Flag_CSCTightHaloTrkMuUnvetoFilter,
	  Flag_CSCTightHalo2015Filter,
	  Flag_globalTightHalo2016Filter,
	  Flag_globalSuperTightHalo2016Filter,
	  Flag_HcalStripHaloFilter,
	  Flag_hcalLaserEventFilter,
	  Flag_EcalDeadCellTriggerPrimitiveFilter,
	  Flag_EcalDeadCellBoundaryEnergyFilter,
	  Flag_goodVertices,
	  Flag_eeBadScFilter,
	  Flag_ecalLaserCorrFilter,
	  Flag_trkPOGFilters,
	  Flag_chargedHadronTrackResolutionFilter,
	  Flag_muonBadTrackFilter,
	  Flag_BadChargedCandidateFilter,
	  Flag_BadPFMuonFilter,
	  Flag_BadChargedCandidateSummer16Filter,
	  Flag_BadPFMuonSummer16Filter,
	  Flag_trkPOG_manystripclus53X,
	  Flag_trkPOG_toomanystripclus53X,
	  Flag_trkPOG_logErrorTooManyClusters,
	  Flag_METFilters
      };
      
      //-----------------------------
      //Assigning GenEvent Properties
      //-----------------------------
      if (isMC()){
	vhtm::GenEvent gev;
	gev.evtWeight    = genWeight; //generator weight
	gev.genWeight    = Generator_weight; //MC generator weight
	gev.qScalePdf    = Generator_scalePDF; //Q2 scale for PDF
	gev.x1           = Generator_x1; //x1 fraction of proton momentum carried by the first parton
	gev.x2           = Generator_x2; //x2 fraction of proton momentum carried by the second parton
	gev.x1Pdf        = Generator_xpdf1; //x*pdf(x) for the first parton
	gev.x2Pdf        = Generator_xpdf2; //x*pdf(x) for the second parton

	listGenEv_->push_back(gev);
      }
      
      //--------------------------------
      //Assigning GenParticle Properties
      //--------------------------------
      if (isMC()){
	for (size_t i = 0; i < nGenPart; ++i){
	  vhtm::GenParticle gp;
	  gp.eta            = GenPart_eta[i];
	  gp.phi            = GenPart_phi[i];
	  gp.pt             = GenPart_pt[i];
	  gp.mass           = GenPart_mass[i]; // Mass stored for all particles with mass > 10 GeV and photons with mass > 1 GeV. 
	                                             // For other particles you can lookup from PDGID
	  gp.pdgId          = GenPart_pdgId[i];
	  gp.status         = GenPart_status[i];
	  gp.motherIndex    = GenPart_genPartIdxMother[i]; //index of the mother particle
	  gp.statusFlag     = GenPart_statusFlags[i]; //gen status flags stored bitwise, bits are: 0 : isPrompt, 
	                                                    // 1 : isDecayedLeptonHadron, 2 : isTauDecayProduct, 
	                                                    // 3 : isPromptTauDecayProduct, 4 : isDirectTauDecayProduct, 
	                                                    // 5 : isDirectPromptTauDecayProduct, 6 : isDirectHadronDecayProduct, 
	                                                    // 7 : isHardProcess, 8 : fromHardProcess, 9 : isHardProcessTauDecayProduct, 
	                                                    // 10 : isDirectHardProcessTauDecayProduct, 11 : fromHardProcessBeforeFSR, 
	                                                    // 12 : isFirstCopy, 13 : isLastCopy, 14 : isLastCopyBeforeFSR, 
	  listGenP_-.push_back(gp);
	}
      }


      //---------------------------
      //Assigning GenJet Properties
      //---------------------------
      if (isMC()){
	//***nGenJet : slimmedGenJets, i.e. ak4 Jets made with visible genparticles***
	for (size_t i = 0; i < nGenJet; ++i){
	  vhtm::GenJet gj;
	  gj.eta            = GenJet_eta[i]; 
	  gj.mass           = GenJet_mass[i];
	  gj.phi            = GenJet_phi[i];
	  gj.pt             = GenJet_pt[i];
	  gj.hadFlv         = GenJet_hadronFlavour[i];
	  gj.parFlv         = GenJet_partonFlavour[i];
	  
	  listGenJet_->push_back(gj);
	}
      }

      //------------------------
      //Assigning MET Properties
      //------------------------
      vhtm::MET m;
      //corrmet......I think so
      m.met                    = MET_pt;
      m.metphi                 = MET_phi;
      m.sumet                  = MET_sumEt;
      m.metSignificance        = MET_significance;
      m.covXX                  = MET_covXX;
      m.covXY                  = MET_covXY;
      m.covYY                  = MET_covYY;
      m.MetUnclustEnUpDeltaX   = MET_MetUnclustEnUpDeltaX; //Delta (METx_mod-METx) Unclustered Energy Up
      m.MetUnclustEnUpDeltaY   = MET_MetUnclustEnUpDeltaY; //Delta (METy_mod-METy) Unclustered Energy Up
      listMET_->push_back(m);

      vhtm::MET mp;
      //puppi
      mp.met             = PuppiMET_pt;
      mp.phi             = PuppiMET_phi;
      mp.sumet           = PuppiMET_sumEt;
      listPMET_->push_back(mp);

      vhtm::MET mr;
      //raw
      mr.met             = RawMET_pt;
      mr.phi             = RawMET_phi;
      mr.sumet           = RawMET_sumEt;
      listRMET_->push_back(mr);


      //------------------------
      //Assigning Tau Properties
      //------------------------
      //***nTau : slimmedTaus after basic selection (pt > 18 && tauID('decayModeFindingNewDMs') && (tauID('byLooseCombinedIsolationDeltaBetaCorr3Hits') || tauID('byVLooseIsolationMVArun2v1DBoldDMwLT') || tauID('byVLooseIsolationMVArun2v1DBnewDMwLT') || tauID('byVLooseIsolationMVArun2v1DBdR03oldDMwLT') || tauID('byVVLooseIsolationMVArun2v1DBoldDMwLT2017v1') || tauID('byVVLooseIsolationMVArun2v1DBoldDMwLT2017v2') || tauID('byVVLooseIsolationMVArun2v1DBnewDMwLT2017v2') || tauID('byVVLooseIsolationMVArun2v1DBdR03oldDMwLT2017v2')))
      for (size_t i = 0; i < nTau; ++i){
	vhtm::Tau ta;
	//p4 related
	ta.eta    = Tau_eta[i];
	ta.phi    = Tau_phi[i];
	ta.pt     = Tau_pt[i];
	ta.charge = Tau_charge[i];
	ta.mass   = Tau_mass[i];

	//impact parameters
	ta.dxyPV     = Tau_dxy[i]; //d_{xy} of lead track with respect to PV, in cm (with sign)
	ta.dzPV      = Tau_dz[i];  //d_{z} of lead track with respect to PV, in cm (with sign)

	ta.chargedIso        = Tau_chargedIso[i];
	ta.neutralIso        = Tau_neutralIso[i];
	ta.rawIso            = Tau_rawIso[i];
	ta.rawIsodR03        = Tau_rawIsodR03[i];
	ta.rawMVAnewDM       = Tau_rawMVAnewDM[i];
	ta.rawMVAoldDM       = Tau_rawMVAoldDM[i];
	ta.rawMVAoldDMdR03   = Tau_rawMVAoldDMdR03[i];
	

	ta.idAntiEle         = Tau_idAntiEle[i];
	ta.idAntiMu          = Tau_idAntiMu[i];

	ta.decayModeFinding       = Tau_idDecayMode[i];
	ta.decayModeFindingNewDMs = Tau_idDecayModeNewDMs[i];

	ta.idMVAnew       = Tau_idMVAnew[i];
	ta.idMVAoldDM     = Tau_idMVAoldDM[i];
	ta.idMVAoldDMdR03 = Tau_idMVAoldDMdR03[i];

	ta.leadTkDeltaEta      = Tau_leadTkDeltaEta[i];
	ta.leadTkDeltaPhi      = Tau_leadTkDeltaPhi[i];
	ta.leadTkPtOverTauPt   = Tau_leadTkPtOverTauPt[i];

	ta.puCorr    = Tau_puCorr[i];
	ta.decayMode = Tau_decayMode[i];
	//** Incomplete
	if (isMC()) {
	  ta.genIdx   = Tau_genPartIdx[i]; // Index into genParticle list for MC matching to status==2 taus
	  ta.genFlv   = Tau_genPartFlav[i]; // Flavour of genParticle for MC matching to status==2 taus: 1 = prompt electron, 2 = prompt muon, 3 = tau->e decay, 4 = tau->mu decay, 5 = hadronic tau decay, 0 = unknown or unmatched
	}
	else {
	  ta.genIdx   = -999;
	  ta.genFlv   = -999;
	}
	ta.jetIdx   = Tau_jetIdx[i];

	listTau_ -> push_back(ta);
	  
      }

      //-------------------------
      //Assigning Muon Properties
      //-------------------------
      
      //***nMuon :: slimmedMuons after basic selection (pt > 3 && track.isNonnull && isLooseMuon)***  
      for (size_t i = 0; i < nMuon; ++i){
	vhtm::Muon mu;
	//IDs
	mu.isPFMuon             = Muon_isPFcand[i];
	mu.isGlobalMuon         = false; //Muon_isGlobal[i] 
	mu.isTrackerMuon        = false; //Muon_isTracker[i] 
	mu.isLooseMuon          = Muon_softId[i];
	mu.isMediumMuon         = Muon_mediumId[i]; //cut-based ID, medium WP
	mu.isTightMuon          = Muon_tightId[i];  //cut-based ID, tight WP
	mu.highPtId             = Muon_highPtId[i]; //high-pT cut-based ID (1 = tracker high pT, 2 = global high pT, which includes tracker high pT)
	mu.tightCharge          = Muon_tightCharge[i]; //Tight charge criterion using pterr/pt of muonBestTrack (0:fail, 2:pass)
	mu.pdgId                = Muon_pdgId[i]; //PDG code assigned by the event reconstruction (not by MC truth)
	mu.mvaTTH               = Muon_mvaTTH[i]; //TTH MVA lepton ID score
	
	//P4 related
	mu.mass   = Muon_mass[i];
	mu.pt     = Muon_pt[i]; 
	mu.ptErr  = Muon_ptErr[i]; //ptError of the muon track
	mu.eta    = Muon_eta[i]; 
	mu.phi    = Muon_phi[i];
	mu.charge = Muon_charge[i];

	//impact parameters
	mu.dxyPV    = Muon_dxy[i]; //dxy (with sign) wrt first PV, in cm
	mu.dxyPVerr = Muon_dxyErr[i];
	mu.dzPV     = Muon_dz[i]; //dz (with sign) wrt first PV, in cm
	mu.dzPVerr  = Muon_dzErr[i];
	mu.dB3D     = Muon_ip3d[i]; //3D impact parameter wrt first PV, in cm
	mu.SIP3D    = Muon_sip3d[i]; //3D impact parameter significance wrt first PV
	
	//isolation variables
	mu.pfRelIso03    = Muon_pfRelIso03_all[i]; //PF relative isolation dR=0.3, total (deltaBeta corrections)
	mu.pfRelChrIso03 = Muon_pfRelIso03_chg[i]; //PF relative isolation dR=0.3, charged component
	mu.pfRelIso04    = Muon_pfRelIso04_all[i]; //PF relative isolation dR=0.4, total (deltaBeta corrections)

	//References
	mu.jetIdx     = Muon_jetIdx[i]; //(index to Jet)index of the associated jet (-1 if none)
	if (isMC()) {
	  mu.genIdx   = Muon_genPartIdx[i]; // Index into genParticle list for MC matching to status==1 muons
	  mu.genFlv   = Muon_genPartFlav[i]; // Flavour of genParticle for MC matching to status==1 
	                                         // muons: 1 = prompt muon (including gamma*->mu mu), 15 = muon from prompt tau, 
	                                         // 5 = muon from b, 4 = muon from c, 3 = muon from light or unknown, 0 = unmatched
	}
	else {
	  mu.genIdx   = -999;
	  mu.genFlv   = -999;
	}

	//Track History
	mu.nMatchedStations      = Muon_nStations[i]; //number of matched stations with default arbitration (segment & track)
	mu.nMatchedTrkLayers     = Muon_nTrackerLayers[i]; //number of layers in the tracker

	listMu_->push_back(mu);
      }

      //-----------------------------
      //Assigning Electron properties
      //-----------------------------

      //***nElectron :slimmedElectrons after basic selection (pt > 5 )***
      for (size_t i = 0; i < nElectron; ++i){

	vhtm::Electron el;

	//p4 related
	el.pt     = Electron_pt[i]; 
	el.eta    = Electron_eta[i];
	el.phi    = Electron_phi[i];
	el.charge = Electron_charge[i];
	el.mass   = Electron_mass[i];	

	//impact parameters
	el.dxyPV    = Electron_dxy[i]; //dxy (with sign) wrt first PV, in cm
	el.dxyPVerr = Electron_dxyErr[i];
	el.dzPV     = Electron_dz[i]; //dz (with sign) wrt first PV, in cm
	el.dzPVerr  = Electron_dzErr[i];
	el.dB3D     = Electron_ip3d[i]; //3D impact parameter wrt first PV, in cm
	el.SIP3D    = Electron_sip3d[i]; //3D impact parameter significance wrt first PV

	//IDs
	el.cutBased           = Electron_cutBased[i]; //cut-based ID (0:fail, 1:veto, 2:loose, 3:medium, 4:tight)
	el.cutBasedHEEP       = Electron_cutBased_HEEP[i]; //cut-based HEEP ID
	el.cutBasedHLTPreSel  = Electron_cutBased_HLTPreSel[i]; // cut-based HLT pre-selection ID
	el.vidNestedWPBitmap  = Electron_vidNestedWPBitmap[i];  // VID compressed bitmap 
	                                                        // (MinPtCut,GsfEleSCEtaMultiRangeCut,GsfEleDEtaInSeedCut,GsfEleDPhiInCut,
	                                                        // GsfEleFull5x5SigmaIEtaIEtaCut,GsfEleHadronicOverEMCut,GsfEleEInverseMinusPInverseCut,
	                                                        // GsfEleEffAreaPFIsoCut,GsfEleConversionVetoCut,GsfEleMissingHitsCut), 2 bits per cut
	el.tightCharge        = Electron_tightCharge[i]; //Tight charge criteria (0:none, 1:isGsfScPixChargeConsistent, 2:isGsfCtfScPixChargeConsistent)
	el.pdgId              = Electron_pdgId[i]; //PDG code assigned by the event reconstruction (not by MC truth)
	el.mvaTTH             = Electron_mvaTTH[i]; //TTH MVA lepton ID score
	el.mvaSpring16GP      = Electron_mvaSpring16GP[i]; //MVA general-purpose ID score //For 2016_CMSSW_80X
	el.mvaSpring16GP_WP80 = Electron_mvaSpring16GP_WP80[i]; //MVA general-purpose ID WP80
	el.mvaSpring16GP_WP90 = Electron_mvaSpring16GP_WP90[i]; //MVA general-purpose ID WP90
	el.mvaSpring16HZZ     = Electron_mvaSpring16HZZ[i]; //MVA HZZ ID score
	el.mvaSpring16HZZ_WPL = Electron_mvaSpring16HZZ_WPL[i]; //MVA HZZ ID loose WP

	
	//References
	el.photonIdx    = Electron_photonIdx[i]; //index of the associated photon (-1 if none)
	if (isMC()) {
	  el.genIdx   = Electron_genPartIdx[i];  // Index into genParticle list for MC matching to status==1 muons
	  el.genFlv   = Electron_genPartFlav[i]; // Flavour of genParticle for MC matching to status==1 
	                                             // electrons or photons: 1 = prompt electron (including gamma*->mu mu), 
	                                             // 15 = electron from prompt tau, 22 = prompt photon (likely conversion), 
	                                             // 5 = electron from b, 4 = electron from c, 3 = electron from light or unknown, 0 = unmatched 
	}
	else {
	  el.genIdx   = -999;
	  el.genFlv   = -999;
	}
	

	//Isolation variables
	el.pfRelIso03    = Electron_pfRelIso03_all[i]; //PF relative isolation dR=0.3, total (with rho*EA PU corrections)
	el.pfRelChrIso03 = Electron_pfRelIso03_chg[i]; //PF relative isolation dR=0.3, charged component

	
	el.deltaEtaTrkSC     = Electron_deltaEtaSC[i]; //delta eta (SC,ele) with sign
	el.sigmaIEtaIEta     = Electron_sieie[i]; //sigma_IetaIeta of the supercluster, calculated with full 5x5 region
	el.hoe               = Electron_hoe[i]; //H over E

	listEle_->push_back(el);
      } //Electron loop ends
      
      //------------------------
      //Assinging Jet Properties
      //------------------------

      //***nJet : slimmedJets, i.e. ak4 PFJets CHS with JECs applied, after basic selection (pt > 15)***
      for (size_t i = 0; i < nJet; ++i){

	vhtm::Jet j;
	//p4 related
	j.eta    = Jet_eta[i];
	j.phi    = Jet_phi[i];
	j.pt     = Jet_pt[i];
	j.mass   = Jet_mass[i];

	j.area                         = Jet_area[i]; //jet catchment area, for JECs
	j.chargedEmEnergyFraction      = Jet_chEmEF[i]; //charged Electromagnetic Energy Fraction
	j.chargedHadronEnergyFraction  = Jet_chHEF[i];  //charged Hadron Energy Fraction
	j.neutralEmEnergyFraction      = Jet_neEmEF[i]; //neutral Electromagnetic Energy Fraction
	j.neutralHadronEnergyFraction  = Jet_neHEF[i];  //neutral Hadron Energy Fraction
	j.muonEnergyFraction           = Jet_muEF[i];
	j.electronEnergyFraction       = -999.9; //?????????????????????????????????

	j.electronMultiplicity    = Jet_nElectrons[i]; //number of electrons in the jet
	j.muonMultiplicity        = Jet_nMuons[i];     //number of muons in the jet
	j.nConstituents           = Jet_nConstituents[i]; //Number of particles in the jet
	
	//b-tagging vars
	j.btagCMVA       = Jet_btagCMVA[i];  // CMVA V2 btag discriminator
	j.btagCSVV2      = Jet_btagCSVV2[i]; // pfCombinedInclusiveSecondaryVertexV2 b-tag discriminator (aka CSVV2)
	                                                                             // oldName (inOldTreeMaker): pfCombinedInclusiveSecondaryVertexV2BJetTags
	j.btagDeepB      = Jet_btagDeepB[i]; //DeepCSV b+bb tag discriminator
	j.btagDeepC      = Jet_btagDeepC[i]; //DeepCSV charm btag discriminator
	j.btagDeepFlavB  = Jet_btagDeepFlavB[i]; //DeepFlavour b+bb+lepb tag discriminator

	//references
	if (isMC()) {
	  j.partonFlavour    = Jet_partonFlavour[i]; //flavour from parton matching
	  j.genJetIdx        = Jet_genJetIdx[i];     //index of matched gen jet
	}
	else {
	  j.partonFlavour    = -999;
	  j.genJetIdx        = -999;
	}
	j.electronIdx1    = Jet_electronIdx1[i]; //index of first matching electron
	j.electronIdx2    = Jet_electronIdx2[i]; //index of second matching electron
	j.muonIdx1        = Jet_muonIdx1[i]; //index of first matching muon
	j.muonIdx2        = Jet_muonIdx2[i]; //index of second matching muon

	j.puId    = Jet_puId[i]; //Pilup ID flags
	j.qgl   = Jet_qgl[i];  //Quark vs Gluon likelihood discriminator
	
	listJet_->push_back(j);
      }//Jet Loop finished
   
      //-----------------------------------
      //Assigning Primary Vertex Properties
      //-----------------------------------
      vhtm::Vertex v;
      v.x        = PV_x;
      v.y        = PV_y;
      v.z        = PV_z;
      v.ndf      = PV_ndof;
      v.chi2     = PV_chi2;
      v.npvs     = PV_npvs; //total number of reconstructed primary vertices
      v.npvsGood = PV_npvsGood; //number of good reconstructed primary vertices. selection:!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2
      
      listPVtx_ -> push_back(v);
      
      //-------------------------------------
      //Assigning Secondary Vertex Properties
      //-------------------------------------

      vhtm::Vertex sv;
      for (size_t i = 0; i < nSV; ++i){
	sv.dlen      = SV_dlen[i]; // decay length in cm
	sv.dlenSig   = SV_dlenSig[i]; // decay length significance
	sv.pAngle    = SV_pAngle[i]; // pointing angle, i.e. acos(p_SV * (SV - PV)) 
	sv.chi2      = SV_chi2[i]; //reduced chi2, i.e. chi/ndof 
	sv.eta       = SV_eta[i];
	sv.ndf       = SV_ndof[i];
	sv.phi       = SV_phi[i];
	sv.pt        = SV_pt[i];
	sv.x         = SV_x[i];
	sv.y         = SV_y[i];
	sv.z         = SV_z[i];
	sv.mass      = SV_mass[i];

	listSVtx_->push_back(sv);
      }

      //---------------------------
      //Assigning GenMET Properties
      //---------------------------
      if (isMC()){
	vhtm::GenMET gm;
	
	gm.met    = GenMET_pt;
	gm.metphi = GenMET_phi;

	listGenMet_ -> push_back(gm);
      }

      //----------------------------------
      //Assigning TriggerObject Properties
      //----------------------------------
      for (size_t i = 0; i < nTrigObj; ++i){
	vhtm::TriggerObject to;
	to.pt          = TrigObj_pt[i];
	to.eta         = TrigObj_eta[i];
	to.phi         = TrigObj_phi[i];
	to.filterBits  = TrigObj_filterBits[i]; //extra bits of associated information: 1 = CaloIdL_TrackIdL_IsoVL, 2 = 1e (WPTight), 4 = 1e (WPLoose), 8 = OverlapFilter PFTau, 16 = 2e, 32 = 1e-1mu, 64 = 1e-1tau, 128 = 3e, 256 = 2e-1mu, 512 = 1e-2mu for Electron (PixelMatched e/gamma); 1 = TrkIsoVVL, 2 = Iso, 4 = OverlapFilter PFTau, 8 = 1mu, 16 = 2mu, 32 = 1mu-1e, 64 = 1mu-1tau, 128 = 3mu, 256 = 2mu-1e, 512 =1mu-2e for Muon; 1 = LooseChargedIso, 2 = MediumChargedIso, 4 = TightChargedIso, 8 = TightID OOSC photons, 16 = HPS, 32 = single-tau + tau+MET, 64 = di-tau, 128 = e-tau, 256 = mu-tau, 512 = VBF+di-tau for Tau; 1 = VBF cross-cleaned from loose iso PFTau for Jet; 
	to.Id          = TrigObj_id[i]; //ID of the object: 11 = Electron (PixelMatched e/gamma), 22 = Photon (PixelMatch-vetoed e/gamma), 13 = Muon, 15 = Tau, 1 = Jet, 6 = FatJet, 2 = MET, 3 = HT, 4 = MHT
	to.l1pt        = TrigObj_l1pt[i]; // pt of associated L1 seed
	to.l1pt_2      = TrigObj_l1pt_2[i]; // pt of associated secondary L1 seed
	to.l2pt        = TrigObj_l2pt[i]; // pt of associated 'L2' seed (i.e. HLT before tracking/PF)
	to.l1iso       = TrigObj_l1iso[i]; // iso of associated L1 seed
	to.l1charge    = TrigObj_l1charge[i]; // charge of associated L1 seed

	listTrgObj_->push_back(to);
      }


      //--------------------------------
      //Assigning Photons
      //--------------------------------
      for (size_t i = 0; i < nPhoton; ++i){
	vhtm::Photon ph;
	
	ph.et           = Photon_pt[i];
	ph.eta          = Photon_eta[i];
	ph.phi          = Photon_phi[i];
	ph.energyErr    = Photon_energyErr[i]; //energy error of the cluster from regression
	
	ph.r9            = Photon_r9[i];
	ph.hoe           = Photon_hoe[i];
	ph.sigmaIEtaIEta = Photon_sieie[i];
	
	ph.electronIdx   = Photon_electronIdx[i];
	ph.jetIdx        = Photon_jetIdx[i];
	
	ph.cutBased      = Photon_cutBased[i]; //cut-based ID (0:fail, 1::loose, 2:medium, 3:tight)
	ph.mvaID         = Photon_mvaID[i];
	ph.pfRelIso03    = Photon_pfRelIso03_all[i];
	ph.pfRelChrIso03 = Photon_pfRelIso03_chg[i];
 
	ph.pdgId        = Photon_pdgId[i];
	ph.vidNestedWPBitmap   = Photon_vidNestedWPBitmap[i]; //(MinPtCut,PhoSCEtaMultiRangeCut,PhoSingleTowerHadOverEmCut,PhoFull5x5SigmaIEtaIEtaCut,PhoAnyPFIsoWithEACut,PhoAnyPFIsoWithEAAndQuadScalingCut,PhoAnyPFIsoWithEACut), 2 bits per cut
	ph.electronVeto     = Photon_electronVeto[i]; //pass electron veto
	ph.isScEtaEB        = Photon_isScEtaEB[i]; //is supercluster eta within barrel acceptance
	ph.isScEtaEE        = Photon_isScEtaEE[i]; //is supercluster eta within endcap acceptance
	ph.mvaID_WP80       = Photon_mvaID_WP80[i];
	ph.mvaID_WP90       = Photon_mvaID_WP90[i];
	ph.pixelSeed        = Photon_pixelSeed[i]; //has pixel seed

	if (isMC()) {
	  ph.genIdx   = -999; //Photon_genPartIdx[i]; 
	  ph.genFlv   = -999; //Photon_genPartFlav[i]; 
	}
	else {
	  ph.genIdx   = -999;
	  ph.genFlv   = -999;
	}
	

      }

      //------------
      //HLT
      //------------
      hltpaths_ {
	"HLTriggerFirstPath","HLT_AK8PFJet360_TrimMass30","HLT_AK8PFHT700_TrimR0p1PT0p03Mass50","HLT_AK8PFHT650_TrimR0p1PT0p03Mass50","HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20","HLT_CaloJet500_NoJetID","HLT_Dimuon13_PsiPrime","HLT_Dimuon13_Upsilon","HLT_Dimuon20_Jpsi","HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf","HLT_DoubleEle33_CaloIdL","HLT_DoubleEle33_CaloIdL_MW","HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW","HLT_DoubleEle33_CaloIdL_GsfTrkIdVL","HLT_DoubleEle37_Ele27_CaloIdL_GsfTrkIdVL","HLT_DoubleMediumIsoPFTau32_Trk1_eta2p1_Reg","HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg","HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg","HLT_DoubleMu33NoFiltersNoVtx","HLT_DoubleMu38NoFiltersNoVtx","HLT_DoubleMu23NoFiltersNoVtxDisplaced","HLT_DoubleMu28NoFiltersNoVtxDisplaced","HLT_DoubleMu4_3_Bs","HLT_DoubleMu4_3_Jpsi_Displaced","HLT_DoubleMu4_JpsiTrk_Displaced","HLT_DoubleMu4_LowMassNonResonantTrk_Displaced","HLT_DoubleMu3_Trk_Tau3mu","HLT_DoubleMu4_PsiPrimeTrk_Displaced","HLT_Mu7p5_L2Mu2_Jpsi","HLT_Mu7p5_L2Mu2_Upsilon","HLT_Mu7p5_Track2_Jpsi","HLT_Mu7p5_Track3p5_Jpsi","HLT_Mu7p5_Track7_Jpsi","HLT_Mu7p5_Track2_Upsilon","HLT_Mu7p5_Track3p5_Upsilon","HLT_Mu7p5_Track7_Upsilon",	"HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing","HLT_Dimuon0er16_Jpsi_NoVertexing","HLT_Dimuon6_Jpsi_NoVertexing","HLT_DoublePhoton60",	"HLT_DoublePhoton85","HLT_Ele22_eta2p1_WPLoose_Gsf","HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1","HLT_Ele23_WPLoose_Gsf","HLT_Ele23_WPLoose_Gsf_WHbbBoost","HLT_Ele24_eta2p1_WPLoose_Gsf","HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20","HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1","HLT_Ele25_WPTight_Gsf","HLT_Ele25_eta2p1_WPLoose_Gsf","HLT_Ele25_eta2p1_WPTight_Gsf",
	"HLT_Ele27_WPLoose_Gsf",
	"HLT_Ele27_WPLoose_Gsf_WHbbBoost",
	"HLT_Ele27_WPTight_Gsf",
	"HLT_Ele27_eta2p1_WPLoose_Gsf",
	"HLT_Ele27_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1",
	"HLT_Ele27_eta2p1_WPLoose_Gsf_DoubleMediumIsoPFTau32_Trk1_eta2p1_Reg",
	"HLT_Ele27_eta2p1_WPLoose_Gsf_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg",
	"HLT_Ele27_eta2p1_WPLoose_Gsf_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg",
	"HLT_Ele27_eta2p1_WPTight_Gsf",
	"HLT_Ele32_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1",
	"HLT_Ele32_eta2p1_WPTight_Gsf",
	"HLT_Ele35_WPLoose_Gsf",
	"HLT_Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50",
	"HLT_Ele45_WPLoose_Gsf",
	"HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50",
	"HLT_Ele105_CaloIdVT_GsfTrkIdT",
	"HLT_Ele30WP60_SC4_Mass55",
	"HLT_Ele30WP60_Ele8_Mass55",
	"HLT_HT200",
	"HLT_HT275",
	"HLT_HT325",
	"HLT_HT425",
	"HLT_HT575",
	"HLT_HT410to430",
	"HLT_HT430to450",
	"HLT_HT450to470",
	"HLT_HT470to500",
	"HLT_HT500to550",
	"HLT_HT550to650",
	"HLT_HT650",
	"HLT_Mu16_eta2p1_MET30",
	"HLT_IsoMu16_eta2p1_MET30",
	"HLT_IsoMu16_eta2p1_MET30_LooseIsoPFTau50_Trk30_eta2p1",
	"HLT_IsoMu17_eta2p1_LooseIsoPFTau20",
	"HLT_IsoMu17_eta2p1_LooseIsoPFTau20_SingleL1",
	"HLT_DoubleIsoMu17_eta2p1",
	"HLT_DoubleIsoMu17_eta2p1_noDzCut",
	"HLT_IsoMu18",
	"HLT_IsoMu19_eta2p1_LooseIsoPFTau20",
	"HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1",
	"HLT_IsoMu19_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg",
	"HLT_IsoMu20",
	"HLT_IsoMu21_eta2p1_LooseIsoPFTau20_SingleL1",
	"HLT_IsoMu21_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg",
	"HLT_IsoMu22",
	"HLT_IsoMu22_eta2p1",
	"HLT_IsoMu24",
	"HLT_IsoMu27",
	"HLT_IsoTkMu18",
	"HLT_IsoTkMu20",
	"HLT_IsoTkMu22",
	"HLT_IsoTkMu22_eta2p1",
	"HLT_IsoTkMu24",
	"HLT_IsoTkMu27",
	"HLT_JetE30_NoBPTX3BX",
	"HLT_JetE30_NoBPTX",
	"HLT_JetE50_NoBPTX3BX",
	"HLT_JetE70_NoBPTX3BX",
	"HLT_L1SingleMu18",
	"HLT_L2Mu10",
	"HLT_L1SingleMuOpen",
	"HLT_L1SingleMuOpen_DT",
	"HLT_L2DoubleMu23_NoVertex",
	"HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10",
	"HLT_L2DoubleMu38_NoVertex_2Cha_Angle2p5_Mass10",
	"HLT_L2Mu10_NoVertex_NoBPTX3BX",
	"HLT_L2Mu10_NoVertex_NoBPTX",
	"HLT_L2Mu35_NoVertex_3Sta_NoBPTX3BX",
	"HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX",
	"HLT_LooseIsoPFTau50_Trk30_eta2p1",
	"HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80",
	"HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90",
	"HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110",
	"HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120",
	"HLT_PFTau120_eta2p1",
	"HLT_VLooseIsoPFTau120_Trk50_eta2p1",
	"HLT_VLooseIsoPFTau140_Trk50_eta2p1",
	"HLT_Mu17_Mu8",
	"HLT_Mu17_Mu8_DZ",
	"HLT_Mu17_Mu8_SameSign",
	"HLT_Mu17_Mu8_SameSign_DZ",
	"HLT_Mu20_Mu10",
	"HLT_Mu20_Mu10_DZ",
	"HLT_Mu20_Mu10_SameSign",
	"HLT_Mu20_Mu10_SameSign_DZ",
	"HLT_Mu17_TkMu8_DZ",
	"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL",
	"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
	"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL",
	"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ",
	"HLT_Mu25_TkMu0_dEta18_Onia",
	"HLT_Mu27_TkMu8",
	"HLT_Mu30_TkMu11",
	"HLT_Mu30_eta2p1_PFJet150_PFJet50",
	"HLT_Mu40_TkMu11",
	"HLT_Mu40_eta2p1_PFJet200_PFJet50",
	"HLT_Mu20",
	"HLT_TkMu20",
	"HLT_Mu24_eta2p1",
	"HLT_TkMu24_eta2p1",
	"HLT_Mu27",
	"HLT_TkMu27",
	"HLT_Mu45_eta2p1",
	"HLT_Mu50",
	"HLT_TkMu50",
	"HLT_Mu38NoFiltersNoVtx_Photon38_CaloIdL",
	"HLT_Mu42NoFiltersNoVtx_Photon42_CaloIdL",
	"HLT_Mu28NoFiltersNoVtxDisplaced_Photon28_CaloIdL",
	"HLT_Mu33NoFiltersNoVtxDisplaced_Photon33_CaloIdL",
	"HLT_Mu23NoFiltersNoVtx_Photon23_CaloIdL",
	"HLT_DoubleMu18NoFiltersNoVtx",
	"HLT_Mu33NoFiltersNoVtxDisplaced_DisplacedJet50_Tight",
	"HLT_Mu33NoFiltersNoVtxDisplaced_DisplacedJet50_Loose",
	"HLT_Mu28NoFiltersNoVtx_DisplacedJet40_Loose",
	"HLT_Mu38NoFiltersNoVtxDisplaced_DisplacedJet60_Tight",
	"HLT_Mu38NoFiltersNoVtxDisplaced_DisplacedJet60_Loose",
	"HLT_Mu38NoFiltersNoVtx_DisplacedJet60_Loose",
	"HLT_Mu28NoFiltersNoVtx_CentralCaloJet40",
	"HLT_PFHT300_PFMET100",
	"HLT_PFHT300_PFMET110",
	"HLT_PFHT550_4JetPt50",
	"HLT_PFHT650_4JetPt50",
	"HLT_PFHT750_4JetPt50",
	"HLT_PFJet15_NoCaloMatched",
	"HLT_PFJet25_NoCaloMatched",
	"HLT_DiPFJet15_NoCaloMatched",
	"HLT_DiPFJet25_NoCaloMatched",
	"HLT_DiPFJet15_FBEta3_NoCaloMatched",
	"HLT_DiPFJet25_FBEta3_NoCaloMatched",
	"HLT_DiPFJetAve15_HFJEC",
	"HLT_DiPFJetAve25_HFJEC",
	"HLT_DiPFJetAve35_HFJEC",
	"HLT_AK8PFJet40",
	"HLT_AK8PFJet60",
	"HLT_AK8PFJet80",
	"HLT_AK8PFJet140",
	"HLT_AK8PFJet200",
	"HLT_AK8PFJet260",
	"HLT_AK8PFJet320",
	"HLT_AK8PFJet400",
	"HLT_AK8PFJet450",
	"HLT_AK8PFJet500",
	"HLT_PFJet40",
	"HLT_PFJet60",
	"HLT_PFJet80",
	"HLT_PFJet140",
	"HLT_PFJet200",
	"HLT_PFJet260",
	"HLT_PFJet320",
	"HLT_PFJet400",
	"HLT_PFJet450",
	"HLT_PFJet500",
	"HLT_DiPFJetAve40",
	"HLT_DiPFJetAve60",
	"HLT_DiPFJetAve80",
	"HLT_DiPFJetAve140",
	"HLT_DiPFJetAve200",
	"HLT_DiPFJetAve260",
	"HLT_DiPFJetAve320",
	"HLT_DiPFJetAve400",
	"HLT_DiPFJetAve500",
	"HLT_DiPFJetAve60_HFJEC",
	"HLT_DiPFJetAve80_HFJEC",
	"HLT_DiPFJetAve100_HFJEC",
	"HLT_DiPFJetAve160_HFJEC",
	"HLT_DiPFJetAve220_HFJEC",
	"HLT_DiPFJetAve300_HFJEC",
	"HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu140",
	"HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu80",
	"HLT_DiCentralPFJet55_PFMET110",
	"HLT_DiCentralPFJet170",
	"HLT_SingleCentralPFJet170_CFMax0p1",
	"HLT_DiCentralPFJet170_CFMax0p1",
	"HLT_DiCentralPFJet220_CFMax0p3",
	"HLT_DiCentralPFJet330_CFMax0p5",
	"HLT_DiCentralPFJet430",
	"HLT_PFHT125",
	"HLT_PFHT200",
	"HLT_PFHT250",
	"HLT_PFHT300",
	"HLT_PFHT350",
	"HLT_PFHT400",
	"HLT_PFHT475",
	"HLT_PFHT600",
	"HLT_PFHT650",
	"HLT_PFHT800",
	"HLT_PFHT900",
	"HLT_PFHT200_PFAlphaT0p51",
	"HLT_PFHT200_DiPFJetAve90_PFAlphaT0p57",
	"HLT_PFHT200_DiPFJetAve90_PFAlphaT0p63",
	"HLT_PFHT250_DiPFJetAve90_PFAlphaT0p55",
	"HLT_PFHT250_DiPFJetAve90_PFAlphaT0p58",
	"HLT_PFHT300_DiPFJetAve90_PFAlphaT0p53",
	"HLT_PFHT300_DiPFJetAve90_PFAlphaT0p54",
	"HLT_PFHT350_DiPFJetAve90_PFAlphaT0p52",
	"HLT_PFHT350_DiPFJetAve90_PFAlphaT0p53",
	"HLT_PFHT400_DiPFJetAve90_PFAlphaT0p51",
	"HLT_PFHT400_DiPFJetAve90_PFAlphaT0p52",
	"HLT_MET60_IsoTrk35_Loose",
	"HLT_MET75_IsoTrk50",
	"HLT_MET90_IsoTrk50",
	"HLT_PFMET120_BTagCSV_p067",
	"HLT_PFMET120_Mu5",
	"HLT_PFMET170_NoiseCleaned",
	"HLT_PFMET170_HBHECleaned",
	"HLT_PFMET170_JetIdCleaned",
	"HLT_PFMET170_NotCleaned",
	"HLT_PFMET170_BeamHaloCleaned",
	"HLT_PFMET90_PFMHT90_IDTight",
	"HLT_PFMET100_PFMHT100_IDTight",
	"HLT_PFMET110_PFMHT110_IDTight",
	"HLT_PFMET120_PFMHT120_IDTight",
	"HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight_BTagCSV_p067",
	"HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight",
	"HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq200",
	"HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq460",
	"HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq240",
	"HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq500",
	"HLT_QuadPFJet_VBF",
	"HLT_L1_TripleJet_VBF",
	"HLT_QuadJet45_TripleBTagCSV_p087",
	"HLT_QuadJet45_DoubleBTagCSV_p087",
	"HLT_DoubleJet90_Double30_TripleBTagCSV_p087",
	"HLT_DoubleJet90_Double30_DoubleBTagCSV_p087",
	"HLT_DoubleJetsC100_DoubleBTagCSV_p026_DoublePFJetsC160",
	"HLT_DoubleJetsC100_DoubleBTagCSV_p014_DoublePFJetsC100MaxDeta1p6",
	"HLT_DoubleJetsC112_DoubleBTagCSV_p026_DoublePFJetsC172",
	"HLT_DoubleJetsC112_DoubleBTagCSV_p014_DoublePFJetsC112MaxDeta1p6",
	"HLT_DoubleJetsC100_SingleBTagCSV_p026",
	"HLT_DoubleJetsC100_SingleBTagCSV_p014",
	"HLT_DoubleJetsC100_SingleBTagCSV_p026_SinglePFJetC350",
	"HLT_DoubleJetsC100_SingleBTagCSV_p014_SinglePFJetC350",
	"HLT_Photon135_PFMET100",
	"HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_PFMET40",
	"HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_VBF",
	"HLT_Photon250_NoHE",
	"HLT_Photon300_NoHE",
	"HLT_Photon26_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon16_AND_HE10_R9Id65_Eta2_Mass60",
	"HLT_Photon36_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon22_AND_HE10_R9Id65_Eta2_Mass15",
	"HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_PFMET40",
	"HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_VBF",
	"HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_PFMET40",
	"HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_VBF",
	"HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_PFMET40",
	"HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_VBF",
	"HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_PFMET40",
	"HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_VBF",
	"HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_PFMET40",
	"HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_VBF",
	"HLT_Mu8_TrkIsoVVL",
	"HLT_Mu17_TrkIsoVVL",
	"HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30",
	"HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30",
	"HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30",
	"HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30",
	"HLT_BTagMu_DiJet20_Mu5",
	"HLT_BTagMu_DiJet40_Mu5",
	"HLT_BTagMu_DiJet70_Mu5",
	"HLT_BTagMu_DiJet110_Mu5",
	"HLT_BTagMu_DiJet170_Mu5",
	"HLT_BTagMu_Jet300_Mu5",
	"HLT_BTagMu_AK8Jet300_Mu5",
	"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
	"HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
	"HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL",
	"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
	"HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL",
	"HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL",
	"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
	"HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
	"HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL",
	"HLT_Mu37_Ele27_CaloIdL_GsfTrkIdVL",
	"HLT_Mu27_Ele37_CaloIdL_GsfTrkIdVL",
	"HLT_Mu8_DiEle12_CaloIdL_TrackIdL",
	"HLT_Mu12_Photon25_CaloIdL",
	"HLT_Mu12_Photon25_CaloIdL_L1ISO",
	"HLT_Mu12_Photon25_CaloIdL_L1OR",
	"HLT_Mu17_Photon22_CaloIdL_L1ISO",
	"HLT_Mu17_Photon30_CaloIdL_L1ISO",
	"HLT_Mu17_Photon35_CaloIdL_L1ISO",
	"HLT_DiMu9_Ele9_CaloIdL_TrackIdL",
	"HLT_TripleMu_5_3_3",
	"HLT_TripleMu_12_10_5",
	"HLT_Mu3er_PFHT140_PFMET125",
	"HLT_Mu6_PFHT200_PFMET80_BTagCSV_p067",
	"HLT_Mu6_PFHT200_PFMET100",
	"HLT_Mu14er_PFMET100",
	"HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL",
	"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
	"HLT_Ele12_CaloIdL_TrackIdL_IsoVL",
	"HLT_Ele17_CaloIdL_GsfTrkIdVL",
	"HLT_Ele17_CaloIdL_TrackIdL_IsoVL",
	"HLT_Ele23_CaloIdL_TrackIdL_IsoVL",
	"HLT_AK8DiPFJet280_200_TrimMass30",
	"HLT_AK8DiPFJet250_200_TrimMass30",
	"HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20",
	"HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20",
	"HLT_PFHT650_WideJetMJJ900DEtaJJ1p5",
	"HLT_PFHT650_WideJetMJJ950DEtaJJ1p5",
	"HLT_Photon22",
	"HLT_Photon30",
	"HLT_Photon36",
	"HLT_Photon50",
	"HLT_Photon75",
	"HLT_Photon90",
	"HLT_Photon120",
	"HLT_Photon175",
	"HLT_Photon165_HE10",
	"HLT_Photon22_R9Id90_HE10_IsoM",
	"HLT_Photon30_R9Id90_HE10_IsoM",
	"HLT_Photon36_R9Id90_HE10_IsoM",
	"HLT_Photon50_R9Id90_HE10_IsoM",
	"HLT_Photon75_R9Id90_HE10_IsoM",
	"HLT_Photon90_R9Id90_HE10_IsoM",
	"HLT_Photon120_R9Id90_HE10_IsoM",
	"HLT_Photon165_R9Id90_HE10_IsoM",
	"HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90",
	"HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelSeedMatch_Mass70",
	"HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55",
	"HLT_Diphoton30_18_Solid_R9Id_AND_IsoCaloId_AND_HE_R9Id_Mass55",
	"HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55",
	"HLT_Dimuon0_Jpsi_Muon",
	"HLT_Dimuon0_Upsilon_Muon",
	"HLT_QuadMuon0_Dimuon0_Jpsi",
	"HLT_QuadMuon0_Dimuon0_Upsilon",
	"HLT_Rsq0p25",
	"HLT_Rsq0p30",
	"HLT_RsqMR240_Rsq0p09_MR200",
	"HLT_RsqMR240_Rsq0p09_MR200_4jet",
	"HLT_RsqMR270_Rsq0p09_MR200",
	"HLT_RsqMR270_Rsq0p09_MR200_4jet",
	"HLT_Rsq0p02_MR300_TriPFJet80_60_40_BTagCSV_p063_p20_Mbb60_200",
	"HLT_Rsq0p02_MR300_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200",
	"HLT_HT200_DisplacedDijet40_DisplacedTrack",
	"HLT_HT250_DisplacedDijet40_DisplacedTrack",
	"HLT_HT350_DisplacedDijet40_DisplacedTrack",
	"HLT_HT350_DisplacedDijet80_DisplacedTrack",
	"HLT_HT350_DisplacedDijet80_Tight_DisplacedTrack",
	"HLT_HT350_DisplacedDijet40_Inclusive",
	"HLT_HT400_DisplacedDijet40_Inclusive",
	"HLT_HT500_DisplacedDijet40_Inclusive",
	"HLT_HT550_DisplacedDijet40_Inclusive",
	"HLT_HT650_DisplacedDijet80_Inclusive",
	"HLT_HT750_DisplacedDijet80_Inclusive",
	"HLT_VBF_DisplacedJet40_DisplacedTrack",
	"HLT_VBF_DisplacedJet40_DisplacedTrack_2TrackIP2DSig5",
	"HLT_VBF_DisplacedJet40_TightID_DisplacedTrack",
	"HLT_VBF_DisplacedJet40_Hadronic",
	"HLT_VBF_DisplacedJet40_Hadronic_2PromptTrack",
	"HLT_VBF_DisplacedJet40_TightID_Hadronic",
	"HLT_VBF_DisplacedJet40_VTightID_Hadronic",
	"HLT_VBF_DisplacedJet40_VVTightID_Hadronic",
	"HLT_VBF_DisplacedJet40_VTightID_DisplacedTrack",
	"HLT_VBF_DisplacedJet40_VVTightID_DisplacedTrack",
	"HLT_PFMETNoMu90_PFMHTNoMu90_IDTight",
	"HLT_PFMETNoMu100_PFMHTNoMu100_IDTight",
	"HLT_PFMETNoMu110_PFMHTNoMu110_IDTight",
	"HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",
	"HLT_MonoCentralPFJet80_PFMETNoMu90_PFMHTNoMu90_IDTight",
	"HLT_MonoCentralPFJet80_PFMETNoMu100_PFMHTNoMu100_IDTight",
	"HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight",
	"HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight",
	"HLT_Ele27_eta2p1_WPLoose_Gsf_HT200",
	"HLT_Photon90_CaloIdL_PFHT500",
	"HLT_DoubleMu8_Mass8_PFHT250",
	"HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT250",
	"HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT250",
	"HLT_DoubleMu8_Mass8_PFHT300",
	"HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300",
	"HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300",
	"HLT_Mu10_CentralPFJet30_BTagCSV_p13",
	"HLT_DoubleMu3_PFMET50",
	"HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV_p13",
	"HLT_Ele15_IsoVVVL_BTagCSV_p067_PFHT400",
	"HLT_Ele15_IsoVVVL_PFHT350_PFMET50",
	"HLT_Ele15_IsoVVVL_PFHT600",
	"HLT_Ele15_IsoVVVL_PFHT350",
	"HLT_Ele15_IsoVVVL_PFHT400_PFMET50",
	"HLT_Ele15_IsoVVVL_PFHT400",
	"HLT_Ele50_IsoVVVL_PFHT400",
	"HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60",
	"HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60",
	"HLT_Mu15_IsoVVVL_BTagCSV_p067_PFHT400",
	"HLT_Mu15_IsoVVVL_PFHT350_PFMET50",
	"HLT_Mu15_IsoVVVL_PFHT600",
	"HLT_Mu15_IsoVVVL_PFHT350",
	"HLT_Mu15_IsoVVVL_PFHT400_PFMET50",
	"HLT_Mu15_IsoVVVL_PFHT400",
	"HLT_Mu50_IsoVVVL_PFHT400",
	"HLT_Dimuon16_Jpsi",
	"HLT_Dimuon10_Jpsi_Barrel",
	"HLT_Dimuon8_PsiPrime_Barrel",
	"HLT_Dimuon8_Upsilon_Barrel",
	"HLT_Dimuon0_Phi_Barrel",
	"HLT_Mu16_TkMu0_dEta18_Onia",
	"HLT_Mu16_TkMu0_dEta18_Phi",
	"HLT_TrkMu15_DoubleTrkMu5NoFiltersNoVtx",
	"HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx",
	"HLT_Mu8",
	"HLT_Mu17",
	"HLT_Mu3_PFJet40",
	"HLT_Ele8_CaloIdM_TrackIdM_PFJet30",
	"HLT_Ele12_CaloIdM_TrackIdM_PFJet30",
	"HLT_Ele17_CaloIdM_TrackIdM_PFJet30",
	"HLT_Ele23_CaloIdM_TrackIdM_PFJet30",
	"HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet140",
	"HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165",
	"HLT_PFHT400_SixJet30_DoubleBTagCSV_p056",
	"HLT_PFHT450_SixJet40_BTagCSV_p056",
	"HLT_PFHT400_SixJet30",
	"HLT_PFHT450_SixJet40",
	"HLT_Ele115_CaloIdVT_GsfTrkIdT",
	"HLT_Mu55",
	"HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15",
	"HLT_Photon90_CaloIdL_PFHT600",
	"HLT_PixelTracks_Multiplicity60ForEndOfFill",
	"HLT_PixelTracks_Multiplicity85ForEndOfFill",
	"HLT_PixelTracks_Multiplicity110ForEndOfFill",
	"HLT_PixelTracks_Multiplicity135ForEndOfFill",
	"HLT_PixelTracks_Multiplicity160ForEndOfFill",
	"HLT_FullTracks_Multiplicity80",
	"HLT_FullTracks_Multiplicity100",
	"HLT_FullTracks_Multiplicity130",
	"HLT_FullTracks_Multiplicity150",
	"HLT_ECALHT800",
	"HLT_DiSC30_18_EIso_AND_HE_Mass70",
	"HLT_MET200",
	"HLT_Ele27_HighEta_Ele20_Mass55",
	"HLT_L1FatEvents",
	"HLT_Physics",
	"HLT_Physics_part0",
	"HLT_Physics_part1",
	"HLT_Physics_part2",
	"HLT_Physics_part3",
	"HLT_Random",
	"HLT_ZeroBias",
	"HLT_AK4CaloJet30",
	"HLT_AK4CaloJet40",
	"HLT_AK4CaloJet50",
	"HLT_AK4CaloJet80",
	"HLT_AK4CaloJet100",
	"HLT_AK4PFJet30",
	"HLT_AK4PFJet50",
	"HLT_AK4PFJet80",
	"HLT_AK4PFJet100",
	"HLT_HISinglePhoton10",
	"HLT_HISinglePhoton15",
	"HLT_HISinglePhoton20",
	"HLT_HISinglePhoton40",
	"HLT_HISinglePhoton60",
	"HLT_EcalCalibration",
	"HLT_HcalCalibration",
	"HLT_GlobalRunHPDNoise",
	"HLT_L1BptxMinus",
	"HLT_L1BptxPlus",
	"HLT_L1NotBptxOR",
	"HLT_L1BeamGasMinus",
	"HLT_L1BeamGasPlus",
	"HLT_L1BptxXOR",
	"HLT_L1MinimumBiasHF_OR",
	"HLT_L1MinimumBiasHF_AND",
	"HLT_HcalNZS",
	"HLT_HcalPhiSym",
	"HLT_ZeroBias_FirstCollisionAfterAbortGap",
	"HLT_ZeroBias_FirstCollisionAfterAbortGap_TCDS",
	"HLT_ZeroBias_IsolatedBunches",
	"HLT_Photon500",
	"HLT_Photon600",
	"HLT_Mu300",
	"HLT_Mu350",
	"HLT_MET250",
	"HLT_MET300",
	"HLT_MET600",
	"HLT_MET700",
	"HLT_PFMET300",
	"HLT_PFMET400",
	"HLT_PFMET500",
	"HLT_PFMET600",
	"HLT_HT2000",
	"HLT_HT2500",
	  "HLT_IsoTrackHE",
	"HLT_IsoTrackHB"
	  };
   



      hltresults_ {
	HLTriggerFirstPath,
	HLT_AK8PFJet360_TrimMass30,
	HLT_AK8PFHT700_TrimR0p1PT0p03Mass50,
	HLT_AK8PFHT650_TrimR0p1PT0p03Mass50,
	HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20,
	HLT_CaloJet500_NoJetID,
	HLT_Dimuon13_PsiPrime,
	HLT_Dimuon13_Upsilon,
	HLT_Dimuon20_Jpsi,
	HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf,
	HLT_DoubleEle33_CaloIdL,
	HLT_DoubleEle33_CaloIdL_MW,
	HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW,
	HLT_DoubleEle33_CaloIdL_GsfTrkIdVL,
	HLT_DoubleEle37_Ele27_CaloIdL_GsfTrkIdVL,
	HLT_DoubleMediumIsoPFTau32_Trk1_eta2p1_Reg,
	HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg,
	HLT_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg,
	HLT_DoubleMu33NoFiltersNoVtx,
	HLT_DoubleMu38NoFiltersNoVtx,
	HLT_DoubleMu23NoFiltersNoVtxDisplaced,
	HLT_DoubleMu28NoFiltersNoVtxDisplaced,
	HLT_DoubleMu4_3_Bs,
	HLT_DoubleMu4_3_Jpsi_Displaced,
	HLT_DoubleMu4_JpsiTrk_Displaced,
	HLT_DoubleMu4_LowMassNonResonantTrk_Displaced,
	HLT_DoubleMu3_Trk_Tau3mu,
	HLT_DoubleMu4_PsiPrimeTrk_Displaced,
	HLT_Mu7p5_L2Mu2_Jpsi,
	HLT_Mu7p5_L2Mu2_Upsilon,
	HLT_Mu7p5_Track2_Jpsi,
	HLT_Mu7p5_Track3p5_Jpsi,
	HLT_Mu7p5_Track7_Jpsi,
	HLT_Mu7p5_Track2_Upsilon,
	HLT_Mu7p5_Track3p5_Upsilon,
	HLT_Mu7p5_Track7_Upsilon,
	HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing,
	HLT_Dimuon0er16_Jpsi_NoVertexing,
	HLT_Dimuon6_Jpsi_NoVertexing,
	HLT_DoublePhoton60,
	HLT_DoublePhoton85,
	HLT_Ele22_eta2p1_WPLoose_Gsf,
	HLT_Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1,
	HLT_Ele23_WPLoose_Gsf,
	HLT_Ele23_WPLoose_Gsf_WHbbBoost,
	HLT_Ele24_eta2p1_WPLoose_Gsf,
	HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20,
	HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1,
	HLT_Ele25_WPTight_Gsf,
	HLT_Ele25_eta2p1_WPLoose_Gsf,
	HLT_Ele25_eta2p1_WPTight_Gsf,
	HLT_Ele27_WPLoose_Gsf,
	HLT_Ele27_WPLoose_Gsf_WHbbBoost,
	HLT_Ele27_WPTight_Gsf,
	HLT_Ele27_eta2p1_WPLoose_Gsf,
	HLT_Ele27_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1,
	HLT_Ele27_eta2p1_WPLoose_Gsf_DoubleMediumIsoPFTau32_Trk1_eta2p1_Reg,
	HLT_Ele27_eta2p1_WPLoose_Gsf_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg,
	HLT_Ele27_eta2p1_WPLoose_Gsf_DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg,
	HLT_Ele27_eta2p1_WPTight_Gsf,
	HLT_Ele32_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1,
	HLT_Ele32_eta2p1_WPTight_Gsf,
	HLT_Ele35_WPLoose_Gsf,
	HLT_Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50,
	HLT_Ele45_WPLoose_Gsf,
	HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50,
	HLT_Ele105_CaloIdVT_GsfTrkIdT,
	HLT_Ele30WP60_SC4_Mass55,
	HLT_Ele30WP60_Ele8_Mass55,
	HLT_HT200,
	HLT_HT275,
	HLT_HT325,
	HLT_HT425,
	HLT_HT575,
	HLT_HT410to430,
	HLT_HT430to450,
	HLT_HT450to470,
	HLT_HT470to500,
	HLT_HT500to550,
	HLT_HT550to650,
	HLT_HT650,
	HLT_Mu16_eta2p1_MET30,
	HLT_IsoMu16_eta2p1_MET30,
	HLT_IsoMu16_eta2p1_MET30_LooseIsoPFTau50_Trk30_eta2p1,
	HLT_IsoMu17_eta2p1_LooseIsoPFTau20,
	HLT_IsoMu17_eta2p1_LooseIsoPFTau20_SingleL1,
	HLT_DoubleIsoMu17_eta2p1,
	HLT_DoubleIsoMu17_eta2p1_noDzCut,
	HLT_IsoMu18,
	HLT_IsoMu19_eta2p1_LooseIsoPFTau20,
	HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1,
	HLT_IsoMu19_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg,
	HLT_IsoMu20,
	HLT_IsoMu21_eta2p1_LooseIsoPFTau20_SingleL1,
	HLT_IsoMu21_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg,
	HLT_IsoMu22,
	HLT_IsoMu22_eta2p1,
	HLT_IsoMu24,
	HLT_IsoMu27,
	HLT_IsoTkMu18,
	HLT_IsoTkMu20,
	HLT_IsoTkMu22,
	HLT_IsoTkMu22_eta2p1,
	HLT_IsoTkMu24,
	HLT_IsoTkMu27,
	HLT_JetE30_NoBPTX3BX,
	HLT_JetE30_NoBPTX,
	HLT_JetE50_NoBPTX3BX,
	HLT_JetE70_NoBPTX3BX,
	HLT_L1SingleMu18,
	HLT_L2Mu10,
	HLT_L1SingleMuOpen,
	HLT_L1SingleMuOpen_DT,
	HLT_L2DoubleMu23_NoVertex,
	HLT_L2DoubleMu28_NoVertex_2Cha_Angle2p5_Mass10,
	HLT_L2DoubleMu38_NoVertex_2Cha_Angle2p5_Mass10,
	HLT_L2Mu10_NoVertex_NoBPTX3BX,
	HLT_L2Mu10_NoVertex_NoBPTX,
	HLT_L2Mu35_NoVertex_3Sta_NoBPTX3BX,
	HLT_L2Mu40_NoVertex_3Sta_NoBPTX3BX,
	HLT_LooseIsoPFTau50_Trk30_eta2p1,
	HLT_LooseIsoPFTau50_Trk30_eta2p1_MET80,
	HLT_LooseIsoPFTau50_Trk30_eta2p1_MET90,
	HLT_LooseIsoPFTau50_Trk30_eta2p1_MET110,
	HLT_LooseIsoPFTau50_Trk30_eta2p1_MET120,
	HLT_PFTau120_eta2p1,
	HLT_VLooseIsoPFTau120_Trk50_eta2p1,
	HLT_VLooseIsoPFTau140_Trk50_eta2p1,
	HLT_Mu17_Mu8,
	HLT_Mu17_Mu8_DZ,
	HLT_Mu17_Mu8_SameSign,
	HLT_Mu17_Mu8_SameSign_DZ,
	HLT_Mu20_Mu10,
	HLT_Mu20_Mu10_DZ,
	HLT_Mu20_Mu10_SameSign,
	HLT_Mu20_Mu10_SameSign_DZ,
	HLT_Mu17_TkMu8_DZ,
	HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL,
	HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ,
	HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL,
	HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ,
	HLT_Mu25_TkMu0_dEta18_Onia,
	HLT_Mu27_TkMu8,
	HLT_Mu30_TkMu11,
	HLT_Mu30_eta2p1_PFJet150_PFJet50,
	HLT_Mu40_TkMu11,
	HLT_Mu40_eta2p1_PFJet200_PFJet50,
	HLT_Mu20,
	HLT_TkMu20,
	HLT_Mu24_eta2p1,
	HLT_TkMu24_eta2p1,
	HLT_Mu27,
	HLT_TkMu27,
	HLT_Mu45_eta2p1,
	HLT_Mu50,
	HLT_TkMu50,
	HLT_Mu38NoFiltersNoVtx_Photon38_CaloIdL,
	HLT_Mu42NoFiltersNoVtx_Photon42_CaloIdL,
	HLT_Mu28NoFiltersNoVtxDisplaced_Photon28_CaloIdL,
	HLT_Mu33NoFiltersNoVtxDisplaced_Photon33_CaloIdL,
	HLT_Mu23NoFiltersNoVtx_Photon23_CaloIdL,
	HLT_DoubleMu18NoFiltersNoVtx,
	HLT_Mu33NoFiltersNoVtxDisplaced_DisplacedJet50_Tight,
	HLT_Mu33NoFiltersNoVtxDisplaced_DisplacedJet50_Loose,
	HLT_Mu28NoFiltersNoVtx_DisplacedJet40_Loose,
	HLT_Mu38NoFiltersNoVtxDisplaced_DisplacedJet60_Tight,
	HLT_Mu38NoFiltersNoVtxDisplaced_DisplacedJet60_Loose,
	HLT_Mu38NoFiltersNoVtx_DisplacedJet60_Loose,
	HLT_Mu28NoFiltersNoVtx_CentralCaloJet40,
	HLT_PFHT300_PFMET100,
	HLT_PFHT300_PFMET110,
	HLT_PFHT550_4JetPt50,
	HLT_PFHT650_4JetPt50,
	HLT_PFHT750_4JetPt50,
	HLT_PFJet15_NoCaloMatched,
	HLT_PFJet25_NoCaloMatched,
	HLT_DiPFJet15_NoCaloMatched,
	HLT_DiPFJet25_NoCaloMatched,
	HLT_DiPFJet15_FBEta3_NoCaloMatched,
	HLT_DiPFJet25_FBEta3_NoCaloMatched,
	HLT_DiPFJetAve15_HFJEC,
	HLT_DiPFJetAve25_HFJEC,
	HLT_DiPFJetAve35_HFJEC,
	HLT_AK8PFJet40,
	HLT_AK8PFJet60,
	HLT_AK8PFJet80,
	HLT_AK8PFJet140,
	HLT_AK8PFJet200,
	HLT_AK8PFJet260,
	HLT_AK8PFJet320,
	HLT_AK8PFJet400,
	HLT_AK8PFJet450,
	HLT_AK8PFJet500,
	HLT_PFJet40,
	HLT_PFJet60,
	HLT_PFJet80,
	HLT_PFJet140,
	HLT_PFJet200,
	HLT_PFJet260,
	HLT_PFJet320,
	HLT_PFJet400,
	HLT_PFJet450,
	HLT_PFJet500,
	HLT_DiPFJetAve40,
	HLT_DiPFJetAve60,
	HLT_DiPFJetAve80,
	HLT_DiPFJetAve140,
	HLT_DiPFJetAve200,
	HLT_DiPFJetAve260,
	HLT_DiPFJetAve320,
	HLT_DiPFJetAve400,
	HLT_DiPFJetAve500,
	HLT_DiPFJetAve60_HFJEC,
	HLT_DiPFJetAve80_HFJEC,
	HLT_DiPFJetAve100_HFJEC,
	HLT_DiPFJetAve160_HFJEC,
	HLT_DiPFJetAve220_HFJEC,
	HLT_DiPFJetAve300_HFJEC,
	HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu140,
	HLT_DiPFJet40_DEta3p5_MJJ600_PFMETNoMu80,
	HLT_DiCentralPFJet55_PFMET110,
	HLT_DiCentralPFJet170,
	HLT_SingleCentralPFJet170_CFMax0p1,
	HLT_DiCentralPFJet170_CFMax0p1,
	HLT_DiCentralPFJet220_CFMax0p3,
	HLT_DiCentralPFJet330_CFMax0p5,
	HLT_DiCentralPFJet430,
	HLT_PFHT125,
	HLT_PFHT200,
	HLT_PFHT250,
	HLT_PFHT300,
	HLT_PFHT350,
	HLT_PFHT400,
	HLT_PFHT475,
	HLT_PFHT600,
	HLT_PFHT650,
	HLT_PFHT800,
	HLT_PFHT900,
	HLT_PFHT200_PFAlphaT0p51,
	HLT_PFHT200_DiPFJetAve90_PFAlphaT0p57,
	HLT_PFHT200_DiPFJetAve90_PFAlphaT0p63,
	HLT_PFHT250_DiPFJetAve90_PFAlphaT0p55,
	HLT_PFHT250_DiPFJetAve90_PFAlphaT0p58,
	HLT_PFHT300_DiPFJetAve90_PFAlphaT0p53,
	HLT_PFHT300_DiPFJetAve90_PFAlphaT0p54,
	HLT_PFHT350_DiPFJetAve90_PFAlphaT0p52,
	HLT_PFHT350_DiPFJetAve90_PFAlphaT0p53,
	HLT_PFHT400_DiPFJetAve90_PFAlphaT0p51,
	HLT_PFHT400_DiPFJetAve90_PFAlphaT0p52,
	HLT_MET60_IsoTrk35_Loose,
	HLT_MET75_IsoTrk50,
	HLT_MET90_IsoTrk50,
	HLT_PFMET120_BTagCSV_p067,
	HLT_PFMET120_Mu5,
	HLT_PFMET170_NoiseCleaned,
	HLT_PFMET170_HBHECleaned,
	HLT_PFMET170_JetIdCleaned,
	HLT_PFMET170_NotCleaned,
	HLT_PFMET170_BeamHaloCleaned,
	HLT_PFMET90_PFMHT90_IDTight,
	HLT_PFMET100_PFMHT100_IDTight,
	HLT_PFMET110_PFMHT110_IDTight,
	HLT_PFMET120_PFMHT120_IDTight,
	HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight_BTagCSV_p067,
	HLT_CaloMHTNoPU90_PFMET90_PFMHT90_IDTight,
	HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq200,
	HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq460,
	HLT_QuadPFJet_BTagCSV_p016_p11_VBF_Mqq240,
	HLT_QuadPFJet_BTagCSV_p016_VBF_Mqq500,
	HLT_QuadPFJet_VBF,
	HLT_L1_TripleJet_VBF,
	HLT_QuadJet45_TripleBTagCSV_p087,
	HLT_QuadJet45_DoubleBTagCSV_p087,
	HLT_DoubleJet90_Double30_TripleBTagCSV_p087,
	HLT_DoubleJet90_Double30_DoubleBTagCSV_p087,
	HLT_DoubleJetsC100_DoubleBTagCSV_p026_DoublePFJetsC160,
	HLT_DoubleJetsC100_DoubleBTagCSV_p014_DoublePFJetsC100MaxDeta1p6,
	HLT_DoubleJetsC112_DoubleBTagCSV_p026_DoublePFJetsC172,
	HLT_DoubleJetsC112_DoubleBTagCSV_p014_DoublePFJetsC112MaxDeta1p6,
	HLT_DoubleJetsC100_SingleBTagCSV_p026,
	HLT_DoubleJetsC100_SingleBTagCSV_p014,
	HLT_DoubleJetsC100_SingleBTagCSV_p026_SinglePFJetC350,
	HLT_DoubleJetsC100_SingleBTagCSV_p014_SinglePFJetC350,
	HLT_Photon135_PFMET100,
	HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_PFMET40,
	HLT_Photon22_R9Id90_HE10_Iso40_EBOnly_VBF,
	HLT_Photon250_NoHE,
	HLT_Photon300_NoHE,
	HLT_Photon26_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon16_AND_HE10_R9Id65_Eta2_Mass60,
	HLT_Photon36_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon22_AND_HE10_R9Id65_Eta2_Mass15,
	HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_PFMET40,
	HLT_Photon36_R9Id90_HE10_Iso40_EBOnly_VBF,
	HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_PFMET40,
	HLT_Photon50_R9Id90_HE10_Iso40_EBOnly_VBF,
	HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_PFMET40,
	HLT_Photon75_R9Id90_HE10_Iso40_EBOnly_VBF,
	HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_PFMET40,
	HLT_Photon90_R9Id90_HE10_Iso40_EBOnly_VBF,
	HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_PFMET40,
	HLT_Photon120_R9Id90_HE10_Iso40_EBOnly_VBF,
	HLT_Mu8_TrkIsoVVL,
	HLT_Mu17_TrkIsoVVL,
	HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30,
	HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30,
	HLT_Ele17_CaloIdL_TrackIdL_IsoVL_PFJet30,
	HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30,
	HLT_BTagMu_DiJet20_Mu5,
	HLT_BTagMu_DiJet40_Mu5,
	HLT_BTagMu_DiJet70_Mu5,
	HLT_BTagMu_DiJet110_Mu5,
	HLT_BTagMu_DiJet170_Mu5,
	HLT_BTagMu_Jet300_Mu5,
	HLT_BTagMu_AK8Jet300_Mu5,
	HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ,
	HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ,
	HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL,
	HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL,
	HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL,
	HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL,
	HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL,
	HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL,
	HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL,
	HLT_Mu37_Ele27_CaloIdL_GsfTrkIdVL,
	HLT_Mu27_Ele37_CaloIdL_GsfTrkIdVL,
	HLT_Mu8_DiEle12_CaloIdL_TrackIdL,
	HLT_Mu12_Photon25_CaloIdL,
	HLT_Mu12_Photon25_CaloIdL_L1ISO,
	HLT_Mu12_Photon25_CaloIdL_L1OR,
	HLT_Mu17_Photon22_CaloIdL_L1ISO,
	HLT_Mu17_Photon30_CaloIdL_L1ISO,
	HLT_Mu17_Photon35_CaloIdL_L1ISO,
	HLT_DiMu9_Ele9_CaloIdL_TrackIdL,
	HLT_TripleMu_5_3_3,
	HLT_TripleMu_12_10_5,
	HLT_Mu3er_PFHT140_PFMET125,
	HLT_Mu6_PFHT200_PFMET80_BTagCSV_p067,
	HLT_Mu6_PFHT200_PFMET100,
	HLT_Mu14er_PFMET100,
	HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL,
	HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL,
	HLT_Ele12_CaloIdL_TrackIdL_IsoVL,
	HLT_Ele17_CaloIdL_GsfTrkIdVL,
	HLT_Ele17_CaloIdL_TrackIdL_IsoVL,
	HLT_Ele23_CaloIdL_TrackIdL_IsoVL,
	HLT_AK8DiPFJet280_200_TrimMass30,
	HLT_AK8DiPFJet250_200_TrimMass30,
	HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20,
	HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV_p20,
	HLT_PFHT650_WideJetMJJ900DEtaJJ1p5,
	HLT_PFHT650_WideJetMJJ950DEtaJJ1p5,
	HLT_Photon22,
	HLT_Photon30,
	HLT_Photon36,
	HLT_Photon50,
	HLT_Photon75,
	HLT_Photon90,
	HLT_Photon120,
	HLT_Photon175,
	HLT_Photon165_HE10,
	HLT_Photon22_R9Id90_HE10_IsoM,
	HLT_Photon30_R9Id90_HE10_IsoM,
	HLT_Photon36_R9Id90_HE10_IsoM,
	HLT_Photon50_R9Id90_HE10_IsoM,
	HLT_Photon75_R9Id90_HE10_IsoM,
	HLT_Photon90_R9Id90_HE10_IsoM,
	HLT_Photon120_R9Id90_HE10_IsoM,
	HLT_Photon165_R9Id90_HE10_IsoM,
	HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90,
	HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelSeedMatch_Mass70,
	HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55,
	HLT_Diphoton30_18_Solid_R9Id_AND_IsoCaloId_AND_HE_R9Id_Mass55,
	HLT_Diphoton30EB_18EB_R9Id_OR_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55,
	HLT_Dimuon0_Jpsi_Muon,
	HLT_Dimuon0_Upsilon_Muon,
	HLT_QuadMuon0_Dimuon0_Jpsi,
	HLT_QuadMuon0_Dimuon0_Upsilon,
	HLT_Rsq0p25,
	HLT_Rsq0p30,
	HLT_RsqMR240_Rsq0p09_MR200,
	HLT_RsqMR240_Rsq0p09_MR200_4jet,
	HLT_RsqMR270_Rsq0p09_MR200,
	HLT_RsqMR270_Rsq0p09_MR200_4jet,
	HLT_Rsq0p02_MR300_TriPFJet80_60_40_BTagCSV_p063_p20_Mbb60_200,
	HLT_Rsq0p02_MR300_TriPFJet80_60_40_DoubleBTagCSV_p063_Mbb60_200,
	HLT_HT200_DisplacedDijet40_DisplacedTrack,
	HLT_HT250_DisplacedDijet40_DisplacedTrack,
	HLT_HT350_DisplacedDijet40_DisplacedTrack,
	HLT_HT350_DisplacedDijet80_DisplacedTrack,
	HLT_HT350_DisplacedDijet80_Tight_DisplacedTrack,
	HLT_HT350_DisplacedDijet40_Inclusive,
	HLT_HT400_DisplacedDijet40_Inclusive,
	HLT_HT500_DisplacedDijet40_Inclusive,
	HLT_HT550_DisplacedDijet40_Inclusive,
	HLT_HT650_DisplacedDijet80_Inclusive,
	HLT_HT750_DisplacedDijet80_Inclusive,
	HLT_VBF_DisplacedJet40_DisplacedTrack,
	HLT_VBF_DisplacedJet40_DisplacedTrack_2TrackIP2DSig5,
	HLT_VBF_DisplacedJet40_TightID_DisplacedTrack,
	HLT_VBF_DisplacedJet40_Hadronic,
	HLT_VBF_DisplacedJet40_Hadronic_2PromptTrack,
	HLT_VBF_DisplacedJet40_TightID_Hadronic,
	HLT_VBF_DisplacedJet40_VTightID_Hadronic,
	HLT_VBF_DisplacedJet40_VVTightID_Hadronic,
	HLT_VBF_DisplacedJet40_VTightID_DisplacedTrack,
	HLT_VBF_DisplacedJet40_VVTightID_DisplacedTrack,
	HLT_PFMETNoMu90_PFMHTNoMu90_IDTight,
	HLT_PFMETNoMu100_PFMHTNoMu100_IDTight,
	HLT_PFMETNoMu110_PFMHTNoMu110_IDTight,
	HLT_PFMETNoMu120_PFMHTNoMu120_IDTight,
	HLT_MonoCentralPFJet80_PFMETNoMu90_PFMHTNoMu90_IDTight,
	HLT_MonoCentralPFJet80_PFMETNoMu100_PFMHTNoMu100_IDTight,
	HLT_MonoCentralPFJet80_PFMETNoMu110_PFMHTNoMu110_IDTight,
	HLT_MonoCentralPFJet80_PFMETNoMu120_PFMHTNoMu120_IDTight,
	HLT_Ele27_eta2p1_WPLoose_Gsf_HT200,
	HLT_Photon90_CaloIdL_PFHT500,
	HLT_DoubleMu8_Mass8_PFHT250,
	HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT250,
	HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT250,
	HLT_DoubleMu8_Mass8_PFHT300,
	HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300,
	HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300,
	HLT_Mu10_CentralPFJet30_BTagCSV_p13,
	HLT_DoubleMu3_PFMET50,
	HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV_p13,
	HLT_Ele15_IsoVVVL_BTagCSV_p067_PFHT400,
	HLT_Ele15_IsoVVVL_PFHT350_PFMET50,
	HLT_Ele15_IsoVVVL_PFHT600,
	HLT_Ele15_IsoVVVL_PFHT350,
	HLT_Ele15_IsoVVVL_PFHT400_PFMET50,
	HLT_Ele15_IsoVVVL_PFHT400,
	HLT_Ele50_IsoVVVL_PFHT400,
	HLT_Mu8_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT300_PFMETNoMu60,
	HLT_Mu10_TrkIsoVVL_DiPFJet40_DEta3p5_MJJ750_HTT350_PFMETNoMu60,
	HLT_Mu15_IsoVVVL_BTagCSV_p067_PFHT400,
	HLT_Mu15_IsoVVVL_PFHT350_PFMET50,
	HLT_Mu15_IsoVVVL_PFHT600,
	HLT_Mu15_IsoVVVL_PFHT350,
	HLT_Mu15_IsoVVVL_PFHT400_PFMET50,
	HLT_Mu15_IsoVVVL_PFHT400,
	HLT_Mu50_IsoVVVL_PFHT400,
	HLT_Dimuon16_Jpsi,
	HLT_Dimuon10_Jpsi_Barrel,
	HLT_Dimuon8_PsiPrime_Barrel,
	HLT_Dimuon8_Upsilon_Barrel,
	HLT_Dimuon0_Phi_Barrel,
	HLT_Mu16_TkMu0_dEta18_Onia,
	HLT_Mu16_TkMu0_dEta18_Phi,
	HLT_TrkMu15_DoubleTrkMu5NoFiltersNoVtx,
	HLT_TrkMu17_DoubleTrkMu8NoFiltersNoVtx,
	HLT_Mu8,
	HLT_Mu17,
	HLT_Mu3_PFJet40,
	HLT_Ele8_CaloIdM_TrackIdM_PFJet30,
	HLT_Ele12_CaloIdM_TrackIdM_PFJet30,
	HLT_Ele17_CaloIdM_TrackIdM_PFJet30,
	HLT_Ele23_CaloIdM_TrackIdM_PFJet30,
	HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet140,
	HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165,
	HLT_PFHT400_SixJet30_DoubleBTagCSV_p056,
	HLT_PFHT450_SixJet40_BTagCSV_p056,
	HLT_PFHT400_SixJet30,
	HLT_PFHT450_SixJet40,
	HLT_Ele115_CaloIdVT_GsfTrkIdT,
	HLT_Mu55,
	HLT_Photon42_R9Id85_OR_CaloId24b40e_Iso50T80L_Photon25_AND_HE10_R9Id65_Eta2_Mass15,
	HLT_Photon90_CaloIdL_PFHT600,
	HLT_PixelTracks_Multiplicity60ForEndOfFill,
	HLT_PixelTracks_Multiplicity85ForEndOfFill,
	HLT_PixelTracks_Multiplicity110ForEndOfFill,
	HLT_PixelTracks_Multiplicity135ForEndOfFill,
	HLT_PixelTracks_Multiplicity160ForEndOfFill,
	HLT_FullTracks_Multiplicity80,
	HLT_FullTracks_Multiplicity100,
	HLT_FullTracks_Multiplicity130,
	HLT_FullTracks_Multiplicity150,
	HLT_ECALHT800,
	HLT_DiSC30_18_EIso_AND_HE_Mass70,
	HLT_MET200,
	HLT_Ele27_HighEta_Ele20_Mass55,
	HLT_L1FatEvents,
	HLT_Physics,
	HLT_Physics_part0,
	HLT_Physics_part1,
	HLT_Physics_part2,
	HLT_Physics_part3,
	HLT_Random,
	HLT_ZeroBias,
	HLT_AK4CaloJet30,
	HLT_AK4CaloJet40,
	HLT_AK4CaloJet50,
	HLT_AK4CaloJet80,
	HLT_AK4CaloJet100,
	HLT_AK4PFJet30,
	HLT_AK4PFJet50,
	HLT_AK4PFJet80,
	HLT_AK4PFJet100,
	HLT_HISinglePhoton10,
	HLT_HISinglePhoton15,
	HLT_HISinglePhoton20,
	HLT_HISinglePhoton40,
	HLT_HISinglePhoton60,
	HLT_EcalCalibration,
	HLT_HcalCalibration,
	HLT_GlobalRunHPDNoise,
	HLT_L1BptxMinus,
	HLT_L1BptxPlus,
	HLT_L1NotBptxOR,
	HLT_L1BeamGasMinus,
	HLT_L1BeamGasPlus,
	HLT_L1BptxXOR,
	HLT_L1MinimumBiasHF_OR,
	HLT_L1MinimumBiasHF_AND,
	HLT_HcalNZS,
	HLT_HcalPhiSym,
	HLT_ZeroBias_FirstCollisionAfterAbortGap,
	HLT_ZeroBias_FirstCollisionAfterAbortGap_TCDS,
	HLT_ZeroBias_IsolatedBunches,
	HLT_Photon500,
	HLT_Photon600,
	HLT_Mu300,
	HLT_Mu350,
	HLT_MET250,
	HLT_MET300,
	HLT_MET600,
	HLT_MET700,
	HLT_PFMET300,
	HLT_PFMET400,
	HLT_PFMET500,
	HLT_PFMET600,
	HLT_HT2000,
	HLT_HT2500,
	HLT_IsoTrackHE,
	HLT_IsoTrackHB


	  };

      hltprescales_ (532, 1.0);   
   }
} 
