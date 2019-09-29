#include "interface/PhysicsObjects.h"

#define NEL(x) (sizeof((x))/sizeof((x)[0]))

ClassImp(vhtm::Event)
ClassImp(vhtm::GenEvent)
ClassImp(vhtm::Electron)
ClassImp(vhtm::GenParticle)
ClassImp(vhtm::GenJet)
ClassImp(vhtm::GenMET)
ClassImp(vhtm::MET)
ClassImp(vhtm::Tau)
ClassImp(vhtm::Muon)
ClassImp(vhtm::Jet)
ClassImp(vhtm::Vertex)
ClassImp(vhtm::TriggerObject)
ClassImp(vhtm::Candidate)
ClassImp(vhtm::Photon)
ClassImp(vhtm::PackedPFCandidate)

vhtm::Candidate::Candidate():
  pt(-999), eta(-999), phi(-999) 
{} 

vhtm::Candidate::Candidate(float _pt, float _eta, float _phi):
  pt(_pt), eta(_eta), phi(_phi) {} 

vhtm::PackedPFCandidate::PackedPFCandidate():
  pt(-999.),
  eta(-999.), 
  phi(-999.),
  energy(-999.),
  trackHighPurity(false),
  pdgId(0),
  charge(-999),
  vx(-999.),
  vy(-999.),
  vz(-999.),
  fromPV(-999),
  dxy(999.),
  dz(999.),
  dxyError(999.),
  dzError(999.),
  dxyEV(999),
  dzEV(999),
  dzAssociatedPV(999),
  numberOfHits(-1),
  numberOfPixelHits(-1),
  pixelLayersWithMeasurement(-1),
  stripLayersWithMeasurement(-1),
  lostInnerHits(-99)
{
  isolationMap.clear();
}

vhtm::Event::Event():
  run(0),
  event(0),
  lumis(0),
  btagWeight_CSVV2(999.0),
  btagWeight_CMVA(999.0),
  fixedGridRhoFastjetAll(-1.0),
  fixedGridRhoFastjetCentralCalo(-1.0),
  fixedGridRhoFastjetCentralNeutral(-1.0),

  HBHENoiseFilter(false),
  HBHENoiseIsoFilter(false),
  CSCTightHaloFilter(false),
  CSCTightHaloTrkMuUnvetoFilter(false),
  CSCTightHalo2015Filter(false),
  globalTightHalo2016Filter(false),
  globalSuperTightHalo2016Filter(false),
  HcalStripHaloFilter(false),
  hcalLaserEventFilter(false),
  EcalDeadCellTriggerPrimitiveFilter(false),
  EcalDeadCellBoundaryEnergyFilter(false),
  goodVertices(false),
  eeBadScFilter(false),
  ecalLaserCorrFilter(false),
  trkPOGFilters(false),
  chargedHadronTrackResolutionFilter(false),
  muonBadTrackFilter(false),
  BadChargedCandidateFilter(false),
  BadPFMuonFilter(false),
  BadChargedCandidateSummer16Filter(false),
  BadPFMuonSummer16Filter(false),
  trkPOG_manystripclus53X(false),
  trkPOG_toomanystripclus53X(false),
  trkPOG_logErrorTooManyClusters(false),
  METFilters(false),

  nTrueInt(-1.0),
  nPU(-1),
  sumEOOT(-1),
  sumLOOT(-1),

  nvtx(-1)
{}

vhtm::GenEvent::GenEvent():
  evtWeight(1),
  genWeight(1),
  qScalePdf(0),
  x1(0), 
  x2(0),
  x1Pdf(0),
  x2Pdf(0)
{}

vhtm::Electron::Electron():
  eta(-999),
  phi(-999),
  pt(-999),
  charge(-999),
  mass(-999),

  dxyPV(-999),
  dxyPVerr(-999),
  dzPV(-999),
  dzPVerr(-999),
  dB3D(-999),
  SIP3D(-999),
  
  cutBased(-999),
  cutBasedHEEP(false),
  cutBasedHLTPreSel(-999),
  vidNestedWPBitmap(-999),
  tightCharge(-999),
  pdgId(-999),
  mvaTTH(-999),
  mvaSpring16GP(-999),
  mvaSpring16HZZ(-999),
  mvaSpring16GP_WP80(false),
  mvaSpring16GP_WP90(false),
  mvaSpring16HZZ_WPL(false),

  photonIdx(-999),
  genIdx(-999),
  genFlv(-999)

  pfRelIso03(-999),
  pfRelChrIso03(-999),
  hoe(-999),
  sigmaIEtaIEta(-999),
  deltaEtaTrkSC(-999)
{}
vhtm::GenParticle::GenParticle():
  eta(-999),
  phi(-999),
  pt(-999),
  mass(-999),
  pdgId(-999),
  status(-999),
  motherIndex(-1),
  statusFlag(-1)
{}

vhtm::GenJet::GenJet():
  eta(-999),
  mass(-999),
  phi(-999),
  pt(-999),
  hadFlv(255),
  parFlv(-999)
{}

vhtm::MET::MET():
  met(-999),
  metphi(-999),
  sumet(-999),
  covXX(-999),
  covXY(-999),
  covYY(-999),
  MetUnclustEnUpDeltaX(-999),
  MetUnclustEnUpDeltaY(-999)
{}

vhtm::GenMET::GenMET():
  met(-999.),
  metphi(-999.)
{}

vhtm::Tau::Tau():
  eta(-999),
  phi(-999),
  pt(-999),
  charge(-999),
  mass(-999),
  dxyPV(-999),
  dzPV(-999),
  chargedIso(-1),
  neutralIso(-1),
  rawIso(-1),
  rawIsodR03(-999),
  rawMVAnewDM(-999),
  rawMVAoldDM(-999),
  rawMVAoldDMdR03(-999),
  
  idAntiEle(-999),
  idAntiMu(-999),

  decayModeFinding(false),
  decayModeFindingNewDMs(false),
  idMVAnew(false),
  idMVAoldDM(false),
  idMVAoldDMdR03(false),

  leadTkDeltaEta(-999),
  leadTkDeltaPhi(-999),
  leadTkPtOverTauPt(-999),
  puCorr(-999),
  decayMode(-999),
  genFlv(255),
  genIdx(-999),
  jetIdx(-999)
{}

vhtm::Muon::Muon():
  isPFMuon(false),
  isGlobalMuon(false),
  isTrackerMuon(false),
  isLooseMuon(false),
  isMediumMuon(false),
  isTightMuon(false),
  highPtId(-1),
  tightCharge(-1),
  mvaTTH(-1),
  mass(-999),
  pt(-999),
  ptErr(-999),
  eta(-999),
  phi(-999),
  charge(-9999),

  dxyPV(-999),
  dxyPVerr(-999),
  dzPV(-999),
  dzPVerr(-999),
  dB3D(-999),
  SIP3D(-999),
  
  pfRelIso03(-999),   
  pfRelIso03_chg(-999), 
  pfRelIso04(-999),
  
  jetIdx(-999),
  genFlv(255),
  genIdx(-999),
  
  nMatchedStations(-999),
  nMatchedTrkLayers(-999)

{}

vhtm::Jet::Jet():
  ta(-999),
  phi(-999),
  pt(-999),
  mass(-999),
  area(-999),
  chargedEmEnergyFraction(-999),
  chargedHadronEnergyFraction(-999),
  electronEnergyFraction(-999),
  muonEnergyFraction(-999),
  neutralEmEnergyFraction(-999),
  neutralHadronEnergyFraction(-999),
  electronMultiplicity(-999),
  muonMultiplicity(-999),
  nConstituents(-999),

  btagCMVA(-999),
  btagCSVV2(-999),
  btagDeepB(-999),
  btagDeepC(-999),
  btagDeepFlavB(-999),

  partonFlavour(-1),
  enJetIdx(-1),

  electronIdx1(-1),
  electronIdx2(-1),
  muonIdx1(-1),
  muonIdx2(-1),

  puId(-999),
  qgl(-999)

{}

vhtm::Vertex::Vertex():
  x(-999),
  y(-999),
  z(-999),
  ndf(-1.),
  chi2(999.),
  npvs(-999),
  npvsGood(-999),
  //SecondaryVtx
  mass(-999.),
  pt(-999.),
  eta(-999.),
  phi(-999.),
  pAngle(-999.),
  dlenSig(-999.),
  dlen(-999.)
{}

vhtm::TriggerObject::TriggerObject():
  pt(-999),
  eta(-999),
  phi(-999),
  filterBits(-1),
  Id(-1),
  l1pt(-999.),
  l1pt_2(-999.),
  l2pt(-999.),
  l1iso(-999),
  l1charge(-999)
{}

vhtm::Photon::Photon():
  et(-999),
  eta(-999),
  phi(-999),
  energyErr(-999),
  
  r9(-999),
  hoe(-999),
  sigmaIEtaIEta(-999),
  
  electronIdx(-999),
  jetIdx(-999),
  
  cutBased(-999), 
  mvaID(-999),
  pfRelIso03(-999),
  pfRelChrIso03(-999),
  
  pdgId(-999),
  vidNestedWPBitmap(-999),
  
  electronVeto(false),
  isScEtaEB(false),   
  isScEtaEE(false),   
  mvaID_WP80(false),  
  mvaID_WP90(false),  
  pixelSeed(false),   
 
  genFlv(255),
  genIdx(-999)

{}
