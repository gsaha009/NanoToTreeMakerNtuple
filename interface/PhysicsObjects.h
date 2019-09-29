#ifndef __AnalysisSpace_TreeMaker_PhysicsObjects_h
#define __AnalysisSpace_TreeMaker_PhysicsObjects_h

#include <vector>
#include <map>
#include <string>
#include "TLorentzVector.h"
#include "TObject.h"

namespace vhtm {
  class Candidate: public TObject {
  public:
    Candidate();
    Candidate(float pt, float eta, float phi);
    virtual ~Candidate() {}

    float pt;
    float eta;
    float phi;
   
    ClassDef(Candidate, 1)
  };

  class PackedPFCandidate: public TObject {
  public:
    PackedPFCandidate();
    virtual ~PackedPFCandidate() {}

    float pt;
    float eta;
    float phi;
    float phiAtVtx;
    float energy;
    bool trackHighPurity;
    
    int pdgId;
    int charge;
    
    double vx;
    double vy;
    double vz;
   
    int fromPV;
    // w.r.t PV (essentially pv[0])
    float dxy;
    float dz;
    float dxyError;
    float dzError;

    // wrt PV with highest sumPtTracks
    float dxyEV;
    float dzEV;
    float dzAssociatedPV;

    int numberOfHits;
    int numberOfPixelHits;
    int pixelLayersWithMeasurement;
    int stripLayersWithMeasurement;
    int lostInnerHits;

    std::map<std::string, std::vector<double>> isolationMap;

    ClassDef(PackedPFCandidate, 1)
  };
 
  class Event: public TObject {
  public:
    Event();
    virtual ~Event() {}
  
    UInt_t run;
    ULong64_t  event;
    UInt_t lumis;
  
    Float_t btagWeight_CSVV2;
    Float_t btagWeight_CMVA;

    Float_t         fixedGridRhoFastjetAll;
    Float_t         fixedGridRhoFastjetCentralCalo;
    Float_t         fixedGridRhoFastjetCentralNeutral;
    //PU
    Float_t         nTrueInt;
    Int_t           nPU;
    Int_t           sumEOOT;
    Int_t           sumLOOT;
    
    UInt_t          nvtx; 

    ClassDef(Event, 1)
  };
  class GenEvent: public TObject {
  public:
    GenEvent();
    virtual ~GenEvent() {}
  
    Float_t evtWeight;
    Float_t genWeight;
    Float_t qScalePdf;
    Float_t x1;
    Float_t x2;
    Float_t x1Pdf;
    Float_t x2Pdf;
  
    ClassDef(GenEvent, 1)
  };
  
  class Electron: public TObject {
  public:
    Electron();
    ~Electron() {}

    // Vertex association variables
    Float_t dxyPV;
    Float_t dxyPVerr;
    Float_t dzPV;
    Float_t dzPVerr;
    Float_t dB3D;
    Float_t SIP3D;

    //p4 related
    Float_t eta;
    Float_t phi;
    Float_t pt;
    Int_t   charge;
    Float_t mass;

    //References
    Int_t photonIdx;
    Int_t genIdx;
    UInt8_t genFlv;

    // Isolation variable
    Float_t pfRelIso03;
    Float_t pfRelChrIso03;
    Float_t hoe;
    Float_t sigmaIEtaIEta;
    Float_t deltaEtaTrkSC;

    // ID variables
    Int_t cutBased ;
    Bool_t cutBasedHEEP;
    Int_t cutBasedHLTPreSel;
    Int_t vidNestedWPBitmap;
    Int_t tightCharge;
    Int_t pdgId;
    Float_t mvaTTH;
    Float_t mvaSpring16GP;
    Float_t mvaSpring16HZZ;
    Bool_t mvaSpring16GP_WP80;
    Bool_t mvaSpring16GP_WP90;
    Bool_t mvaSpring16HZZ_WPL;


    ClassDef(Electron, 1)
  };
  class GenParticle: public TObject {
  public:
    GenParticle();
    ~GenParticle() {}
  
    Float_t eta;
    Float_t phi;
    Float_t pt;
    Flat_t mass;
    Int_t pdgId;
    Int_t status;
    Int_t motherIndex;
    Int_t statusFlag;
  
    ClassDef(GenParticle, 1)
  };
  class GenJet: public TObject {
  public:
    GenJet();
    ~GenJet() {}
  
    Float_t eta;
    Flaot_t mass;
    Float_t phi;
    Float_t pt;
    UInt8_t hadFlv;
    Int_t parFlv;
  
    ClassDef(GenJet, 1)
  };
  class MET: public TObject {
  public:
    MET();
    ~MET() {}
  
    Float_t met;
    Float_t metphi;
    Float_t sumet;
    Float_t covXX;
    Float_t covXY;
    Float_t covYY;
    Float_t MetUnclustEnUpDeltaX;
    Float_t MetUnclustEnUpDeltaY;

    ClassDef(MET, 1)
  };
  class Tau: public TObject {
  public:
    Tau();
    ~Tau() {}
  
    Float_t eta;
    Float_t phi;
    Float_t pt;
    Int_t charge;
    Float_t mass;
  
    Float_t dxyPV;
    Float_t dzPV;

    //Isolation
    Float_t chargedIso;
    Float_t neutralIso;
    Float_t rawIso;
    Float_t rawIsodR03;
    Float_t rawMVAnewDM;
    Float_t rawMVAoldDM;
    Float_t rawMVAoldDMdR03;

    UInt8_t idAntiEle;
    UInt8_t idAntiMu;

    Bool_t decayModeFinding;
    Bool_t decayModeFindingNewDMs;
    Bool_t idMVAnew;
    Bool_t idMVAoldDM;
    Bool_t idMVAoldDMdR03;

    Float_t leadTkDeltaEta;
    Float_t leadTkDeltaPhi;
    Float_t leadTkPtOverTauPt;

    Float_t puCorr;

    Int_t decayMode;
    
    UInt8_t genFlv;
    Int_t genIdx;
    Int_t jetIdx;
  
    ClassDef(Tau, 1)
  };
  class Muon: public TObject {
  public:
    Muon();
    ~Muon() {}
    Bool_t isPFMuon;
    Bool_t isGlobalMuon;
    Bool_t isTrackerMuon;
    Bool_t isLooseMuon;
    Bool_t isMediumMuon;
    Bool_t isTightMuon;
    UInt8_t highPtId;
    Int_t tightCharge;
    Float_t mvaTTH;
    Float_t mass;
    Float_t pt;
    Float_t ptErr;
    Float_t eta;
    Float_t phi;
    Float_t charge;

    Float_t dxyPV;
    Float_t dxyPVerr;
    Float_t dzPV;
    Float_t dzPVerr;
    Float_t dB3D;
    Float_t SIP3D;

    Float_t pfRelIso03;   
    Float_t pfRelIso03_chg; 
    Float_t pfRelIso04;
    
    Int_t jetIdx;
    UInt8_t genFlv;
    Int_t genIdx;

    Int_t nMatchedStations;
    Int_t nMatchedTrkLayers;

    ClassDef(Muon, 1)
  };
  class Jet: public TObject {
  public:
    Jet();
    ~Jet() {}
    Float_t eta;
    Float_t phi;
    Float_t pt;
    Float_t mass;
    Float_t area;
    Float_t chargedEmEnergyFraction;
    Float_t chargedHadronEnergyFraction;
    Float_t electronEnergyFraction;
    Float_t muonEnergyFraction;
    Float_t neutralEmEnergyFraction;
    Float_t neutralHadronEnergyFraction;
    Int_t electronMultiplicity;
    Int_t muonMultiplicity;
    Int_t nConstituents;

    Float_t btagCMVA;
    Float_t btagCSVV2;
    Float_t btagDeepB;
    Float_t btagDeepC;
    Float_t btagDeepFlavB;

    Int_t partonFlavour;
    Int_t genJetIdx;

    Int_t electronIdx1;
    Int_t electronIdx2;
    Int_t muonIdx1;
    Int_t muonIdx2;

    Int_t puId;
    Float_t qgl;
  
    ClassDef(Jet, 1)
  };
  class Vertex: public TObject {
  public:
    Vertex();
    virtual ~Vertex() {}
  
    Float_t x;
    Float_t y;
    Float_t z;
    Float_t ndf;
    Float_t chi2;
    Int_t npvs;
    Int_t npvsGood;
    //SecondaryVtx & decayLength Properties
    Float_t mass;
    Float_t pt;
    Float_t phi;
    Float_t eta;
    Float_t pAngle;
    Float_t dlenSig;
    Float_t dlen;

    ClassDef(Vertex, 1)
  };
  class GenMET: public TObject {
  public:
    GenMET();
    virtual ~GenMET() {}
  
    double met;
    double metphi;
    double sumet;
  
    ClassDef(GenMET, 1)
  };
  class TriggerObject: public TObject {
  public:
    TriggerObject();
    virtual ~TriggerObject() {}
  
    Float_t pt;
    Float_t eta;
    Float_t phi;
    Int_t filterBits;
    Int_t Id;
    Float_t l1pt;
    Float_t l1pt_2;
    Float_t l2pt;

    Int_t l1iso;
    Int_t l1charge;
  
    ClassDef(TriggerObject, 1)
  };
  class Photon : public TObject {
  public:
    Photon();
    virtual ~Photon() {}
    
    Float_t et;
    Float_t eta;
    Float_ phi;
    Float_t energyErr;

    Float_t r9;
    Float_t hoe;
    Float_t sigmaIEtaIEta;

    Int_t electronIdx;
    Int_t jetIdx;
    
    Int_t cutBased; 
    Float_t mvaID;
    Float_t pfRelIso03;
    Float_t pfRelChrIso03;

    Int_t pdgId;
    Int_t vidNestedWPBitmap;

    Bool_t electronVeto;
    Bool_t isScEtaEB;   
    Bool_t isScEtaEE;   
    Bool_t mvaID_WP80;  
    Bool_t mvaID_WP90;  
    Bool_t pixelSeed;   
 
    UInt8_t genFlv;
    Int_t genIdx;

    

    ClassDef(Photon, 1)
  };
}
#endif
