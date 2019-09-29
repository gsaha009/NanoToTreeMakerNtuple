#ifndef __NanoToTMntuple__hh
#define __NanoToTMntuple__hh

#define NEL(x) (sizeof((x))/sizeof((x)[0]))

#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <bitset>

#include "TLorentzVector.h"
#include "TVector.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"

#include "PhysicsObjects.h"
#include "NanoUtil.h"
#include "NanoBase.h"


class NanoToTMntuple: public NanoBase {
    
public:
  NanoToTMntuple();
  virtual ~NanoToTMntuple();
  virtual void eventLoop() final;  // the main analysis 
  virtual bool beginJob() override;
  virtual void closeFiles() override;
  void clearLists();

private:
  std::vector<vhtm::Event>* listEv_;
  std::vector<vhtm::Vertex>* listPVtx_;
  std::vector<vhtm::Vertex>* listSVtx_;

  std::vector<vhtm::Electron>* listEle_;
  std::vector<vhtm::Muon>* listMu_;
  std::vector<vhtm::Jet>* listJet_;
  std::vector<vhtm::GenParticle>* listGenP_;
  std::vector<vhtm::PackedPFCandidate>* listPPF_;
  std::vector<vhtm::MET>* listMET_;
  std::vector<vhtm::MET>* listPMET_;
  std::vector<vhtm::MET>* listRMET_;
  std::vector<vhtm::Tau>* listTau_;
  std::vector<vhtm::Photon>* listPhoton_;
  std::vector<std::string>* hltpaths_;
  std::vector<bool>* hltresults_;
  std::vector<int>* hltprescales_;
  std::vector<std::string>* flags_;
  std::vector<bool>* flagValues_;
  //new
  if (isMC()) {
    std::vector<vhtm::GenEvent> listGenEv_;
    std::vector<vhtm::GenParticle> listGenP_;
    std::vector<vhtm::GenJet> listGenJet_;
    std::vector<vhtm::GenMET> listGenMet_;
  }

  std::vector<std::string> eventFilelist_;
  std::unordered_map<std::string, int> eventIdStore_;
  TTree *outTree;
};
#endif
