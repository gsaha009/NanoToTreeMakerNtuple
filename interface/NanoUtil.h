#ifndef __NANOUTIL__HH
#define __NANOUTIL__HH

#include <string>
#include <vector>
#include <map>
#include <unordered_map>

#include "TMath.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TH1K.h"

using std::ostream;

namespace NanoUtil {
  // Templated functioned must be defined in the header itself
  void tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters=" ");
  void bit_print(int value, int pos=32, ostream& os=std::cout);
  template <typename T>
  T deltaPhiT(T phi1, T phi2) {
    T result = phi1 - phi2;
    while (result > TMath::Pi()) result -= 2*TMath::Pi();
    while (result <= -TMath::Pi()) result += 2*TMath::Pi();
    return result;
  }  
  double deltaPhi(double phia, double phib);
  double deltaPhi(const TLorentzVector& a, const TLorentzVector& b);
  double deltaR(const TLorentzVector& a, const TLorentzVector& b);
  bool sameObject(const TLorentzVector& lv1, const TLorentzVector& lv2);
  double cutValue(const std::map<std::string, double>& m, const std::string& cname);
  const std::map<std::string, double>& cutMap(const std::map<std::string, std::map<std::string, double>>& hmap, const std::string& mkey);
  void buildList(const std::vector<std::string>& tokens, std::vector<std::string>& list);
  void buildMap(const std::vector<std::string>& tokens, std::map<std::string, int>& hmap);
  void buildMap(const std::vector<std::string>& tokens, std::unordered_map<std::string, int>& hmap);
  void storeCuts(const std::vector<std::string>& tokens, std::map<std::string, std::map<std::string, double>>& hmap);
  void showCuts(const std::map<std::string, std::map<std::string, double> >& hmap, ostream& os=std::cout);

  template <class T>
  void showList(const T& coll, const char* optcstr="", std::ostream& os=std::cout) {
    os << optcstr << ", Total # = " << coll.size() << ":" << std::endl;
    for (auto const& v: coll)
      os << v << std::endl;
  }  
  template <class T1, class T2>
  void showMap(const std::map<T1,T2>& m, const char* optcstr="", std::ostream& os=std::cout) {
    os << optcstr << std::endl;
    for (auto const& k: m)
      os << k.first << std::endl;
  }
  /*
  template <class T1, class T2>
  void showMap(const std::unordered_map<T1,T2>& m, const char* optcstr="", std::ostream& os=std::cout) {
    os << optcstr << std::endl;
    for (auto const& k: m)
      os << k.first << std::endl;
      }*/
  /* copyList */
  template <class T>
  void copyList (const T& sourceColl, T& destColl) {
    destColl.clear();
    for (auto const& v: sourceColl)   
      destColl.push_back(v); 
  }
}
#endif
