#ifndef FatJetSelection_h
#define FatJetSelection_h

#include <iostream>
using namespace std;


#include "TLorentzVector.h"
#include <vector>
#include "LQEvent.h"
#include "KFatJet.h"
#include "KMuon.h"
#include "KElectron.h"
#include "BaseSelection.h"

class FatJetSelection : public BaseSelection {

 public:
  FatJetSelection(LQEvent ev);
  ~FatJetSelection();

  FatJetSelection& operator= (const FatJetSelection& obj);
  FatJetSelection(const FatJetSelection& bs);

 
  void Selection (std::vector<snu::KFatJet>& jetColl);
  void Selection (std::vector<snu::KFatJet>& jetColl, bool LepVeto, std::vector<snu::KMuon>& muonColl, std::vector<snu::KElectron>& electronColl);
  void BasicSelection (std::vector<snu::KFatJet>& jetColl);
  
  bool PassUserID (ID id, snu::KFatJet jet);
  bool PassUserID (snu::KFatJet jet, TString ID);
  bool PassUserID_PFFatJetLoose( snu::KFatJet jet);
  bool PassUserID_PFFatJetMedium( snu::KFatJet jet);
  bool PassUserID_PFFatJetTight( snu::KFatJet jet);
  

  void SelectFatJets(std::vector<snu::KFatJet>& jetColl, std::vector<snu::KMuon> muonColl, std::vector<snu::KElectron> electronColl, TString ID,float ptcut=-999., float etacut=-999.);
  void SelectFatJets(std::vector<snu::KFatJet>& jetColl, TString ID,float ptcut=-999., float etacut=-999.);



};

#endif
