#ifndef DiElectron_SS_h
#define DiElectron_SS_h

#include "AnalyzerCore.h"
class DiElectron_SS : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  DiElectron_SS();
  ~DiElectron_SS();

  /// Functions from core
  virtual void BeginCycle() throw( LQError );
  virtual void BeginEvent()throw( LQError );
  virtual void ExecuteEvents()throw( LQError );
  virtual void EndCycle()throw( LQError );
  virtual void ClearOutputVectors()throw( LQError );
  
  void InitialiseAnalysis() throw( LQError );
  void MakeHistograms();
  void FillCutFlow(TString cut, float w);
 private:
  
  //
  // The output variables 
  //
  /// Vectors for output objetcs
  std::vector<snu::KMuon> out_muons;
  std::vector<snu::KElectron> out_electrons;


  ClassDef ( DiElectron_SS, 1);
};
#endif
