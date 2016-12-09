#ifndef DiElectron_OS_h
#define DiElectron_OS_h

#include "AnalyzerCore.h"
class DiElectron_OS : public AnalyzerCore {

 public:
  //// constructors                                                                                                                                                             
  DiElectron_OS();
  ~DiElectron_OS();

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


  ClassDef ( DiElectron_OS, 1);
};
#endif
