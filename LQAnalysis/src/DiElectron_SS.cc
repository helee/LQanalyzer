// $Id: DiElectron_SS.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQDiElectron_SS Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "DiElectron_SS.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (DiElectron_SS);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
DiElectron_SS::DiElectron_SS() :  AnalyzerCore(), out_muons(0)  {
  
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("DiElectron_SS");
  
  Message("In DiElectron_SS constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();
//  MakeCleverHistograms(sighist_mm,"DiMuon");


}


void DiElectron_SS::InitialiseAnalysis() throw( LQError ) {
  
  /// Initialise histograms
  MakeHistograms();  
  //
  // You can out put messages simply with Message function. Message( "comment", output_level)   output_level can be VERBOSE/INFO/DEBUG/WARNING 
  // You can also use m_logger << level << "comment" << int/double  << LQLogger::endmsg;
  //
  
  Message("Making clever hists for Z ->ll test code", INFO);
  
  /// only available in v7-6-X branch and newer
  //// default lumimask is silver ////
  //// In v7-6-2-(current) the default is changed to gold (since METNoHF bug)
  ///When METNoHF isfixed the default will be back to silver
  /// set to gold if you want to use gold json in analysis
  /// To set uncomment the line below:

  return;
}


void DiElectron_SS::ExecuteEvents()throw( LQError ){

  /// Apply the gen weight 
  if(!isData) weight*=MCweight;
    
  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;
   
  FillCutFlow("NoCut", weight);
  
  if(isData) FillHist("Nvtx_nocut_data",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);
  else  FillHist("Nvtx_nocut_mc",  eventbase->GetEvent().nVertices() ,weight, 0. , 50., 50);


  if(!PassMETFilter()) return;     /// Initial event cuts : 
  FillCutFlow("EventCut", weight);

   /// #### CAT::: triggers stored are all HLT_Ele/HLT_DoubleEle/HLT_Mu/HLT_TkMu/HLT_Photon/HLT_DoublePhoton
   
  if(!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex    
  FillCutFlow("VertexCut", weight);
 
  float pileup_reweight=(1.0);
  if(!isData) { pileup_reweight = TempPileupWeight();}
   
  TString analysis_trigger="HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";   
  vector<TString> trignames;
  trignames.push_back(analysis_trigger);

//  std::vector<snu::KElectron> electrons =  GetElectrons(false,false,"ELECTRON_NOCUT");

  TString elid = "ELECTRON_POG_TIGHT";
  if(k_running_nonprompt) elid = "ELECTRON_HN_FAKELOOSE_NOD0";  
  std::vector<snu::KElectron> electrons =  GetElectrons(true, false, elid);   
  std::vector<snu::KElectron> electronHNVeto = GetElectrons(BaseSelection::ELECTRON_HN_VETO);

/*  std::vector<snu::KElectron> electrons =  GetElectrons(BaseSelection::ELECTRON_NOCUT);  ... WONT WORK
  std::vector<snu::KElectron> electrons =  GetElectrons("ELECTRON_NOCUT");               ... WILL WORK  
   
  std::vector<snu::KElectron> electrons =  GetElectrons(BaseSelection::ELECTRON_POG_TIGHT);  ... WILL WORK  
  std::vector<snu::KElectron> electrons =  GetElectrons("ELECTRON_POG_TIGHT");                ... WILL WORK  
   
   */

   //  std::vector<snu::KElectron> electrons2 =  GetElectrons(BaseSelection::ELECTRON_HN_FAKELOOSE_NOD0);

  std::vector<snu::KJet> jets = GetJets("JET_HN");
  int njet = jets.size();
  int nbjet = NBJet(GetJets("JET_HN"));
  std::vector<snu::KMuon> muons = GetMuons("MUON_HN_TIGHT"); 
  std::vector<snu::KMuon> muonHNVeto = GetMuons(BaseSelection::MUON_HN_VETO);

//   bool trig_pass= true;//PassTrigger("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v", muons, prescale);
//   CorrectMuonMomentum(muons);

  float trig_pass= TriggerEff(analysis_trigger, electrons);
  if(isData) trig_pass = PassTrigger(analysis_trigger);
   
  float trigger_sf(1.);
  float id_iso_sf(1.);
  float trigger_ps(1.);
  float reco_weight=1.;

  float ev_weight = 1;
  ev_weight *= trig_pass;
  if(!isData && !k_running_nonprompt && !k_running_chargeflip){
    trigger_sf  = TriggerScaleFactor(electrons, muons, analysis_trigger);
//     id_iso_sf   = ElectronScaleFactor("ELECTRON_POG_TIGHT", electrons, 0); 
//     reco_weight = ElectronRecoScaleFactor(electrons);
    trigger_ps  = WeightByTrigger(analysis_trigger, TargetLumi);
    for(unsigned int iel=0; iel < electrons.size(); iel++){
      id_iso_sf *= ElectronScaleFactor("ELECTRON_POG_TIGHT", electrons);
      reco_weight *= ElectronRecoScaleFactor(electrons);
    }
    ev_weight = weight * pileup_reweight * trigger_sf * id_iso_sf * reco_weight * trigger_ps;
  }
  if(k_running_nonprompt){
    ev_weight=1.; /// In case... should not be needed
    ev_weight *= Get_DataDrivenWeight_EE(electrons);
  }
  if(isData && k_running_chargeflip){
    ev_weight=1.;
    ev_weight *= WeightCFEvent(electrons, k_running_chargeflip);
  }

  if(trig_pass > 0.){
  weight *= trig_pass;
  FillCutFlow("TriggerCut", weight);
  if(electrons.size() == 2){
    if(muonHNVeto.size() > 0) return;
    if(electronHNVeto.size() > 2) return; 
    FillCutFlow("DiEl", ev_weight);
    if(!k_running_chargeflip){ 
      if(electrons.at(0).Charge()*electrons.at(1).Charge() == -1) return; //require SS
    }
    FillCutFlow("SameSign", ev_weight);
    snu::KParticle X = electrons.at(0) + electrons.at(1);
    if(X.M() < 10.) return;
    if((X.M() > 80.) && (X.M() < 100.)) return;
    if((electrons.at(0).Pt() > 20.) && (electrons.at(1).Pt() > 20.)){
      FillCutFlow("Presel", ev_weight);
      FillHist("mass_ee_alljet", X.M(), ev_weight, 0., 500., 250);
      FillHist("MET_alljet", eventbase->GetEvent().PFMET(), ev_weight, 0., 500., 250);
      FillHist("pT_ee_alljet", X.Pt(), ev_weight, 0., 500., 250);
      FillHist("Njets_alljet", jets.size(), ev_weight, 0., 7., 7);
      FillHist("Nbjets_alljet", nbjet, ev_weight, 0., 5., 5);
      FillHist("Nvtx_alljet", eventbase->GetEvent().nVertices(), ev_weight, 0., 45., 45);
      FillHist("pT_e1_alljet", electrons.at(0).Pt(), ev_weight, 0., 500., 250);
      FillHist("pT_e2_alljet", electrons.at(1).Pt(), ev_weight, 0., 500., 250);
      FillHist("eta_e1_alljet", electrons.at(0).Eta(), ev_weight, -3., 3., 60);
      FillHist("eta_e2_alljet", electrons.at(1).Eta(), ev_weight, -3., 3., 60);

      if(njet == 1){
        FillHist("mass_ee_1jet", X.M(), ev_weight, 0., 500., 250);
        FillHist("MET_1jet", eventbase->GetEvent().PFMET(), ev_weight, 0., 500., 250);
        FillHist("pT_ee_1jet", X.Pt(), ev_weight, 0., 500., 250);
        FillHist("Nvtx_1jet", eventbase->GetEvent().nVertices(), ev_weight, 0., 45., 45);
        FillHist("pT_e1_1jet", electrons.at(0).Pt(), ev_weight, 0., 500., 250);
        FillHist("pT_e2_1jet", electrons.at(1).Pt(), ev_weight, 0., 500., 250);
        FillHist("eta_e1_1jet", electrons.at(0).Eta(), ev_weight, -3., 3., 60);
        FillHist("eta_e2_1jet", electrons.at(1).Eta(), ev_weight, -3., 3., 60);
      }
      if(njet > 1){
        FillHist("mass_ee_2jet", X.M(), ev_weight, 0., 500., 250);
        FillHist("MET_2jet", eventbase->GetEvent().PFMET(), ev_weight, 0., 500., 250);
        FillHist("pT_ee_2jet", X.Pt(), ev_weight, 0., 500., 250);
        FillHist("Njets_2jet", jets.size(), ev_weight, 0., 7., 7);
        FillHist("Nvtx_2jet", eventbase->GetEvent().nVertices(), ev_weight, 0., 45., 45);
        FillHist("pT_e1_2jet", electrons.at(0).Pt(), ev_weight, 0., 500., 250);
        FillHist("pT_e2_2jet", electrons.at(1).Pt(), ev_weight, 0., 500., 250);
        FillHist("eta_e1_2jet", electrons.at(0).Eta(), ev_weight, -3., 3., 60);
        FillHist("eta_e2_2jet", electrons.at(1).Eta(), ev_weight, -3., 3., 60);
      }    
    }
  }
  }

/*  if(isData && k_running_chargeflip){
    ev_weight=1.;
    ev_weight *= WeightCFEvent(electrons, k_running_chargeflip);
  
    if(electrons.size() == 2){ 
    if(muonHNVeto.size() > 0) return;
    if(electronHNVeto.size() > 2) return;
//    if(electrons.at(0).Charge()*electrons.at(1).Charge() == 1) return;
    snu::KParticle X = electrons.at(0) + electrons.at(1);
    if(X.M() < 10.) return;
    if((X.M() > 80.) && (X.M() < 100.)) return;
    if((electrons.at(0).Pt() > 20.) && (electrons.at(1).Pt() > 15.)){
      if(njet == 1){
        FillHist("mass_ee_1jet_CF", X.M(), ev_weight, 0., 500., 250);
        FillHist("MET_1jet_CF", eventbase->GetEvent().PFMET(), ev_weight, 0., 500., 250);
        FillHist("pT_ee_1jet_CF", X.Pt(), ev_weight, 0., 500., 250);
        FillHist("Nvtx_1jet_CF", eventbase->GetEvent().nVertices(), ev_weight, 0., 45., 45);
        FillHist("pT_e1_1jet_CF", electrons.at(0).Pt(), ev_weight, 0., 500., 250);
        FillHist("pT_e2_1jet_CF", electrons.at(1).Pt(), ev_weight, 0., 500., 250);
        FillHist("eta_e1_1jet_CF", electrons.at(0).Eta(), ev_weight, -3., 3., 60);
        FillHist("eta_e2_1jet_CF", electrons.at(1).Eta(), ev_weight, -3., 3., 60);
      }
      if(njet > 1){
        FillHist("mass_ee_2jet_CF", X.M(), ev_weight, 0., 500., 250);
        FillHist("MET_2jet_CF", eventbase->GetEvent().PFMET(), ev_weight, 0., 500., 250);
        FillHist("pT_ee_2jet_CF", X.Pt(), ev_weight, 0., 500., 250);
        FillHist("Njets_2jet_CF", jets.size(), ev_weight, 0., 7., 7);
        FillHist("Nvtx_2jet_CF", eventbase->GetEvent().nVertices(), ev_weight, 0., 45., 45);
        FillHist("pT_e1_2jet_CF", electrons.at(0).Pt(), ev_weight, 0., 500., 250);
        FillHist("pT_e2_2jet_CF", electrons.at(1).Pt(), ev_weight, 0., 500., 250);
        FillHist("eta_e1_2jet_CF", electrons.at(0).Eta(), ev_weight, -3., 3., 60);
        FillHist("eta_e2_2jet_CF", electrons.at(1).Eta(), ev_weight, -3., 3., 60);
      }
    }
    }
  }*/
/*   if(jets.size() > 3){
     if(nbjet > 0){
       if(muons.size() ==2) {
	 if(electrons.size() >= 1){
	   cout << "electrons is tight = " << electrons.at(0).PassTight() << endl;
	   if(!SameCharge(muons)){
	     if(muons.at(0).Pt() > 20. && muons.at(1).Pt() > 10.){
	       if(eventbase->GetEvent().PFMET() > 30){
		 if(trig_pass){
		   FillHist("Massmumu", GetDiLepMass(muons), ev_weight, 0., 200.,400);
		   FillHist("Massmumu_zoomed", GetDiLepMass(muons), ev_weight, 0.,50.,200);
		   FillCLHist(sighist_mm, "DiMuon", eventbase->GetEvent(), muons,electrons,jets, ev_weight);
		 }
	       }
	     }
	   }
	 }
       }
     }
   }*/
  return;
}// End of execute event loop
  


void DiElectron_SS::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void DiElectron_SS::BeginCycle() throw( LQError ){
  
  Message("In begin Cycle", INFO);
  
  //
  //If you wish to output variables to output file use DeclareVariable
  // clear these variables in ::ClearOutputVectors function
  //DeclareVariable(obj, label, treename );
  //DeclareVariable(obj, label ); //-> will use default treename: LQTree
  //  DeclareVariable(out_electrons, "Signal_Electrons", "LQTree");
  //  DeclareVariable(out_muons, "Signal_Muons");

  
  return;
  
}

DiElectron_SS::~DiElectron_SS() {
  
  Message("In DiElectron_SS Destructor" , INFO);
  
}

void DiElectron_SS::FillCutFlow(TString cut, float weight){

  if(GetHist("cutflow")) {
    GetHist("cutflow")->Fill(cut,weight);
  }
  else{
    AnalyzerCore::MakeHistograms("cutflow", 10, 0., 10.);

    GetHist("cutflow")->GetXaxis()->SetBinLabel(1,"NoCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(2,"EventCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(3,"VertexCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(4,"TriggerCut");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(5,"DiEl");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(6,"SameSign");
    GetHist("cutflow")->GetXaxis()->SetBinLabel(7,"Presel");
  }
}



void DiElectron_SS::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void DiElectron_SS::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this DiElectron_SSCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void DiElectron_SS::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}



