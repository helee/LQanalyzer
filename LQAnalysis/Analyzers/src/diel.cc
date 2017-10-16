// $Id: diel.cc 1 2013-11-26 10:23:10Z jalmond $
/***************************************************************************
 * @Project: LQdiel Frame - ROOT-based analysis framework for Korea SNU
 * @Package: LQCycles
 *
 * @author John Almond       <jalmond@cern.ch>           - SNU
 *
 ***************************************************************************/

/// Local includes
#include "diel.h"

//Core includes
#include "EventBase.h"                                                                                                                           
#include "BaseSelection.h"


//// Needed to allow inheritance for use in LQCore/core classes
ClassImp (diel);


 /**
  *   This is an Example Cycle. It inherits from AnalyzerCore. The code contains all the base class functions to run the analysis.
  *
  */
diel::diel() :  AnalyzerCore(), out_muons(0)  {
  
  
  // To have the correct name in the log:                                                                                                                            
  SetLogName("diel");
  
  Message("In diel constructor", INFO);
  //
  // This function sets up Root files and histograms Needed in ExecuteEvents
  InitialiseAnalysis();

}


void diel::InitialiseAnalysis() throw( LQError ) {
  
  /// Initialise histograms
  MakeHistograms();  
  // You can out put messages simply with Message function. Message( "comment", output_level)   output_level can be VERBOSE/INFO/DEBUG/WARNING 
  // You can also use m_logger << level << "comment" << int/double  << LQLogger::endmsg;
  
  Message("Making clever hists for Z ->ll test code", INFO);
  
  /// only available in v7-6-X branch and newer
  //// default lumimask is silver ////
  //// In v7-6-2-(current) the default is changed to gold (since METNoHF bug)
  ///When METNoHF isfixed the default will be back to silver
  /// set to gold if you want to use gold json in analysis
  /// To set uncomment the line below:

  return;
}


void diel::ExecuteEvents()throw( LQError ){

  if(k_running_nonprompt){
    while(!fake_configured){
      std::map<TString, std::pair<std::pair<TString,TString>  ,std::pair<float,TString> > >fake_hists;
      fake_hists["fr_muon_central"] = std::make_pair(std::make_pair("Muon_Data_v7_SIP3_FR.root","Muon_Data_v7_SIP3_FR_Awayjet40"), std::make_pair(70., "TH2D"));
      fake_hists["fr_electron_central"] = std::make_pair(std::make_pair("Electron_Data_v7_FR.root","Electron_Data_v7_FR_Awayjet40") , std::make_pair(70., "TH2D"));

      ConfigureFakeHists("/data1/LQAnalyzer_rootfiles_for_analysis/CATAnalysis2016/Fake/DiLep/", fake_hists);
    }
  }
  /// Set weight for MC (nevents / xsec) 
  if(!isData) weight *= MCweight;
  FillHist("GenWeight", weight, 1., -2., 2., 20);
    
  m_logger << DEBUG << "RunNumber/Event Number = "  << eventbase->GetEvent().RunNumber() << " : " << eventbase->GetEvent().EventNumber() << LQLogger::endmsg;
  m_logger << DEBUG << "isData = " << isData << LQLogger::endmsg;
   
  if(isData) FillHist("Nvtx_nocut_data", eventbase->GetEvent().nVertices(), weight, 0., 50., 50);
  else FillHist("Nvtx_nocut_mc", eventbase->GetEvent().nVertices(), weight, 0., 50., 50);

  FillHist("cutflow", 0.5, 1., 0., 10., 10);

  if(!PassMETFilter()) return;     /// Initial event cuts : 
  if (!eventbase->GetEvent().HasGoodPrimaryVertex()) return; //// Make cut on event wrt vertex                                                                                 
  TString trig = "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v";

  TString elid = "ELECTRON_HN_TIGHTv4";
  if(k_running_nonprompt){ elid = "ELECTRON_HN_FAKELOOSEv7"; }
  std::vector<snu::KElectron> electrons = GetElectrons(false, false, elid);
  std::vector<snu::KElectron> electrons_veto = GetElectrons("ELECTRON_HN_VETO");
  std::vector<snu::KMuon> muons = GetMuons("MUON_HN_TIGHT"); 
  std::vector<snu::KMuon> muons_veto = GetMuons("MUON_HN_VETO", true);

  std::vector<snu::KJet> jets = GetJets("JET_HN");
  std::vector<snu::KJet> jets5 = GetJets("JET_HN", 10., 5.0);
  int njet = jets.size();
  int nbjet = NBJet(GetJets("JET_NOLEPTONVETO", 20., 2.5));
  int njet5 = jets5.size();
  int njetopt[5][6] = {{0}};

  snu::KParticle jetopt[5][6][njet5];
  double x[5] = {10., 15., 20., 25., 30.};
  double y[6] = {4.5, 4.6, 4.7, 4.8, 4.9, 5.0};
  TString SS1jet_pT_fwd[5][6]; TString SS1jet_eta_fwd[5][6];
  TString Tch_Nevent[5][6]; TString Tch_2jet_Nevent[5][6]; TString Tch_pT_fwd[5][6]; TString Tch_eta_fwd[5][6];
  TString Tch_Nevent_unweight[5][6]; TString Tch_2jet_Nevent_unweight[5][6];
  int ptcut[5] = {10, 15, 20, 25, 30};
  int etacut[6] = {45, 46, 47, 48, 49, 50};

  for(unsigned int ij1=0; ij1<5; ij1++){
    for(unsigned int ij2=0; ij2<6; ij2++){
      int ijet = 0;
      for(int ij0=0; ij0<njet5; ij0++){
        if((jets5[ij0].Pt() > x[ij1]) && (fabs(jets5[ij0].Eta()) < y[ij2])){ jetopt[ij1][ij2][ijet] = jets5[ij0]; ijet++; }
      }
      njetopt[ij1][ij2] = ijet;
    }
  }

  double puweight = 1.;
  double reco_sf = 1.;
  double id_iso_sf = 1.;
  double trigger_sf = 1.;
  double trigger_ps = 1.;
  double ev_weight = 1.;

  if(!isData){
    puweight = eventbase->GetEvent().PileUpWeight(snu::KEvent::down);
    reco_sf = mcdata_correction->ElectronRecoScaleFactor(electrons);
    id_iso_sf = mcdata_correction->ElectronScaleFactor(elid, electrons, 0);
    trigger_ps = WeightByTrigger(trig, TargetLumi);
    std::vector<snu::KMuon> mu;
    double trigger_eff_Data = mcdata_correction->TriggerEfficiencyLegByLeg(electrons, elid, mu, "muid", 1, 0, 0);
    double trigger_eff_MC = mcdata_correction->TriggerEfficiencyLegByLeg(electrons, elid, mu, "muid", 1, 1, 0);
    trigger_sf = trigger_eff_Data/trigger_eff_MC;
    ev_weight = weight*puweight*reco_sf*id_iso_sf*trigger_sf*trigger_ps*GetKFactor()*MC_CR_Correction(0);
  }
  
  CorrectedMETRochester(muons);
  double met = eventbase->GetEvent().PFMET();
  double ST = eventbase->GetEvent().PFMET();
  double met2st = met*met/ST;  

  bool charge = false;

  FillHist("cutflow", 1.5, 1., 0., 10., 10);

  if(PassTrigger(trig)){
    FillHist("cutflow", 2.5, 1., 0., 10., 10);
    if(electrons.size() == 2){      
    FillHist("cutflow", 3.5, 1., 0., 10., 10);
      if((!k_running_chargeflip) && (electrons[0].Charge() == electrons[1].Charge())){ charge = true; }
      if(k_running_chargeflip && (electrons[0].Charge()*electrons[1].Charge() < 0)){ charge = true; }
      FillHist("charge_lep1", electrons[0].Charge()*electrons[1].Charge(), 1., -2., 2., 4);
      if(charge){
        FillHist("charge_lep2", electrons[0].Charge()*electrons[1].Charge(), 1., -2., 2., 4);
        FillHist("cutflow", 4.5, 1., 0., 10., 10);
        if(k_running_nonprompt){
          ev_weight = m_datadriven_bkg->Get_DataDrivenWeight_EE(false, electrons, "ELECTRON_HN_TIGHTv4", "ptcone", "fr_electron_central",0);
        }
        if(k_running_chargeflip){
          ShiftElectronEnergy(electrons, "ELECTRON_HN_TIGHTv4", true);
          ev_weight = GetCFweight(0, electrons, true, "ELECTRON_HN_TIGHTv4");
        }
        if((electrons[0].Pt() > 25.) && (electrons[1].Pt() > 15.)){
          FillHist("cutflow", 5.5, 1., 0., 10., 10);
          if((muons_veto.size() == 0) && (electrons_veto.size() == 2)){
            FillHist("cutflow", 6.5, 1., 0., 10., 10);
            snu::KParticle X = electrons[0] + electrons[1];
            if((X.M() > 10.) && (fabs(X.M()-91.1876) > 10.)){
              FillHist("cutflow", 7.5, 1., 0., 10., 10);
              ST = ST +  electrons[0].Pt() + electrons[1].Pt();
              for(unsigned int j=0; j<njet; j++){ ST += jets[j].Pt(); }
              met2st = met*met/ST;

              if(njet == 1){
                for(unsigned int il1=0; il1<5; il1++){
                  for(unsigned int il2=0; il2<6; il2++){
                    if(njetopt[il1][il2] > 0){
                      SS1jet_pT_fwd[il1][il2] = Form("SS1jet_pT_fwd_pT_%d_eta_%d", ptcut[il1], etacut[il2]);
                      SS1jet_eta_fwd[il1][il2] = Form("SS1jet_eta_fwd_pT_%d_eta_%d", ptcut[il1], etacut[il2]);
                      FillHist("SS1jet_pT_lep1", electrons[0].Pt(), ev_weight, 0., 500., 500);
                      FillHist("SS1jet_pT_lep2", electrons[1].Pt(), ev_weight, 0., 500., 500);
                      FillHist("SS1jet_eta_lep1", electrons[0].Eta(), ev_weight, -2.5, 2.5, 50);
                      FillHist("SS1jet_eta_lep2", electrons[1].Eta(), ev_weight, -2.5, 2.5, 50);
                      for(unsigned int il3=0; il3<njetopt[il1][il2]; il3++){
                        if(fabs(jetopt[il1][il2][il3].Eta()) > 2.5){
                          FillHist(SS1jet_pT_fwd[il1][il2], jetopt[il1][il2][il3].Pt(), ev_weight, 0., 300., 300);
                          FillHist(SS1jet_eta_fwd[il1][il2], jetopt[il1][il2][il3].Eta(), ev_weight, -5., 5., 100);
                        }
                      }
                    }
                  }
                }
              }

              if(njet > 1){
                FillHist("Presel_Nevent", 0.5, ev_weight, 0., 2., 2);
                FillHist("cutflow", 8.5, 1., 0., 10., 10);
                if(nbjet > 0) return;
                double wmass = 10000.; int j1 = 0; int j2 = 0;
                for(unsigned int k=0; k<njet; k++){
                  for(unsigned int l=k+1; l<njet; l++){
                    snu::KParticle JJ = jets[k] + jets[l];
                    if(fabs(JJ.M()-80.4) < wmass){ wmass = fabs(JJ.M()-80.4); j1=k; j2=l; }
                  }
                }

                snu::KParticle ll = electrons[0] + electrons[1];
                snu::KParticle l1jj = electrons[0] + jets[j1] + jets[j2];
                snu::KParticle l2jj = electrons[1] + jets[j1] + jets[j2];
                snu::KParticle lljj = electrons[0] + electrons[1] + jets[j1] + jets[j2];
                snu::KParticle W = jets[j1] + jets[j2];

                if((met2st < 15.) && (W.M() < 150.)){
                  FillHist("cutflow", 9.5, 1., 0., 10., 10);
                  FillHist("Sch_Nevent", 0.5, ev_weight, 0., 2., 2);
                  FillHist("Sch_pT_lep1", electrons[0].Pt(), ev_weight, 0., 500., 500);
                  FillHist("Sch_pT_lep2", electrons[1].Pt(), ev_weight, 0., 500., 500);
                  FillHist("Sch_pT_jet1", jets5[0].Pt(), ev_weight, 0., 500., 500);
                  FillHist("Sch_pT_jet2", jets5[1].Pt(), ev_weight, 0., 500., 500);
                  FillHist("Sch_pT_ll", ll.Pt(), ev_weight, 0., 500., 500);
                  FillHist("Sch_eta_lep1", electrons[0].Eta(), ev_weight, -2.5, 2.5, 50);
                  FillHist("Sch_eta_lep2", electrons[1].Eta(), ev_weight, -2.5, 2.5, 50);
                  FillHist("Sch_mass_ll", ll.M(), ev_weight, 0., 1500., 1500);
                  FillHist("Sch_mass_l1jj", l1jj.M(), ev_weight, 0., 1500., 1500);
                  FillHist("Sch_mass_l2jj", l2jj.M(), ev_weight, 0., 1500., 1500);
                  FillHist("Sch_mass_lljj", lljj.M(), ev_weight, 0., 1500., 1500);

                  bool fwd1 = false, fwd2 = false;
                  for(unsigned int ik1=0; ik1<5; ik1++){
                    for(unsigned int ik2=0; ik2<6; ik2++){
                      if(njetopt[ik1][ik2] > 0){
                        Tch_Nevent[ik1][ik2] = Form("Tch_Nevent_pT_%d_eta_%d", ptcut[ik1], etacut[ik2]);
                        Tch_Nevent_unweight[ik1][ik2] = Form("Tch_Nevent_unweight_pT_%d_eta_%d", ptcut[ik1], etacut[ik2]);
                        Tch_pT_fwd[ik1][ik2] = Form("Tch_pT_fwd_pT_%d_eta_%d", ptcut[ik1], etacut[ik2]);
                        Tch_eta_fwd[ik1][ik2] = Form("Tch_eta_fwd_pT_%d_eta_%d", ptcut[ik1], etacut[ik2]);
                        for(unsigned int ik4=0; ik4<njetopt[ik1][ik2]; ik4++){
                          if(fabs(jetopt[ik1][ik2][ik4].Eta()) > 2.5){ 
                            FillHist(Tch_Nevent[ik1][ik2], 0.5, ev_weight, 0., 2., 2);
                            FillHist(Tch_Nevent_unweight[ik1][ik2], 0.5, 1., 0., 2., 2);
                            FillHist(Tch_pT_fwd[ik1][ik2], jetopt[ik1][ik2][ik4].Pt(), ev_weight, 0., 300., 300);
                            FillHist(Tch_eta_fwd[ik1][ik2], jetopt[ik1][ik2][ik4].Eta(), ev_weight, -5., 5., 100);
                          }
                        }
                      }
                      if(njetopt[ik1][ik2] > 1){
                        Tch_2jet_Nevent[ik1][ik2] = Form("Tch_2jet_Nevent_pT_%d_eta_%d", ptcut[ik1], etacut[ik2]);
                        Tch_2jet_Nevent_unweight[ik1][ik2] = Form("Tch_2jet_Nevent_unweight_pT_%d_eta_%d", ptcut[ik1], etacut[ik2]);
                        for(unsigned int ik3=0; ik3<njetopt[ik1][ik2]; ik3++){
                          if(jetopt[ik1][ik2][ik3].Eta() > 2.5){ fwd1 = true; }
                          if(jetopt[ik1][ik2][ik3].Eta() < -2.5){ fwd2 = true; }
                        }
                        if(fwd1 && fwd2){
                          FillHist(Tch_2jet_Nevent[ik1][ik2], 0.5, ev_weight, 0., 2., 2);
                          FillHist(Tch_2jet_Nevent_unweight[ik1][ik2], 0.5, 1., 0., 2., 2);
                        }
                      }
                    }
                  }
                }
              }
            } 
          }
        }
      }
    }
  }
  return;
}// End of execute event loop
  


void diel::EndCycle()throw( LQError ){
  
  Message("In EndCycle" , INFO);

}


void diel::BeginCycle() throw( LQError ){
  
  Message("In begin Cycle", INFO);
   
  //
  //If you wish to output variables to output file use DeclareVariable
  // clear these variables in ::ClearOutputVectors function
  //DeclareVariable(obj, label, treename );
  //DeclareVariable(obj, label ); //-> will use default treename: LQTree
  //  DeclareVariable(out_electrons, "Signal_Electrons", "LQTree");
  //  DeclareVariable(out_muons, "Signal_Muons");
  ConfigureFake();
  
  return;
  
}

diel::~diel() {
  
  Message("In diel Destructor" , INFO);
  
}


void diel::BeginEvent( )throw( LQError ){

  Message("In BeginEvent() " , DEBUG);

  return;
}



void diel::MakeHistograms(){
  //// Additional plots to make
    
  maphist.clear();
  AnalyzerCore::MakeHistograms();
  Message("Made histograms", INFO);
  /**
   *  Remove//Overide this dielCore::MakeHistograms() to make new hists for your analysis
   **/
  
}


void diel::ClearOutputVectors() throw(LQError) {

  // This function is called before every execute event (NO need to call this yourself.
  
  // Add any new output vector you create to this list. 
  // if you do not the vector will keep keep getting larger when it is filled in ExecuteEvents and will use excessive amoun of memory
  //
  // Reset all variables declared in Declare Variable
  //
  out_muons.clear();
  out_electrons.clear();
}



