
//This is the latest Jan 7, 2022 version with no hgc, but all eff corrected, livetime, h(p)tracking are done run by run. Lets see the diff in Yield as compared to last time when the tracking, TLT were also slope based. 

//my simc yield is 25% higher, so relaxing some cuts
//note most of the cuts are already applied in the rootfile that is used here:
//updated code to apply geom cuts July 27,2022
#include <TSystem.h>
#include <TString.h>
#include "TFile.h"
#include "TTree.h"
#include <TNtuple.h>
#include "TCanvas.h"
#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TH1F.h"
#include <TH2.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TMath.h>
#include <TLegend.h>
#include <TPaveLabel.h>
#include <TProfile.h>
#include <TPolyLine.h>
#include <TObjArray.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
using namespace std;
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <istream>
#include <vector>
#include <cmath>
#include <ios>
#include <iosfwd>
#include <iomanip>
#include <streambuf>

#include "nlohmann/json.hpp"
using json = nlohmann::json;
#include "ROOT/RDataFrame.hxx"
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/VectorUtil.h"
#include "ROOT/RVec.hxx"
#include "TVector3.h"
#include "ROOT/RSnapshotOptions.hxx"


void ptsidis_skim_coutRootfile_geom_cuts(){
  gStyle->SetPalette(1,0);
  gStyle->SetOptStat(1000011);
  gStyle->SetOptFit(11);
  gStyle->SetTitleOffset(1.,"Y");
  gStyle->SetTitleOffset(.7,"X");
  gStyle->SetLabelSize(0.04,"XY");
  gStyle->SetTitleSize(0.06,"XY");
  gStyle->SetPadLeftMargin(0.12);
  gROOT->SetBatch("kTRUE");

  ofstream txtfile("txtfile_current_cut/yield_jin_jout.txt",ios::app);

  json j_info;//input
  {
    std::ifstream ifs("/u/group/c-csv/hdbhatt/my_analysis/ptsidis_runsinfo/db2/ptsidis_csv_eff_april2.json");
    //std::ifstream ifs("/u/group/c-csv/hdbhatt/yielddec21/json_files/all_eff_corrJan2022.json");

    //std::ifstream ifs("/u/group/c-csv/hdbhatt/yielddec21/json_files/run_charge_curr_eff_run7642.json");
    // std::ifstream ifs("/u/group/c-csv/hdbhatt/yielddec21/json_files/short_all_eff_corrJan2022.json");
    ifs>>j_info;
  }
  json j_out1;//for output json at the end of each run
  json j_out2;//for output json at the end of whole script
  for (auto it = j_info.begin();it!=j_info.end();it++){

    int RunNumber = std::stoi(it.key());
    double efficiency = j_info[std::to_string(RunNumber).c_str()]["ALL_CORR"].get<double>();//=======ok??
    double charge = j_info[std::to_string(RunNumber).c_str()]["ChargeCorr"].get<double>();//=======ok??
    double charge_times_eff = efficiency*charge;
    cout<<RunNumber<<"\t\t"<<charge<<" * "<<efficiency<<" = "<<charge_times_eff<<endl;
    //peter ptc.f, july 27, 2022
    Double_t hcal_1pr_zpos = 338.69;
    Double_t hcal_2ta_zpos = 349.69;
    Double_t hcal_3ta_zpos = 360.69;
    Double_t hcal_4ta_zpos = 371.69;
    Double_t hcal_left     =  35.00;//	!actual size of calorimeter
    Double_t hcal_right    = -35.00;
    Double_t hcal_top      = -69.66;
    Double_t hcal_bottom   =  60.34;
    //SHMS Calorimeter position

    Double_t scal_1pr_zpos = 292.6;
    Double_t scal_2ta_zpos = 306.4;
    Double_t scal_3ta_zpos = 323.7;
    Double_t scal_4ta_zpos = 341.0;
    Double_t scal_left     =  63.00;
    Double_t scal_right    = -63.00;
    Double_t scal_top      = -70.00;
    Double_t scal_bottom   =  70.00;


    
    // TString filename =  Form("/lustre19/expphy/volatile/hallc/c-csv/hdbhatt/ROOTfiles/skimdataROOTfiles/skimROOTfilesjuly27/run_%d.root",RunNumber);
    TString filename =  Form("/lustre19/expphy/volatile/hallc/c-csv/hdbhatt/ROOTfiles/skimROOTfilesjuly27/run_%d.root",RunNumber);
    TFile *f = new TFile(filename,"READ");
    if(f->IsZombie()){
      std::cout << "Error opening file "<< std::endl;
      continue;
    }
   
    //TFile *f1 = new TFile(Form("yield_root_20bins/yield_hist_%d.root",RunNumber),"RECREATE");

    
    cout<< "Opening Run Num = "<< RunNumber <<endl; 
    if(!it.value()["ChargeCorr"].empty()){
   
      /* if(f->IsZombie()){
	 std::cout << "Error opening file "<< std::endl;
	 exit(-1);
	 }
      */

      TH1F* h_pcal_clean = new TH1F("h_pcal_clean","hist_pcal_clean",100, 0,2);
      TH1F* h_pcal_coin = new TH1F("h_pcal_coin","h_pcal_coin",100, 0,2);
      TH1F* h_pcal_acc = new TH1F("h_pcal_acc","h_pcal_acc",100, 0,2);
      TH1F* h_ctime = new TH1F("h_ctime","h_ctime",300,0,100);
      TH1F* h_ctime_coin = new TH1F("h_ctime_coin","h_ctime_coin",300,0,100);
      TH1F* h_ctime_acc = new TH1F("h_ctime_acc","h_ctime_acc",300,0,100);
      TH1F* h_ctime_clean = new TH1F("h_ctime_clean"," h_ctime_clean",300,0,100);
      TH1F* h_mom = new TH1F("h_mom"," h_mom",70, 0, 7);
      TH2F* h_aero_y_vs_x = new TH2F("h_aero_y_vs_x","Aero_y_vs_x;yAtAero;xAtAero",100,-100,100,100,-100,100);
      TH2F* h_showerback_y_vs_x = new TH2F("h_showerback_y_vs_x","Showerback_y_vs_x;yAtCal;xAtCal",100,-70,70,100,-70,70);



      TH1F* h_hStartTime = new TH1F("h_hStartTime",";HMS Starttime",200,0,200.);
      TH1F* h_hStartTime_track = new TH1F("h_hStartTime_track",";HMS Starttime",200,0,200.);
      TH1F* h_hFpTime = new TH1F("h_hFpTime",";HMS Fptime",200,0,200.);
      TH1F* h_etotnorm = new TH1F("h_etotnorm"," ; Etot norm",100,0.,2.);
      TH1F* h_npeSum = new TH1F("h_npeSum"," ; Cer Npe",100,0.,20.);
      TH1F* h_ev2 = new TH1F("h_ev2"," ; hEv2",10,0.,10.);
      TH1F* h_goodev2 = new TH1F("h_goodev2"," ; hgoodEv2",10,0.,10.);
      TH1F* h_goodev = new TH1F("h_goodev"," ; Good ev2 with PID cuts",10,0.,10.);
      TH1F* h_goodevTrack = new TH1F("h_goodevTrack"," ; Good ev2 with PID/track cuts",10,0.,10.);
     
      TH1F *hist_zhadprime = new TH1F("hist_zhadprime","zhadron after targetmass corr;zhadron;Entries",20,0,1);
      TH1F *hist_xbjprime = new TH1F("hist_xbjprime","xBj after targetmass corr;xBj;Entries",20,0,1);
      TH1F *hist_zhadprime_acc = new TH1F("hist_zhadprime_acc","zhadron after targetmass corr;zhadron;Entries",20,0,1);
      TH1F *hist_xbjprime_acc = new TH1F("hist_xbjprime_acc","xBj after targetmass corr;xBj;Entries",20,0,1);
      

      TH1F* hist_pt = new TH1F("hist_pt"," Pt Plot;Pt(GeV);Entries",20, 0,1);
      TH1F* hist_pt_mom = new TH1F("hist_pt_mom"," Pt Plot;Pt(GeV);Entries",20, 0,1);
      TH1F* hist_pt_acc = new TH1F("hist_pt_acc"," Pt Plot;Pt(GeV);Entries",20, 0,1);
      TH1F* hist_pt_mom_acc = new TH1F("hist_pt_mom_acc"," Pt Plot;Pt(GeV);Entries",20, 0,1);
      TH1F *h_pt_clean_mom_counts =  new TH1F("h_pt_clean_mom_counts","Zhad Clean Data;Pt;Entries",20,0,1);


      TH1F* hist_phi = new TH1F("hist_phi"," Phi Plot;Phi(Rad);Entries",20, -1,7);
      TH1F* hist_phi_mom = new TH1F("hist_phi_mom"," Phi Plot;Phi(Rad);Entries",20, -1,7);
      TH1F* hist_phi_acc = new TH1F("hist_phi_acc"," Phi Plot;Phi(Rad);Entries",20, -1,7);
      TH1F* hist_phi_mom_acc = new TH1F("hist_phi_mom_acc"," Phi Plot;Phi(Rad);Entries",20, -1,7);
      TH1F *h_phi_clean_mom_counts =  new TH1F("h_phi_clean_mom_counts","Zhad Clean Data;Phi (Rad);Entries",20,-1,7);



      TH1F *hist_zhad = new TH1F("hist_zhad","zhadron before targetmass corr;zhadron;Entries",20,0,1);
      TH1F *hist_xbj = new TH1F("hist_xbj","xBj before targetmass corr;xBj;Entries",20,0,1);
      TH1F *hist_zhad_acc = new TH1F("hist_zhad_acc","zhadron before targetmass corr;zhadron;Entries",20,0,1);
      TH1F *hist_xbj_acc = new TH1F("hist_xbj_acc","xBj before targetmass corr;xBj;Entries",20,0,1);
      
      TH1F *hist_zhad_mom = new TH1F("hist_zhad_mom","zhad before targetmass corr;zhadron;Entries",20,0,1);
      TH1F *hist_xbj_mom = new TH1F("hist_xbj_mom","xBj before targetmass corr;xBj;Entries",20,0,1);
      TH1F *hist_zhad_mom_acc = new TH1F("hist_zhad_mom_acc","zhad before targetmass corr;zhadron;Entries",20,0,1);
      TH1F *hist_xbj_mom_acc = new TH1F("hist_xbj_mom_acc","xBj before targetmass corr;xBj;Entries",20,0,1);
      TH1F *h_xbj_clean_mom_counts =  new TH1F("h_xbj_clean_mom_counts","Xbj Clean Data;Xbj;Entries)",20,0,1);
      TH1F *h_zhad_clean_mom_counts =  new TH1F("h_zhad_clean_mom_counts","Zhad Clean Data;Zhad;Entries",20,0,1);
      TH1F *h_zhadprime_clean_mom_counts =  new TH1F("h_zhadprime_clean_mom_counts","Zhad Clean Data;Zhad;Entries",20,0,1);



  
  
      Double_t pionmass   = 0.1395701835; 
      Double_t protonmass = 0.93827231;//GeV/c2
      Double_t pi = 3.1416;

      TTree *tt = (TTree*) f->Get("T");
      TTree *ss = (TTree*) f->Get("TSP");
      //
      Double_t  Scal_evNumber;
      ss->SetBranchAddress("evNumber",&Scal_evNumber);
      Double_t  Scal_BCM1_charge;
      ss->SetBranchAddress("P_BCM1_scalerCharge",&Scal_BCM1_charge);
      Double_t  Scal_BCM1_current;
      ss->SetBranchAddress("P_BCM1_scalerCurrent",&Scal_BCM1_current);
      Double_t  Scal_time;
      ss->SetBranchAddress("P_1MHz_scalerTime",&Scal_time);
      Double_t  Scal_EDTM;
      ss->SetBranchAddress("P_EDTM_scaler",&Scal_EDTM);
      Double_t  Scal_TRIG2;
      ss->SetBranchAddress("P_pTRIG2_scaler",&Scal_TRIG2);
      Double_t  Scal_TRIG3;
      ss->SetBranchAddress("P_pTRIG3_scaler",&Scal_TRIG3);
      Double_t  Scal_TRIG1;
      ss->SetBranchAddress("P_pTRIG1_scaler",&Scal_TRIG1);
      Double_t  Scal_TRIG4;
      ss->SetBranchAddress("P_pTRIG4_scaler",&Scal_TRIG4);
      Double_t  Scal_TRIG5;
      ss->SetBranchAddress("P_pTRIG5_scaler",&Scal_TRIG5);
      Double_t  Scal_TRIG6;
      ss->SetBranchAddress("P_pTRIG6_scaler",&Scal_TRIG6);
  
  
      //loop through scalers
      Int_t nscal_reads=0;
      Int_t nscal_reads_cut=0;
      Double_t prev_read=-1;
      Double_t charge_sum=0;
      Double_t charge_sum_cut=0;
      Double_t prev_charge=0;
      Double_t charge_sum_corr=0;
      Double_t charge_sum_cut_corr=0;
      Double_t prev_charge_corr=0;
      Double_t event_flag[100000];
      Double_t scal_event_number[100000];
      Double_t tot_scal_EDTM=0;
      Double_t tot_scal_cut_EDTM=0;
      Double_t prev_EDTM=0;
      Double_t tot_scal_hEL_CLEAN=0;
      Double_t tot_scal_cut_hEL_CLEAN=0;
      Double_t prev_hEL_CLEAN=0;
      Double_t tot_scal_hEL_REAL=0;
      Double_t tot_scal_cut_hEL_REAL=0;
      Double_t prev_hEL_REAL=0;
      Double_t tot_scal_TRIG2=0;
      Double_t tot_scal_TRIG3=0;
      Double_t prev_TRIG2=0;
      Double_t prev_TRIG3=0;
      Double_t tot_scal_cut_TRIG2=0;
      Double_t tot_scal_cut_TRIG3=0;
      Double_t tot_scal_TRIG1=0;
      Double_t tot_scal_TRIG4=0;
      Double_t prev_TRIG1=0;
      Double_t prev_TRIG4=0;
      Double_t tot_scal_cut_TRIG1=0;
      Double_t tot_scal_cut_TRIG4=0;
      Double_t tot_scal_cut_time=0;
      Double_t tot_scal_TRIG5=0;
      Double_t prev_TRIG5=0;
      Double_t tot_scal_cut_TRIG5=0;
      Double_t tot_scal_TRIG6=0;
      Double_t prev_TRIG6=0;
      Double_t tot_scal_cut_TRIG6=0;
      Double_t threshold_cut=0.;
      Double_t ctmin = 2;
      Double_t ctmax = 2;
      //
      Double_t tot_scal_time=0;
      Double_t prev_time=0;
      //
      Long64_t scal_entries = ss->GetEntries();
      //cout << " scal ent = " << scal_entries << endl;
      Double_t nlast=float(scal_entries);
      TH1F *h_cur_entry = new TH1F("h_cur_entry","; ENtry;current",nlast,0,nlast);
      TH1F *h_cur = new TH1F("h_cur","; Current ;",200,0,100);
      Long64_t data_entries = tt->GetEntries();
      for (int i = 0; i < scal_entries; i++) {
	ss->GetEntry(i);
	h_cur_entry->Fill(float(i),Scal_BCM1_current);
	if (Scal_BCM1_current > 3) h_cur->Fill(Scal_BCM1_current);
      }
   
      Double_t peak_current = h_cur->GetBinCenter(h_cur->GetMaximumBin());
      //cout << " Peak current = " << peak_current  <<" " <<  h_cur->GetMaximumBin() << endl;
      Double_t Scal_BCM1_charge_corr=0;
      for (int i = 0; i < scal_entries; i++) {
	ss->GetEntry(i);
	event_flag[nscal_reads] = 0;
	scal_event_number[nscal_reads] = Scal_evNumber;
	Double_t BCM1_correction=1.;
	if (Scal_BCM1_current >2.) {
	  if (Scal_BCM1_current <= 60) {
	    BCM1_correction =1.0 + 0.045* ( log(60.)-log(Scal_BCM1_current))/( log(60.)-log(2.) );
	  } else {
	    BCM1_correction =1.0 + 0.010*(Scal_BCM1_current-60)/25.;
	  } 
	}
	Scal_BCM1_charge_corr+=Scal_BCM1_current*(Scal_time-prev_time)*BCM1_correction;
	///cout << Scal_BCM1_charge << " "  << Scal_BCM1_charge_corr << " " << BCM1_correction << endl;
	//if (TMath::Abs(Scal_BCM1_current-peak_current) < threshold_cut) {
	if (TMath::Abs(Scal_BCM1_current) > threshold_cut) {


	  event_flag[nscal_reads] = 1;
	  tot_scal_cut_time+=(Scal_time-prev_time);
	  tot_scal_cut_TRIG2+=(Scal_TRIG2-prev_TRIG2);
	  tot_scal_cut_TRIG3+=(Scal_TRIG3-prev_TRIG3);
	  tot_scal_cut_TRIG1+=(Scal_TRIG1-prev_TRIG1);
	  tot_scal_cut_TRIG4+=(Scal_TRIG4-prev_TRIG4);
	  tot_scal_cut_TRIG5+=(Scal_TRIG5-prev_TRIG5);
	  tot_scal_cut_TRIG6+=(Scal_TRIG6-prev_TRIG6);

	  charge_sum_cut+=(Scal_BCM1_charge-prev_charge);
	  charge_sum_cut_corr+=(Scal_BCM1_charge_corr-prev_charge_corr);
	  nscal_reads_cut++;
	}
	prev_charge = Scal_BCM1_charge;
	prev_charge_corr = Scal_BCM1_charge_corr;
	prev_time = Scal_time;
	prev_TRIG2 = Scal_TRIG2;
	prev_TRIG3 = Scal_TRIG3;
	prev_TRIG1 = Scal_TRIG1;
	prev_TRIG4 = Scal_TRIG4;
	prev_TRIG5 = Scal_TRIG5;
	prev_TRIG6 = Scal_TRIG6;

	nscal_reads++;

	charge_sum=Scal_BCM1_charge;
	charge_sum_corr=Scal_BCM1_charge_corr;
	tot_scal_TRIG2=Scal_TRIG2;
	tot_scal_TRIG3=Scal_TRIG3;
	tot_scal_TRIG1=Scal_TRIG1;
	tot_scal_TRIG4=Scal_TRIG4;
	tot_scal_TRIG5=Scal_TRIG5;
	tot_scal_TRIG6=Scal_TRIG6;
	tot_scal_time=Scal_time;
      }
      //cout << "nscal_reads_cut " << nscal_reads_cut<< endl;
  

      //
      Double_t gevtyp;
      tt->SetBranchAddress("g_evtyp",&gevtyp);
      Double_t gevnum;
      tt->SetBranchAddress("g_evnum",&gevnum);
      Double_t mom;
      tt->SetBranchAddress("P_gtr_p",&mom);
      Double_t pbeta;
      tt->SetBranchAddress("P_gtr_beta",&pbeta);
      Double_t pdelta;
      tt->SetBranchAddress("P_gtr_dp",&pdelta);
      Double_t pcaletottrack;
      tt->SetBranchAddress("P_cal_etottracknorm",&pcaletottrack);
      Double_t pcaletrack;
      tt->SetBranchAddress("P_cal_etracknorm",&pcaletrack); 
      Double_t hgcer;
      tt->SetBranchAddress("P_hgcer_npeSum",&hgcer);
      Double_t hcaletottrack;
      tt->SetBranchAddress("H_cal_etottracknorm",&hcaletottrack);
      Double_t ctime;
      tt->SetBranchAddress("CTime_ePiCoinTime_ROC2",&ctime);
      double_t hgcer_xAtcer;
      tt->SetBranchAddress("P_hgcer_xAtCer", &hgcer_xAtcer);
      double_t hgcer_yAtcer;
      tt->SetBranchAddress("P_hgcer_yAtCer", &hgcer_yAtcer);
      Double_t ptrack_x_fp;
      tt->SetBranchAddress("P_dc_x_fp", &ptrack_x_fp);
      Double_t ptrack_y_fp;
      tt->SetBranchAddress("P_dc_y_fp", &ptrack_y_fp);
      Double_t ptrack_xp_fp;
      tt->SetBranchAddress("P_dc_xp_fp", &ptrack_xp_fp);
      Double_t ptrack_yp_fp;
      tt->SetBranchAddress("P_dc_yp_fp", &ptrack_yp_fp);
      Double_t dW2;
      tt->SetBranchAddress("H_kin_primary_W2",&dW2);
      Double_t dnu;
      tt->SetBranchAddress("H_kin_primary_omega",&dnu);
      Double_t dpmiss;
      tt->SetBranchAddress("P_kin_secondary_pmiss",&dpmiss);
      Double_t theta;
      tt->SetBranchAddress("P_kin_secondary_th_xq", &theta);
      Double_t dataphi;
      tt->SetBranchAddress("P_kin_secondary_ph_xq", &dataphi);
      Double_t xbj;  
      tt->SetBranchAddress("H_kin_primary_x_bj",&xbj);
      Double_t Q2;
      tt->SetBranchAddress("H_kin_primary_Q2",&Q2);
      Double_t htrack_x_fp;
      tt->SetBranchAddress("H_dc_x_fp", &htrack_x_fp);
      Double_t htrack_y_fp;
      tt->SetBranchAddress("H_dc_y_fp", &htrack_y_fp);
      Double_t htrack_xp_fp;
      tt->SetBranchAddress("H_dc_xp_fp", &htrack_xp_fp);
      Double_t htrack_yp_fp;
      tt->SetBranchAddress("H_dc_yp_fp", &htrack_yp_fp);


      

      //cout<<"good"<<endl;
      Int_t nscal_reads_2=0;
      for (int i = 0; i < data_entries; i++) {
	tt->GetEntry(i);
	//
	if (event_flag[nscal_reads_2]==1) {
	  Bool_t pidcut = pcaletottrack<0.8 && hcaletottrack > 0.8;
	  if (pidcut){
	    h_ctime->Fill(ctime);
	    h_mom->Fill(mom);
	  }
	}//evt flag
	if (gevnum > scal_event_number[nscal_reads_2])  nscal_reads_2++;
      }//entries
      //
    
      //accidental cut
      Double_t momentum = h_mom->GetMean();
      Int_t bin_max = h_ctime->GetMaximumBin();
      Double_t max_value = h_ctime->GetBinCenter(bin_max); 


      Bool_t accL_cut, accR_cut, accidental_cut, ctime_cut, coin_cut,acc_cut, central_hole, right_or_left_of_hole, remove_region, position_cut;
      Double_t xathgcer,yathgcer,mhp_square, xbjprime, zhadprime, pt, xatpreshower, yatpreshower, xatshower, yatshower,xatpreshowerback , yatpreshowerback, xatshowerback, yatshowerback, Epi, zhad, Mx2,xataero,yataero;;
      //geom cut related
      Bool_t peter_all_geom_cut, common_cuts,  pcal_geom_cut_hem,  paero_geom_cut,hcal_geom_cut,hms_dc_fp_cut;
      Bool_t  hms_hodo_cut, pcal_geom_cut, shms_dc_fp_cut, hourglass_cut, geom_cut_and_hgcer_cut;
      Double_t xathhodo,xatpcal,yatpcal,xathcal,yathcal;
      
      //cout<<"max_val = "<<max_value<<endl;
      Int_t nscal_reads_3=0;
      for (int i = 0; i < data_entries; i++) {
	tt->GetEntry(i);
	if (event_flag[nscal_reads_3]==1) {

	  accL_cut        =   ctime > (max_value - 18 ) && ctime < (max_value - 10);
	  accR_cut        =   ctime > (max_value + 14 ) && ctime < (max_value + 22);
	  accidental_cut = (accR_cut || accL_cut);
	  ctime_cut = ctime > (max_value - ctmax) && ctime < (max_value + ctmin);
	  coin_cut = ctime_cut && pcaletottrack<0.8 && hcaletottrack > 0.7;
	  acc_cut = accidental_cut && pcaletottrack<0.8 && hcaletottrack > 0.7;




	  xathgcer = ptrack_x_fp + 156.27  * ptrack_xp_fp;
	  yathgcer = ptrack_y_fp + 156.27  * ptrack_yp_fp;

	  central_hole = (pow(yathgcer-1.33, 2) +  pow(xathgcer-0.83, 2) >= pow(6.0, 2));//radius = 6
	 
	  right_or_left_of_hole = ((yathgcer > 1.33 && (xathgcer < 0. ||  xathgcer > 3.0) ) ||  (yathgcer < 1.33 && (xathgcer < 0.||  xathgcer > 3)  ));
	 
	  remove_region = central_hole && right_or_left_of_hole;//ok


	 
	  xatpreshower = ptrack_x_fp + 292.64 * (ptrack_xp_fp);
	  yatpreshower = ptrack_y_fp + 292.64 * (ptrack_yp_fp);

	  xatshower = ptrack_x_fp + 306.44 * (ptrack_xp_fp);
	  yatshower = ptrack_y_fp + 306.44 * (ptrack_yp_fp);

	  xatpreshowerback = ptrack_x_fp + 302.64 * (ptrack_xp_fp);
	  yatpreshowerback = ptrack_y_fp + 302.64 * (ptrack_yp_fp);

	  xatshowerback = ptrack_x_fp + 356.44 * (ptrack_xp_fp);
	  yatshowerback = ptrack_y_fp + 356.44 * (ptrack_yp_fp);
		
	  xataero = ptrack_x_fp + 231.0 * (ptrack_xp_fp);
	  yataero = ptrack_y_fp + 231.0 * (ptrack_yp_fp);

	  // my original pcal_geom_cut Not used HERE
	  pcal_geom_cut_hem =  xatshowerback>-60  &&  xatshowerback < 60  && yatshowerback>-60  && yatshowerback < 60 && xatpreshowerback  > -60 && xatpreshowerback  < 60  && yatpreshowerback >-60 && yatpreshowerback < 60;
	  //aero geom cut
	  paero_geom_cut =  xataero >-50 && xataero < 50 && yataero >-50 && yataero<50;

	  
	  //position_cut =  xatshowerback>-60  &&  xatshowerback < 60  && yatshowerback>-60  && yatshowerback < 60 && xatpreshowerback  > -60 && xatpreshowerback  < 60  && yatpreshowerback >-60 && yatpreshowerback < 60 && xataero >-50 && xataero < 50 && yataero >-50 && yataero<50;
	  ///geom cuts added July 27 START!!!
	  //peter line 1933
	  //hcal geom
	  xathcal = htrack_x_fp + hcal_4ta_zpos * htrack_xp_fp;
	  yathcal = htrack_y_fp + hcal_4ta_zpos * htrack_yp_fp;
	  hcal_geom_cut =  (yathcal <= (hcal_left-2.0)  &&  yathcal >= (hcal_right+2.0)   &&  xathcal <= (hcal_bottom-2.0) && xathcal >= (hcal_top+2.0));
	  //hms focal plane cuts
	  hms_dc_fp_cut = abs(htrack_x_fp)<=58;
	  //hms hodo cut
	  xathhodo = htrack_x_fp + 318 *  htrack_xp_fp;
	  hms_hodo_cut = abs(xathhodo)<=59;
	  //shms cal cut (peters)
	  xatpcal = ptrack_x_fp + scal_4ta_zpos *ptrack_xp_fp;
	  yatpcal = ptrack_y_fp + scal_4ta_zpos *ptrack_yp_fp;

	  pcal_geom_cut = (yatpcal <= (scal_left-2.0) &&  yatpcal >= (scal_right+2.0)) &&  (xatpcal <= (scal_bottom-2.0) &&   xatpcal>= (scal_top+2.0));

	  //shms focal plane cuts
	  shms_dc_fp_cut = abs(ptrack_x_fp) <=38 && abs(ptrack_y_fp) <=38;
	  //! Hourglass cut: 
	  hourglass_cut = (ptrack_y_fp <= (10 + abs(ptrack_x_fp)))  && (ptrack_y_fp >= (-10 - abs(ptrack_x_fp)));

      
	  //all geom cuts 
	  peter_all_geom_cut =  pcal_geom_cut && hcal_geom_cut && paero_geom_cut &&  shms_dc_fp_cut && hms_dc_fp_cut && hms_hodo_cut &&   hourglass_cut;

	  geom_cut_and_hgcer_cut =   peter_all_geom_cut  && remove_region;
      
	  //all geom cuts DONE!!!




	  ///geom cuts added July 27 DONE!!!


	  
	  Epi = sqrt(pow(pionmass,2) + pow(mom,2));
	  zhad = Epi/dnu;
	  Mx2 = (protonmass + dnu - zhad*dnu) *  (protonmass + dnu - zhad*dnu) - abs(dpmiss)*abs(dpmiss);
	  pt=mom*sin(theta);
	  xbjprime = 2*xbj/(1+sqrt(1+4*xbj*protonmass*protonmass/Q2));//target mass corrected x
	  mhp_square = pionmass*pionmass+pt*pt;
	  zhadprime = (zhad/2)*(xbjprime/xbj)*(1+sqrt(1-(4*xbj*xbj*protonmass*protonmass*mhp_square)/(zhad*zhad*Q2*Q2)));//target mass corrected zhad
	 
	  if(dataphi<0)dataphi=2*pi+dataphi;
	  if(dataphi>=0)dataphi=dataphi;

	  Bool_t common_cuts = Mx2 > 2.6 && dW2 > 4 &&  peter_all_geom_cut;//================This has all geomcuts
	  if (coin_cut && common_cuts)
	    {
	      h_aero_y_vs_x->Fill(yataero,xataero);
	      h_showerback_y_vs_x->Fill(yatpreshowerback,xatpreshowerback);
	      hist_zhad->Fill(zhad);
	      hist_pt->Fill(pt);
	     
	      hist_phi->Fill(dataphi);

	      hist_xbj->Fill(xbj);
	      hist_xbjprime->Fill(xbjprime);
	      //hist_zhadprime->Fill(zhadprime);
	      // hist_pt->Fill(pt);
	    }

	  if(mom>=3){
	    if (coin_cut && common_cuts && remove_region && hgcer > 0.5){
	      h_pcal_coin->Fill(pcaletottrack); 
	      h_ctime_coin->Fill(ctime);
	      hist_zhad_mom->Fill(zhad);
	      hist_pt_mom->Fill(pt);
	      hist_xbj_mom->Fill(xbj);
	      hist_zhadprime->Fill(zhadprime);
	      hist_phi_mom->Fill(dataphi);



	    }
	    if (acc_cut && common_cuts &&  remove_region && hgcer > 0.5){
	      h_pcal_acc->Fill(pcaletottrack); 
	      h_ctime_acc->Fill(ctime);
	      hist_zhad_mom_acc->Fill(zhad);
	      hist_pt_mom_acc->Fill(pt);
	      hist_phi_mom_acc->Fill(dataphi);
	      hist_xbj_mom_acc->Fill(xbj);
	      hist_zhadprime_acc->Fill(zhadprime);
	    }
	  }



	  if(mom < 3){
	    if (coin_cut && common_cuts ){
	      h_pcal_coin->Fill(pcaletottrack); 
	      h_ctime_coin->Fill(ctime); 
	      hist_zhad_mom->Fill(zhad);
	      hist_pt_mom->Fill(pt);
	      hist_phi_mom->Fill(dataphi);
	      hist_xbj_mom->Fill(xbj);	   
	      hist_zhadprime->Fill(zhadprime);

	    }
	    if (acc_cut && common_cuts){
	      h_pcal_acc->Fill(pcaletottrack); 
	      h_ctime_acc->Fill(ctime);
	      hist_zhad_mom_acc->Fill(zhad);
	      hist_pt_mom_acc->Fill(pt);
	      hist_phi_mom_acc->Fill(dataphi);
	      hist_xbj_mom_acc->Fill(xbj);
	      hist_zhadprime_acc->Fill(zhadprime);

	    }
	  }
	    
	}//evt flag
	if (gevnum > scal_event_number[nscal_reads_3])  nscal_reads_3++;
      }//entries
 
      //
   
      //acc subtraction just do math
      Double_t coin_integral = h_pcal_coin->Integral(); 
      Double_t acc_integral_4 = h_pcal_acc->Integral(); 
      Double_t acc_integral = acc_integral_4 *(ctmin+ctmax)/16; 
      Double_t clean_integral = coin_integral-acc_integral;

      //cout<<runs[r]<<"\t\t"<<"clean_integral =====>   "<<clean_integral<<endl;
      //cout<<runs[r]<<"\t\t"<<"coin_integral  =====>   "<<coin_integral<<endl;
      //cout<<runs[r]<<"\t\t"<<"acc_integral   =====>   "<<acc_integral<<endl;
   
      //verify these yields by actual subtraction of coin and acc hists
      TH1D *  h_pcal_acc_copy = (TH1D*)h_pcal_acc->Clone();
      h_pcal_acc_copy->Scale(0.25); 
      TH1D *  h_pcal_coin_copy = (TH1D*)h_pcal_coin->Clone();
      TH1D *  h_pcal_acc_only = (TH1D*)h_pcal_acc_copy->Clone();
      h_pcal_clean->Add(h_pcal_coin_copy,h_pcal_acc_only,1.0,-1.0);
    

      TH1D *  h_ctime_acc_copy = (TH1D*)h_ctime_acc->Clone();
      h_ctime_acc_copy->Scale(0.25); 
      TH1D *  h_ctime_coin_copy = (TH1D*)h_ctime_coin->Clone();
      TH1D *  h_ctime_acc_only = (TH1D*)h_ctime_acc_copy->Clone();
      h_ctime_clean->Add(h_ctime_coin_copy,h_ctime_acc_only,1.0,-1.0);

      //acc for xbj
      TH1D *  h_xbj_acc_copy = (TH1D*)hist_xbj_mom_acc->Clone();
      h_xbj_acc_copy->Scale(0.25); 
      TH1D *  h_xbj_coin_copy = (TH1D*)hist_xbj_mom->Clone();
      TH1D *  h_xbj_acc_only = (TH1D*)h_xbj_acc_copy->Clone();
      h_xbj_clean_mom_counts->Add(h_xbj_coin_copy,h_xbj_acc_only,1.0,-1.0);
      // acc sub for z


      TH1D *  h_zhad_acc_copy = (TH1D*)hist_zhad_mom_acc->Clone();
      h_zhad_acc_copy->Scale(0.25); 
      TH1D *  h_zhad_coin_copy = (TH1D*)hist_zhad_mom->Clone();
      TH1D *  h_zhad_acc_only = (TH1D*)h_zhad_acc_copy->Clone();
      h_zhad_clean_mom_counts->Add(h_zhad_coin_copy,h_zhad_acc_only,1.0,-1.0);


      TH1D *  h_pt_acc_copy = (TH1D*)hist_pt_mom_acc->Clone();
      h_pt_acc_copy->Scale(0.25); 
      TH1D *  h_pt_coin_copy = (TH1D*)hist_pt_mom->Clone();
      TH1D *  h_pt_acc_only = (TH1D*)h_pt_acc_copy->Clone();
      h_pt_clean_mom_counts->Add(h_pt_coin_copy,h_pt_acc_only,1.0,-1.0);


      TH1D *  h_phi_acc_copy = (TH1D*)hist_phi_mom_acc->Clone();
      h_phi_acc_copy->Scale(0.25); 
      TH1D *  h_phi_coin_copy = (TH1D*)hist_phi_mom->Clone();
      TH1D *  h_phi_acc_only = (TH1D*)h_phi_acc_copy->Clone();
      h_phi_clean_mom_counts->Add(h_phi_coin_copy,h_phi_acc_only,1.0,-1.0);


      //acc sub for z prime
      TH1D *  h_zhadprime_acc_copy = (TH1D*)hist_zhadprime_acc->Clone();
      h_zhadprime_acc_copy->Scale(0.25); 
      TH1D *  h_zhadprime_coin_copy = (TH1D*)hist_zhadprime->Clone();
      TH1D *  h_zhadprime_acc_only = (TH1D*)h_zhadprime_acc_copy->Clone();
      h_zhadprime_clean_mom_counts->Add(h_zhadprime_coin_copy,h_zhadprime_acc_only,1.0,-1.0);




      //acc sub done
      double hgcer_eff;
      if(momentum>=3.3){hgcer_eff=0.986;}
      if(momentum>=3 && momentum<3.3){hgcer_eff=0.982;}
      if(momentum<3 && momentum>= 2.9){hgcer_eff=0.98;}
      if(momentum<2.9){hgcer_eff=0.995;}

      TH1D *  h_xbj_clean_yield = (TH1D*)h_xbj_clean_mom_counts->Clone();
      h_xbj_clean_yield->Scale(1/(charge_times_eff*hgcer_eff));
      TH1D *  h_zhad_clean_yield = (TH1D*)h_zhad_clean_mom_counts->Clone();
      h_zhad_clean_yield->Scale(1/(charge_times_eff*hgcer_eff));

      //plot run by run yield and save as ROOTFILE Run By Run
      // TFile *f1 = new TFile(Form("yield_root_20bins/yield_hist_%d.root",RunNumber),"RECREATE");
      std::string rootfile_yield_name = Form("yield_root_20bins_hgcer_pt/yield_histo_geom_cuts_%d.root",RunNumber);
      TFile *rootfile_yield = new TFile(rootfile_yield_name.c_str(),"RECREATE");
      h_xbj_clean_yield->Write();
      h_zhad_clean_yield->Write();
      // h_p_clean_yield->Write();

      rootfile_yield->Close();

      std::string rootfile_count_name = Form("yield_root_20bins_hgcer_mx22p4_pt/count_histo_geom_cuts_%d.root",RunNumber);
      TFile *rootfile_count = new TFile(rootfile_count_name.c_str(),"RECREATE");
      if(momentum>=3.3)
	{ 
	  h_xbj_clean_mom_counts->Scale(1/0.986);
	  h_zhad_clean_mom_counts->Scale(1/0.986);
	  h_pt_clean_mom_counts->Scale(1/0.986);
	  h_zhadprime_clean_mom_counts->Scale(1/0.986);
	  h_phi_clean_mom_counts->Scale(1/0.986);

	}
      if(momentum>=3.0 && momentum < 3.3)
	{ 
	  h_xbj_clean_mom_counts->Scale(1/0.982);
	  h_zhad_clean_mom_counts->Scale(1/0.982);
	  h_pt_clean_mom_counts->Scale(1/0.982);
	  h_phi_clean_mom_counts->Scale(1/0.982);
	  h_zhadprime_clean_mom_counts->Scale(1/0.982);
	}
      if(momentum<3.0 && momentum >= 2.9)
	{ 
	  h_xbj_clean_mom_counts->Scale(1/0.98);
	  h_zhad_clean_mom_counts->Scale(1/0.98);
	  h_pt_clean_mom_counts->Scale(1/0.98);
	  h_phi_clean_mom_counts->Scale(1/0.98);
	  h_zhadprime_clean_mom_counts->Scale(1/0.98);
	}
      if(momentum<3.0)
	{ 
	  h_xbj_clean_mom_counts->Scale(1/0.995);
	  h_zhad_clean_mom_counts->Scale(1/0.995);
	  h_pt_clean_mom_counts->Scale(1/0.995);
	  h_phi_clean_mom_counts->Scale(1/0.995);
	  h_zhadprime_clean_mom_counts->Scale(1/0.995);
	}


      h_xbj_clean_mom_counts->Write();
      h_zhad_clean_mom_counts->Write();
      h_pt_clean_mom_counts->Write();
      h_phi_clean_mom_counts->Write();
      h_zhadprime_clean_mom_counts->Write();
 
      rootfile_count->Close();

      TCanvas * d = new TCanvas ("d", "d", 1200, 800);
      d->Divide(2,3);

      d->cd(1);
      h_xbj_coin_copy->Draw("e,p1");
      h_xbj_clean_mom_counts->Draw("e,p1 same");
      h_xbj_acc_only->Draw("e,p1 same");
      h_xbj_clean_mom_counts->SetLineColor(kBlack);
      h_xbj_acc_only->SetLineColor(kRed);
      h_xbj_coin_copy->SetTitle(Form("%d Xbj Counts",RunNumber));
      
      d->cd(2);
      h_zhad_coin_copy->Draw("e,p1");
      h_zhad_clean_mom_counts->Draw("e,p1 same");
      h_zhad_acc_only->Draw("e,p1 same");
      h_zhad_clean_mom_counts->SetLineColor(kBlack);
      h_zhad_acc_only->SetLineColor(kRed);
      h_zhad_coin_copy->SetTitle(Form("%d Zhad Counts",RunNumber));


      d->cd(3);
      h_xbj_clean_yield->Draw("e, p1");
      h_xbj_clean_yield->SetTitle(Form("%d Xbj Yield",RunNumber));

      d->cd(4);
      h_zhad_clean_yield->Draw("e, p1");
      h_zhad_clean_yield->SetTitle(Form("%d Zhad Yield",RunNumber));

      d->cd(5);
      h_pt_coin_copy->Draw("e,p1");
      h_pt_clean_mom_counts->Draw("e,p1 same");
      h_pt_acc_only->Draw("e,p1 same");
      h_pt_clean_mom_counts->SetLineColor(kBlack);
      h_pt_acc_only->SetLineColor(kRed);
      h_pt_coin_copy->SetTitle(Form("%d Pt Counts",RunNumber));

      d->cd(6);
      h_phi_coin_copy->Draw("e,p1");
      h_phi_clean_mom_counts->Draw("e,p1 same");
      h_phi_acc_only->Draw("e,p1 same");
      h_phi_clean_mom_counts->SetLineColor(kBlack);
      h_phi_acc_only->SetLineColor(kRed);
      h_phi_coin_copy->SetTitle(Form("%d Phi Counts",RunNumber));


      d->SaveAs(Form("pdfjan13/counts_yield_hist_%d_geom_cuts.pdf", RunNumber));




      cout<<"  "<<endl;
      cout<<"  "<<endl;
      cout<<"  "<<endl;
      cout<<RunNumber<<" "<<momentum<<" "<<"HGCER EFF ="<< hgcer_eff<<endl;

      //yield histograms

      Double_t coin_integral_hist = h_pcal_coin_copy->Integral(); 
      Double_t acc_integral_4_hist = h_pcal_acc_copy->Integral(); 
      Double_t acc_integral_hist = h_pcal_acc_only->Integral();
      Double_t clean_integral_hist =  h_pcal_clean->Integral();


      //Double_t clean_integral_hist =  h_xbj_clean_yield->Integral();

      
      // cout<<runs[r]<<"\t\t"<<"clean_integral =====>   "<<clean_integral_hist<<endl;
      // cout<<runs[r]<<"\t\t"<<"coin_integral  =====>   "<<coin_integral_hist<<endl;
      // cout<<runs[r]<<"\t\t"<<"acc_integral   =====>   "<<acc_integral_hist<<endl;
      // verified that these yields are exactly same  by actual subtraction of coin and acc hists

      Double_t charge_cut5uA_bcmcorr  = charge_sum_cut_corr/1000.0;//mC
      Double_t yield;
      Double_t yield_err;
      //getting correct yield for each run:
      if(momentum >= 3 && momentum < 3.3 ){
	yield = clean_integral/(charge_cut5uA_bcmcorr* efficiency * 0.982);
	yield_err = sqrt(clean_integral)/(charge_cut5uA_bcmcorr*efficiency *0.982 );
      }
      if(momentum >=3.3){
	yield = clean_integral/(charge_cut5uA_bcmcorr*efficiency * 0.986);
	yield_err = sqrt(clean_integral)/(charge_cut5uA_bcmcorr*efficiency * 0.986);
      }
      if(momentum < 3.0 && momentum >=2.9 ){
	yield = clean_integral/(charge_cut5uA_bcmcorr* efficiency * 0.98);
	yield_err = sqrt(clean_integral)/(charge_cut5uA_bcmcorr*efficiency *0.98);
      }
      if(momentum<2.9){
	yield = clean_integral/(charge_cut5uA_bcmcorr*efficiency * 0.995 );
	yield_err = sqrt(clean_integral)/(charge_cut5uA_bcmcorr*efficiency *0.995);
      }
      ///cout<<runs[r]<<setprecision(4)<<fixed<<"\t\t"<<"curr, charge, yield, yield_err ="<<peak_current<<"\t\t"<<charge_cut5uA_bcmcorr<<"\t\t"<<yield<<"\t\t"<<yield_err<<endl;
      cout<<"  "<<endl;
      cout<<"  "<<endl;
      cout<<"Final Yield = "<<RunNumber<<"  "<<yield<<"      "<<h_xbj_clean_yield->Integral()<<endl;
      cout<<"  "<<endl;
      cout<<"  "<<endl;
      txtfile<<RunNumber<<setprecision(6)<<fixed<<"\t\t"<<efficiency<<"\t\t"<<peak_current<<"\t\t"<<charge_cut5uA_bcmcorr<<"\t\t"<<yield<<"\t\t"<<yield_err<<endl;

      TCanvas * cc = new TCanvas ("cc", "cc", 1200, 800);
      cc->Divide(3,3);

      cc->cd(1);
      h_cur_entry->Draw();
      h_cur_entry->SetTitle(Form("current %d",RunNumber));
      h_cur_entry->SetMarkerStyle(20);
      h_cur_entry->SetLineWidth(2);


      cc->cd(2);gPad->SetLogy();
      h_cur->Draw();
      h_cur->SetLineWidth(2);

      cc->cd(3);
      h_pcal_coin->Draw();
      h_pcal_coin->SetLineWidth(2);
      h_pcal_acc->Draw("sames");
      h_pcal_acc->SetLineWidth(2);
      h_pcal_acc->SetLineColor(kRed);
      h_pcal_coin->SetTitle(Form("shms cal %d",RunNumber));

   
      cc->cd(4);
      h_ctime->Draw();
      h_ctime->SetLineWidth(2);
      h_ctime_coin->Draw("sames");
      h_ctime_coin->SetLineWidth(2);
      h_ctime_coin->SetLineColor(kGreen);
      h_ctime_acc->Draw("sames");
      h_ctime_acc->SetLineWidth(2);
      h_ctime_acc->SetLineColor(kRed);

      cc->cd(5);
      h_ctime_clean->Draw("hist");
      h_ctime_clean->SetLineWidth(2);


      cc->cd(6);
      h_pcal_clean->Draw("hist");
      h_pcal_clean->SetLineWidth(2);

      cc->cd(7);
      h_mom->Draw("hist");
      h_mom->SetLineWidth(2);

      cc->cd(8);gPad->SetLogz();
      h_aero_y_vs_x->Draw("colz");

      cc->cd(9);gPad->SetLogz();
      h_showerback_y_vs_x->Draw("colz");
      cc->SaveAs(Form("pdfjan/base_json_%d_geom_cuts.pdf", RunNumber));

      TCanvas * c = new TCanvas ("c", "c", 1200, 800);
      c->Divide(3,2);

      c->cd(1);
      hist_zhad->Draw("hist");
      hist_zhad->SetLineWidth(2);
      hist_zhadprime->Draw("hist same");
      hist_zhadprime->SetLineWidth(2);
      hist_zhadprime->SetLineColor(kRed);

      c->cd(2);
      hist_xbj->Draw("hist");
      hist_xbj->SetLineWidth(2);
      hist_xbjprime->Draw("hist same");
      hist_xbjprime->SetLineWidth(2);
      hist_xbjprime->SetLineColor(kRed);
      c->cd(3);
      hist_pt->Draw("hist");
      hist_pt->SetLineWidth(2);

      c->cd(4);
      hist_phi->Draw("hist");
      hist_phi->SetLineWidth(2);


      c->SaveAs(Form("pdfjan/x_z_pt_%d_geom_cuts.pdf", RunNumber));


      //you should fill the output json here to get filled after each run
      j_out1[it.key()]["yield"]=yield;
      j_out1[it.key()]["yield_err"]=yield_err;
      j_info[it.key()]["yield"]=yield;
      j_info[it.key()]["yield_err"]=yield_err;
      // std::ofstream ofs("/home/hdbhatt/yield/json_files/yield2022_effgt5_nohgc1.json");

      std::ofstream ofs("/home/hdbhatt/yield/json_files/yield_june27.json");
      ofs<<j_info.dump(4)<<std::endl;
    }//if charge value is not empty
    
    // std::ofstream ofs("json_files/yield_each_run_extra.json");
    //ofs<<j_info.dump(4)<<std::endl;

    //you should fill the output json here to get filled after whole script is processed
 
  } //j_info.begin...
      
 
 
}


