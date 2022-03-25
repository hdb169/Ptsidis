 #define ROOT_TMath
    
#include "TMathBase.h"
    
#include "TError.h"
#include <algorithm>
#include <limits>
#include <cmath>
#define NRUNS 1
void pt_vs_phi_csvruns(){
  // Int_t runs[NRUNS]     ={5696};//{5799};
  Int_t runs[NRUNS]     ={1};//{5799};
 TFile *f;
  TTree *tt;
  for(Int_t i=0; i<NRUNS; i++){
    // f=new TFile(Form("/lustre19/expphy/volatile/hallc/c-csv/hdbhatt/ROOTfiles/ROOTfiles_ptsidis/pt_vs_phi_run5700_%d.root", runs[i]),"read");
 

    f=new TFile(Form("/lustre19/expphy/volatile/hallc/c-csv/hdbhatt/ROOTfiles/skimROOTfilesfeb11/all_csv_runs_pip_pim_%d.root", runs[i]),"read");

    // f=new TFile(Form("/lustre19/expphy/volatile/hallc/c-csv/hdbhatt/ROOTfiles/skimROOTfilesfeb11/ptsidis_skim_all_%d.root", runs[i]),"read");
    //  f=new TFile(Form("/lustre19/expphy/volatile/hallc/c-csv/hdbhatt/ROOTfiles/skimROOTfilesfeb11/run_%d.root", runs[i]),"read");

    //f=new TFile(Form("/lustre19/expphy/volatile/hallc/c-csv/hdbhatt/ROOTfiles/ROOTfiles_ptsidis/ptsidis_giant_%d.root", runs[i]),"read");
  
 

     //  f=new TFile(Form("/lustre19/expphy/volatile/hallc/c-csv/hdbhatt/ROOTfiles/ROOTfiles_ptsidis/run_%d_-1.root", runs[i]),"read");
     gStyle->SetPalette(1,0);
     // gStyle->SetOptStat(1000011);
     gStyle->SetOptStat(0);
     gStyle->SetOptFit(11);
     // gStyle->SetTitleOffset(1.,"Y");
     // gStyle->SetTitleOffset(.5,"X");
     gStyle->SetLabelSize(0.04,"XY");
     gStyle->SetTitleSize(0.04,"XY");
     gStyle->SetPadLeftMargin(0.12);

     if ((!f) || (f->IsZombie())){
       cout << "No file found or Zombie file" << endl;
       continue;
    }

     Double_t rad_to_deg = 180./3.14159265358979323846;
    tt=(TTree *)f->Get("T");
    Long64_t nentriesD = tt->GetEntries();
    TProfile *hprof = new TProfile("hprof","Profile of pt versus phi",100,-10,370,0,1);
    //   TProfile *hprof = new TProfile("hprof","Profile of pt versus phi",100,-1,7,0,1);
//TH2D *ptvsphi = new TH2D("ptvsphi"," pt versus phi",100,-1,7,100,0,1);
    TH2D *ptvsphi = new TH2D("ptvsphi"," pt versus phi",100,-10,370,100,0,1);
// double Epi,zhad;


    // Double_t dnu;
    // tt->SetBranchAddress("H.kin.primary.omega",&dnu);
    // Double_t hdelta;
    // tt->SetBranchAddress("H.gtr.dp", &hdelta);                                                                            
    // Double_t pdelta;
    // tt->SetBranchAddress("P.gtr.dp", &pdelta);                                                                            
    // Double_t shmsP;
    // tt->SetBranchAddress("P.gtr.p", &shmsP); 
    // Double_t hmsP;
    // tt->SetBranchAddress("H.gtr.p", &hmsP); 
    // Double_t xbj;
    // tt->SetBranchAddress("H.kin.primary.x_bj",&xbj);

 Double_t xathgcer,yathgcer,mhp_square, xbjprime, zhadprime, pt, xatpreshower, yatpreshower, xatshower, yatshower,xatpreshowerback , yatpreshowerback, xatshowerback, yatshowerback, Epi, zhad, Mx2,xataero,yataero;;
 
Double_t aero;
      tt->SetBranchAddress("P_aero_npeSum",&aero);
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
      Double_t hdelta;
      tt->SetBranchAddress("H_gtr_dp",&hdelta); 
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
      // double_t hgcer_xAtcer;
      // tt->SetBranchAddress("P_hgcer_xAtCer", &hgcer_xAtcer);
      // double_t hgcer_yAtcer;
      // tt->SetBranchAddress("P_hgcer_yAtCer", &hgcer_yAtcer);
      // Double_t track_x_fp;
      // tt->SetBranchAddress("P_dc_x_fp", &track_x_fp);
      // Double_t track_y_fp;
      // tt->SetBranchAddress("P_dc_y_fp", &track_y_fp);
      // Double_t track_xp_fp;
      // tt->SetBranchAddress("P_dc_xp_fp", &track_xp_fp);
      // Double_t track_yp_fp;
      // tt->SetBranchAddress("P_dc_yp_fp", &track_yp_fp);
      Double_t dW2;
      tt->SetBranchAddress("H_kin_primary_W2",&dW2);
      Double_t dnu;
      tt->SetBranchAddress("H_kin_primary_omega",&dnu);
      Double_t dpmiss;
      tt->SetBranchAddress("P_kin_secondary_pmiss",&dpmiss);
      Double_t theta;
      tt->SetBranchAddress("P_kin_secondary_th_xq", &theta);
      Double_t dataphi,phi;
      tt->SetBranchAddress("P_kin_secondary_ph_xq", &dataphi);
      Double_t xbj;  
      tt->SetBranchAddress("H_kin_primary_x_bj",&xbj);
      Double_t Q2;
      tt->SetBranchAddress("H_kin_primary_Q2",&Q2);




    Double_t pionmass   = 0.1395701835; 
    Double_t protonmass = 0.93827231;//GeV/c2
    Double_t pi = 3.1416;
    Double_t charge = 10;//say

    for (int kk=0; kk<nentriesD;  kk++){
	tt->GetEntry(kk);

       if(dataphi<0)dataphi=2*pi+dataphi;
	  if(dataphi>=0)dataphi=dataphi;
	  double phi=dataphi*rad_to_deg;
	  pt=mom*sin(theta);


	Epi = sqrt(pow(pionmass,2) + pow(mom,2));
	zhad = Epi/dnu;
	  Mx2 = (protonmass + dnu - zhad*dnu) *  (protonmass + dnu - zhad*dnu) - abs(dpmiss)*abs(dpmiss);




	  if(pdelta>-10&&pdelta<20 && hdelta>-8 && hdelta<8&&Mx2>2.6&&dW2>4&&ctime>42&&ctime<46&&aero>4&&pcaletottrack<0.6&&hcaletottrack>0.8){
	  //hprof->Fill(zhad,xbj,1);
	  hprof->Fill(phi,pt,1);
	  ptvsphi->Fill(phi,pt);


      }
    }//event loop done
  
   

    TCanvas *c1 = new TCanvas("c1", "Tracking Eff",1200,800); 
    c1->Divide(1,1);
    c1->cd(1);
    hprof->Draw("");
    hprof->SetLineWidth(2);
    hprof->SetLineColor(kRed);
    hprof->GetXaxis()->CenterTitle();
    hprof->GetXaxis()->SetTitleOffset(1.2);
    hprof->GetXaxis()->SetTitle("Phi (Degree)");
    hprof->GetYaxis()->SetTitle("Pt (GeV)");
    hprof->SetTitle("Pt vs Phi");
    hprof->GetXaxis()->CenterTitle();
    hprof->GetYaxis()->CenterTitle();
    
    c1->cd(1);gPad->SetGrid();
    ptvsphi->Draw("SCAT");
    ptvsphi->SetLineWidth(2);
    ptvsphi->SetLineColor(kRed);
    ptvsphi->SetMarkerColor(kRed);
    //ptvsphi->SetMarkerStyle(8);
    ptvsphi->GetXaxis()->CenterTitle();
    ptvsphi->GetXaxis()->SetTitleOffset(1.2);
    ptvsphi->GetXaxis()->SetTitle("Phi (Degree)");
    ptvsphi->GetYaxis()->SetTitle("Pt (GeV)");
    ptvsphi->SetTitle("Pt vs Phi");
    ptvsphi->GetXaxis()->CenterTitle();
    ptvsphi->GetYaxis()->CenterTitle();

    // c1->SaveAs("pdf/pt_vs_phi_ptsidis.pdf");   
     c1->SaveAs("pdf/pt_vs_phi_csv.pdf");

  }
}
