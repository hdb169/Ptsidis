
#define NRUNS 1
void pt_vs_phi_csvruns(){
  // Int_t runs[NRUNS]     ={5600};//{5799};
  Int_t runs[NRUNS]     ={1};//{5799};
 TFile *f;
  TTree *tt;
  for(Int_t i=0; i<NRUNS; i++){
    
     f=new TFile(Form("/lustre19/expphy/volatile/hallc/c-csv/hdbhatt/ROOTfiles/skimROOTfilesfeb11/ptsidis_skim_all_%d.root", runs[i]),"read");


 if ((!f) || (f->IsZombie())){
      cout << "No file found or Zombie file" << endl;
      continue;
    }

    tt=(TTree *)f->Get("T");
    Long64_t nentriesD = tt->GetEntries();
    TProfile *hprof = new TProfile("hprof","Profile of pt versus phi",100,0,7,0,1);
    // double Epi,zhad;


   

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
      Double_t dataphi;
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
	  pt=mom*sin(theta);


	Epi = sqrt(pow(pionmass,2) + pow(mom,2));
	zhad = Epi/dnu;
	  Mx2 = (protonmass + dnu - zhad*dnu) *  (protonmass + dnu - zhad*dnu) - abs(dpmiss)*abs(dpmiss);




	  if(pdelta>-10&&pdelta<20 && hdelta>-8 && hdelta<8&&Mx2>2.6&&dW2>4&&ctime>42&&ctime<46&&aero>4&&pcaletottrack<0.6&&hcaletottrack>0.8){
	  //hprof->Fill(zhad,xbj,1);
	  hprof->Fill(dataphi,pt,1);


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
       hprof->GetXaxis()->SetTitle("phi");
       hprof->GetYaxis()->SetTitle("pt");
       hprof->SetTitle("pt vs phi");
       c1->SaveAs("pt_vs_phi_ptsidis.pdf");
  }
}
