To check at local computer: Full working code is ptsidis_skim_coutRootfile_geom_cuts.cxx

#define NRUNS 1//2
void cuts_study(){

  //  ofstream txtfile("textfiles/feynman_x.txt");
  ofstream txt("feynmanx_csv.txt",ios::app);
  //gROOT->SetBatch(kTRUE);
  for(Int_t r = 0; r<NRUNS;r++){
    Int_t runs[NRUNS]     = {6094};//{6542};//{6542,6543}
    //TString filename  = Form("ROOTfiles/run_%d.root", runs[r]);
    // TString filename  = Form("ROOTfiles/added_6543_and_%d.root", runs[r]);
    TString filename  = Form("ROOTfiles/run_%d.root", runs[r]);

    TFile *f = new TFile(filename, "READ");
    if(f->IsZombie()){
      std::cout << "Error opening file "<< std::endl;
      continue;
    }
    //hists
    TH1F* hfyn_x = new TH1F("hfyn_x","Feynman X",100,0,10);
    TH1F* hx = new TH1F("hx","XBj",100,0,1);
    TH1F* hmsmom = new TH1F("hmsmom","hms P",100,0,7);
    TH1F* pmom = new TH1F("pmom","shms P",100,0,7);
    TH1F* hW2  = new TH1F("hW2","W2",100,0,12);
    TH1F* hQ2  = new TH1F("hQ2","Q2",100,0,12);
    TH1F* hpt = new TH1F("hpt","Pt",100,0, 1);
    TH1F* hpt2 = new TH1F("hpt2","Pt2",100,0, 1);
    TH1F* hz = new TH1F("hz","Zhad",100,0, 1);

    TH2F* h_aero_y_vs_x = new TH2F("h_aero_y_vs_x","Aero_y_vs_x(After Cuts);yAtAero;xAtAero",100,-100,100,100,-100,100);
    TH2F* h_showerback_y_vs_x = new TH2F("h_showerback_y_vs_x","Showerback_y_vs_x;yAtCal;xAtCal",100,-70,70,100,-70,70);
    TH2F* h_aero_y_vs_x_bc = new TH2F("h_aero_y_vs_x_bc","Aero_y_vs_x (before cuts);yAtAero;xAtAero",100,-100,100,100,-100,100);
    TH2F* h_showerback_y_vs_x_bc = new TH2F("h_showerback_y_vs_x_bc","Showerback_y_vs_x (before cuts);yAtCal;xAtCal",100,-70,70,100,-70,70);
    TH2F* h_hgcer_y_vs_x = new TH2F("h_hgcer_y_vs_x","Hgcer_y_vs_x (After Cuts);yAtHgcer;xAtHgcer",100,-100,100,100,-100,100);
    TH2F* h_hgcer_y_vs_x_bc = new TH2F("h_hgcer_y_vs_x_bc","Hgcer_y_vs_x (Before Cuts);yAtHgcer;xAtHgcer",100,-100,100,100,-100,100);

    TH1F* hmsdelta = new TH1F("hmsdelta","hms Delta (after cuts);hms Delra; Counts",100,-10,10);
    TH1F* shmsdelta = new TH1F("shmsdelta","shms Delta (after cuts);shms Delra; Counts",100,-12,22);
    TH1F* hmsdelta_bc = new TH1F("hmsdelta_bc","hms Delta (before cuts);hms Delta; Counts",100,-10,10);//bc = before cuts
    TH1F* shmsdelta_bc = new TH1F("shmsdelta_bc","shms Delta (before cuts);shms Delta; Counts",100,-12,22);

    ////hists

   
    Double_t pionmass   = 0.1395701835; 
    Double_t protonmass = 0.93827231;//GeV/c2
    Double_t pi = 3.1416;
    TTree * tt = (TTree *)f->Get("T");
    Long64_t data_entries  = tt->GetEntries();
   
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

    
    
    double pt2,inv_mass,fyn_x,Epi,zhad;
    Bool_t accL_cut, accR_cut, accidental_cut, ctime_cut, coin_cut,acc_cut, central_hole, right_or_left_of_hole, remove_region, peter_all_geom_cut, common_cuts,  pcal_geom_cut_hem,  paero_geom_cut,hcal_geom_cut,hms_dc_fp_cut;
    Bool_t  hms_hodo_cut, pcal_geom_cut, shms_dc_fp_cut, hourglass_cut, geom_cut_and_hgcer_cut;
    Double_t xathgcer,yathgcer,mhp_square, xbjprime, zhadprime, pt, xatpreshower, yatpreshower, xatshower, yatshower,xatpreshowerback , yatpreshowerback, xatshowerback, yatshowerback,   Mx2,xataero,yataero,xathcal,yathcal;
    Double_t xathhodo,xatpcal,yatpcal;

    //cout<<runs[r]<<endl;
    Double_t gevtyp;
    tt->SetBranchAddress("g_evtyp",&gevtyp);
    Double_t gevnum;
    tt->SetBranchAddress("g_evnum",&gevnum);
    Double_t mom;
    tt->SetBranchAddress("P_gtr_p",&mom);
    Double_t hmom;
    tt->SetBranchAddress("H_gtr_p",&hmom); 
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

    Double_t htrack_x_fp;
    tt->SetBranchAddress("H_dc_x_fp", &htrack_x_fp);
    Double_t htrack_y_fp;
    tt->SetBranchAddress("H_dc_y_fp", &htrack_y_fp);
    Double_t htrack_xp_fp;
    tt->SetBranchAddress("H_dc_xp_fp", &htrack_xp_fp);
    Double_t htrack_yp_fp;
    tt->SetBranchAddress("H_dc_yp_fp", &htrack_yp_fp);
    
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


    for (int i = 0; i < data_entries; i++) {
      tt->GetEntry(i);

     
      Epi = sqrt(pow(pionmass,2) + pow(mom,2));
      zhad = Epi/dnu;
      pt=mom*sin(theta);
      pt2=pt*pt;
      inv_mass=sqrt(dW2);
      fyn_x=2*sqrt(fabs(mom*mom-pt2))/inv_mass;
      //cuts study
      
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
      
      
      Epi = sqrt(pow(pionmass,2) + pow(mom,2));
      zhad = Epi/dnu;
      Mx2 = (protonmass + dnu - zhad*dnu) *  (protonmass + dnu - zhad*dnu) - abs(dpmiss)*abs(dpmiss);
      pt=mom*sin(theta);
      xbjprime = 2*xbj/(1+sqrt(1+4*xbj*protonmass*protonmass/Q2));//target mass corrected x
      mhp_square = pionmass*pionmass+pt*pt;
      zhadprime = (zhad/2)*(xbjprime/xbj)*(1+sqrt(1-(4*xbj*xbj*protonmass*protonmass*mhp_square)/(zhad*zhad*Q2*Q2)));//target mass corrected zhad
	 
      if(dataphi<0)dataphi=2*pi+dataphi;
      if(dataphi>=0)dataphi=dataphi;

      if(pdelta>-10&&pdelta<20&&hdelta>-8&&hdelta<8&&gevtyp>3){


	hmsdelta_bc->Fill(hdelta);
	shmsdelta_bc->Fill(pdelta);

	
	h_aero_y_vs_x_bc->Fill(yataero,xataero);
	h_showerback_y_vs_x_bc->Fill(yatpreshowerback,xatpreshowerback);
	h_hgcer_y_vs_x_bc->Fill(yathgcer,xathgcer);

	if(geom_cut_and_hgcer_cut){
	  h_aero_y_vs_x->Fill(yataero,xataero);
	  h_showerback_y_vs_x->Fill(yatpreshowerback,xatpreshowerback);
	  hmsdelta->Fill(hdelta);
	  shmsdelta->Fill(pdelta);
	  h_hgcer_y_vs_x->Fill(yathgcer,xathgcer);
	}
	
	hQ2->Fill(Q2);
      }//pid

    }//entries



    TCanvas *c1 = new TCanvas("c1","c1", 1200,800);
    c1->Divide(4,3);
    c1->cd(1);gPad->SetLogz();
    h_aero_y_vs_x_bc->Draw("colz");

    c1->cd(2);gPad->SetLogz();
    h_aero_y_vs_x->Draw("colz");


    

    c1->cd(3);gPad->SetLogz();
    h_showerback_y_vs_x_bc->Draw("colz");

    
    c1->cd(4);gPad->SetLogz();
    h_showerback_y_vs_x->Draw("colz");

    
 
    c1->cd(5);gPad->SetLogz();
    h_hgcer_y_vs_x_bc->Draw("colz");

  
    c1->cd(6);gPad->SetLogz();

    h_hgcer_y_vs_x->Draw("colz");

    c1->cd(7);
    hmsdelta_bc->Draw();
    hmsdelta_bc->SetLineColor(kBlue);

    c1->cd(8);
    hmsdelta->Draw("");
    hmsdelta->SetLineColor(kRed);


    c1->cd(9);
    shmsdelta_bc->Draw();
    shmsdelta_bc->SetLineColor(kBlue);

    c1->cd(10);

    shmsdelta->Draw("");
    shmsdelta->SetLineColor(kRed);


    
  
    
    c1->SaveAs(Form("pdf/cuts_study_%d.pdf", runs[r]));

  }//runs
}//void
