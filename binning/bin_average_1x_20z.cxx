//For 1x, 20 z bins, bin avearging example:
//this is an example of 1 x setting has 20 z bins, and grid of 20*20 in (x,z) for one setting.
#define NRUNS 1
void bin_average_1x_20z(){
  Int_t runs[NRUNS]     ={5839};
  TFile *f;
  TTree *tt;
  for(Int_t i=0; i<NRUNS; i++){
    f=new TFile(Form("/lustre19/expphy/volatile/hallc/c-csv/hdbhatt/ROOTfiles/ROOTfiles_ptsidis/run_%d_-1.root", runs[i]),"read");
    if ((!f) || (f->IsZombie())){
      cout << "No file found or Zombie file" << endl;
      continue;
    }

    tt=(TTree *)f->Get("T");
    Long64_t nentriesD = tt->GetEntries();

    double x_sum[20]={0.};
    double x_av[20]={0.};
    Int_t x_num[20]={0};
    double yield_x[20]={0.};
    double z_sum[20][20]={0.};
    double z_av[20][20]={0.};
    Int_t z_num[20][20]={0};
    double yield_z[20][20]={0.};
    Double_t pionMass = 0.1395701835;   

    Double_t Epi,zhad ;
    TH2D *h_xz     = new TH2D("h_xz","", 20,0,1,20,0,1);                                      
    TH1D *h_x     = new TH1D("h_x","", 20,0,1);                                      
    TH1D *h_z[20];
    for(int k=0;k<20;k++){
      h_z[k]     = new TH1D(Form("h_z[%d]",k),"", 20,0,1);                                      
    }

    Double_t dnu;
    tt->SetBranchAddress("H.kin.primary.omega",&dnu);
    Double_t hdelta;
    tt->SetBranchAddress("H.gtr.dp", &hdelta);                                                                            
    Double_t pdelta;
    tt->SetBranchAddress("P.gtr.dp", &pdelta);                                                                            
    Double_t shmsP;
    tt->SetBranchAddress("P.gtr.p", &shmsP); 
    Double_t hmsP;
    tt->SetBranchAddress("H.gtr.p", &hmsP); 
    Double_t xbj;
    tt->SetBranchAddress("H.kin.primary.x_bj",&xbj);


    Double_t pionmass   = 0.1395701835; 
    Double_t protonmass = 0.93827231;//GeV/c2
    Double_t pi = 3.1416;
    Double_t charge = 10;//say

    for (int kk=0; kk<nentriesD;  kk++){
      

	tt->GetEntry(kk);
	Epi = sqrt(pow(pionmass,2) + pow(shmsP,2));
	zhad = Epi/dnu;
	if(pdelta>-10&&pdelta<20 && hdelta>-8 && hdelta<8){

	  for(int k=0;k<20;k++){//20 bins in x
	    if(xbj>0+0.05*k  && xbj< 0.05*0.05*k)
	      {
	 	x_sum[k]+=xbj;
		x_num[k]++;
	      }

	    for(int j=0;j<20;j++){
	      if(zhad>=0+0.05*j  &&  zhad<0.05+0.05*j){
		z_sum[k][j]+=zhad;
		z_num[k][j]++;
	      }//z
	    }
	      
	    h_xz->Fill(xbj,zhad);
	    h_x->Fill(xbj);
	    for(int j=0;j<20;j++){
	      h_z[j]->Fill(zhad);
	    }

	  }
	}
      }//event loop done
  
 
    for(int k=0;k<20;k++){
      x_av[k]=x_sum[k]/x_num[k];
      yield_x[k]=x_num[k]/charge;
      if(x_num[k]==0)yield_x[k]=0;

      for(int j=0;j<20;j++){
	z_av[k][j]=z_sum[k][j]/z_num[k][j];
	yield_z[k][j]=z_num[k][j]/charge;
	if(z_num[k][j]==0)yield_z[k][j]=0;

	cout<<h_x->GetBinCenter(k+1)<<"\t\t"<<x_av[k]<<"\t\t"<<yield_x[k]<<"\t\t"<<" z had "<< h_z[k]->GetBinCenter(j+1)<<"\t\t"<<z_av[k][j]<<"\t\t"<<yield_z[k][j]<<endl;
      }
    }


     TCanvas *c1 = new TCanvas("c1", "Tracking Eff",1200,800); 
     c1->Divide(1,1);
       c1->cd(1);
       h_xz->Draw("colz");
       h_xz->SetLineWidth(2);
       h_xz->SetLineColor(kRed);
       h_xz->GetXaxis()->CenterTitle();
       h_xz->GetXaxis()->SetTitleOffset(1.2);
       h_xz->GetXaxis()->SetTitle("x");
       h_xz->GetYaxis()->SetTitle("z");
       h_xz->SetTitle("z vs x");

       /* c1->cd(2);
       h_x->Draw("colz");
       h_x->SetLineWidth(2);
       h_x->SetLineColor(kRed);
       h_x->GetXaxis()->CenterTitle();
       h_x->GetXaxis()->SetTitleOffset(1.2);
       h_x->GetXaxis()->SetTitle("x");
       h_x->GetYaxis()->SetTitle("z");
       h_x->SetTitle("xBj");

       c1->cd();
       h_z->Draw("colz");
       h_z->SetLineWidth(2);
       h_z->SetLineColor(kRed);
       h_z->GetXaxis()->CenterTitle();
       h_z->GetXaxis()->SetTitleOffset(1.2);
       h_z->GetXaxis()->SetTitle("x");
       h_z->GetYaxis()->SetTitle("z");
       h_z->SetTitle("zhad");
       */
  }
}

