void weighted_av(){
  gStyle->SetOptFit(1);
  //  const double kin[] = {3.98,4.75,4.775,5.5,4.0};
  //   const double xbj[] = {0.5,0.6,0.55,0.6,0.35};
  const double kin[] = {4.0,3.98,4.75,4.764,5.5};
  const double xbj[] = {0.35,0.45,0.45,0.55,0.6};


   
 const int nkin = sizeof(kin)/sizeof(*kin);
  const double xval[] = {0.325,0.375,0.425,0.475,0.525,0.575,0.625,0.675,0.725,0.775};
  const int nx = sizeof(xval)/sizeof(*xval);
  for(int ikin=0;ikin<nkin;ikin++){
   ofstream outfile(Form("./TextFilesSum/weighted_average_%.2f.csv",kin[ikin]));
   for(int ix=0;ix<nx;ix++){
    double kinematics, x, y, err;
    vector<double>Kinematics, X, Y, Err;
    ifstream infile("./TextFilesSum/weighted_avg.csv");
    if(!infile){
      cout<<"File doesn't exist. Quitting!!!"<<endl;
      exit(0);
    }
    infile.ignore(10000,'\n');

    while(infile>>kinematics>>x>>y>>err){
     if(kinematics==kin[ikin] && x==xval[ix]){
       Kinematics.push_back(kinematics);
       X.push_back(x);
       Y.push_back(y);
       Err.push_back(err);
     }
    }
    
    int npt = X.size();
    double ent[npt];

    for(int i=0;i<X.size();i++){
      ent[i] = i;
      cout<<Kinematics[i]<<"\t"<<X[i]<<"\t"<<Y[i]<<"\t"<<Err[i]<<endl;
    }
    if(npt!=0){
     TCanvas* c1 = new TCanvas();
     TGraphErrors* gr = new TGraphErrors(npt,ent,&Y[0],0,&Err[0]);
     gr->SetMarkerStyle(20);
     gr->SetTitle(Form("Weighted Average for kin = %.2f and x = %.3f;npt;y",kin[ikin],xval[ix]));
     gr->Draw("AP");
     gr->Fit("pol0");
     TF1* f1= (TF1*)gr->GetFunction("pol0");
     gr->GetXaxis()->Set(npt,-0.5,npt-0.5);
     outfile<<kin[ikin]<<"\t"<<xval[ix]<<"\t"<<f1->GetParameter(0)<<"\t"<<f1->GetParError(0)<<endl;
     c1->SaveAs(Form("./plots/temp/plot_%.2f_%.3f.pdf",kin[ikin],xval[ix]));
     infile.close();
    }
  }
  outfile.close();
 }
 gSystem->Exec(Form("pdfunite ./plots/temp/plot_*.pdf ./plots/plot_all.pdf"));
 gSystem->Exec(Form("rm -rf ./plots/temp/plot_*.pdf"));

 TCanvas* c_multi = new TCanvas("c_multi");

 TGraphErrors* gr1[nkin];

 for(int ikin=0;ikin<nkin;ikin++){
  ifstream inf(Form("./TextFilesSum/weighted_average_%.2f.csv",kin[ikin]));
  double kk, xx, yy, ee;
  vector<double> KK, XX, YY, EE;
  while(inf>>kk>>xx>>yy>>ee){
    KK.push_back(kk);
    XX.push_back(xx);
    YY.push_back(yy);
    EE.push_back(ee);
  }
  int nn = KK.size();
  gr1[ikin] = new TGraphErrors(nn,&XX[0],&YY[0],0,&EE[0]);
  gr1[ikin]->SetMarkerStyle(20);
  // gr1[ikin]->SetMarkerColor(ikin+1);
  // gr1[ikin]->SetLineColor(ikin+1);
  gr1[ikin]->SetMarkerColor(ikin+1);
  gr1[ikin]->SetLineColor(ikin+1);
 }
 TLegend* lg = new TLegend(0.15,0.70,0.30,0.90);
 // TMultiGraph* mg = new TMultiGraph("mg","Weighted Y vs X;X;Y (weighted)");
  TMultiGraph* mg = new TMultiGraph("mg","Sum Ratio of H2 +/- and D2 +/- (CSV Data);Z; Sum Ratio");
for(int ikin=0;ikin<nkin;ikin++){
   mg->Add(gr1[ikin]);
   // lg->AddEntry(gr1[ikin],Form("kin = %.2f GeV2",kin[ikin]),"lp");
   lg->AddEntry(gr1[ikin],Form("Q2 = %.3f GeV2, x =  %.2f",kin[ikin],xbj[ikin]),"lp");

 }
 mg->Draw("AP");//gPad->SetGrid();
 mg->GetYaxis()->SetRangeUser(0.3,0.8);
 mg->GetXaxis()->SetRangeUser(0,0.77);

 
 lg->Draw();
 c_multi->SaveAs("./plots/weighted_sumratio_vs_z.pdf");
}
