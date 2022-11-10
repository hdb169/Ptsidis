This script will include 5 textfiles, that are to be compared and plotted on top of each other.
  these textfilesare under textfile directory. This file compares my multiplicity vs Peter's multiplicity. 

void plot_multiplicity(){
  const int nfile = 5;
  const int nmult = 2;
  const int mult[2] = {18,25};
  const int color[nmult] = {1,2};
  TGraphErrors* gPpnr[nfile][nmult];
  TGraphErrors* gPpwr[nfile][nmult];
  TGraphErrors* gDpnr[nfile][nmult];
  TGraphErrors* gDpwr[nfile][nmult];
  TGraphErrors* gPmnr[nfile][nmult];
  TGraphErrors* gPmwr[nfile][nmult];
  TGraphErrors* gDmnr[nfile][nmult];
  TGraphErrors* gDmwr[nfile][nmult];
  TMultiGraph* mg[nfile][5];
  TCanvas* c1[nfile];
  double myx[nfile];
  double myq2[nfile];
 
  for(int ifile=0;ifile<nfile;ifile++){
    for(int imult=0;imult<nmult;imult++){
      ifstream infile(Form("textfiles/mult_%d_%d.txt",mult[imult],ifile+1));
      if(ifile==0){
        infile.ignore(10000,'\n');
        infile.ignore(10000,'\n');
      }

      vector<double> Kin, X, Q2, Z, Z_off, Ppnr, Ppwr, Pper, Dpnr, Dpwr, Dper, Pmnr, Pmwr, Pmer, Dmnr, Dmwr, Dmer;
      double kin, x, q2, z, ppnr, ppwr, pper, dpnr, dpwr, dper, pmnr, pmwr, pmer, dmnr, dmwr, dmer;

      while(infile>>kin>>x>>q2>>z>>ppnr>>ppwr>>pper>>dpnr>>dpwr>>dper>>pmnr>>pmwr>>pmer>>dmnr>>dmwr>>dmer){
        Kin.push_back(kin);
        X.push_back(x);
        Q2.push_back(q2);
        Z.push_back(z);
        Z_off.push_back(z+0.01);
        Ppnr.push_back(ppnr);
        Ppwr.push_back(ppwr);
        Pper.push_back(pper);
        Dpnr.push_back(dpnr);
        Dpwr.push_back(dpwr);
        Dper.push_back(dper);
        Pmnr.push_back(pmnr);
        Pmwr.push_back(pmwr);
        Pmer.push_back(pmer);
        Dmnr.push_back(dmnr);
        Dmwr.push_back(dmwr);
        Dmer.push_back(dmer);
      }
      int npt = Kin.size();
      cout<<"Data points: "<<npt<<endl;
      myx[ifile] = X[0];
      myq2[ifile] = Q2[0];
      for(int ipt=0;ipt<npt;ipt++){
       cout<<Kin[ipt]<<"\t"<<X[ipt]<<"\t"<<Q2[ipt]<<"\t"<<Z[ipt]<<"\t"<<Z_off[ipt]<<"\t"<<Ppnr[ipt]<<"\t"<<Ppwr[ipt]<<"\t"<<Pper[ipt]<<"\t"<<Dpnr[ipt]<<"\t"<<Dpwr[ipt]<<"\t"<<Dper[ipt]<<"\t"<<Pmnr[ipt]<<"\t"<<Pmwr[ipt]<<"\t"<<Pmer[ipt]<<"\t"<<Dmnr[ipt]<<"\t"<<Dmwr[ipt]<<"\t"<<Dmer[ipt]<<endl;
      }
      gPpnr[ifile][imult] = new TGraphErrors(npt,&Z[0],&Ppnr[0],0,&Pper[0]);
      gPpwr[ifile][imult] = new TGraphErrors(npt,&Z_off[0],&Ppwr[0],0,&Pper[0]);
      gDpnr[ifile][imult] = new TGraphErrors(npt,&Z[0],&Dpnr[0],0,&Pper[0]);
      gDpwr[ifile][imult] = new TGraphErrors(npt,&Z_off[0],&Dpwr[0],0,&Pper[0]);
      gPmnr[ifile][imult] = new TGraphErrors(npt,&Z[0],&Pmnr[0],0,&Pmer[0]);
      gPmwr[ifile][imult] = new TGraphErrors(npt,&Z_off[0],&Pmwr[0],0,&Pmer[0]);
      gDmnr[ifile][imult] = new TGraphErrors(npt,&Z[0],&Dmnr[0],0,&Pmer[0]);
      gDmwr[ifile][imult] = new TGraphErrors(npt,&Z_off[0],&Dmwr[0],0,&Pmer[0]);

      gPpnr[ifile][imult]->SetMarkerStyle(20);
      gPpwr[ifile][imult]->SetMarkerStyle(24);
      gDpnr[ifile][imult]->SetMarkerStyle(20);
      gDpwr[ifile][imult]->SetMarkerStyle(24);
      gPmnr[ifile][imult]->SetMarkerStyle(20);
      gPmwr[ifile][imult]->SetMarkerStyle(24);
      gDmnr[ifile][imult]->SetMarkerStyle(20);
      gDmwr[ifile][imult]->SetMarkerStyle(24);

      gPpnr[ifile][imult]->SetMarkerSize(0.5);
      gPpwr[ifile][imult]->SetMarkerSize(0.5);
      gDpnr[ifile][imult]->SetMarkerSize(0.5);
      gDpwr[ifile][imult]->SetMarkerSize(0.5);
      gPmnr[ifile][imult]->SetMarkerSize(0.5);
      gPmwr[ifile][imult]->SetMarkerSize(0.5);
      gDmnr[ifile][imult]->SetMarkerSize(0.5);
      gDmwr[ifile][imult]->SetMarkerSize(0.5);

      gPpnr[ifile][imult]->SetMarkerColor(color[imult]);
      gPpwr[ifile][imult]->SetMarkerColor(color[imult]);
      gDpnr[ifile][imult]->SetMarkerColor(color[imult]);
      gDpwr[ifile][imult]->SetMarkerColor(color[imult]);
      gPmnr[ifile][imult]->SetMarkerColor(color[imult]);
      gPmwr[ifile][imult]->SetMarkerColor(color[imult]);
      gDmnr[ifile][imult]->SetMarkerColor(color[imult]);
      gDmwr[ifile][imult]->SetMarkerColor(color[imult]);
    }
  }
  for(int ifile=0;ifile<nfile;ifile++){
    c1[ifile] = new TCanvas(Form("c1[%d]",ifile),"canvas",1200,600);
    c1[ifile]->Divide(2,2);
    c1[ifile]->cd(1);gPad->SetLogy();
    mg[ifile][0] = new TMultiGraph(Form("mg[%d][0]",ifile+1),Form("p pi+ x=%.3f Q^{2}=%.3f GeV2;Z;multiplicity",myx[0],myq2[0]));
    mg[ifile][0]->Add(gPpnr[ifile][0],"P");
    mg[ifile][0]->Add(gPpnr[ifile][1],"P");
    mg[ifile][0]->Add(gPpwr[ifile][0],"P");
    mg[ifile][0]->Add(gPpwr[ifile][1],"P");
    mg[ifile][0]->Draw("AP");
    mg[ifile][0]->GetYaxis()->SetRangeUser(0.05,1.5);;

    TLegend *leg0 = new TLegend(0.7,0.6,0.9,0.9);
    leg0->AddEntry(gPpnr[ifile][0],"multi 18 no rho","ep");
    leg0->AddEntry(gPpwr[ifile][0],"multi 18 with rho","ep");
    leg0->AddEntry(gPpnr[ifile][1],"multi 25 no rho","ep");
    leg0->AddEntry(gPpwr[ifile][1],"multi 25 with rho","ep");
    leg0->SetTextSize(0.04);
    leg0->SetBorderSize(0);
    leg0->SetFillStyle(0);
    leg0->Draw();
    c1[ifile]->cd(2);gPad->SetLogy();
    mg[ifile][1] = new TMultiGraph(Form("mg[%d][1]",ifile+1),Form("d pi+ x=%.3f Q^{2}=%.3f GeV2;Z;multiplicity",myx[1],myq2[1]));
    mg[ifile][1]->Add(gDpnr[ifile][0],"P");
    mg[ifile][1]->Add(gDpnr[ifile][1],"P");
    mg[ifile][1]->Add(gDpwr[ifile][0],"P");
    mg[ifile][1]->Add(gDpwr[ifile][1],"P");
    mg[ifile][1]->Draw("AP");
    mg[ifile][1]->GetYaxis()->SetRangeUser(0.05,1.5);;

    TLegend *leg1 = new TLegend(0.7,0.6,0.9,0.9);
    leg1->AddEntry(gDpnr[ifile][0],"multi 18 no rho","ep");
    leg1->AddEntry(gDpwr[ifile][0],"multi 18 with rho","ep");
    leg1->AddEntry(gDpnr[ifile][1],"multi 25 no rho","ep");
    leg1->AddEntry(gDpwr[ifile][1],"multi 25 with rho","ep");
    leg1->SetTextSize(0.04);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    leg1->Draw();
    c1[ifile]->cd(3);gPad->SetLogy();
    mg[ifile][2] = new TMultiGraph(Form("mg[%d][2]",ifile+1),Form("p pi- x=%.3f Q^{2}=%.3f GeV2 ;Z;multiplicity",myx[2],myq2[2]));
    mg[ifile][2]->Add(gPmnr[ifile][0],"P");
    mg[ifile][2]->Add(gPmnr[ifile][1],"P");
    mg[ifile][2]->Add(gPmwr[ifile][0],"P");
    mg[ifile][2]->Add(gPmwr[ifile][1],"P");
    mg[ifile][2]->Draw("AP");
    mg[ifile][2]->GetYaxis()->SetRangeUser(0.05,1.5);;
 TLegend *leg2 = new TLegend(0.7,0.6,0.9,0.9);
    leg2->AddEntry(gPmnr[ifile][0],"multi 18 no rho","ep");
    leg2->AddEntry(gPmwr[ifile][0],"multi 18 with rho","ep");
    leg2->AddEntry(gPmnr[ifile][1],"multi 25 no rho","ep");
    leg2->AddEntry(gPmwr[ifile][1],"multi 25 with rho","ep");
    leg2->SetTextSize(0.04);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->Draw();
    c1[ifile]->cd(4);gPad->SetLogy();
    mg[ifile][3] = new TMultiGraph(Form("mg[%d][3]",ifile+1),Form("d pi- x=%.3f Q^{2}=%.3f GeV2;Z;multiplicity",myx[3],myq2[3]));
    mg[ifile][3]->Add(gDmnr[ifile][0],"P");
    mg[ifile][3]->Add(gDmnr[ifile][1],"P");
    mg[ifile][3]->Add(gDmwr[ifile][0],"P");
    mg[ifile][3]->Add(gDmwr[ifile][1],"P");
    mg[ifile][3]->Draw("AP");
    mg[ifile][3]->GetYaxis()->SetRangeUser(0.05,1.5);;

    TLegend *leg3 = new TLegend(0.7,0.6,0.9,0.9);
    leg3->AddEntry(gDmnr[ifile][0],"multi 18 no rho","ep");
    leg3->AddEntry(gDmwr[ifile][0],"multi 18 with rho","ep");
    leg3->AddEntry(gDmnr[ifile][1],"multi 25 no rho","ep");
    leg3->AddEntry(gDmwr[ifile][1],"multi 25 with rho","ep");
    leg3->SetTextSize(0.04);
    leg3->SetBorderSize(0);
    leg3->SetFillStyle(0);
    leg3->Draw();
    c1[ifile]->SaveAs(Form("plots/plot%d.pdf",ifile));
  }
  gSystem->Exec(Form("pdfunite plots/plot*.pdf plots/multiplicity_plots.pdf"));
  gSystem->Exec(Form("rm -rf plots/plot*.pdf"));
}
