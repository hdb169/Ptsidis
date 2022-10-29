//Ydata = rad corr
//Ysimc = inc Norad
//Needs to rewrittten-Oct18-2020
#include "ROOT/RDataFrame.hxx"
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TLine.h"
#include "TH1.h"
#include "THStack.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "nlohmann/json.hpp"
using json = nlohmann::json;
#include "Get_all_eff.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>

double Get_weighted_mean(std::vector<double> value,std::vector<double> err){//weighted av of yields for each bin
  double weighted_mean = 0;
  double weighted_sigma = 0;
  for(int i = 0;i<err.size();++i){
    if(err[i]!=0){
      weighted_mean += value[i]/(err[i]*err[i]);
      weighted_sigma += 1/(err[i]*err[i]);
    }
  }
  double mean = 0;
  if(weighted_mean!=0){
    mean = weighted_mean/weighted_sigma;
  }
  return mean;
};
double Get_weighted_sigma(std::vector<double> value,std::vector<double> err){//error of yield of each bin
  double weighted_sigma = 0;
  for(int i = 0;i<err.size();++i){
    if(err[i]!=0){
      weighted_sigma += 1/(err[i]*err[i]);
    }
  }
  double sig = 0;
  if(weighted_sigma!=0){
    sig = std::sqrt(1/weighted_sigma);
  }
  return sig;
};
std::vector<double> Get_weighted_averages(std::vector<int> RunNumbers,int bins){
  std::vector<double> weighted_yield;
  std::vector<double> weighted_yield_err;
  json j_info;
  {
    //  std::ifstream ifs("/home/hdbhatt/yieldnew/db2/runs_info_6073_6067_6022_trash.json");//all run by run eff etc
    // std::ifstream ifs("/home/hdbhatt/yieldnew/db2/ptsidis_csv_eff_april2_no6073.json");//all run by run eff etc for csv
    std::ifstream ifs("/home/hdbhatt/yieldnew/db2/ptsidis_csv_eff_april2.json");//all run by run eff etc for csv

    ifs>>j_info;
  }
  std::cout<<"Hello"<<std::endl;
  for(int i = 1;i<=bins;++i){
    std::vector<double> yields;
    std::vector<double> yield_errs;
    double bin_center = 0;
    for(auto it = RunNumbers.begin();it!=RunNumbers.end();++it){
      int RunNumber = *it;
      double charge = j_info[(std::to_string(RunNumber)).c_str()]["ChargeCorr"].get<double>();
      //can get charge
      TFile *root_file = new TFile(("/home/hdbhatt/yield/yield_skim_oct20/count_"+std::to_string(RunNumber)+".root").c_str());//all the data count hists
      TH1D *h_z = new TH1D("","",bins,0,1);
      h_z = (TH1D*)root_file->Get("h_zhad_clean_mom_counts");
      // TH1D *h_z_bg = new TH1D("","",bins,0,1);
      // h_z_bg = (TH1D*)root_file->Get("z_bg");
      double EFF = Get_all_eff(RunNumber);
      // h_z->Add(h_z_bg,-1.0/6);
      yields.push_back(h_z->GetBinContent(i)/(charge*EFF));
      yield_errs.push_back(sqrt(h_z->GetBinContent(i))/(charge*EFF));
      bin_center = h_z->GetBinCenter(i);
    }
    double yield_i = Get_weighted_mean(yields,yield_errs);
    //double yield_err_i = Get_weighted_sigma(yields,yield_errs);
    weighted_yield.push_back(yield_i);
    //weighted_yield_err.push_back(yield_err_i);
  }
  return weighted_yield;//one for one bin, many bins => vect
};
std::vector<double> Get_weighted_errs(std::vector<int> RunNumbers,int bins){
  std::vector<double> weighted_yield;
  std::vector<double> weighted_yield_err;
  json j_info;
  {
    //  std::ifstream ifs("/home/hdbhatt/yieldnew/db2/runs_info_6073_6067_6022_trash.json");;//all run by run eff etc
    // std::ifstream ifs("/home/hdbhatt/yieldnew/db2/ptsidis_csv_eff_april2_no6073.json");
    std::ifstream ifs("/home/hdbhatt/yieldnew/db2/ptsidis_csv_eff_april2.json");

    ifs>>j_info;
  }
  
  for(int i = 1;i<=bins;++i){
    std::vector<double> yields;
    std::vector<double> yield_errs;
    double bin_center = 0;
    for(auto it = RunNumbers.begin();it!=RunNumbers.end();++it){
      int RunNumber = *it;
      double charge = j_info[(std::to_string(RunNumber)).c_str()]["ChargeCorr"].get<double>();

      TFile *root_file =new TFile(("/home/hdbhatt/yield/yield_skim_oct20/count_"+std::to_string(RunNumber)+".root").c_str());//counts hists
      //cout<<"charge = "<<charge<<endl;


      TH1D *h_z = new TH1D("","",bins,0,1);
      h_z = (TH1D*)root_file->Get("h_zhad_clean_mom_counts");
      // TH1D *h_z_bg = new TH1D("","",bins,0,1);
      // h_z_bg = (TH1D*)root_file->Get("z_bg");
      double EFF = Get_all_eff(RunNumber);
      // h_z->Add(h_z_bg,-1.0/6);
      yields.push_back(h_z->GetBinContent(i)/(charge*EFF));
      yield_errs.push_back(sqrt(h_z->GetBinContent(i))/(charge*EFF));
      bin_center = h_z->GetBinCenter(i);
    }
    // double yield_i = Get_weighted_mean(yields,yield_errs);
    double yield_err_i = Get_weighted_sigma(yields,yield_errs);
    //weighted_yield.push_back(yield_i);
    weighted_yield_err.push_back(yield_err_i);
  }
  return weighted_yield_err;//earror of each bin = vect
};
std::vector<double> Get_x_axis(std::vector<int> RunNumbers,int bins){
  std::vector<double> x_axis;
  
  for(int i = 1;i<=bins;++i){
    double bin_center = 0;
    int RunNumber = RunNumbers[0];
    //TFile *root_file =new TFile(("/home/hdbhatt/yield/yield_root_20bins_hgcer/count_histo_"+std::to_string(RunNumber)+".root").c_str());//counts hists
    TFile *root_file =new TFile(("/home/hdbhatt/yield/yield_skim_oct20/count_"+std::to_string(RunNumber)+".root").c_str());//counts hists
    TH1D *h_z = new TH1D("","",bins,0,1);
    h_z = (TH1D*)root_file->Get("h_zhad_clean_mom_counts");
    bin_center = h_z->GetBinCenter(i);
    x_axis.push_back(bin_center);
  }
  return x_axis;
};


//int plot_Q2x_ratio_yieldcorr(){
//int rd_meas_20bins_hgcerindata_notinsimc(){//this file copied to below
//int  ptsidis_factorizationKin500510520_D2overH2(){
//int Finalsimc_csv_h2_d2_factorization_plot_april29(){
//int testH2simc_kin1314_ptsidis_errorJune3(){

//int original_yieldBack_testjuly1(){
//int test_cs_by_printing(){
//  int test_cs_by_printing_Ynoradcorr_simcincrad(){
int data_simc_comp_CS(){
  //int yield_and_yield_err_poster(){
  //int test_cs_by_printing_Ynoradcorr_simcincrad(){
  // ofstream txtsum("new500_txt_sum_500_d2overh2_yeshgcer_july24.txt",ios::app);
  // ofstream txtdiff("new500_txt_diff_500_d2overh2_yeshgcer_july24.txt",ios::app);

  // ofstream txtsum1("new500_check_sum_yeshgcer_july24.txt",ios::app);
  // ofstream txtdiff1("new500_check_diff_yeshgcer_july24.txt",ios::app);

  ofstream txtsum("textfiles/may24_pt_d2overh2_sum_yeshgcer_july24.txt",ios::app);
  ofstream txtdiff("textfiles/may24_pt_d2overh2_diff_yeshgcer_july24.txt",ios::app);
  //july8
  ofstream txtd2pos("textfilesjuly/ypos_data_over_simc_rad_yeshgcer_july24.txt",ios::app);
  ofstream txtd2neg("textfilesjuly/yneg_data_over_simc_rad_yeshgcer_july24.txt",ios::app);
  ofstream txth2pos("textfilesjuly/yh2pos_data_over_simc_rad_yeshgcer_july24.txt",ios::app);
  ofstream txth2neg("textfilesjuly/yh2neg_data_over_simc_rad_yeshgcer_july24.txt",ios::app);

  ofstream txtd2pos_hgc("textfilesjuly_yeshgc/ypos_data_over_simc_rad_yeshgcer_july24.txt",ios::app);
  ofstream txtd2neg_hgc("textfilesjuly_yeshgc/yneg_data_over_simc_rad_yeshgcer_july24.txt",ios::app);
  ofstream txth2pos_hgc("textfilesjuly_yeshgc/yh2pos_data_over_simc_rad_yeshgcer_july24.txt",ios::app);
  ofstream txth2neg_hgc("textfilesjuly_yeshgc/yh2neg_data_over_simc_rad_yeshgcer_july24.txt",ios::app);

  //july8

  ofstream txtsum_peter("textfiles/peter_sum_yeshgcer_july24.txt",ios::app);
  ofstream txtdiff_peter("textfiles/peter_diff_yeshgcer_july24.txt",ios::app);

  ofstream txtd2diff_oversum("textfiles/d2diff_over_d2sum_yeshgcer_july24.txt",ios::app);
  ofstream txth2diff_oversum("textfiles/h2diff_over_h2sum_yeshgcer_july24.txt",ios::app);

  //Aug1
  /*
  ofstream txtyieldd2pos("textfilesposter/d2pos_yield.txt",ios::app);
  ofstream txtyieldd2neg("textfilesposter/d2neg_yield.txt",ios::app);
  ofstream txtyieldh2pos("textfilesposter/h2pos_yield.txt",ios::app);
  ofstream txtyieldh2neg("textfilesposter/h2neg_yield.txt",ios::app);
*/
  
  ofstream txtyieldd2pos("textfilesposter/yd2pos_data_over_simc_rad_poster.txt",ios::app);
  ofstream txtyieldd2neg("textfilesposter/yd2neg_data_over_simc_rad_poster.txt",ios::app);
  ofstream txtyieldh2pos("textfilesposter/yh2pos_data_over_simc_rad_poster.txt",ios::app);
  ofstream txtyieldh2neg("textfilesposter/yh2neg_data_over_simc_rad_poster.txt",ios::app);
  
  ofstream txtyieldd2pos_aug15("textfiles_datasimc_comp/yd2pos_data_simc_rad.txt",ios::app);
  ofstream txtyieldd2neg_aug15("textfiles_datasimc_comp/yd2neg_data_simc_rad.txt",ios::app);
  ofstream txtyieldh2pos_aug15("textfiles_datasimc_comp/yh2pos_data_simc_rad.txt",ios::app);
  ofstream txtyieldh2neg_aug15("textfiles_datasimc_comp/yh2neg_data_simc_rad.txt",ios::app);
  //oct20, D2-/D2+, H2-/D2+, H2+/D2+
 ofstream txtfile_d2neg_to_d2pos("textfiles_ratio_tod2pos/yd2neg_to_d2pos1.txt",ios::app);
 ofstream txtfile_h2neg_to_d2pos("textfiles_ratio_tod2pos/yh2neg_to_d2pos1.txt",ios::app);
 ofstream txtfile_h2pos_to_d2pos("textfiles_ratio_tod2pos/yh2pos_to_d2pos1.txt",ios::app);
 ofstream txtfile_all("textfiles_ratio_tod2pos/h2p_d2p_h2n_d2n.txt",ios::app);

  //oct20 done


  //  ofstream txtsum_sim("may11_csv_d2overh2_sum_sim_yeshgcer_july24.txt",ios::app);
  //ofstream txtdiff_sim("may11_csv_d2overh2_diff_sim_yeshgcer_july24.txt",ios::app);
 
  //ofstream txtsum1("relative_err_sum_yeshgcer_july24.txt",ios::app);
  // ofstream txtdiff1("relative_err_diff_yeshgcer_july24.txt",ios::app);

  json j_Q2x;
  {
    //   std::ifstream runs("db2/kin_group_xQ2z_6073_trash.json");
    //std::ifstream runs("/home/hdbhatt/yieldnew/db2/kin_group_xQ2z_6073_trash_trash.json");
    //GOOD std::ifstream runs("/home/hdbhatt/yieldnew/db2/kin_group_H2D2.json");//took only H2 D2 rund
    //std::ifstream runs("/home/hdbhatt/yieldnew/db2/kin_group_H2D2_trash.json");//took only H2 D2 rund

    // std::ifstream runs("/home/hdbhatt/yieldnew/db2/kin_group_H2D2_trash_no180RunGr.json");//took only H2 D2 rund
    // std::ifstream runs("/home/hdbhatt/yieldnew/db2/kin_group_H2D2_no100_180_550_540_RunGr.json");//took only H2 D2 rund
    // std::ifstream runs("/home/hdbhatt/yieldnew/db2/kin_group_H2D2_70_80_90_100.json");//took only H2 D2 rund
    //======april29   std::ifstream runs("/home/hdbhatt/yieldnew/db2/kin_group_H2D2_500_510_520.json");//took only H2 D2 rund
    //The following is the best one for CSV H2 D2 data 
   

     std::ifstream runs("/home/hdbhatt/yieldnew/db2/kin_group_H2D2_csv_only.json");
    //std::ifstream runs("/home/hdbhatt/yieldnew/db2/kin_group_xQ2z_no6073.json");



//has on D2-H2 rootfile json for Fact check
    //  std::ifstream runs("/home/hdbhatt/yieldnew/db2/kin_group_xQ2z_ptsidikin1314.json");//has on D2-H2 rootfile json for Fact check
    //std::ifstream runs("db2/kin_group_xQ2z_ptsidikin1314_590_to_770_not_710_no5404no590.json");//has on D2-H2 rootfile json for Fact check

    //     std::ifstream runs("db2/trash620.json");//has on D2-H2 rootfile json for Fact check
    //     std::ifstream runs("db2/trash600.json");//has on D2-H2 rootfile json for Fact check
    // std::ifstream runs("db2/kin_group_runs_missing.json");

    // std::ifstream runs("db2/kin_group_no200.json");
    //std::ifstream runs("db2/kin_group_beyond200.json");
    // std::ifstream runs("db2/kin_group_200_only1.json");
    runs>>j_Q2x;
  }
  json j_info;
  {
    //std::ifstream ifs("/home/hdbhatt/yieldnew/db2/runs_info_6073_6067_6022_trash.json");
    std::ifstream ifs("/home/hdbhatt/yieldnew/db2/ptsidis_csv_eff_april2.json");
    ifs>>j_info;
  }
  json j_cuts;
  {
    std::ifstream ifs("/home/hdbhatt/yieldnew/db2/all_cut.json");
    ifs>>j_cuts;
  }
  int bins = j_cuts["bins"].get<int>();

  int coolcolor[11] = {4,3,7,39,38,37,36,35,34,33,32};
  int warmcolor[11] = {2,6,46,45,44,43,42,41,40,47,48};
  auto mg_all = new TMultiGraph(); 
  std::vector<int> i_which_x = {0,1,2,0,3,1,0,2,1,2,3,3};
  std::vector<int> i_which_Q2= {0,0,0,1,0,1,2,1,2,2,1,2};
  int i_whichx = 0; 
  json jout;

  for(json::iterator it = j_Q2x.begin();it!=j_Q2x.end();++it){
    double xbj = std::stod(it.key());
    auto j_Q2z = it.value();
    auto mg_x_all = new TMultiGraph();
    for(json::iterator it  = j_Q2z.begin();it!=j_Q2z.end();++it){
      double Q2 = std::stod(it.key());
      auto j_z = it.value();//h2, d2, al
      // std::string canvas_name = "x:"+std::to_string(xbj).substr(0,4)+",Q2:"+std::to_string(Q2).substr(0,5);
      std::string canvas_name = "x =   "+std::to_string(xbj).substr(0,4)+",Q2 =   "+std::to_string(Q2).substr(0,5);

      std::string canvas_filename = "x_Q2_"+std::to_string(100*xbj).substr(0,2)+"_"+std::to_string(1000*Q2).substr(0,4);
      std::string q2x_name = "x:"+std::to_string(xbj).substr(0,4)+",Q2:"+std::to_string(Q2).substr(0,5)+"_yieldratio";
      std::string q2x_filename = "x_Q2_"+std::to_string(100*xbj).substr(0,2)+"_"+std::to_string(1000*Q2).substr(0,4)+"_yieldratio";
      int i_color = 1;
      auto mg = new TMultiGraph();
      auto mg_frag = new TMultiGraph();
      auto mg_RD = new TMultiGraph();
      //THStack* hs = new THStack("yield_ratio","yield ratio");


      TGraphErrors* sum_ratio = new TGraphErrors();
      // std::string sum_string = "sum z setting "+(std::to_string(z)).substr(0,4);
      // sum_ratio->SetName(sum_string.c_str());

      TGraphErrors* diff_ratio = new TGraphErrors();
      // std::string diff_string = "diff z setting "+(std::to_string(z)).substr(0,4);
      // diff_ratio->SetName(diff_string.c_str());


      TGraphErrors* sum_ratio_sim = new TGraphErrors();//May11
      TGraphErrors* diff_ratio_sim = new TGraphErrors();//May11



      if(xbj !=0 && Q2!=0){  //z loop

        for(json::iterator it = j_z.begin();it!=j_z.end();++it){
          double z = std::stod(it.key());
          std::string q2xz_str = "x:"+std::to_string(xbj).substr(0,4)+",Q2:"+std::to_string(Q2).substr(0,5)+",z:"+std::to_string(z).substr(0,4);
          std::string q2xz_str_filename = "x_Q2_z_"+std::to_string(100*xbj).substr(0,2)+"_"+std::to_string(1000*Q2).substr(0,4)+"_"+std::to_string(100*z).substr(0,2);
	  // std::string q2xz_str_filename = "x_Q2_z_"+std::to_string(100*xbj).substr(0,2)+"_"+std::to_string(1000*Q2).substr(0,4)+"_"+std::to_string(100*z).substr(0,2);


          std::cout<<"xbj "<<xbj<<" Q2 "<<Q2<<" z "<<z<<std::endl;
          std::vector<int> neg_D2_runs, pos_D2_runs,neg_Dummy_runs,pos_Dummy_runs;
          std::vector<double> neg_D2_yields,neg_D2_errs;
          std::vector<double> pos_D2_yields,pos_D2_errs;
          std::vector<double> neg_Dummy_yields,neg_Dummy_errs;
          std::vector<double> pos_Dummy_yields,pos_Dummy_errs;
          //added 1
	  std::vector<int> neg_H2_runs, pos_H2_runs;
          std::vector<double> neg_H2_yields,neg_H2_errs;
          std::vector<double> pos_H2_yields,pos_H2_errs;
          
          ///add1 done



	  std::vector<double> x_axis;



	  std::cout<< "working fine 1..."<<std::endl;//----------cout 1

          auto runjs = it.value();

          int RunGroup = runjs["group_num"].get<int>();
          std::cout<<"RunGroup "<<RunGroup<<std::endl;

	  // TFile *rootfile_neg_simH2 = new TFile(("/home/hdbhatt/yieldsimc2022/results/yield20binsptt/kinematics_yield_"+std::to_string(RunGroup)+"_simc.root").c_str());
          //TFile *rootfile_neg_sim = new TFile(("/home/hdbhatt/yieldsimc2022/results/yield20binsptt/kinematics_D2_yield_"+std::to_string(RunGroup)+"_simc.root").c_str());

	  TFile *rootfile_neg_simH2 = new TFile(("/home/hdbhatt/yieldsimc2022/results/yield20binsptt_oct/kinematics_yield_"+std::to_string(RunGroup)+"_simc_geomcuts_nomx2now2.root").c_str());
          TFile *rootfile_neg_sim = new TFile(("/home/hdbhatt/yieldsimc2022/results/yield20binsptt_oct/kinematics_D2_yield_"+std::to_string(RunGroup)+"_simc_geomcuts_nomx2now2.root").c_str());


          TH1D *h_z_neg_sim = new TH1D("","pi- sim sum",bins,0,1);
          TH1D *h_z_neg_sim_incrad = new TH1D("","pi- sim SIDIS",bins,0,1);
	  
	  //doing one by one simc hist cout D2 


	  h_z_neg_sim_incrad = (TH1D*)rootfile_neg_sim->Get("z_neg_inc_rad");


	  TH1D *h_z_neg_sim_incnorad = new TH1D("","pi- sim SIDIS",bins,0,1);
          h_z_neg_sim_incnorad = (TH1D*)rootfile_neg_sim->Get("z_neg_inc_norad");

	  //sig
	  TH1D *h_z_neg_sim_incnorad_sig = new TH1D("","pi- sim SIDIS",bins,0,1);
          h_z_neg_sim_incnorad_sig = (TH1D*)rootfile_neg_sim->Get("z_neg_inc_norad_sighad");

	  //sigdone


          TH1D *h_z_neg_sim_excrad = new TH1D("","pi- sim exc",bins,0,1);
	 
	  h_z_neg_sim_excrad = (TH1D*)rootfile_neg_sim->Get("z_neg_exc_rad");
        

	  TH1D *h_z_neg_sim_rho = new TH1D("","pi- sim rho",bins,0,1);
          h_z_neg_sim_rho = (TH1D*)rootfile_neg_sim->Get("z_neg_rho");
       

	  TH1D *h_z_neg_sim_delta = new TH1D("","pi- sim delta",bins,0,1);
          h_z_neg_sim_delta = (TH1D*)rootfile_neg_sim->Get("z_neg_delta");
       	 


	  TH2D *h_xs_z_neg_sim = new TH2D("","pi- sim xs",bins,0,1,bins,0,1);
          h_xs_z_neg_sim = (TH2D*)rootfile_neg_sim->Get("xs_z_neg_inc_rad");


        

          h_z_neg_sim->Add(h_z_neg_sim_incrad,1);
	  h_z_neg_sim->Add(h_z_neg_sim_excrad,1);
          h_z_neg_sim->Add(h_z_neg_sim_rho,1);
	  h_z_neg_sim->Add(h_z_neg_sim_delta,1);

	  TCanvas *m = new TCanvas ("m","m",1200,800);
	  h_z_neg_sim->Draw();
	  m->SaveAs("plot.pdf");

	  //
	  
	  std::cout<<"RunGroup "<<RunGroup<<"   "<<"h inc rad neg D2 = "<< h_z_neg_sim_incrad<< std::endl;/////////////////////////Check simc file exits
	  std::cout<<"RunGroup "<<RunGroup<<"   "<<"h inc no rad neg D2 = "<< h_z_neg_sim_incnorad<< std::endl;
	  std::cout<<"RunGroup "<<RunGroup<<"   "<<"h inc exc rad neg D2 = "<<h_z_neg_sim_excrad<< std::endl;
	  std::cout<<"RunGroup "<<RunGroup<<"   "<<"h neg rho neg  D2 = "<<h_z_neg_sim_rho<< std::endl;
	  std::cout<<"RunGroup "<<RunGroup<<"   "<<"h neg delta neg D2 = "<< h_z_neg_sim_delta<< std::endl;
	  std::cout<<"RunGroup "<<RunGroup<<"   "<<"h neg D2 Sum   = "<< h_z_neg_sim<< std::endl;

	  //H2===================================================================================H2H222222222222222222222H22222222222222222H2222222222222H222222222HHHHHHHHHHHHHHHHHH2HHHHHHHHHHHHHH2HHHHHHH2HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222HHHHHHHHHHHHHHH222222222222

          TH1D *h2_z_neg_sim = new TH1D("","pi- sim sum",bins,0,1);
          TH1D *h2_z_neg_sim_incrad = new TH1D("","pi- sim SIDIS",bins,0,1);
	  h2_z_neg_sim_incrad = (TH1D*)rootfile_neg_simH2->Get("z_neg_inc_rad");
	  
	  TH1D *h2_z_neg_sim_incnorad = new TH1D("","pi- sim SIDIS",bins,0,1);
          h2_z_neg_sim_incnorad = (TH1D*)rootfile_neg_simH2->Get("z_neg_inc_norad");
         
	  //sig
 
	  TH1D *h2_z_neg_sim_incnorad_sig = new TH1D("","pi- sim SIDIS",bins,0,1);
          h2_z_neg_sim_incnorad_sig = (TH1D*)rootfile_neg_simH2->Get("z_neg_inc_norad_sighad");
         
	  //sig

	  TH1D *h2_z_neg_sim_rho = new TH1D("","pi- sim rho",bins,0,1);
          h2_z_neg_sim_rho = (TH1D*)rootfile_neg_simH2->Get("z_neg_rho");
	  TH1D *h2_z_neg_sim_delta = new TH1D("","pi- sim delta",bins,0,1);
          h2_z_neg_sim_delta= (TH1D*)rootfile_neg_simH2->Get("z_neg_delta");

	  //cs H2-
	  TH1D *h2_z_neg_sim_sig_incrad = new TH1D("","pi- sim sigma",bins,0,1);
          h2_z_neg_sim_sig_incrad = (TH1D*)rootfile_neg_simH2->Get("sighad_neg_inc_rad");

	  std::cout<<"RunGroup "<<RunGroup<<"   "<<"h inc rad neg h2 = "<< h2_z_neg_sim_incrad<< std::endl;/////////////////////////Check simc file exits
	  std::cout<<"RunGroup "<<RunGroup<<"   "<<"h inc no rad neg h2 = "<< h2_z_neg_sim_incnorad<< std::endl;
	  std::cout<<"RunGroup "<<RunGroup<<"   "<<"h neg rho neg  h2 = "<<h2_z_neg_sim_rho<< std::endl;
	  std::cout<<"RunGroup "<<RunGroup<<"   "<<"h neg delta neg h2 = "<< h2_z_neg_sim_delta<< std::endl;



	  h2_z_neg_sim->Add(h2_z_neg_sim_incrad,1);
          h2_z_neg_sim->Add(h2_z_neg_sim_rho,1);
	  h2_z_neg_sim->Add(h2_z_neg_sim_delta,1);
	  //these are hist mamnes 

	  
	  // // TH1D *h2_z_neg_sim = new TH1D("","pi- sim sum",bins,0,1);
	  // // TH1D *h2_z_neg_sim_incrad = new TH1D("","pi- sim SIDIS",bins,0,1);
       
	  // //h2_z_neg_sim_incrad = (TH1D*)rootfile_neg_simH2->Get("z_neg_inc_rad");
          
	  // //TH1D *h2_z_neg_sim_incnorad = new TH1D("","pi- sim SIDIS",bins,0,1);
	  // //h2_z_neg_sim_incnorad = (TH1D*)rootfile_neg_simH2->Get("z_neg_inc_norad");
	  //   //TH1D *h2_z_neg_sim_excrad = new TH1D("","pi- sim exc",bins,0,1);
	  //   // h2_z_neg_sim_excrad = (TH1D*)rootfile_neg_sim->Get("z_neg_exc_rad");
	  // // TH1D *h2_z_neg_sim_rho = new TH1D("","pi- sim rho",bins,0,1);
	  // // h2_z_neg_sim_rho = (TH1D*)rootfile_neg_simH2->Get("z_neg_rho");
	   
	  //  TH2D *h2_xs_z_neg_sim = new TH2D("","pi- sim xs",bins,0,1,bins,0,1);
	  //  h2_xs_z_neg_sim = (TH2D*)rootfile_neg_sim->Get("xs_z_neg_inc_rad");
	  //  h2_z_neg_sim->Add(h2_z_neg_sim_incrad,1);
	  // h2_z_neg_sim->Add(h2_z_neg_sim_excrad,1);
	  //  h2_z_neg_sim->Add(h2_z_neg_sim_rho,1);
	  //h2 don
	    

	  TFile *rootfile_pos_simH2 = new TFile(("/home/hdbhatt/yieldsimc2022/results/yield20binsptt_oct/kinematics_yield_"+std::to_string(RunGroup)+"_simc_geomcuts_nomx2now2.root").c_str());//H2
	  TFile *rootfile_pos_sim   = new TFile(("/home/hdbhatt/yieldsimc2022/results/yield20binsptt_oct/kinematics_D2_yield_"+std::to_string(RunGroup)+"_simc_geomcuts_nomx2now2.root").c_str());//D2



	  TH1D *h_z_pos_sim = new TH1D("","pi+ sim sum",bins,0,1);
          TH1D *h_z_pos_sim_incrad = new TH1D("","pi+ sim SIDIS",bins,0,1);
          h_z_pos_sim_incrad = (TH1D*)rootfile_pos_sim->Get("z_pos_inc_rad");
         
	  TH1D *h_z_pos_sim_incnorad = new TH1D("","pi+ sim SIDIS",bins,0,1);
	  h_z_pos_sim_incnorad = (TH1D*)rootfile_pos_sim->Get("z_pos_inc_norad");

	  //sig
	  TH1D *h_z_pos_sim_incnorad_sig = new TH1D("","pi+ sim SIDIS",bins,0,1);
	  h_z_pos_sim_incnorad_sig = (TH1D*)rootfile_pos_sim->Get("z_pos_inc_norad_sighad");
	  //sig

	  

         
	  TH1D *h_z_pos_sim_excrad = new TH1D("","pi+ sim exc",bins,0,1);
          h_z_pos_sim_excrad = (TH1D*)rootfile_pos_sim->Get("z_pos_exc_rad");
          TH1D *h_z_pos_sim_rho = new TH1D("","pi+ sim rho",bins,0,1);
          h_z_pos_sim_rho = (TH1D*)rootfile_pos_sim->Get("z_pos_rho");
          TH1D *h_z_pos_sim_delta = new TH1D("","pi+ sim delta",bins,0,1);
          h_z_pos_sim_delta = (TH1D*)rootfile_pos_sim->Get("z_pos_delta");
          TH2D *h_xs_z_pos_sim = new TH2D("","pi+ sim xs",bins,0,1,bins,0,1);
          h_xs_z_pos_sim = (TH2D*)rootfile_pos_sim->Get("xs_z_pos_inc_rad");

	 
	  std::cout<<"RunGroup "<<RunGroup<<"   "<<"h inc rad pos D2 = "<< h_z_pos_sim_incrad<< std::endl;/////////////////////////Check simc file exits
	  std::cout<<"RunGroup "<<RunGroup<<"   "<<"h inc no rad pos D2 = "<< h_z_pos_sim_incnorad<< std::endl;
	  std::cout<<"RunGroup "<<RunGroup<<"   "<<"h inc exc rad pos D2 = "<<h_z_pos_sim_excrad<< std::endl;
	  std::cout<<"RunGroup "<<RunGroup<<"   "<<"h pos rho pos  D2 = "<<h_z_pos_sim_rho<< std::endl;
	  std::cout<<"RunGroup "<<RunGroup<<"   "<<"h pos delta pos D2 = "<< h_z_pos_sim_delta<< std::endl;


          h_z_pos_sim->Add(h_z_pos_sim_incrad,1);
          h_z_pos_sim->Add(h_z_pos_sim_excrad,1);
          h_z_pos_sim->Add(h_z_pos_sim_rho,1);
	  h_z_pos_sim->Add(h_z_pos_sim_delta,1);
	  


	  //H2
	 

	  // // TH1D *h2_z_pos_sim = new TH1D("","pi+ sim sum",bins,0,1);
          // // TH1D *h2_z_pos_sim_incrad = new TH1D("","pi+ sim SIDIS",bins,0,1);
          // // h2_z_pos_sim_incrad = (TH1D*)rootfile_pos_simH2->Get("z_pos_inc_rad");
          // // TH1D *h2_z_pos_sim_incnorad = new TH1D("","pi+ sim SIDIS",bins,0,1);
          // // h2_z_pos_sim_incnorad = (TH1D*)rootfile_pos_simH2->Get("z_pos_inc_norad");
          // // TH1D *h2_z_pos_sim_rho = new TH1D("","pi+ sim rho",bins,0,1);
          // // h2_z_pos_sim_rho = (TH1D*)rootfile_pos_sim->Get("z_pos_rho");
          // // h2_z_pos_sim->Add(h2_z_pos_sim_incrad,1);
          // // h2_z_pos_sim->Add(h2_z_pos_sim_rho,1);


	  TH1D *h2_z_pos_sim = new TH1D("","pi+ sim sum",bins,0,1);
	  TH1D *h2_z_pos_sim_incrad = new TH1D("","pi+ sim SIDIS",bins,0,1);
	  h2_z_pos_sim_incrad = (TH1D*)rootfile_pos_simH2->Get("z_pos_inc_rad");
	  TH1D *h2_z_pos_sim_incnorad = new TH1D("","pi+ sim SIDIS",bins,0,1);
	  h2_z_pos_sim_incnorad = (TH1D*)rootfile_pos_simH2->Get("z_pos_inc_norad");
        
	  //sig
	  TH1D *h2_z_pos_sim_incnorad_sig = new TH1D("","pi+ sim SIDIS",bins,0,1);
	  h2_z_pos_sim_incnorad_sig = (TH1D*)rootfile_pos_simH2->Get("z_pos_inc_norad_sighad");
        
	  //sig


	  TH1D *h2_z_pos_sim_excrad = new TH1D("","pi+ sim exc",bins,0,1);
	  h2_z_pos_sim_excrad = (TH1D*)rootfile_pos_simH2->Get("z_pos_exc_rad");
          TH1D *h2_z_pos_sim_rho = new TH1D("","pi+ sim rho",bins,0,1);
          h2_z_pos_sim_rho = (TH1D*)rootfile_pos_simH2->Get("z_pos_rho");
	  TH1D *h2_z_pos_sim_delta = new TH1D("","pi+ sim delta",bins,0,1);
	  h2_z_pos_sim_delta = (TH1D*)rootfile_pos_simH2->Get("z_pos_delta");
          TH2D *h2_xs_z_pos_sim = new TH2D("","pi+ sim xs",bins,0,1,bins,0,1);
          h2_xs_z_pos_sim = (TH2D*)rootfile_pos_simH2->Get("xs_z_pos_inc_rad");


	  //CS H2+

	  TH1D *h2_z_pos_sim_sig_incrad = new TH1D("","pi+ sim sigma",bins,0,1);
          h2_z_pos_sim_sig_incrad = (TH1D*)rootfile_pos_simH2->Get("sighad_pos_inc_rad");



	  std::cout<<"RunGroup "<<RunGroup<<"   "<<"h inc rad pos H2 = "<< h2_z_pos_sim_incrad<< std::endl;/////////////////////////Check simc file exits
	  std::cout<<"RunGroup "<<RunGroup<<"   "<<"h inc no rad pos H2 = "<< h2_z_pos_sim_incnorad<< std::endl;
	  std::cout<<"RunGroup "<<RunGroup<<"   "<<"h inc exc rad pos H2 = "<<h2_z_pos_sim_excrad<< std::endl;
	  std::cout<<"RunGroup "<<RunGroup<<"   "<<"h pos rho pos  H2 = "<<h2_z_pos_sim_rho<< std::endl;
	  std::cout<<"RunGroup "<<RunGroup<<"   "<<"h pos delta pos H2 = "<< h2_z_pos_sim_delta<< std::endl;




	  h2_z_pos_sim->Add(h2_z_pos_sim_incrad,1);
          h2_z_pos_sim->Add(h2_z_pos_sim_excrad,1);
          h2_z_pos_sim->Add(h2_z_pos_sim_rho,1);
	  h2_z_pos_sim->Add(h2_z_pos_sim_delta,1);

	  // //H2





	  std::cout<< "working fine 2..."<<std::endl;//----------cout

	  if(z!=0 && !runjs["neg"]["H2"].empty()){// data yield //=======================================================>
	    // if(z!=0){// data yield //=======================================================>


            neg_D2_runs = runjs["neg"]["D2"].get<std::vector<int>>();
            pos_D2_runs = runjs["pos"]["D2"].get<std::vector<int>>();
	    //add2
	    neg_H2_runs = runjs["neg"]["H2"].get<std::vector<int>>();
            pos_H2_runs = runjs["pos"]["H2"].get<std::vector<int>>();
	    ///add2 done


            neg_Dummy_runs = runjs["neg"]["Dummy"].get<std::vector<int>>();
            pos_Dummy_runs = runjs["pos"]["Dummy"].get<std::vector<int>>();

	    //added
	    // if(!neg_H2_runs.empty() && !pos_H2_runs.empty()){
	   
	    std::cout<<"THIS IS IMP"<<std::endl;
	    std::cout<<"check neg first run = "<<neg_D2_runs[0]<<std::endl;//working 
            std::cout<<"check neg second run = "<<neg_D2_runs[1]<<std::endl;//working 
	    std::cout<<"check neg third run = "<<neg_D2_runs[2]<<std::endl;//working 
	    std::cout<<"check pos first run = "<<pos_D2_runs[0]<<std::endl;//working 
            std::cout<<"check pos second run = "<<pos_D2_runs[1]<<std::endl;//working 
	    std::cout<<"check pos third run = "<<pos_D2_runs[2]<<std::endl;//working 
	   
	    std::cout<<"check neg first run = "<<neg_H2_runs[0]<<std::endl;//working 
            std::cout<<"check neg second run = "<<neg_H2_runs[1]<<std::endl;//working 
	    std::cout<<"check neg third run = "<<neg_H2_runs[2]<<std::endl;//working 
	    std::cout<<"check pos first run = "<<pos_H2_runs[0]<<std::endl;//working 
            std::cout<<"check pos second run = "<<pos_H2_runs[1]<<std::endl;//working 
	    std::cout<<"check pos third run = "<<pos_H2_runs[2]<<std::endl;//working 
	  


	    std::cout<<"THIS IS IMP"<<std::endl;

	    std::cout<<"  "<<std::endl;

	    neg_D2_yields = Get_weighted_averages(neg_D2_runs,bins);//not working?old message
            
	    std::cout<<"check neg sec run = "<<neg_D2_runs[1]<<std::endl;

	    std::cout<<neg_D2_yields[9]<<std::endl;
            neg_D2_errs = Get_weighted_errs(neg_D2_runs,bins);
            pos_D2_yields = Get_weighted_averages(pos_D2_runs,bins);
	    std::cout<<"checking yield = "<<pos_D2_yields[9]<<std::endl;

	    std::cout<<pos_D2_yields[9]<<std::endl;
            pos_D2_errs = Get_weighted_errs(pos_D2_runs,bins);

	    //add3
	    neg_H2_yields = Get_weighted_averages(neg_H2_runs,bins);//not working?old message
            
	    std::cout<<"check neg sec run = "<<neg_H2_runs[1]<<std::endl;

	    std::cout<<neg_H2_yields[9]<<std::endl;
            neg_H2_errs = Get_weighted_errs(neg_H2_runs,bins);
            pos_H2_yields = Get_weighted_averages(pos_H2_runs,bins);
	    std::cout<<"checking yield = "<<pos_H2_yields[9]<<std::endl;

	    std::cout<<pos_H2_yields[9]<<std::endl;
            pos_H2_errs = Get_weighted_errs(pos_H2_runs,bins);
	    ///add3 done

            neg_Dummy_yields = Get_weighted_averages(neg_Dummy_runs,bins);
            neg_Dummy_errs = Get_weighted_errs(neg_Dummy_runs,bins);
            pos_Dummy_yields = Get_weighted_averages(pos_Dummy_runs,bins);
            pos_Dummy_errs = Get_weighted_errs(pos_Dummy_runs,bins);
            std::cout<<neg_Dummy_yields[9]<<std::endl;
            std::cout<<pos_Dummy_yields[9]<<std::endl;
            x_axis = Get_x_axis(neg_D2_runs,bins);
            std::cout<<x_axis[9]<<std::endl;
          }//if z not 0



	  std::cout<<"Now Working 1"<<std::endl;
          auto h_z_neg_incnorad = h_z_neg_sim_incnorad;
	  auto h_z_neg_incnorad_sig = h_z_neg_sim_incnorad_sig;
	  auto h_z_neg_incrad = h_z_neg_sim_incrad;
          auto h_z_pos_incnorad = h_z_pos_sim_incnorad;
	  auto h_z_pos_incnorad_sig = h_z_pos_sim_incnorad_sig;
	  auto h_z_pos_incrad = h_z_pos_sim_incrad;
	  std::cout<<"Now Working  2"<<std::endl;

	  //H2

	  std::cout<<"Now Working 1"<<std::endl;
          auto h2_z_neg_incnorad = h2_z_neg_sim_incnorad;
	  auto h2_z_neg_incnorad_sig = h2_z_neg_sim_incnorad_sig;
	  auto h2_z_neg_incrad = h2_z_neg_sim_incrad;
          auto h2_z_pos_incnorad = h2_z_pos_sim_incnorad;
	  auto h2_z_pos_incnorad_sig = h2_z_pos_sim_incnorad_sig;
	  auto h2_z_pos_incrad = h2_z_pos_sim_incrad;
	  std::cout<<"Now Working  2"<<std::endl;


	  //H2


	  //TH1* rp_radia_corr_neg = (TH1*)h_z_neg_incnorad->Clone("rp_radia_corr");
          //TH1* rp_radia_corr_pos = (TH1*)h_z_pos_incnorad->Clone("rp_radia_corr");
	 
	  // doing rad corr/norad at a time  = norad/rad/norad
	  
	  //In fact this is not rad correction (oct 27, because it is 1/rad)
	  TH1D* rp_radia_corr_neg_pre = new TH1D("rp_radia_corr_neg_pre","",bins,0,1);
          TH1D* rp_radia_corr_pos_pre = new TH1D("rp_radia_corr_pos_pre","",bins,0,1);
	  rp_radia_corr_neg_pre ->Divide(h_z_neg_incnorad,h_z_neg_incrad);          
	  rp_radia_corr_pos_pre ->Divide(h_z_pos_incnorad,h_z_pos_incrad);   
	  
	  TH1D* rp_radia_corr_neg = new TH1D("rp_radia_corr_neg","",bins,0,1);
          TH1D* rp_radia_corr_pos = new TH1D("rp_radia_corr_pos","",bins,0,1);
	  rp_radia_corr_neg ->Divide(rp_radia_corr_neg_pre,h_z_neg_incnorad);          
	  rp_radia_corr_pos ->Divide( rp_radia_corr_pos_pre,h_z_pos_incnorad);   

	  
	  TH1D* rph2_radia_corr_neg_pre = new TH1D("rph2_radia_corr_neg_pre","",bins,0,1);
          TH1D* rph2_radia_corr_pos_pre = new TH1D("rph2_radia_corr_pos_pre","",bins,0,1);
	  rph2_radia_corr_neg_pre ->Divide(h2_z_neg_incnorad,h2_z_neg_incrad);          
	  rph2_radia_corr_pos_pre ->Divide(h2_z_pos_incnorad,h2_z_pos_incrad);
	  
	  TH1D* rph2_radia_corr_neg = new TH1D("rph2_radia_corr_neg","",bins,0,1);
          TH1D* rph2_radia_corr_pos = new TH1D("rph2_radia_corr_pos","",bins,0,1);
	  rph2_radia_corr_neg ->Divide(rph2_radia_corr_neg_pre,h2_z_neg_incnorad);          
	  rph2_radia_corr_pos ->Divide(rph2_radia_corr_pos_pre,h2_z_pos_incnorad);



	  //done
	  /*
        
	  TH1D* rp_radia_corr_neg = new TH1D("rp_radia_corr_neg","",bins,0,1);
          TH1D* rp_radia_corr_pos = new TH1D("rp_radia_corr_pos","",bins,0,1);
	  rp_radia_corr_neg ->Divide(h_z_neg_incnorad,h_z_neg_incrad);          
	  rp_radia_corr_pos ->Divide(h_z_pos_incnorad,h_z_pos_incrad);   


	  TH1D* rph2_radia_corr_neg = new TH1D("rph2_radia_corr_neg","",bins,0,1);
          TH1D* rph2_radia_corr_pos = new TH1D("rph2_radia_corr_pos","",bins,0,1);
	  rph2_radia_corr_neg ->Divide(h2_z_neg_incnorad,h2_z_neg_incrad);          
	  rph2_radia_corr_pos ->Divide(h2_z_pos_incnorad,h2_z_pos_incrad);
*/

	  //sighad cross -section is the ratio of "weight*sighad / weight" hist ratios in z (one is weighted by weight*sighad, and other by weight only)
	  TH1D* rp_sighad_corr_neg = new TH1D("rp_sighad_corr_neg","",bins,0,1);
          TH1D* rp_sighad_corr_pos = new TH1D("rp_sighad_corr_pos","",bins,0,1);
	  rp_sighad_corr_neg ->Divide(h_z_neg_incnorad_sig,h_z_neg_incrad);          
	  rp_sighad_corr_pos ->Divide(h_z_pos_incnorad_sig,h_z_pos_incrad);   


	  TH1D* rph2_sighad_corr_neg = new TH1D("rph2_sighad_corr_neg","",bins,0,1);
          TH1D* rph2_sighad_corr_pos = new TH1D("rph2_sighad_corr_pos","",bins,0,1);
	  rph2_sighad_corr_neg ->Divide(h2_z_neg_incnorad_sig,h2_z_neg_incrad);          
	  rph2_sighad_corr_pos ->Divide(h2_z_pos_incnorad_sig,h2_z_pos_incrad);


	  /* ORIGINAL

	     TH1D* rp_radia_corr_neg = new TH1D("rp_radia_corr_neg","",bins,0,1);
	     TH1D* rp_radia_corr_pos = new TH1D("rp_radia_corr_pos","",bins,0,1);
	     rp_radia_corr_neg ->Divide(h_z_neg_incnorad,h_z_neg_incrad);          
	     rp_radia_corr_pos ->Divide(h_z_pos_incnorad,h_z_pos_incrad);   


	     TH1D* rph2_radia_corr_neg = new TH1D("rph2_radia_corr_neg","",bins,0,1);
	     TH1D* rph2_radia_corr_pos = new TH1D("rph2_radia_corr_pos","",bins,0,1);
	     rph2_radia_corr_neg ->Divide(h2_z_neg_incnorad,h2_z_neg_incrad);          
	     rph2_radia_corr_pos ->Divide(h2_z_pos_incnorad,h2_z_pos_incrad);

	     //sighad cross -section is the ratio of "weight*sighad / weight" hist ratios in z (one is weighted by weight*sighad, and other by weight only)
	     TH1D* rp_sighad_corr_neg = new TH1D("rp_sighad_corr_neg","",bins,0,1);
	     TH1D* rp_sighad_corr_pos = new TH1D("rp_sighad_corr_pos","",bins,0,1);
	     rp_sighad_corr_neg ->Divide(h_z_neg_incnorad_sig,h_z_neg_incrad);          
	     rp_sighad_corr_pos ->Divide(h_z_pos_incnorad_sig,h_z_pos_incrad);   


	     TH1D* rph2_sighad_corr_neg = new TH1D("rph2_sighad_corr_neg","",bins,0,1);
	     TH1D* rph2_sighad_corr_pos = new TH1D("rph2_sighad_corr_pos","",bins,0,1);
	     rph2_sighad_corr_neg ->Divide(h2_z_neg_incnorad_sig,h2_z_neg_incrad);          
	     rph2_sighad_corr_pos ->Divide(h2_z_pos_incnorad_sig,h2_z_pos_incrad);


	     ORIGINAL

	  */


	  //july8 note I am not doing rad corr in data. so without chaning script, I will divide by same quantity or not multiply by rad corr factor
	  
	  //H2

 
     	  std::cout<<"Now Working  4"<<std::endl;

	  //H2


	  TCanvas* ch2_radia_pos = new TCanvas();
          rph2_radia_corr_pos->Draw();
	  std::cout<<"Now Working  5"<<std::endl;
          std::string ch2_radia_pos_name = "/u/group/c-csv/hdbhatt/yielddec21/yieldJan_2022/results/yield/statistics_corr/radiah2_corr_neg"+q2xz_str_filename+"_pos_hgcer.pdf";
          ch2_radia_pos->SaveAs(ch2_radia_pos_name.c_str());
          TCanvas* ch2_radia_neg = new TCanvas();
          rph2_radia_corr_neg->Draw();
          std::string ch2_radia_neg_name = "/u/group/c-csv/hdbhatt/yielddec21/yieldJan_2022/results/yield/statistics_corr/radiah2_corr_neg"+q2xz_str_filename+"_neg_hgcer.pdf";
          ch2_radia_neg->SaveAs(ch2_radia_neg_name.c_str());

     	  std::cout<<"Now Working  5"<<std::endl;

	  //D2

          TCanvas* c_radia_pos = new TCanvas();
          rp_radia_corr_pos->Draw();
	  std::cout<<"Now Working  5"<<std::endl;
          std::string c_radia_pos_name = "/u/group/c-csv/hdbhatt/yielddec21/yieldJan_2022/results/yield/statistics_corr/radia_corr_neg"+q2xz_str_filename+"_pos_hgcer.pdf";
          c_radia_pos->SaveAs(c_radia_pos_name.c_str());
          TCanvas* c_radia_neg = new TCanvas();
          rp_radia_corr_neg->Draw();
          std::string c_radia_neg_name = "/u/group/c-csv/hdbhatt/yielddec21/yieldJan_2022/results/yield/statistics_corr/radia_corr_neg"+q2xz_str_filename+"_neg_hgcer.pdf";
          c_radia_neg->SaveAs(c_radia_neg_name.c_str());


	 
     	  std::cout<<"Now Working  6"<<std::endl;


	  //D2
          TGraphErrors* g_z_pos_yield = new TGraphErrors();
          g_z_pos_yield->SetTitle("pi+ data yield");
          TGraphErrors* g_z_pos_Dummy_yield = new TGraphErrors();
          g_z_pos_Dummy_yield->SetTitle("pi+ Dummy yield");
          TGraphErrors* g_z_pos_yield_corr = new TGraphErrors();
          g_z_pos_yield_corr->SetTitle("pi+ data yield");
          int nbins_yield_pos = bins;
	 


	  //add4
	  TGraphErrors* g_z_pos_yieldH2 = new TGraphErrors();
          g_z_pos_yieldH2->SetTitle("pi+ data yield H2");
	  // TGraphErrors* g_z_pos_Dummy_yield = new TGraphErrors();
          //g_z_pos_Dummy_yield->SetTitle("pi+ Dummy yield");
          TGraphErrors* g_z_pos_yieldH2_corr = new TGraphErrors();
          g_z_pos_yieldH2_corr->SetTitle("pi+ data yield H2");
          //int nbins_yield_pos = bins;
	  //add4 done


	  std::cout<<"Now Working  6"<<std::endl;




	  //just to see the yield

	  int ii_yield_pos = 0;
	  //add?????????????????

          for(int i = 0;i<nbins_yield_pos;i++){
            double x = x_axis[i]; 
            double y = pos_D2_yields[i];
            double y_dummy = pos_Dummy_yields[i]; 
            double error = pos_D2_errs[i]; 
            double error_dummy = pos_Dummy_errs[i]; 
            double y_corr = y-0.245*y_dummy;
            double error_corr = sqrt(error*error+0.245*0.245*error_dummy*error_dummy);
	  
	    ///add5 begin
            double yH2 = pos_H2_yields[i];
            double errorH2 = pos_H2_errs[i]; 
	    double yH2_corr = yH2-0.262*y_dummy;//==================================================0.262
            double errorH2_corr = sqrt(errorH2*errorH2+0.262*0.262*error_dummy*error_dummy);

	    ///add5 done


            //std::cout<<i<<" x "<<x<<" y "<<y<<" err_corr "<<error_corr<<std::endl;
            if(y!=0){
              g_z_pos_yield->SetPoint(ii_yield_pos,x,y);
              g_z_pos_yield->SetPointError(ii_yield_pos,0,error);
              g_z_pos_Dummy_yield->SetPoint(ii_yield_pos,x,y_dummy);
              g_z_pos_Dummy_yield->SetPointError(ii_yield_pos,0,error_dummy);
              g_z_pos_yield_corr->SetPoint(ii_yield_pos,x,y_corr);
              g_z_pos_yield_corr->SetPointError(ii_yield_pos,0,error_corr);
	      //add 6 begin
	     
	      g_z_pos_yieldH2->SetPoint(ii_yield_pos,x,yH2);
              g_z_pos_yieldH2->SetPointError(ii_yield_pos,0,errorH2);
              //g_z_pos_Dummy_yield->SetPoint(ii_yield_pos,x,y_dummy);
              //g_z_pos_Dummy_yield->SetPointError(ii_yield_pos,0,error_dummy);
              g_z_pos_yieldH2_corr->SetPoint(ii_yield_pos,x,yH2_corr);
              g_z_pos_yieldH2_corr->SetPointError(ii_yield_pos,0,errorH2_corr);

	      //add 6 done

              ii_yield_pos++;
            }
          }

          TCanvas *c_yield_pos = new TCanvas();
          gStyle->SetOptTitle(0);
          TMultiGraph* mg_z_yield_pos = new TMultiGraph(); 
          g_z_pos_yield->GetYaxis()->SetTitle("yield");
          g_z_pos_yield->GetYaxis()->SetTitleSize(0.53);
          g_z_pos_yield->GetXaxis()->SetRangeUser(0.3,1);
          g_z_pos_yield->SetMarkerStyle(4);
         
	  mg_z_yield_pos->Add(g_z_pos_yield,"P");

          g_z_pos_Dummy_yield->SetMarkerStyle(4);
          g_z_pos_Dummy_yield->SetMarkerColor(kRed);
          g_z_pos_Dummy_yield->SetTitle("data Dummy_yield");
          mg_z_yield_pos->Add(g_z_pos_Dummy_yield,"P");
          g_z_pos_yield_corr->SetMarkerStyle(4);
          g_z_pos_yield_corr->SetMarkerColor(kRed);
          g_z_pos_yield_corr->SetTitle("data yield final");
          mg_z_yield_pos->Add(g_z_pos_yield_corr,"P");


          TGraph* g_z_pos_sim = new TGraph(h_z_pos_sim);
          TGraph* g_z_pos_sim_incrad = new TGraph(h_z_pos_sim_incrad);
          TGraph* g_z_pos_sim_incnorad = new TGraph(h_z_pos_sim_incnorad);
          TGraph* g_z_pos_sim_excrad = new TGraph(h_z_pos_sim_excrad);
          TGraph* g_z_pos_sim_rho = new TGraph(h_z_pos_sim_rho);
          TGraph* g_z_pos_sim_delta = new TGraph(h_z_pos_sim_delta);

          mg_z_yield_pos->Add(g_z_pos_sim); 
          g_z_pos_sim_incrad->SetLineColor(kOrange);
          g_z_pos_sim_incrad->SetTitle("pi+ sim SIDIS");
          mg_z_yield_pos->Add(g_z_pos_sim_incrad); 
	  g_z_pos_sim_excrad->SetLineColor(kBlue);
	  g_z_pos_sim_excrad->SetTitle("pi+ sim exc");
          mg_z_yield_pos->Add(g_z_pos_sim_excrad); 
          g_z_pos_sim_rho->SetLineColor(kRed);
          g_z_pos_sim_rho->SetTitle("pi+ sim rho");
          mg_z_yield_pos->Add(g_z_pos_sim_rho); 
          g_z_pos_sim_delta->SetLineColor(46);
          g_z_pos_sim_delta->SetTitle("pi+ sim delta");
	  mg_z_yield_pos->Add(g_z_pos_sim_delta); 
          g_z_pos_sim_incnorad->SetLineColor(40);
          g_z_pos_sim_incnorad->SetTitle("pi+ sim inc norad");
          mg_z_yield_pos->Add(g_z_pos_sim_incnorad,"L");
          mg_z_yield_pos->SetTitle(q2xz_str.c_str());
          mg_z_yield_pos->GetYaxis()->SetTitle("yield");
          mg_z_yield_pos->GetXaxis()->SetTitle("z (Pi+)");
          mg_z_yield_pos->Draw("AL");
          c_yield_pos->BuildLegend(0.75,0.75,1,1);
          std::string c_yield_pos_name = "/u/group/c-csv/hdbhatt/yielddec21/yieldJan_2022/results/yield/statistics_corr/yield_"+q2xz_str_filename+"_pos_hgcer.pdf";
          c_yield_pos->SaveAs(c_yield_pos_name.c_str());
	  




	  ///////////// ADD H2 Pos


	  TCanvas *c_yieldH2_pos = new TCanvas();
          gStyle->SetOptTitle(0);
          TMultiGraph* mg_z_yieldH2_pos = new TMultiGraph(); 
          g_z_pos_yieldH2->GetYaxis()->SetTitle("yield");
          g_z_pos_yieldH2->GetYaxis()->SetTitleSize(0.53);
          g_z_pos_yieldH2->GetXaxis()->SetRangeUser(0.3,1);
          g_z_pos_yieldH2->SetMarkerStyle(4);
         
	  mg_z_yieldH2_pos->Add(g_z_pos_yieldH2,"P");

	  // g_z_pos_Dummy_yield->SetMarkerStyle(4);
          //g_z_pos_Dummy_yield->SetMarkerColor(kRed);
	  // g_z_pos_Dummy_yield->SetTitle("data Dummy_yield");
          mg_z_yieldH2_pos->Add(g_z_pos_Dummy_yield,"P");
          g_z_pos_yieldH2_corr->SetMarkerStyle(4);
          g_z_pos_yieldH2_corr->SetMarkerColor(kRed);
          g_z_pos_yieldH2_corr->SetTitle("data yield final");
          mg_z_yieldH2_pos->Add(g_z_pos_yieldH2_corr,"P");


          TGraph* gh2_z_pos_sim = new TGraph(h2_z_pos_sim);
          TGraph* gh2_z_pos_sim_incrad = new TGraph(h2_z_pos_sim_incrad);
          TGraph* gh2_z_pos_sim_incnorad = new TGraph(h2_z_pos_sim_incnorad);
          TGraph* gh2_z_pos_sim_excrad = new TGraph(h2_z_pos_sim_excrad);
          TGraph* gh2_z_pos_sim_rho = new TGraph(h2_z_pos_sim_rho);
          TGraph* gh2_z_pos_sim_delta = new TGraph(h2_z_pos_sim_delta);

          mg_z_yieldH2_pos->Add(g_z_pos_sim); 
	  g_z_pos_sim_incrad->SetLineColor(kOrange);
	  g_z_pos_sim_incrad->SetTitle("pi+ sim SIDIS");
          mg_z_yieldH2_pos->Add(g_z_pos_sim_incrad); 
          g_z_pos_sim_excrad->SetLineColor(kBlue);
          g_z_pos_sim_excrad->SetTitle("pi+ sim exc");
	  mg_z_yieldH2_pos->Add(g_z_pos_sim_excrad); 
          g_z_pos_sim_rho->SetLineColor(kRed);
          g_z_pos_sim_rho->SetTitle("pi+ sim rho");
          mg_z_yieldH2_pos->Add(g_z_pos_sim_rho); 
	  g_z_pos_sim_delta->SetLineColor(46);
          g_z_pos_sim_delta->SetTitle("pi+ sim delta");
	  mg_z_yieldH2_pos->Add(g_z_pos_sim_delta); 
          gh2_z_pos_sim_incnorad->SetLineColor(40);
          gh2_z_pos_sim_incnorad->SetTitle("pi+ sim inc norad");


          mg_z_yieldH2_pos->Add(g_z_pos_sim_incnorad,"L");
          mg_z_yieldH2_pos->SetTitle(q2xz_str.c_str());
          mg_z_yieldH2_pos->GetYaxis()->SetTitle("yield");
          mg_z_yieldH2_pos->GetXaxis()->SetTitle("z (Pi+)");
          mg_z_yieldH2_pos->Draw("AL");
          c_yieldH2_pos->BuildLegend(0.75,0.75,1,1);
          std::string c_yieldH2_pos_name = "/u/group/c-csv/hdbhatt/yielddec21/yieldJan_2022/results/yield/statistics_corr/yieldH2_"+q2xz_str_filename+"_pos_hgcer.pdf";
          c_yieldH2_pos->SaveAs(c_yieldH2_pos_name.c_str());
	  






	  //////ADD H2 Pos Plot done










  
          TGraphErrors* g_z_neg_yield = new TGraphErrors();
          g_z_neg_yield->SetTitle("pi- data yield");
          TGraphErrors* g_z_neg_Dummy_yield = new TGraphErrors();
          g_z_neg_Dummy_yield->SetTitle("pi- Dummy yield");
          TGraphErrors* g_z_neg_yield_corr = new TGraphErrors();
          g_z_neg_yield_corr->SetTitle("pi- data yield");
	
	  ////add 7 start
	  TGraphErrors* g_z_neg_yieldH2 = new TGraphErrors();
          g_z_neg_yieldH2->SetTitle("pi- data yieldH2");
          //TGraphErrors* g_z_neg_Dummy_yield = new TGraphErrors();
	  // g_z_neg_Dummy_yield->SetTitle("pi- Dummy yield");
          TGraphErrors* g_z_neg_yieldH2_corr = new TGraphErrors();
          g_z_neg_yieldH2_corr->SetTitle("pi- data yieldH2");

	  ////add 7 done





	  int nbins_yield_neg = bins;




          int ii_yield_neg = 0;//////////////////////////add not required for this
          for(int i = 0;i<nbins_yield_neg;i++){
            double x = x_axis[i]; 
            double y = neg_D2_yields[i];
            double y_dummy = neg_Dummy_yields[i]; 
            double error = neg_D2_errs[i]; 
            double error_dummy = neg_Dummy_errs[i]; 
            double y_corr = y-0.245*y_dummy;
            double error_corr = sqrt(error*error+0.245*0.245*error_dummy*error_dummy);

	    ///////add 8

	    double yH2 = neg_H2_yields[i];
	    // double y_dummy = neg_Dummy_yields[i]; 
            double errorH2 = neg_H2_errs[i]; 
            //double error_dummy = neg_Dummy_errs[i]; 
            double yH2_corr = yH2-0.262*y_dummy;
            double errorH2_corr = sqrt(errorH2*errorH2+0.262*0.262*error_dummy*error_dummy);

	    //////add 8 done

            //std::cout<<i<<" x "<<x<<" y "<<y<<std::endl;
            if(y!=0){
              g_z_neg_yield->SetPoint(ii_yield_neg,x,y);
              g_z_neg_yield->SetPointError(ii_yield_neg,0,error);
              g_z_neg_Dummy_yield->SetPoint(ii_yield_neg,x,y_dummy);
              g_z_neg_Dummy_yield->SetPointError(ii_yield_neg,0,error_dummy);
              g_z_neg_yield_corr->SetPoint(ii_yield_neg,x,y_corr);
              g_z_neg_yield_corr->SetPointError(ii_yield_neg,0,error_corr);

	      //////add 9 begin
	      g_z_neg_yieldH2->SetPoint(ii_yield_neg,x,yH2);
              g_z_neg_yieldH2->SetPointError(ii_yield_neg,0,errorH2);
	      // g_z_neg_Dummy_yield->SetPoint(ii_yield_neg,x,y_dummy);
	      // g_z_neg_Dummy_yield->SetPointError(ii_yield_neg,0,error_dummy);
              g_z_neg_yieldH2_corr->SetPoint(ii_yield_neg,x,yH2_corr);
              g_z_neg_yieldH2_corr->SetPointError(ii_yield_neg,0,errorH2_corr);


	      /////add9 done





              ii_yield_neg++;
            }
          }





          TCanvas *c_yield_neg = new TCanvas();
          gStyle->SetOptTitle(0);
          TMultiGraph* mg_z_yield_neg = new TMultiGraph(); 
          g_z_neg_yield->GetYaxis()->SetTitle("yield");
          g_z_neg_yield->GetYaxis()->SetTitleSize(0.53);
          g_z_neg_yield->GetXaxis()->SetRangeUser(0.3,1);
          g_z_neg_yield->SetMarkerStyle(4);
          mg_z_yield_neg->Add(g_z_neg_yield,"P");
          g_z_neg_Dummy_yield->SetMarkerStyle(4);
          g_z_neg_Dummy_yield->SetMarkerColor(kRed);
          g_z_neg_Dummy_yield->SetTitle("data Dummy_yield");
          mg_z_yield_neg->Add(g_z_neg_Dummy_yield,"P");
          g_z_neg_yield_corr->SetMarkerStyle(4);
          g_z_neg_yield_corr->SetMarkerColor(kRed);
          g_z_neg_yield_corr->SetTitle("data yield final");
          mg_z_yield_neg->Add(g_z_neg_yield_corr,"P");
          TGraph* g_z_neg_sim = new TGraph(h_z_neg_sim);
          TGraph* g_z_neg_sim_incrad = new TGraph(h_z_neg_sim_incrad);
          TGraph* g_z_neg_sim_incnorad = new TGraph(h_z_neg_sim_incnorad);
          TGraph* g_z_neg_sim_excrad = new TGraph(h_z_neg_sim_excrad);
          TGraph* g_z_neg_sim_rho = new TGraph(h_z_neg_sim_rho);
          TGraph* g_z_neg_sim_delta = new TGraph(h_z_neg_sim_delta);
          mg_z_yield_neg->Add(g_z_neg_sim); 
          g_z_neg_sim_incrad->SetLineColor(kOrange);
          g_z_neg_sim_incrad->SetTitle("pi- sim SIDIS");
          mg_z_yield_neg->Add(g_z_neg_sim_incrad); 
          g_z_neg_sim_excrad->SetLineColor(kBlue);
          g_z_neg_sim_excrad->SetTitle("pi- sim exc");
          mg_z_yield_neg->Add(g_z_neg_sim_excrad); 
          g_z_neg_sim_rho->SetLineColor(kRed);
          g_z_neg_sim_rho->SetTitle("pi- sim rho");
          mg_z_yield_neg->Add(g_z_neg_sim_rho); 
          g_z_neg_sim_delta->SetLineColor(46);
          g_z_neg_sim_delta->SetTitle("pi- sim delta");
          mg_z_yield_neg->Add(g_z_neg_sim_delta); 
          g_z_neg_sim_incnorad->SetLineColor(40);
          g_z_neg_sim_incnorad->SetTitle("pi- sim inc norad");
          mg_z_yield_neg->Add(g_z_neg_sim_incnorad,"L");


          mg_z_yield_neg->SetTitle(q2xz_str.c_str());
          mg_z_yield_neg->GetYaxis()->SetTitle("yield");
          mg_z_yield_neg->GetXaxis()->SetTitle("z (Pi-)");
          mg_z_yield_neg->Draw("AL");
          c_yield_neg->BuildLegend(0.75,0.75,1,1);
          std::string c_yield_neg_name = "/u/group/c-csv/hdbhatt/yielddec21/yieldJan_2022/results/yield/statistics_corr/yield_"+q2xz_str_filename+"_neg_hgcer.pdf";
          c_yield_neg->SaveAs(c_yield_neg_name.c_str());
	 

	  ////add10  started here

	  TCanvas *c_yieldH2_neg = new TCanvas();
          gStyle->SetOptTitle(0);
          TMultiGraph* mg_z_yieldH2_neg = new TMultiGraph(); 
          g_z_neg_yieldH2->GetYaxis()->SetTitle("yield");
          g_z_neg_yieldH2->GetYaxis()->SetTitleSize(0.53);
          g_z_neg_yieldH2->GetXaxis()->SetRangeUser(0.3,1);
          g_z_neg_yieldH2->SetMarkerStyle(4);
          mg_z_yieldH2_neg->Add(g_z_neg_yieldH2,"P");
          g_z_neg_Dummy_yield->SetMarkerStyle(4);
          g_z_neg_Dummy_yield->SetMarkerColor(kRed);
	  g_z_neg_Dummy_yield->SetTitle("data Dummy_yield");
          mg_z_yieldH2_neg->Add(g_z_neg_Dummy_yield,"P");
          g_z_neg_yieldH2_corr->SetMarkerStyle(4);
          g_z_neg_yieldH2_corr->SetMarkerColor(kRed);
          g_z_neg_yieldH2_corr->SetTitle("data yield final");
          mg_z_yieldH2_neg->Add(g_z_neg_yieldH2_corr,"P");
         

	  TGraph* gh2_z_neg_sim = new TGraph(h2_z_neg_sim);
          TGraph* gh2_z_neg_sim_incrad = new TGraph(h2_z_neg_sim_incrad);
          TGraph* gh2_z_neg_sim_incnorad = new TGraph(h2_z_neg_sim_incnorad);
          //TGraph* gh2_z_neg_sim_excrad = new TGraph(h2_z_neg_sim_excrad);
          TGraph* gh2_z_neg_sim_rho = new TGraph(h2_z_neg_sim_rho);
          TGraph* gh2_z_neg_sim_delta = new TGraph(h2_z_neg_sim_delta);
          
	  // // mg_z_yieldH2_neg->Add(g_z_neg_sim); ///============================================================Mistake JH
	  // // gh2_z_neg_sim_incrad->SetLineColor(kOrange);
          // // gh2_z_neg_sim_incrad->SetTitle("pi- sim SIDIS");
          // // mg_z_yieldH2_neg->Add(g_z_neg_sim_incrad); 
	 

	  mg_z_yieldH2_neg->Add(gh2_z_neg_sim); ///============================================================Corrected
	  gh2_z_neg_sim_incrad->SetLineColor(kOrange);
          gh2_z_neg_sim_incrad->SetTitle("pi- sim SIDIS");
          mg_z_yieldH2_neg->Add(gh2_z_neg_sim_incrad); ///============================================================Corrected
	  // gh2_z_neg_sim_excrad->SetLineColor(kBlue);
	  // gh2_z_neg_sim_excrad->SetTitle("pi- sim exc");
	  // mg_z_yieldH2_neg->Add(gh2_z_neg_sim_excrad); 


          gh2_z_neg_sim_rho->SetLineColor(kRed);
          gh2_z_neg_sim_rho->SetTitle("pi- sim rho");
          mg_z_yieldH2_neg->Add(gh2_z_neg_sim_rho); 
          gh2_z_neg_sim_delta->SetLineColor(46);
          gh2_z_neg_sim_delta->SetTitle("pi- sim delta");
	  mg_z_yieldH2_neg->Add(gh2_z_neg_sim_delta); 
          gh2_z_neg_sim_incnorad->SetLineColor(40);
          gh2_z_neg_sim_incnorad->SetTitle("pi- sim inc norad");
          mg_z_yieldH2_neg->Add(gh2_z_neg_sim_incnorad,"L");


          mg_z_yieldH2_neg->SetTitle(q2xz_str.c_str());
          mg_z_yieldH2_neg->GetYaxis()->SetTitle("yield");
          mg_z_yieldH2_neg->GetXaxis()->SetTitle("z (Pi-)");
          mg_z_yieldH2_neg->Draw("AL");
          c_yieldH2_neg->BuildLegend(0.75,0.75,1,1);
          std::string c_yieldH2_neg_name = "/u/group/c-csv/hdbhatt/yielddec21/yieldJan_2022/results/yield/statistics_corr/yieldH2_"+q2xz_str_filename+"_neg_hgcer.pdf";
          c_yieldH2_neg->SaveAs(c_yieldH2_neg_name.c_str());


	  ///add10  done






          h_z_neg_sim_incnorad->Divide(h_z_pos_sim_incnorad);//D2

          int nbins = bins; 

          // TGraphErrors* g_yield_ratio = new TGraphErrors(h_z_neg_all);
          TGraphErrors* g_yield_ratio = new TGraphErrors();
          std::string z_string = "R_Y z setting "+(std::to_string(z)).substr(0,4);
          g_yield_ratio->SetName(z_string.c_str());
          //TGraphErrors* g_frag_ratio = new TGraphErrors();
          //std::string frag_z_string = "frag z setting "+(std::to_string(z)).substr(0,4);
          //g_frag_ratio->SetName(frag_z_string.c_str());
          TGraphErrors* g_RDmeas_ratio = new TGraphErrors();
          std::string RDmeas_z_string = "R^D_meas z setting "+(std::to_string(z)).substr(0,4);
          g_RDmeas_ratio->SetName(RDmeas_z_string.c_str());


	  //
	  //
	  //
	  //
	  //add 10+1===========================================================================================================
	  TGraphErrors* gh2_yield_ratio = new TGraphErrors();
          std::string zh2_string = "R_YH2 z setting "+(std::to_string(z)).substr(0,4);
          gh2_yield_ratio->SetName(zh2_string.c_str());
          //TGraphErrors* g_frag_ratio = new TGraphErrors();
          //std::string frag_z_string = "frag z setting "+(std::to_string(z)).substr(0,4);
          //g_frag_ratio->SetName(frag_z_string.c_str());
          TGraphErrors* gh2_RDmeas_ratio = new TGraphErrors();
          std::string RDmeasH2_z_string = "R^D_meas H2 z setting "+(std::to_string(z)).substr(0,4);
          gh2_RDmeas_ratio->SetName(RDmeasH2_z_string.c_str());

	  //sum diff ratio
	  

	  //add 10+1===========================================================================================================


          int ii = 0;
          for(int i = 0;i<nbins;i++){
            //std::cout<<i<<std::endl;
            double x = x_axis[i];
            double y_neg_D2 = neg_D2_yields[i]; 
            double y_neg_Dummy = neg_Dummy_yields[i]; 
            double y_neg_delta = h_z_neg_sim_delta->GetBinContent(i+1);
            double y_neg_exc = h_z_neg_sim_excrad->GetBinContent(i+1);
            double y_neg_rho = h_z_neg_sim_rho->GetBinContent(i+1);//=====================A
	    double y_neg_incnorad =  h_z_neg_sim_incnorad->GetBinContent(i+1);//May11
	    //sig	   
	    double y_neg_incnorad_sig =  h_z_neg_sim_incnorad_sig->GetBinContent(i+1);//May11
	    //sig
	    //cs
	    ///	double y_neg_sigrad = h_z_neg_sim_sig_incrad->GetBinContent(i+1);//D2
	    // double y_neg = y_neg_D2-0.245*y_neg_Dummy-y_neg_delta-y_neg_exc;//============================================================SIMC SUBTRACTION EXCL,DELTA========D2
	    //working:  double y_neg = y_neg_D2-0.245*y_neg_Dummy-y_neg_exc;//============================================================SIMC SUBTRACTION EXCL,DELTA========D2
	    double y_neg = y_neg_D2-0.245*y_neg_Dummy-y_neg_exc-y_neg_delta-y_neg_rho;//============================================================SIMC SUBTRACTION EXCL,DELTA========D2



	    double error_neg_D2 = neg_D2_errs[i]; 
            double error_neg_Dummy = neg_Dummy_errs[i]; 
            double error_neg = std::sqrt(error_neg_D2*error_neg_D2+0.245*0.245*error_neg_Dummy*error_neg_Dummy);
            double y_pos_D2 = pos_D2_yields[i]; 
            double y_pos_Dummy = pos_Dummy_yields[i]; 
	    double y_pos_delta = h_z_pos_sim_delta->GetBinContent(i+1);
            double y_pos_exc = h_z_pos_sim_excrad->GetBinContent(i+1);
            double y_pos_rho = h_z_pos_sim_rho->GetBinContent(i+1);
	    double y_pos_incnorad =  h_z_pos_sim_incnorad->GetBinContent(i+1);//May11
	    //sig
	    	    double y_pos_incnorad_sig =  h_z_pos_sim_incnorad_sig->GetBinContent(i+1);//May11
//double y_pos_sigrad = h_z_pos_sim_sig_incrad->GetBinContent(i+1);//D2

	    // double y_pos = y_pos_D2-0.245*y_pos_Dummy-y_pos_delta-y_pos_exc;//============================================================SIMC SUBTRACTION EXCL,DELTA============D2
	    //working  double y_pos = y_pos_D2-0.245*y_pos_Dummy-y_pos_exc;//============================================================SIMC SUBTRACTION EXCL,DELTA============D2
	    double y_pos = y_pos_D2-0.245*y_pos_Dummy-y_pos_exc-y_pos_delta-y_pos_rho;//============================================================SIMC SUBTRACTION EXCL,DELTA============D2


	    double error_pos_D2 = pos_Dummy_errs[i]; 
            double error_pos_Dummy = pos_Dummy_errs[i]; 
            double error_pos = std::sqrt(error_pos_D2*error_pos_D2+0.245*0.245*error_pos_Dummy*error_pos_Dummy);
            double radia_corr_neg = rp_radia_corr_neg->GetBinContent(i+1);//vector starts from 0, bins start from 1, 
            double radia_corr_pos = rp_radia_corr_pos->GetBinContent(i+1);
	    //sig
	    double sighad_corr_neg = rp_sighad_corr_neg->GetBinContent(i+1);
            double sighad_corr_pos = rp_sighad_corr_pos->GetBinContent(i+1);
	   

	    

	    //july8

	    //this is when data is with radiatin
	    /*	    
		    double pos_incrad   = h_z_pos_incrad->GetBinContent(i+1);
		    double neg_incrad   = h_z_neg_incrad->GetBinContent(i+1);
		    double posH2_incrad = h2_z_pos_incrad->GetBinContent(i+1);
		    double negH2_incrad = h2_z_neg_incrad->GetBinContent(i+1);

	    */	   

	    //
	    //
	    //
	    //
	    //This is when rad corr is in data and 
	    double pos_incrad   = h_z_pos_incrad->GetBinContent(i+1);
	    double neg_incrad   = h_z_neg_incrad->GetBinContent(i+1);
	    double posH2_incrad = h2_z_pos_incrad->GetBinContent(i+1);
	    double negH2_incrad = h2_z_neg_incrad->GetBinContent(i+1);
	   
	    

	    //july8
	    /*
            y_neg = y_neg*radia_corr_neg;
            y_pos = y_pos*radia_corr_pos;
            error_pos = error_pos*radia_corr_pos;
            error_neg = error_neg*radia_corr_neg;
	    */
	    //No rad corr in data 

	    //oct 27, This is (data-exc-rho-delta/norad)*(norad/rad)*(signorad/norad): this is dsigdzdpT2 
	    /*
	      y_neg = y_neg*radia_corr_neg*sighad_corr_neg/y_neg_incnorad;
	      y_pos = y_pos*radia_corr_pos*sighad_corr_pos/y_pos_incnorad;
	      error_pos = error_pos*radia_corr_pos*sighad_corr_pos/y_pos_incnorad;
	      error_neg = error_neg*radia_corr_neg*sighad_corr_neg/y_neg_incnorad;
	    */

	    //Not eradia corr is not actual radia corr, it is 1/rad = (norad/rad)/norad
	    y_neg = y_neg*radia_corr_neg*sighad_corr_neg;
            y_pos = y_pos*radia_corr_pos*sighad_corr_pos;
            error_pos = error_pos*radia_corr_pos*sighad_corr_pos;
            error_neg = error_neg*radia_corr_neg*sighad_corr_neg;



	    //july8 I am not doing radiative correction in this script so I commented above 4 lines

	    //CS y_neg_cs means yneg multiplied by CS and same for all targets/polarities:
	    //cs june9
	    //double y_neg_cs = y_neg*y_neg_sigrad;
	    //double y_pos_cs = y_pos*y_pos_sigrad;






	    ////////////add 11 start

	    /*
	    // double x = x_axis[i];
            double y_neg_H2 = neg_H2_yields[i]; 
	    // double y_neg_Dummy = neg_Dummy_yields[i]; 
            //double y_neg_delta = h_z_neg_sim_delta->GetBinContent(i+1);
            //double y_neg_exc = h_z_neg_sim_excrad->GetBinContent(i+1);
            double yH2_neg = y_neg_H2-0.262*y_neg_Dummy;//-y_neg_delta;//-y_neg_exc;
          

	    double errorH2_neg = neg_H2_errs[i]; 
	    // double error_neg_Dummy = neg_Dummy_errs[i]; 
            double error_H2_neg = std::sqrt(error_neg_H2*error_neg_H2+0.262*0.262*error_neg_Dummy*error_neg_Dummy);
           


	    double y_pos_H2 = pos_H2_yields[i]; 
            //double y_pos_Dummy = pos_Dummy_yields[i]; 
            //double y_pos_delta = h_z_pos_sim_delta->GetBinContent(i+1);
	    // double y_pos_exc = h_z_pos_sim_excrad->GetBinContent(i+1);
            double yH2_pos = y_pos_H2-0.262*y_pos_Dummy;//-y_pos_delta;//-y_pos_exc;
            double error_pos_H2 = pos_Dummy_errs[i]; 
            //double error_pos_Dummy = pos_Dummy_errs[i]; 
            double errorH2_pos = std::sqrt(error_pos_H2*error_pos_H2+0.262*0.262*error_pos_Dummy*error_pos_Dummy);
	    double radia_corrH2_neg = rph2_radia_corr_neg->GetBinContent(i+1);//vector starts from 0, bins start from 1, 
	    double radia_corrH2_pos = rph2_radia_corr_pos->GetBinContent(i+1);

	    */

	    double y_neg_H2 = neg_H2_yields[i]; 
            double y_neg_delta_H2 = h2_z_neg_sim_delta->GetBinContent(i+1);
	    // double y_neg_exc_H2 = h2_z_neg_sim_excrad->GetBinContent(i+1);// not exist
            double y_neg_rho_H2 = h2_z_neg_sim_rho->GetBinContent(i+1);//=====================A
	    double y_neg_incnorad_H2 =  h2_z_neg_sim_incnorad->GetBinContent(i+1);//May11
	    //sig
	    double y_neg_incnorad_H2_sig =  h2_z_neg_sim_incnorad_sig->GetBinContent(i+1);//May11

	    // working :double yH2_neg = y_neg_H2-0.245*y_neg_Dummy-y_neg_delta; //-y_neg_exc;//============================================================SIMC SUBTRACTION EXCL,DELTA H2
	    double yH2_neg = y_neg_H2-0.262*y_neg_Dummy-y_neg_delta_H2-y_neg_rho_H2; //-y_neg_exc;//============================================================SIMC SUBTRACTION EXCL,DELTA H2

	    double error_neg_H2 = neg_H2_errs[i]; 
	    double errorH2_neg = std::sqrt(error_neg_H2*error_neg_H2+0.245*0.245*error_neg_Dummy*error_neg_Dummy);
	    double y_pos_H2 = pos_H2_yields[i]; 
	    double y_pos_delta_H2 = h2_z_pos_sim_delta->GetBinContent(i+1);
	    double y_pos_exc_H2 = h2_z_pos_sim_excrad->GetBinContent(i+1);
            double y_pos_rho_H2 = h2_z_pos_sim_rho->GetBinContent(i+1);
	    double y_pos_incnorad_H2 =  h2_z_pos_sim_incnorad->GetBinContent(i+1);//May11
	    //sig	    
	    double y_pos_incnorad_H2_sig =  h2_z_pos_sim_incnorad_sig->GetBinContent(i+1);//May11


            double yH2_pos = y_pos_H2-0.262*y_pos_Dummy-y_pos_delta_H2-y_pos_exc_H2-y_pos_rho_H2;//============================================================SIMC SUBTRACTION EXCL,DELTA H2
            double error_pos_H2 = pos_H2_errs[i];// pos_Dummy_errs[i]; (this was on May 24)//=================May23 see 
            double errorH2_pos = std::sqrt(error_pos_H2*error_pos_H2+0.245*0.245*error_pos_Dummy*error_pos_Dummy);
            double radia_corrH2_neg = rph2_radia_corr_neg->GetBinContent(i+1);//vector starts from 0, bins start from 1, 
            double radia_corrH2_pos = rph2_radia_corr_pos->GetBinContent(i+1);

	    //sig H2
            double sighad_corrH2_neg = rph2_sighad_corr_neg->GetBinContent(i+1);//vector starts from 0, bins start from 1, 
            double sighad_corrH2_pos = rph2_sighad_corr_pos->GetBinContent(i+1);



	    //commented below July8 for no rad corr 
	    
	    //uncomm oct 27
	    //sighad
	    /*
	      yH2_neg = yH2_neg*radia_corrH2_neg*sighad_corrH2_neg/y_neg_incnorad_H2;
	      yH2_pos = yH2_pos*radia_corrH2_pos*sighad_corrH2_pos/y_pos_incnorad_H2;
	      errorH2_pos = errorH2_pos*radia_corrH2_pos*sighad_corrH2_pos/y_pos_incnorad_H2;
	      errorH2_neg = errorH2_neg*radia_corrH2_neg*sighad_corrH2_neg/y_neg_incnorad_H2;
	    */
	    yH2_neg = yH2_neg*radia_corrH2_neg*sighad_corrH2_neg;
	    yH2_pos = yH2_pos*radia_corrH2_pos*sighad_corrH2_pos;
	    errorH2_pos = errorH2_pos*radia_corrH2_pos*sighad_corrH2_pos;
	    errorH2_neg = errorH2_neg*radia_corrH2_neg*sighad_corrH2_neg;
	    
	    //cs
	    //double yH2_neg_cs = yH2_neg*y_neg_sigrad_H2;
	    //double yH2_pos_cs = yH2_pos*y_pos_sigrad_H2;


	    ///////add 11 done

            //This convert yield to cross sectio,choose yield or xsn
            //double D2_xs_z_neg,D2_xs_z_pos;
            //D2_xs_z_neg = h_xs_z_neg_sim->GetYaxis()->GetBinCenter(i);
            //D2_xs_z_pos = h_xs_z_pos_sim->GetYaxis()->GetBinCenter(i);
            //double D2_neg_yield_incrad = h_z_neg_sim_incrad->GetBinContent(i);
            //double D2_pos_yield_incrad = h_z_pos_sim_incrad->GetBinContent(i);
            //y_pos = y_D2_pos*D2_xs_z_pos/D2_pos_yield_incrad;
            //y_neg = y_D2_neg*D2_xs_z_neg/D2_neg_yield_incrad;
            //error_pos = error_pos*D2_xs_z_pos/D2_pos_yield_incrad;
            //error_neg = error_neg*D2_xs_z_neg/D2_neg_yield_incrad;


            double y = (y_neg)/(y_pos);//yield ratio
            double error = (y_neg/y_pos)*std::sqrt((error_neg*error_neg)/(y_neg*y_neg)+(error_pos*error_pos)/(y_pos*y_pos));

            
	    ////addd 12 begin
	    double yH2 = (yH2_neg)/(yH2_pos);//yield ratio
            double errorH2 = (yH2_neg/yH2_pos)*std::sqrt((errorH2_neg*errorH2_neg)/(yH2_neg*yH2_neg)+(errorH2_pos*errorH2_pos)/(yH2_pos*yH2_pos));
	    //july8 coommented because there is norad in data and rad in simc

	    //This is when data has radiation
	    /*

	      double ypos_data_over_simc   = y_pos/pos_incrad; 
	      double yneg_data_over_simc   = y_neg/neg_incrad; 
	      double yH2pos_data_over_simc = yH2_pos/posH2_incrad; 
	      double yH2neg_data_over_simc = yH2_neg/negH2_incrad; 


	    */

	    //taking no rad in simc as data is rad corr
	    //commented: this is when data is rad corr
	    double ypos_data_over_simc   = y_pos/pos_incrad; 
	    double yneg_data_over_simc   = y_neg/neg_incrad; 
	    double yH2pos_data_over_simc = yH2_pos/posH2_incrad; 
	    double yH2neg_data_over_simc = yH2_neg/negH2_incrad;

	    double d2pos_percent_diff = abs(y_pos-pos_incrad)/pos_incrad;
	    double d2neg_percent_diff = abs(y_neg-neg_incrad)/neg_incrad;
	    double h2pos_percent_diff = abs(yH2_pos-posH2_incrad)/posH2_incrad;
	    double h2neg_percent_diff = abs(yH2_neg-negH2_incrad)/negH2_incrad;



	    
	    double ypos_data_over_simc_err   = ypos_data_over_simc * std::sqrt((error_pos*error_pos)/(y_pos*y_pos));
	    double yneg_data_over_simc_err   = yneg_data_over_simc * std::sqrt((error_neg*error_neg)/(y_neg*y_neg));
	    double yH2pos_data_over_simc_err = yH2pos_data_over_simc * std::sqrt((errorH2_pos*errorH2_pos)/(yH2_pos*yH2_pos));
	    double yH2neg_data_over_simc_err = yH2neg_data_over_simc * std::sqrt((errorH2_neg*errorH2_neg)/(yH2_neg*yH2_neg));

	    
	      
	    //july 8
         
	    ////add 12 done

	    ///add 13 start///
	    ////Commented because it is for H2 over D2


	    // double yield_sum_ratio = (yH2_pos + yH2_neg)/(y_pos + y_neg);



	    // double yield_diff_ratio = (yH2_pos - yH2_neg)/(y_pos - y_neg);

	    // double yield_num_err = sqrt(errorH2_neg*errorH2_neg + errorH2_pos*errorH2_pos);
	    // double yield_den_err = sqrt(error_neg*error_neg + error_pos*error_pos);

	    // double sum_ratio_err = yield_sum_ratio * sqrt(((yield_num_err/(yH2_pos + yH2_neg)) * (yield_num_err/(yH2_pos + yH2_neg))) + ((yield_den_err/(y_pos + y_neg)) * (yield_den_err/(y_pos + y_neg)))   );

	    // double diff_ratio_err = yield_diff_ratio * sqrt(((yield_num_err/(yH2_pos - yH2_neg)) * (yield_num_err/(yH2_pos - yH2_neg))) + ((yield_den_err/(y_pos - y_neg)) * (yield_den_err/(y_pos - y_neg)))   );



	    //for quick work, I am just changing D2 over H2 sum diff ratio not the H2 over D2 sum diff, with same names etc...
	    ////it is for D2 over H2

	    double yield_sum_ratio = (y_pos + y_neg)/(yH2_pos + yH2_neg);

	    double yield_diff_ratio =(y_pos - y_neg)/(yH2_pos - yH2_neg);

	    double yield_num_err = sqrt(error_neg*error_neg + error_pos*error_pos);

	    double yield_den_err = sqrt(errorH2_neg*errorH2_neg + errorH2_pos*errorH2_pos);
	 
	    double sum_ratio_err = yield_sum_ratio * sqrt(((yield_den_err/(yH2_pos + yH2_neg)) * (yield_den_err/(yH2_pos + yH2_neg))) + ((yield_num_err/(y_pos + y_neg)) * (yield_num_err/(y_pos + y_neg)))   );
	  

	    double diff_ratio_err = yield_diff_ratio * sqrt(((yield_den_err/(yH2_pos - yH2_neg)) * (yield_den_err/(yH2_pos - yH2_neg))) + ((yield_num_err/(y_pos - y_neg)) * (yield_num_err/(y_pos - y_neg)))   );

	    ///oct 20, 2022, D2-/D2+, H2-/D2+, H2+/D2+ yield 
	    double y_d2neg_over_d2pos =  (y_neg/y_pos);
	    double y_h2neg_over_d2pos =  (yH2_neg/y_pos);
	    double y_h2pos_over_d2pos =  (yH2_pos/y_pos);

	    double ratio_errd2neg_over_d2neg = error_neg/y_neg;
	    double ratio_errd2pos_over_d2pos = error_pos/y_pos;
	    double ratio_errh2pos_over_h2pos = errorH2_pos/yH2_pos;
	    double ratio_errh2neg_over_h2neg = errorH2_neg/yH2_neg;


	    double y_d2neg_over_d2pos_err = y_d2neg_over_d2pos*sqrt((ratio_errd2neg_over_d2neg)*(ratio_errd2neg_over_d2neg)+(ratio_errd2pos_over_d2pos)*(ratio_errd2pos_over_d2pos));
	    double y_h2neg_over_d2pos_err = y_h2neg_over_d2pos*sqrt((ratio_errh2neg_over_h2neg)*(ratio_errh2neg_over_h2neg)+(ratio_errd2pos_over_d2pos)*(ratio_errd2pos_over_d2pos));
	    double y_h2pos_over_d2pos_err = y_h2pos_over_d2pos*sqrt((ratio_errh2pos_over_h2pos)*(ratio_errh2pos_over_h2pos)+(ratio_errd2pos_over_d2pos)*(ratio_errd2pos_over_d2pos));

	    //oct 20, 2022 done



	    //May11
	    double yield_sum_ratio_sim = (y_pos_incnorad + y_neg_incnorad)/(y_pos_incnorad_H2 + y_neg_incnorad_H2);
	    double yield_diff_ratio_sim = (y_pos_incnorad - y_neg_incnorad)/(y_pos_incnorad_H2 - y_neg_incnorad_H2);
	    double yield_sum_ratio_sim_err  = 0;
	    double yield_diff_ratio_sim_err = 0; 

	    //==========June 3============Making D2+ - D2 - / D2+ + D2 - ratio

	    double d2_diff_over_sum = (y_pos - y_neg)/(y_pos + y_neg);
	    double err_d2_num = sqrt(error_neg*error_neg + error_pos*error_pos);
	    double err_d2_den = sqrt(error_neg*error_neg + error_pos*error_pos);
	    double d2_ratio_err = d2_diff_over_sum *sqrt(((err_d2_den/(y_pos + y_neg)) * (err_d2_den/(y_pos + y_neg))) + ((err_d2_num/(y_pos - y_neg)) * (err_d2_num/(y_pos - y_neg)))   );

	    double d2_diff_over_sum_simc = (y_pos_incnorad - y_neg_incnorad)/(y_pos_incnorad + y_neg_incnorad);
	    double err_d2_diff_over_sum_simc = 0;

	    //======================Making H2+ - H2 - / H2+ + D2 - ratio
	    double h2_diff_over_sum = (yH2_pos - yH2_neg)/(yH2_pos + yH2_neg);
	    double err_h2_num = sqrt(errorH2_neg*errorH2_neg + errorH2_pos*errorH2_pos);
	    double err_h2_den = sqrt(errorH2_neg*errorH2_neg + errorH2_pos*errorH2_pos);
	    double h2_ratio_err = h2_diff_over_sum *sqrt(((err_h2_den/(yH2_pos + yH2_neg)) * (err_h2_den/(yH2_pos + yH2_neg))) + ((err_h2_num/(yH2_pos - yH2_neg)) * (err_h2_num/(yH2_pos - yH2_neg)))   );
	    double h2_diff_over_sum_simc = (y_pos_incnorad_H2 - y_neg_incnorad_H2)/(y_pos_incnorad_H2 + y_neg_incnorad_H2);
	    double err_h2_diff_over_sum_simc = 0;





	    //Peter's Formula for error May 23 email;
	    //diff ratio
	    double e1_diffratio = (y_pos + error_pos - y_neg)/(yH2_pos - yH2_neg);
	    double e2_diffratio = (y_pos - (y_neg + error_neg))/(yH2_pos - yH2_neg);
	    double e3_diffratio = (y_pos - y_neg)/(yH2_pos - errorH2_pos - yH2_neg);
	    double e4_diffratio = (y_pos - y_neg) / (yH2_pos - (yH2_neg + errorH2_neg));

	    double error_diffratio = sqrt(  pow((yield_diff_ratio - e1_diffratio),2) + pow((yield_diff_ratio - e2_diffratio),2)+ pow((yield_diff_ratio - e3_diffratio),2)+ pow((yield_diff_ratio - e4_diffratio),2)   );

	    double e1_sumratio = (y_pos + error_pos + y_neg)/(yH2_pos + yH2_neg);
	    double e2_sumratio = (y_pos + (y_neg + error_neg))/(yH2_pos + yH2_neg);
	    double e3_sumratio = (y_pos + y_neg)/(yH2_pos + errorH2_pos + yH2_neg);
	    double e4_sumratio = (y_pos + y_neg) / (yH2_pos + (yH2_neg + errorH2_neg));

	    double error_sumratio = sqrt(  pow((yield_sum_ratio - e1_sumratio),2) + pow((yield_sum_ratio - e2_sumratio),2)+ pow((yield_sum_ratio - e3_sumratio),2)+ pow((yield_sum_ratio - e4_sumratio),2)   );

	    ///peter's error done


	    ////=============i===========================

   

	    //   ///////D2 over H2//////Ratio Flipping  check it again, seems wrong


	    //   double yield_sum_ratio_D2overH2 = (y_pos + y_neg)/(yH2_pos + yH2_neg);

	    //   double yield_diff_ratio_D2overH2 =(y_pos - y_neg)/ (yH2_pos - yH2_neg);
	   
	    //   double yield_num_D2overH2_err = sqrt(error_neg*error_neg + error_pos*error_pos);

	    //   double yield_den_D2overH2_err = sqrt(errorH2_neg*errorH2_neg + errorH2_pos*errorH2_pos);

	    //   double sum_ratio_D2overH2_err = yield_sum_ratio_D2overH2 * sqrt(((yield_num_D2overH2_err/(yH2_pos + yH2_neg)) * (yield_num_D2overH2_err/(yH2_pos + yH2_neg))) + ((yield_den_D2overH2_err/(y_pos + y_neg)) * (yield_den_D2overH2_err/(y_pos + y_neg)))   );

	    //   double diff_ratio_D2overH2_err = yield_diff_ratio_D2overH2 * sqrt(((yield_num_D2overH2_err/(yH2_pos - yH2_neg)) * (yield_num_D2overH2_err/(yH2_pos - yH2_neg))) + ((yield_den_D2overH2_err/(y_pos - y_neg)) * (yield_den_D2overH2_err/(y_pos - y_neg)))   );
   





	    // ///////////////////////////////



	    ///add 13 done////
	    double y_RD = (4*y-1)/(1-y);
          

	    double error_RD = 3*error/((1-y)*(1-y));
	    ///================================================added 13+1
	    double yH2_RD = (4*yH2-1)/(1-yH2);
          

	    double errorH2_RD = 3*errorH2/((1-yH2)*(1-yH2));

	    ///================================================added 13+1

            if(y>0 && error< 0.05){//0.1================================>
              //double y_RD = (4*y-1)/(1-y);
              g_yield_ratio->SetPoint(ii,x,y);
              //g_yield_ratio->SetPoint(ii,x,y_RD);
              g_yield_ratio->SetPointError(ii,0,error);
              g_RDmeas_ratio->SetPoint(ii,x,y_RD);
              g_RDmeas_ratio->SetPointError(ii,0,error_RD);
              ii++;
              std::string z_str = std::to_string(x); 
              if(y_RD>0 && error_RD < 1){
                jout[std::to_string(Q2)][std::to_string(xbj)][z_str]["value"].push_back(y_RD);
                jout[std::to_string(Q2)][std::to_string(xbj)][z_str]["error"].push_back(error_RD);
              }
            }

	    //add 13+1============================================

	    if(yH2>0 && errorH2< 0.05){//0.1================================>
              //double y_RD = (4*y-1)/(1-y);
              gh2_yield_ratio->SetPoint(ii,x,yH2);
              //g_yield_ratio->SetPoint(ii,x,y_RD);
              gh2_yield_ratio->SetPointError(ii,0,errorH2);
              gh2_RDmeas_ratio->SetPoint(ii,x,yH2_RD);
              gh2_RDmeas_ratio->SetPointError(ii,0,errorH2_RD);
              ii++;
              // std::string z_str = std::to_string(x); 
              // if(y_RD>0 && error_RD < 1){
              //   jout[std::to_string(Q2)][std::to_string(xbj)][z_str]["value"].push_back(y_RD);
              //   jout[std::to_string(Q2)][std::to_string(xbj)][z_str]["error"].push_back(error_RD);
              // }
            }
	    //add 13+1===============================================
	    
	    //sum diff plot

	    //if(yield_diff_ratio >0 && sum_ratio_err< 0.05 && yield_sum_ratio > 0){
	    //if(yield_diff_ratio >0.2 && yield_diff_ratio < 2 && x > 0.2 && z > 0.3){
	    if(yield_diff_ratio >0. && yield_sum_ratio > 0.){
	      sum_ratio->SetPoint(ii,x,yield_sum_ratio);
	      sum_ratio->SetPointError(ii,0,sum_ratio_err);
	      sum_ratio_sim->SetPoint(ii,x,yield_sum_ratio_sim);
	      sum_ratio_sim->SetPointError(ii,x,yield_sum_ratio_sim_err);

	    
	      // txtsum<<" "<<Q2<<"   "<<z<<"  "<<x<<"  "<<yield_sum_ratio<< "  "<<sum_ratio_err<<endl;
	      txtsum    <<" "<<Q2<<"     "<<x<<"  "<<yield_sum_ratio<< "  "<<sum_ratio_err<<endl;
	      // txtsum_sim<<" "<<Q2<<"     "<<x<<"  "<<yield_sum_ratio_sim<< "  "<<"0.00"<<endl;
	      txtsum_peter <<" "<<Q2<<"     "<<x<<"  "<<yield_sum_ratio<< "  "<<error_sumratio<<endl;

              //z = central z, x = bin value of z
	      diff_ratio->SetPoint(ii,x,yield_diff_ratio);
              diff_ratio->SetPointError(ii,0,diff_ratio_err);
	      diff_ratio_sim->SetPoint(ii,x,yield_diff_ratio_sim);
	      diff_ratio_sim->SetPointError(ii,x,yield_diff_ratio_sim_err);

	      //  txtdiff<<" "<<Q2<<"   "<<z<< "  "<<x<<"  "<<yield_diff_ratio<< "  "<<diff_ratio_err<<endl;
	      txtdiff    <<" "<<Q2<<"    "<<x<<"  "<<yield_diff_ratio<< "  "<<diff_ratio_err<<endl;
              txtdiff_peter    <<" "<<Q2<<"    "<<x<<"  "<<yield_diff_ratio<< "  "<<error_diffratio<<endl;
	     
	      //May3 D2 diff over sum, H2 diff over sum
	      //NOTE THIS x is Z not the xBJ

              txtd2diff_oversum    <<" "<<Q2<<"     "<<x<<"  "<<d2_diff_over_sum << "  "<<d2_ratio_err<<endl;
              txth2diff_oversum    <<" "<<Q2<<"     "<<x<<"  "<<h2_diff_over_sum << "  "<<h2_ratio_err<<endl;


	      //july x means z=================================================================================================================data/simc july22
	      txtd2pos<<" "<<Q2<<"     "<<x<<"  "<<ypos_data_over_simc << "  "<<ypos_data_over_simc_err<<endl;
	      txtd2neg<<" "<<Q2<<"     "<<x<<"  "<<yneg_data_over_simc << "  "<<yneg_data_over_simc_err<<endl;
	      txth2pos<<" "<<Q2<<"     "<<x<<"  "<<yH2pos_data_over_simc << "  "<<yH2pos_data_over_simc_err<<endl;
	      txth2neg<<" "<<Q2<<"     "<<x<<"  "<<yH2neg_data_over_simc << "  "<<yH2neg_data_over_simc_err<<endl;
	      // july 23
	      txtd2pos_hgc<<" "<<Q2<<"     "<<x<<"  "<<y_pos<<"    "<<pos_incrad<<"    "<< d2pos_percent_diff<<"    "<<ypos_data_over_simc << "  "<<ypos_data_over_simc_err<<endl;
	      txtd2neg_hgc<<" "<<Q2<<"     "<<x<<"  "<<y_neg<<"    "<<neg_incrad<<"    "<< d2neg_percent_diff<<"    "<<yneg_data_over_simc << "  "<<yneg_data_over_simc_err<<endl;
	      txth2pos_hgc<<" "<<Q2<<"     "<<x<<"  "<<yH2_pos<<"    "<<posH2_incrad<<"    "<< h2pos_percent_diff<<"    "<<yH2pos_data_over_simc << "  "<<yH2pos_data_over_simc_err<<endl;
	      txth2neg_hgc<<" "<<Q2<<"     "<<x<<"  "<<yH2_neg<<"    "<<negH2_incrad<<"    "<< h2neg_percent_diff<<"    "<<yH2neg_data_over_simc << "  "<<yH2neg_data_over_simc_err<<endl;



	      txtyieldd2pos<<" "<<Q2<<"     "<<x<<"  "<<y_pos<<  "    "<<error_pos<<endl;	    
	      txtyieldd2neg<<" "<<Q2<<"     "<<x<<"  "<<y_neg<<  "    "<<error_neg<<endl;	    
	      txtyieldh2pos<<" "<<Q2<<"     "<<x<<"  "<<yH2_pos<<"    "<<errorH2_pos<<endl;	    
	      txtyieldh2neg<<" "<<Q2<<"     "<<x<<"  "<<yH2_neg<<"    "<<errorH2_neg<<endl;	    

	      txtyieldd2pos_aug15<<" "<<RunGroup<< "   "<<Q2<<"     "<<x<<"  "<<y_pos<<  "    "<<error_pos<<"   "<<pos_incrad<<endl;
	      txtyieldd2neg_aug15<<" "<<RunGroup<< "   "<<Q2<<"     "<<x<<"  "<<y_neg<<  "    "<<error_neg<<"   "<<neg_incrad<<endl;
	      txtyieldh2pos_aug15<<" "<<RunGroup<< "   "<<Q2<<"     "<<x<<"  "<<yH2_pos<<  "    "<<errorH2_pos<<"   "<<posH2_incrad<<endl;
	      txtyieldh2neg_aug15<<" "<<RunGroup<< "   "<<Q2<<"     "<<x<<"  "<<yH2_neg<<  "    "<<errorH2_neg<<"   "<<negH2_incrad<<endl;

	      //oct20, D2-/D2+, H2-/D2+, H2+/D2+ ratios and Errorrs
	      /*
	      txtfile_d2neg_to_d2pos<<" "<<RunGroup<< "   "<<Q2<<"     "<<x<<"  "<<y_d2neg_over_d2pos<<  "    "<<y_d2neg_over_d2pos_err<<endl;
	      txtfile_h2neg_to_d2pos<<" "<<RunGroup<< "   "<<Q2<<"     "<<x<<"  "<<y_h2neg_over_d2pos<<  "    "<<y_h2neg_over_d2pos_err<<endl;
	      txtfile_h2pos_to_d2pos<<" "<<RunGroup<< "   "<<Q2<<"     "<<x<<"  "<<y_h2pos_over_d2pos<<  "    "<<y_h2pos_over_d2pos_err<<endl;
	      */

	      if(y_d2neg_over_d2pos>0){
		txtfile_d2neg_to_d2pos<<"    "<<Q2<<"     "<<x<<"  "<<y_d2neg_over_d2pos<<  "    "<<y_d2neg_over_d2pos_err<<endl;
	      }
	      if(y_h2neg_over_d2pos>0){
		txtfile_h2neg_to_d2pos<<"    "<<Q2<<"     "<<x<<"  "<<y_h2neg_over_d2pos<<  "    "<<y_h2neg_over_d2pos_err<<endl;
	      }
	      if(y_h2pos_over_d2pos>0){
		txtfile_h2pos_to_d2pos<<"    "<<Q2<<"     "<<x<<"  "<<y_h2pos_over_d2pos<<  "    "<<y_h2pos_over_d2pos_err<<endl;
	      }
	   
	      if(yH2_pos>0 || y_pos>0 || yH2_neg>0 || y_neg){
	      txtfile_all<<setprecision(4)<<fixed<<RunGroup<< "   "<<Q2<<"    "<<x<<"  "<<yH2_pos<<"  "<<errorH2_pos<<"  "<<y_pos<<"  "<<error_pos<<"  "<<yH2_neg<<"   "<<errorH2_neg<<"   "<<y_neg<<"   "<<error_neg<<endl;
	      }


	      //oct 20 done


	      // txtdiff_sim<<" "<<Q2<<"     "<<x<<"  "<<yield_diff_ratio_sim<< "  "<<"0.00"<<endl;

	      // txtsum1<<" "<<Q2<<"   "<<z<< "  "<<x<<"  "<<sum_ratio_err/yield_sum_ratio<< "  "<<yield_num_err/(y_pos + y_neg)<<"   "<<yield_den_err/(yH2_pos + yH2_neg)<<endl;
 
	      //    txtdiff1<<" "<<Q2<<"   "<<z<< "  "<<x<<"  "<<diff_ratio_err/yield_diff_ratio<< "  "<<yield_num_err/(y_pos - y_neg)<<"   "<<yield_den_err/(yH2_pos - yH2_neg)<<endl;

	      //Error in sum ratio/sum ratio
	      //error in numerator/numerator
	      //error in denominator/denominator 

              ii++;
             
	    }//coomm Now
	    //add 13+1===============================================

	  }//nbins

	  //
	 

	  int nbins_sim = h_z_neg_sim_incnorad->GetXaxis()->GetNbins();
	  TGraphErrors* g_yield_ratio_sim = new TGraphErrors();
	  std::string z_string_sim = "z simc "+(std::to_string(z)).substr(0,4);
	  g_yield_ratio_sim->SetName(z_string_sim.c_str());
	  //TGraphErrors* g_yield_ratio_sim = new TGraphErrors(h_z_neg_sim_incnorad);
	  int ii_sim = 0;
	  for(int i = 0;i<nbins_sim;i++){
	    //std::cout<<i<<std::endl;
	    double x = h_z_neg_sim_incnorad->GetBinCenter(i+1);//+0.001*i_color;
	    double y = h_z_neg_sim_incnorad->GetBinContent(i+1);
	    double error = h_z_neg_sim_incnorad->GetBinError(i+1);
	    //std::cout<<i<<" x "<<x<<" y "<<y<<std::endl;
	    if(y!=0){
	      //      double y_RD = (4*y-1)/(1-y);
	      //       g_yield_ratio_sim->SetPoint(ii_sim,x,y_RD);
	      g_yield_ratio_sim->SetPoint(ii_sim,x,y);
	      //g_yield_ratio_sim->SetPointError(ii_sim,0,error);
	      ii_sim++;
	    }//y!=0
	  }//simc bin

	  g_yield_ratio->SetName(q2xz_str_filename.c_str());


	   

	  //c_Q2x_ratio->cd();
	  //h_z_neg_all->SetLineColor(i_color);
	  g_yield_ratio->SetMarkerStyle(4);
	  g_yield_ratio->SetMarkerColor(i_color);
	  g_yield_ratio->SetLineColor(i_color);
	  g_yield_ratio_sim->SetMarkerStyle(5);
	  g_yield_ratio_sim->SetMarkerColorAlpha(i_color,0.35);
	  g_yield_ratio_sim->SetLineColor(i_color);
	  //hs->Add(h_z_neg_all);
	  //h_z_neg_all->Draw("same");
	  mg->Add(g_yield_ratio,"P");
	  mg->Add(g_yield_ratio_sim,"L");
	  g_RDmeas_ratio->SetMarkerStyle(4);
	  g_RDmeas_ratio->SetMarkerColor(i_color);
	  g_RDmeas_ratio->SetLineColor(i_color);
         
	  mg_RD->Add(g_RDmeas_ratio,"P");
	  //c_Q2x_ratio->Update();


	  TCanvas *c_Q2x_z_ratio = new TCanvas(q2x_name.c_str(),q2x_name.c_str(),1900,1000);
	  //TCanvas *c_Q2x_ratio = new TCanvas("",q2x_name.c_str(),1900,1000);
	  //TCanvas *c_Q2x_z_ratio = new TCanvas();
	  // //h_z_neg_all->Draw();
	  g_yield_ratio->GetXaxis()->SetRangeUser(0.1,1);
	  //g_yield_ratio->GetYaxis()->SetRangeUser(0.1,1.2);
	  g_yield_ratio->GetXaxis()->SetTitle("z");
	  g_yield_ratio->GetYaxis()->SetTitle("yield_ratio");
	  g_yield_ratio->GetXaxis()->SetTitleSize(0.053);
	  g_yield_ratio->GetYaxis()->SetTitleSize(0.053);
	  g_yield_ratio->Draw("AP"); 
	  g_yield_ratio_sim->Draw("L");
	  c_Q2x_z_ratio->BuildLegend(0.1,0.1,0.5,0.2,q2xz_str.c_str());
	  std::string zratiopdfname = "/u/group/c-csv/hdbhatt/yielddec21/yieldJan_2022/results/yield/statistics_corr/"+q2xz_str_filename+"_ratio_hgcer.pdf";

	  c_Q2x_z_ratio->SaveAs(zratiopdfname.c_str());

	  TGraphErrors* g_yield_ratio_all = new TGraphErrors();
	  std::string z_string_all = "z setting "+(std::to_string(z)).substr(0,4);
	  g_yield_ratio_all->SetName(z_string_all.c_str());

	  // int ii_all = 0;
	  //for(int i = 0;i<nbins;i++){
	  //  //std::cout<<i<<std::endl;
	  //  double x = h_z_neg_all->GetBinCenter(i)+0.001*i_color+0.6*i_which_x[i_whichx];
	  //  double y = h_z_neg_all->GetBinContent(i)+i_which_Q2[i_whichx];
	  //  double error = h_z_neg_all->GetBinError(i);
	  //  //std::cout<<i<<" x "<<x<<" y "<<y<<std::endl;
	  //  if(y!=0 && y!= 1 && y!=2 && error< 0.2){
	  //    g_yield_ratio_all->SetPoint(ii_all,x,y);
	  //    g_yield_ratio_all->SetPointError(ii_all,0,error);
	  //    ii_all++;
	  //  }
	  //}
	  //  g_yield_ratio_all->SetMarkerStyle(4);
	  //  g_yield_ratio_all->SetMarkerColor(i_color);
	  //  g_yield_ratio_all->SetLineColor(i_color);
	  //  mg_all->Add(g_yield_ratio_all,"P");
	  i_color++;
	}//loop over z
	//	}//===================================

	//plotting outside of z loop so that we have RD vs z plot



	TCanvas *cc = new TCanvas("cc","",1200,800);
	cc->Divide(1,1);
	cc->cd(1);gPad->SetGrid();
	sum_ratio->Draw("AP");
	sum_ratio->SetMarkerStyle(20);;
	sum_ratio->SetMarkerSize(2);
	sum_ratio->GetXaxis()->SetRangeUser(0.,1.0);
	sum_ratio->GetYaxis()->SetRangeUser(1.0,4.0);
	sum_ratio->GetXaxis()->SetTitle("z");
	sum_ratio->GetYaxis()->SetTitle("(YH2+ + YH2-)/(YD2+ + YD2-)");
	sum_ratio->GetXaxis()->CenterTitle();
	sum_ratio->GetYaxis()->CenterTitle();
	//	sum_ratio->Draw("AP");
	gPad->Update();
	TCanvas *dd = new TCanvas("dd","",1200,800);
	dd->Divide(1,1);
	dd->cd(1);gPad->SetGrid();
	diff_ratio->Draw("AP");
	diff_ratio->SetMarkerStyle(20);;
	diff_ratio->SetMarkerSize(2);
	diff_ratio->GetXaxis()->SetRangeUser(0.,1.0);
	diff_ratio->GetYaxis()->SetRangeUser(1.0,4.0);
	diff_ratio->GetXaxis()->SetTitle("z");
	diff_ratio->GetYaxis()->SetTitle("(YH2+ - YH2-)/(YD2+ - YD2-)");
	diff_ratio->GetXaxis()->CenterTitle();
	diff_ratio->GetYaxis()->CenterTitle();
	//diff_ratio->Draw("AP");

	TCanvas *ee = new TCanvas("ee","",1200,800);
	ee->Divide(1,1);
	ee->cd(1);gPad->SetGrid();
	sum_ratio_sim->Draw("AP");
	sum_ratio_sim->SetMarkerStyle(20);;
	sum_ratio_sim->SetMarkerSize(2);
	sum_ratio_sim->GetXaxis()->SetRangeUser(0.,1.0);
	sum_ratio_sim->GetYaxis()->SetRangeUser(1.0,4.0);
	sum_ratio_sim->GetXaxis()->SetTitle("z");
	sum_ratio_sim->GetYaxis()->SetTitle("(YH2+ + YH2-)/(YD2+ + YD2-)");
	sum_ratio_sim->GetXaxis()->CenterTitle();
	sum_ratio_sim->GetYaxis()->CenterTitle();
	//	sum_ratio_sim->Draw("AP");
	gPad->Update();
	TCanvas *ff = new TCanvas("ff","",1200,800);
	ff->Divide(1,1);
	ff->cd(1);gPad->SetGrid();
	diff_ratio_sim->Draw("AP");
	diff_ratio_sim->SetMarkerStyle(20);;
	diff_ratio_sim->SetMarkerSize(2);
	diff_ratio_sim->GetXaxis()->SetRangeUser(0.,1.0);
	diff_ratio_sim->GetYaxis()->SetRangeUser(1.0,4.0);
	diff_ratio_sim->GetXaxis()->SetTitle("z");
	diff_ratio_sim->GetYaxis()->SetTitle("(YH2+ - YH2-)/(YD2+ - YD2-)");
	diff_ratio_sim->GetXaxis()->CenterTitle();
	diff_ratio_sim->GetYaxis()->CenterTitle();
	//diff_ratio_sim->Draw("AP");

	gPad->Update();
	//cc->SaveAs("sum_diff_hgcer.pdf");
	std::string sumname = "/u/group/c-csv/hdbhatt/yielddec21/yieldJan_2022/results/yield/statistics_corr/"+canvas_filename+"_sum_ratio_hgcer.pdf";
	//c_Q2x_ratio->BuildLegend(0.6,0.6,0.95,0.95);
	cc->SaveAs(sumname.c_str());
	std::string diffname = "/u/group/c-csv/hdbhatt/yielddec21/yieldJan_2022/results/yield/statistics_corr/"+canvas_filename+"_diff_ratio_hgcer.pdf";
	//c_Q2x_ratio->BuildLegend(0.6,0.6,0.95,0.95);
	dd->SaveAs(diffname.c_str());
	std::string sumname_sim = "/u/group/c-csv/hdbhatt/yielddec21/yieldJan_2022/results/yield/statistics_corr/"+canvas_filename+"_sum_ratio_sim_hgcer.pdf";
	//c_Q2x_ratio_sim->BuildLegend(0.6,0.6,0.95,0.95);
	ee->SaveAs(sumname_sim.c_str());
	std::string diffname_sim = "/u/group/c-csv/hdbhatt/yielddec21/yieldJan_2022/results/yield/statistics_corr/"+canvas_filename+"_diff_ratio_sim_hgcer.pdf";
	//c_Q2x_ratio_sim->BuildLegend(0.6,0.6,0.95,0.95);
	ff->SaveAs(diffname_sim.c_str());



       
	//i_whichx++;
	TCanvas *c_Q2x_ratio = new TCanvas("",q2x_name.c_str(),1900,1000);
	c_Q2x_ratio->cd();
	c_Q2x_ratio->SetTitle(canvas_name.c_str());
	//hs->Draw();
	//TF1 *f_RD = new TF1("f1","2.5/x-2.5",0,1);
	mg->GetXaxis()->SetTitleSize(0.053);
	mg->GetYaxis()->SetTitleSize(0.053);
	mg->GetXaxis()->SetLabelSize(0.05);
	mg->GetYaxis()->SetLabelSize(0.05);
	mg->SetMinimum(0.3);
	mg->SetMaximum(0.9);
	std::string mg_title = canvas_name+",    Zhadron";
	mg->SetTitle(mg_title.c_str());
	//TPaveText *pt_mg = new TPaveText(0.5,0.8,0.8,1);
	//pt_mg->AddText(canvas_name.c_str());
	//pt_mg->Draw("same");
	mg->Draw("APL");
	//f_RD->Draw("same");
	mg->GetHistogram()->SetTitle(canvas_name.c_str());
	mg->GetXaxis()->SetTitle(mg_title.c_str());
	mg->GetYaxis()->SetTitle("Yield Ratio (pi-/pi+)");
	mg->GetXaxis()->SetLimits(0.3,0.8);
	//  mg->GetYaxis()->SetLimits(0.,1.0);

	//auto hermes_RD = [](double z){return ((1.0-z)*0.083583)/((1.0+z)*1.988);};
	//TF1 *fit = new TF1("HERMES","(1.0-x)**0.083583/(1.0+x)**1.9838",0,1);
	//fit->Draw("same");
	//      std::string ratiopdfname = "results/yield/statistics_corr/"+canvas_filename+"_RDratio_hgcer.pdf";
	std::string ratiopdfname = "/u/group/c-csv/hdbhatt/yielddec21/yieldJan_2022/results/yield/statistics_corr/"+canvas_filename+"_ratio_hgcer.pdf";
	//c_Q2x_ratio->BuildLegend(0.6,0.6,0.95,0.95);
	c_Q2x_ratio->SaveAs(ratiopdfname.c_str());
	TCanvas* c_RD_ratio = new TCanvas();
	mg_RD->GetXaxis()->SetTitleSize(0.053);
	mg_RD->GetYaxis()->SetTitleSize(0.053);
	mg_RD->GetXaxis()->SetLabelSize(0.05);
	mg_RD->GetYaxis()->SetLabelSize(0.05);
	mg_RD->SetMinimum(0);
	mg_RD->SetMaximum(10);
	std::string mg_RD_title = canvas_name+",z";
	mg_RD->SetTitle(mg_RD_title.c_str());
	TPaveText *pt_mg_RD = new TPaveText(0.5,0.8,0.8,1);
	pt_mg_RD->AddText(canvas_name.c_str());
	pt_mg_RD->Draw("same");
	mg_RD->Draw("APL");
	mg_RD->GetHistogram()->SetTitle(canvas_name.c_str());
	mg_RD->GetXaxis()->SetTitle(mg_RD_title.c_str());
	mg_RD->GetYaxis()->SetTitle("RD_Measured");

	mg_RD->GetXaxis()->SetLimits(0.3,0.8);

	std::string RDpdfname = "/u/group/c-csv/hdbhatt/yielddec21/yieldJan_2022/results/yield/statistics_corr/"+canvas_filename+"_RDratio_hgcer.pdf";
	c_RD_ratio->BuildLegend(0.7,0.7,0.9,0.9);
	c_RD_ratio->SaveAs(RDpdfname.c_str());
	//	}//H2!empty
      }//if x,Q2 not 0

      //i_whichq2++;
      //      std::cout<<i_whichx<<std::endl;
    }//loop over Q2
  }//loop over x

  //std::string jout_name = "results/yield_ratio.json";
  //std::ofstream ofs_jout(jout_name.c_str());
  //ofs_jout<<jout.dump(4)<<std::endl;
  //TCanvas* c_yield_ratio_all = new TCanvas();
  //mg_all->GetXaxis()->SetTitle("z");
  //mg_all->GetYaxis()->SetTitle("R^D_meas");
  //mg_all->GetXaxis()->SetTitleSize(0.053);
  //mg_all->GetYaxis()->SetTitleSize(0.053);
  //mg_all->GetXaxis()->SetLabelSize(0.05);
  //mg_all->GetYaxis()->SetLabelSize(0.05);
  //mg_all->Draw("A");
  //std::string c_yield_ratio_all_name = "results/yield/yield_ratio_all_hgcer.pdf";
  //c_yield_ratio_all->SaveAs(c_yield_ratio_all_name.c_str());
  return 0;

}
