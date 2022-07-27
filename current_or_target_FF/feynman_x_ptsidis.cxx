#define NRUNS 700
void feynman_x_ptsidis(){

  //  ofstream txtfile("textfiles/feynman_x.txt");
  ofstream txt("feynmanx_ptsidis_new.txt",ios::app);

  for(Int_t r = 0; r<NRUNS;r++){
    Int_t runs[NRUNS]     = {5401,5402,5403,5404,5405,5406,5407,5408,5409,5410,5411,5412,5413,5414,5415,5416,5417,5418,5419,5420,5421,5422,5423,5424,5425,5426,5427,5428,5429,5430,5431,5432,5433,5434,5435,5436,5437,5438,5439,5440,5441,5442,5443,5444,5445,5446,5447,5448,5449,5450,5451,5452,5453,5454,5455,5456,5457,5458,5459,5460,5461,5462,5463,5464,5465,5466,5467,5468,5469,5470,5471,5472,5473,5474,5475,5476,5477,5478,5479,5480,5481,5482,5483,5484,5485,5486,5487,5488,5489,5490,5491,5492,5493,5494,5495,5496,5497,5498,5499,5500,5501,5502,5503,5504,5505,5506,5507,5508,5509,5510,5511,5512,5513,5514,5515,5516,5517,5518,5519,5520,5521,5522,5523,5524,5525,5526,5527,5528,5529,5530,5531,5532,5533,5534,5535,5536,5537,5538,5539,5540,5541,5542,5543,5544,5545,5546,5547,5548,5549,5550,5551,5552,5553,5554,5555,5556,5557,5558,5559,5560,5561,5562,5563,5564,5565,5566,5567,5568,5569,5570,5571,5572,5573,5574,5575,5576,5577,5578,5579,5580,5581,5582,5583,5584,5585,5586,5587,5588,5589,5590,5591,5592,5593,5594,5595,5596,5597,5598,5599,5600,5601,5602,5603,5604,5605,5606,5607,5608,5609,5610,5611,5612,5613,5614,5615,5616,5617,5618,5619,5620,5621,5622,5623,5624,5625,5626,5627,5628,5629,5630,5631,5632,5633,5634,5635,5636,5637,5638,5639,5640,5641,5642,5643,5644,5645,5646,5647,5648,5649,5650,5651,5652,5653,5654,5655,5656,5657,5658,5659,5660,5661,5662,5663,5664,5665,5666,5667,5668,5669,5670,5671,5672,5673,5674,5675,5676,5677,5678,5679,5680,5681,5682,5683,5684,5685,5686,5687,5688,5689,5690,5691,5692,5693,5694,5695,5696,5697,5698,5699,5700,5701,5702,5703,5704,5705,5706,5707,5708,5709,5710,5711,5712,5713,5714,5715,5716,5717,5718,5719,5720,5721,5722,5723,5724,5725,5726,5727,5728,5729,5730,5731,5732,5733,5734,5735,5736,5737,5738,5739,5740,5741,5742,5743,5744,5745,5746,5747,5748,5749,5750,5751,5752,5753,5754,5755,5756,5757,5758,5759,5760,5761,5762,5763,5764,5765,5766,5767,5768,5769,5770,5771,5772,5773,5774,5775,5776,5777,5778,5779,5780,5781,5782,5783,5784,5785,5786,5787,5788,5789,5790,5791,5792,5793,5794,5795,5796,5797,5798,5799,5800,5801,5802,5803,5804,5805,5806,5807,5808,5809,5810,5811,5812,5813,5814,5815,5816,5817,5818,5819,5820,5821,5822,5823,5824,5825,5826,5827,5828,5829,5830,5831,5832,5833,5834,5835,5836,5837,5838,5839,5840,5841,5842,5843,5844,5845,5846,5847,5848,5849,5850,5851,5852,5853,5854,5855,5856,5857,5858,5859,5860,5861,5862,5863,5864,5865,5866,5867,5868,5869,5870,5871,5872,5873,5874,5875,5876,5877,5878,5879,5880,5881,5882,5883,5884,5885,5886,5887,5888,5889,5890,5891,5892,5893,5894,5895,5896,5897,5898,5899,5900,5901,5902,5903,5904,5905,5906,5907,5908,5909,5910,5911,5912,5913,5914,5915,5916,5917,5918,5919,5920,5921,5922,5923,5924,5925,5926,5927,5928,5929,5930,5931,5932,5933,5934,5935,5936,5937,5938,5939,5940,5941,5942,5943,5944,5945,5946,5947,5948,5949,5950,5951,5952,5953,5954,5955,5956,5957,5958,5959,5960,5961,5962,5963,5964,5965,5966,5967,5968,5969,5970,5971,5972,5973,5974,5975,5976,5977,5978,5979,5980,5981,5982,5983,5984,5985,5986,5987,5988,5989,5990,5991,5992,5993,5994,5995,5996,5997,5998,5999,6000,6001,6002,6003,6004,6005,6006,6007,6008};//

    TString filename  = Form("/lustre19/expphy/volatile/hallc/c-csv/hdbhatt/ROOTfiles/skimROOTfilesfeb11/run_%d.root", runs[r]);
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


    ////hists

    double pt,pt2,inv_mass,fyn_x,Epi,zhad;
  
    Double_t pionmass   = 0.1395701835; 
    Double_t protonmass = 0.93827231;//GeV/c2
    Double_t pi = 3.1416;
    TTree * tt = (TTree *)f->Get("T");
    Long64_t data_entries = tt->GetEntries();

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
    tt->SetBranchAddress("H_gtr_dp",&pdelta);
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
    Double_t track_x_fp;
    tt->SetBranchAddress("P_dc_x_fp", &track_x_fp);
    Double_t track_y_fp;
    tt->SetBranchAddress("P_dc_y_fp", &track_y_fp);
    Double_t track_xp_fp;
    tt->SetBranchAddress("P_dc_xp_fp", &track_xp_fp);
    Double_t track_yp_fp;
    tt->SetBranchAddress("P_dc_yp_fp", &track_yp_fp);
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

      if(pdelta>-10&&pdelta<20&&hdelta>-8&&hdelta<8&&gevtyp>3){
	Epi = sqrt(pow(pionmass,2) + pow(mom,2));
	zhad = Epi/dnu;

	pt=mom*sin(theta);
	pt2=pt*pt;
	inv_mass=sqrt(dW2);
	fyn_x=2*sqrt(fabs(mom*mom-pt2))/inv_mass;

	hfyn_x->Fill(fyn_x);
	hx->Fill(xbj);
	hQ2->Fill(Q2);
	hmsmom->Fill(hmom);
	pmom->Fill(mom);
	hW2->Fill(dW2);
	hpt->Fill(pt);
	hpt2->Fill(pt2);
	hz->Fill(zhad);


      }//pid

    }//entries



    TCanvas *c1 = new TCanvas("c1","c1", 1200,800);
    c1->Divide(3,3);
    c1->cd(1);
    hx->Draw();
    hx->SetLineWidth(2);
    
    c1->cd(2);
    hmsmom->Draw();
    hmsmom->SetLineWidth(2);
    
    c1->cd(3);
    pmom->Draw();
    pmom->SetLineWidth(2);

    c1->cd(4);
    hW2->Draw();
    hW2->SetLineWidth(2);

    c1->cd(5);
    hpt->Draw();
    hpt->SetLineWidth(2);

    c1->cd(6);
    hpt2->Draw();
    hpt2->SetLineWidth(2);

    c1->cd(7);
    hz->Draw();
    hz->SetLineWidth(2);

    c1->cd(8);
    hfyn_x->Draw();
    hfyn_x->SetLineWidth(2);

 c1->cd(9);
    hQ2->Draw();
    hQ2->SetLineWidth(2);


    double  xav = hx->GetMean();
    double  zav = hz->GetMean();
    double  hmsmomav = hmsmom->GetMean();
    double  shmsmomav = pmom->GetMean();
    double  w2av = hW2->GetMean();
    double  ptav = hpt->GetMean();
    double  pt2av = hpt2->GetMean();
    double  feynxav = hfyn_x->GetMean();
    double  Q2av = hQ2->GetMean();

    txt<<runs[r]<<setprecision(3)<<fixed<<"\t\t"<<xav<<"\t\t"<<zav<<"\t\t"<<Q2av<<"\t\t"<<hmsmomav<<"\t\t"<<shmsmomav<<"\t\t"<<w2av<<"\t\t"<<ptav<<"\t\t"<<pt2av<<"\t\t"<<feynxav<<endl;

    c1->SaveAs(Form("pdf/feynman_hists_ptsidis_%d.pdf", runs[r]));

  }//runs
}//void
//Run		xav		zav		Q2av		hmsp		shmsp		W2			        pt               pt2       feynmanx
