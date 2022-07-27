#define NRUNS 700
void feynman_x(){

  //  ofstream txtfile("textfiles/feynman_x.txt");
  ofstream txt("feynmanx_csv.txt",ios::app);

  for(Int_t r = 0; r<NRUNS;r++){
    Int_t runs[NRUNS]     = { 6020, 6021, 6022, 6026, 6027, 6028, 6029, 6030, 6031, 6032, 6043, 6044, 6045, 6046, 6049, 6050, 6051, 6052, 6053, 6054, 6055, 6056, 6057, 6058, 6067, 6068, 6070, 6071, 6073, 6080, 6081, 6087, 6088, 6089, 6090, 6091, 6093, 6094, 6095, 6096, 6097, 6098, 6109, 6110, 6111, 6112, 6114, 6115, 6116, 6119, 6120, 6124, 6125, 6127, 6128, 6129, 6130, 6132, 6133, 6135, 6136, 6137, 6138, 6145, 6146, 6154, 6155, 6166, 6168, 6169, 6170, 6171, 6172, 6179, 6180, 6181, 6182, 6183, 6184, 6185, 6186, 6187, 6194, 6195, 6196, 6197, 6198, 6199, 6205, 6208, 6209, 6210, 6211, 6212, 6213, 6214, 6215, 6216, 6217, 6219, 6220, 6242, 6243, 6244, 6245, 6246, 6247, 6248, 6249, 6250, 6251, 6252, 6261, 6262, 6263, 6264, 6265, 6266, 6267, 6270, 6271, 6272, 6273, 6274, 6275, 6276, 6288, 6289, 6290, 6291, 6292, 6293, 6294, 6295, 6296, 6304, 6305, 6306, 6307, 6308, 6309, 6310, 6311, 6312, 6317, 6318, 6322, 6323, 6325, 6326, 6327, 6328, 6329, 6339, 6340, 6341, 6342, 6343, 6344, 6345, 6346, 6347, 6348, 6349, 6350, 6351, 6352, 6353, 6354, 6359, 6360, 6363, 6364, 6365, 6367, 6368, 6370, 6371, 6372, 6373, 6374, 6375, 6376, 6377, 6378, 6379, 6380, 6381, 6382, 6383, 6385, 6386, 6387, 6388, 6389, 6390, 6391, 6393, 6394, 6395, 6396, 6397, 6398, 6399, 6400, 6401, 6402, 6403, 6404, 6405, 6406, 6407, 6408, 6410, 6411, 6412, 6413, 6415, 6416, 6417, 6419, 6421, 6422, 6423, 6425, 6426, 6427, 6428, 6429, 6430, 6431, 6432, 6433, 6434, 6435, 6436, 6439, 6440, 6441, 6442, 6443, 6444, 6445, 6446, 6447, 6448, 6451, 6452, 6453, 6454, 6455, 6456, 6459, 6460, 6461, 6462, 6463, 6464, 6465, 6466, 6467, 6468, 6469, 6470, 6471, 6472, 6473, 6474, 6475, 6476, 6477, 6478, 6479,  6486, 6487, 6488, 6489, 6490, 6491, 6493, 6494, 6495, 6496, 6497, 6498, 6499, 6500, 6501, 6502, 6503, 6504, 6506, 6507, 6509, 6510, 6512, 6513, 6514, 6515, 6516, 6517, 6518, 6519, 6520, 6521, 6522, 6523, 6524, 6525, 6526, 6527, 6528, 6529, 6531, 6532, 6533, 6534, 6535, 6536, 6537, 6538, 6539, 6541, 6542, 6543, 6544, 6545, 6546, 6547, 6548, 6549, 6550, 6551, 6552, 7593, 7594, 7595, 7596, 7597, 7599, 7600, 7601, 7602, 7603, 7604, 7605, 7606, 7607, 7608, 7609, 7610, 7611, 7612, 7613, 7614, 7615, 7616, 7617, 7618, 7619, 7620, 7621, 7622, 7623, 7625, 7626, 7627, 7628, 7629, 7630, 7631, 7632, 7633, 7634, 7635, 7636, 7637, 7638, 7639, 7640, 7641, 7642, 7643, 7644, 7645, 7646, 7647, 7648, 7649, 7650, 7651, 7652, 7654, 7655, 7657, 7658, 7659, 7660, 7661, 7662, 7663, 7666, 7667, 7668, 7669, 7670, 7671, 7672, 7673, 7674, 7675, 7676, 7677, 7678, 7679, 7680, 7681, 7682, 7684, 7685, 7686, 7687, 7688, 7689, 7690, 7692, 7694, 7695, 7697, 7698, 7699, 7702, 7703, 7706, 7707, 7708, 7709, 7710, 7711, 7712, 7713, 7714, 7715, 7716, 7717, 7718, 7719, 7720, 7721, 7722, 7724, 7725, 7727, 7728, 7735, 7736, 7737, 7738, 7739, 7740, 7747, 7748, 7749, 7750, 7751, 7752, 7757, 7758, 7759, 7760, 7762, 7768, 7769, 7770, 7771, 7772, 7773, 7774, 7775, 7776, 7777, 7783, 7784, 7785, 7786, 7787, 7790, 7791, 7792, 7793, 7794, 7796, 7797, 7798, 7799, 7800, 7801, 7802, 7803, 7804, 7805, 7806, 7807, 7808, 7809, 7811, 7812, 7813, 7814, 7815, 7816, 7817, 7819, 7820 };//{  6021, 6112};

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


    c1->SaveAs(Form("pdf/feynman_hists_%d.pdf", runs[r]));

  }//runs
}//void
