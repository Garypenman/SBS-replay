//#include "gmn_tree.C"
#include "skim_tree.C"
#include "TChain.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TFile.h"
#include <iostream>
#include <fstream>
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TCut.h"
#include "TTreeFormula.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompSVD.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "/home/gpenman/Analysis/SBS-Analysis/GEn/include/gen_ana.h"

//const double Mp=0.938272;

void TOFcal(const char *outputfilename="TOFcal_out.root"){
  gStyle->SetPalette(kRainBow);
  gROOT->SetBatch(kTRUE);
  
  TString currentline;
  
  TChain *C = new TChain("T");
  C->Add("/volatile/halla/sbs/gpenman/data/out_tof_gen2_H2.root");
  
  //Default reference PMTs for aligning all other PMTs:
  int hodorefID = 44; // reference BAR (not PMT); we'll set the LEFT PMT to the reference by convention:
  int hcalrefID = 199; //row 16 column 6
  
  double zhodo = 1.854454; //meters
  double Lbar_hodo = 0.6; //meters
  double etof0 = (1.96+3.0)/0.299792458; // "central" TOF value:
  
  double thetaHCAL = 34.7;
  double dHCAL = 17.0;
  double dSBS = 2.8;
  
  double W2min=0.0, W2max=1.7;
  
  double sbsmaxfield = 1.3; //1.26?
  double sbsfield = 1.0;

  const double Dgap = 48.0*2.54/100.0;
  
  double Ebeam = 4.291;
  
  thetaHCAL *= TMath::Pi()/180.0;
  
  TVector3 zaxis_HCAL(-sin(thetaHCAL),0,cos(thetaHCAL));
  TVector3 xaxis_HCAL(0,-1,0);
  TVector3 yaxis_HCAL = (zaxis_HCAL.Cross(xaxis_HCAL)).Unit();
  
  TVector3 HCALorigin = dHCAL*zaxis_HCAL;
  
  //180 walk correction slopes, 90 propagation speed slopes, 90 paddle rf offsets and 90 zero offsets.
  int nparams_hodo = 180 + 180 + 90;
  //288 zero offsets, 288 rf offsets, 2*288 TW params for hall B like fit
  int nparams_HCAL = 288*4;
  
  //read in previously generated config for iterative calibration
  vector<double> hodoparams, hcalparams, hcalparams_test;
  vector<double> HODOt0, HODOwL, HODOwR, vscint, tpad_rf; 
  vector<double> HCALt0, HCALw0, HCALw1, HCALw2, tpad_rf_hcal;
  //vector<double> HCALt0_test, HCALw0_test, HCALw1_test, HCALw2_test;
  TString hodoconfigfile = "hodoparams.txt";
  TString hcalconfigfile = "hcalparams.txt";
  //TString hcalconfigfile = "hcalparams_good.txt";
  
  //get hodo params if they exist otherwise fill vector of zeros
  ifstream hodoconfigin(hodoconfigfile);
  if(!hodoconfigin){
    cout << "Hodo config file not found, assuming first run of calibration." << endl;
    for (int i=0; i<nparams_hodo; i++)
      hodoparams.push_back(0.0);
  }
  else{
    while(currentline.ReadLine(hodoconfigin)){
      hodoparams.push_back(stod(currentline.Data()));
    }
  }
  hodoconfigin.close();
  
  //get hcal params if they exist otherwise fill vector of zeros
  //HCAL and hodo TDCs both use BigBite trigger as reference time, so both have reference time subtracted:
  //input the known good tw params from offline fit
  ifstream hcalconfigin(hcalconfigfile);
  if(!hcalconfigin){
    cout << "Hcal config file not found, assuming first run of calibration." << endl;
    for (int i=0; i<nparams_HCAL; i++){
      hcalparams.push_back(0.0);
    }
  }else{
    while(currentline.ReadLine(hcalconfigin)){
      hcalparams.push_back(stod(currentline.Data()));
    }
  }
  hcalconfigin.close();

  //now fill the vectors associated with each set of parameters
  for( int i=0; i<90; i++ ){
    HODOt0.push_back( hodoparams[i] );
    vscint.push_back( 1.0/hodoparams[i+90] );
    HODOwL.push_back( hodoparams[i+180] );
    HODOwR.push_back( hodoparams[i+270] );
    tpad_rf.push_back( hodoparams[i+360] );
    }
  
  for (int i=0; i<288; i++ ){
    HCALt0.push_back( hcalparams[i] );
    HCALw0.push_back( hcalparams[i+288] );
    HCALw1.push_back( hcalparams[i+2*288] );
    tpad_rf_hcal.push_back( hcalparams[i+4*288] );
  }
  
  //Ordering of params:
  // t0i, walkileft, walkiright, 1/vi
  // later, rfi
  TMatrixD Mhodo(nparams_hodo,nparams_hodo);
  TMatrixD Mhcal(nparams_HCAL,nparams_HCAL);
  TVectorD bhodo(nparams_hodo);
  TVectorD bhcal(nparams_HCAL);

  for( int ipar=0; ipar<nparams_hodo; ipar++ ){
    for( int jpar=0; jpar<nparams_hodo; jpar++ ){
      Mhodo(ipar,jpar) = 0.0;
    }
    bhodo(ipar) = 0.0;
  }
  
  for( int ipar=0; ipar<nparams_HCAL; ipar++ ){
    for( int jpar=0; jpar<nparams_HCAL; jpar++ ){
      Mhcal(ipar,jpar) = 0.0;
    }
    bhcal(ipar) = 0.0;
  }
  
  
  //declare outfile here incase we want to choose to abort execution later on due to missing input files
  //but not overwrite previously generated outputfiles etc
  TFile *fout = new TFile(outputfilename,"RECREATE");
  fout->cd();
  

  //declare histograms here
  int nbins=250;
  TH2D *htmean_hodoID_old = new TH2D("htmean_hodoID_old","OLD ;bar ID; bar t_{mean} [ns]", 90,-0.5,89.5, nbins,-30,30);
  TH2D *htmean_hodoID_raw = new TH2D("htmean_hodoID_raw","RAW ;bar ID; bar t_{mean} [ns]", 90,-0.5,89.5, nbins,-30,30);
  TH2D *htmean_hodoID_new = new TH2D("htmean_hodoID_new","NEW ;bar ID; bar t_{mean} [ns]", 90,-0.5,89.5, nbins,-30,30);

  TH2D *htdiff_yhodo_old = new TH2D("htdiff_yhodo_old", "OLD; y_{hodo} (m); bar t_{L}-t_{R} [ns]", nbins,-0.3,0.3,nbins,-15,15);
  TH2D *htdiff_yhodo_raw = new TH2D("htdiff_yhodo_raw", "RAW; y_{hodo} (m); bar t_{L}-t_{R} [ns]", nbins,-0.3,0.3,nbins,-15,15);
  TH2D *htdiff_yhodo_new = new TH2D("htdiff_yhodo_new", "NEW; y_{hodo} (m); bar t_{L}-t_{R} [ns]", nbins,-0.3,0.3,nbins,-15,15);

  TH2D *htmean_hcalID_old = new TH2D("htmean_hcalID_old","OLD ;hcal ID; clus t_{mean} [ns]", 288,-0.5,287.5,nbins,-30,30);
  TH2D *htmean_hcalID_new = new TH2D("htmean_hcalID_new","NEW ;hcal ID; clus t_{mean} [ns]", 288,-0.5,287.5,nbins,-30,30);
  
  TH2D *htmean_hcale_old = new TH2D("htmean_hcale_old","OLD ;hcal edep [GeV]; clus t_{mean} [ns]", nbins,0,0.5,nbins,-30,30);
  TH2D *htmean_hcale_new = new TH2D("htmean_hcale_new","NEW ;hcal edep [GeV]; clus t_{mean} [ns]", nbins,0,0.5,nbins,-30,30);
  
  //Really want a large quantity of zero-field LH2 data for this purpose:
  TH2D *htcoin_vs_HCALID_oldold = new TH2D("htcoin_vs_HCALID_oldold", "OLD hodo, OLD HCAL; HCAL ID; t_{HCAL}-t_{HODO,corr} [ns]",288,-0.5,287.5,250,-50,50); 
  TH2D *htcoin_vs_HCALID_oldnew = new TH2D("htcoin_vs_HCALID_oldnew", "NEW hodo, OLD HCAL; HCAL ID; t_{HCAL}-t_{HODO,corr} [ns]",288,-0.5,287.5,250,-50,50); 
  TH2D *htcoin_vs_HCALID_newold = new TH2D("htcoin_vs_HCALID_newold", "OLD hodo, NEW HCAL; HCAL ID; t_{HCAL}-t_{HODO,corr} [ns]",288,-0.5,287.5,250,-50,50); 
  TH2D *htcoin_vs_HCALID_newnew = new TH2D("htcoin_vs_HCALID_newnew", "NEW hodo, NEW HCAL; HCAL ID; t_{HCAL}-t_{HODO,corr} [ns]",288,-0.5,287.5,250,-50,50); 

  TH2D *htcoin_vs_HODOID_oldold = new TH2D("htcoin_vs_HODOID_oldold", "OLD hodo, OLD HCAL; Hodo Bar; t_{HCAL}-t_{HODO,corr} [ns]",90,-0.5,89.5,250,-50,50); 
  TH2D *htcoin_vs_HODOID_oldnew = new TH2D("htcoin_vs_HODOID_oldnew", "NEW hodo, OLD HCAL; Hodo Bar; t_{HCAL}-t_{HODO,corr} [ns]",90,-0.5,89.5,250,-50,50); 
  TH2D *htcoin_vs_HODOID_newold = new TH2D("htcoin_vs_HODOID_newold", "OLD hodo, NEW HCAL; Hodo Bar; t_{HCAL}-t_{HODO,corr} [ns]",90,-0.5,89.5,250,-50,50); 
  TH2D *htcoin_vs_HODOID_newnew = new TH2D("htcoin_vs_HODOID_newnew", "NEW hodo, NEW HCAL; Hodo Bar; t_{HCAL}-t_{HODO,corr} [ns]",90,-0.5,89.5,250,-50,50); 
  
  TH1D *hRes_all = new TH1D("hRes_all","",nbins,-4,4);
  TH1D *hRes_all_corr = new TH1D("hRes_all_corr","",nbins,-4,4);
  TH2D *h2d_ResTW_all = new TH2D("h2d_ResTW_all","",nbins,0,30,nbins,-4,4);
  TH2D *h2d_ResTW_all_corr = new TH2D("h2d_ResTW_all_corr","",nbins,0,30,nbins,-4,4);
  TH1D *hRes[90];
  TH1D *hRes_corr[90];
  fout->mkdir("hRes_bars");
  fout->cd("hRes_bars");
  for( int i=0; i<90; i++){
    hRes[i] = new TH1D(Form("hRes_%i",i),"",nbins,-4,4);
    hRes_corr[i] = new TH1D(Form("hRes_corr_%i",i),"",nbins,-4,4);
  }
  fout->cd();

  TH1D *hdt_vz[90];
  TH1D *hdt_vz_corr[90];
  fout->mkdir("hdt_vz_bars");
  fout->cd("hdt_vz_bars");
  for( int i=0; i<90; i++){
    hdt_vz[i] = new TH1D(Form("hdt_vz_%i",i),"",nbins,200,250);
    hdt_vz_corr[i] = new TH1D(Form("hdt_vz_corr_%i",i),"",nbins,200,250);
  }
  fout->cd();
  TH1D *hdt_vz_all = new TH1D("hdt_vz_all","",nbins,200,250);
  TH1D *hdt_vz_all_corr = new TH1D("hdt_vz_all_corr","",nbins,200,250);
  TH2D *h2d_dt_vz = new TH2D("h2d_dt_vz","",90,-0.5,89.5,nbins,200,250);
  TH2D *h2d_dt_vz_corr = new TH2D("h2d_dt_vz_corr","",90,-0.5,89.5,nbins,200,250);
  TH2D *h2d_dt_vzTW = new TH2D("h2d_dt_vzTW","",nbins,0,30,nbins,200,250);
  TH2D *h2d_dt_vzTW_corr = new TH2D("h2d_dt_vzTW_corr","",nbins,0,30,nbins,200,250);
  //hcal
  TH1D *hResHCAL_all = new TH1D("hResHCAL_all","",nbins,-4,4);
  TH1D *hResHCAL_all_corr = new TH1D("hResHCAL_all_corr","",nbins,-4,4);
  TH1D *hResHCAL[288];
  TH1D *hResHCAL_corr[288];
  fout->mkdir("hResHCAL_bars");
  fout->cd("hResHCAL_bars");
  for( int i=0; i<288; i++){
    hResHCAL[i] = new TH1D(Form("hResHCAL_%i",i),"",nbins,-4,4);
    hResHCAL_corr[i] = new TH1D(Form("hResHCAL_corr_%i",i),"",nbins,-4,4);
  }
  fout->cd();

  TH1D *hdtHCAL_vz[288];
  TH1D *hdtHCAL_vz_corr[288];
  fout->mkdir("hdtHCAL_vz_bars");
  fout->cd("hdtHCAL_vz_bars");
  for( int i=0; i<288; i++){
    hdtHCAL_vz[i] = new TH1D(Form("hdtHCAL_vz_%i",i),"",nbins,125,175);
    hdtHCAL_vz_corr[i] = new TH1D(Form("hdtHCAL_vz_corr_%i",i),"",nbins,125,175);
  }
  fout->cd();
  
  fout->mkdir("hTWHCAL")->cd();
  TH2D *hTWHCAL[288];
  TH2D *hTWHCAL_CORR[288];
  TH2D *hTWHCAL_test[288];
  for( int i=0; i<288; i++){
    hTWHCAL[i] = new TH2D(Form("hTWHCAL_%i",i),"",100,0,0.2,100,-25,25);
    hTWHCAL_CORR[i] = new TH2D(Form("hTWHCAL_CORR_%i",i),"",100,0,0.2,100,-25,25);
  }
  fout->cd();
  
  TH1D *hdtHCAL_vz_all = new TH1D("hdtHCAL_vz_all","",nbins,125,175);
  TH1D *hdtHCAL_vz_all_corr = new TH1D("hdtHCAL_vz_all_corr","",nbins,125,175);
  TH2D *h2d_dtHCAL_vz = new TH2D("h2d_dtHCAL_vz","",288,-0.5,287.5,nbins,125,175);
  TH2D *h2d_dtHCAL_vz_corr = new TH2D("h2d_dtHCAL_vz_corr","",288,-0.5,287.5,nbins,125,175);

  TH2D *htrig_BBvSBS = new TH2D("htrig_BBvSBS","BBTrig as measured in 1190 vs F1; 1190 t_{trig} [ns]; F1 t_{trig} [ns]",nbins,-1,-1,nbins,-1,-1);
  TH2D *hRF_BBvSBS = new TH2D("hRF_BBvSBS","RF as measured in 1190 vs F1; 1190 t_{rf} [ns]; F1 t_{rf} [ns]",nbins,-1,-1,nbins,-1,-1);
  TH2D *hBB_TRIGvRF = new TH2D("hBB_TRIGvRF","Trigger and RF in 1190; t_{trig} [ns]; t_{rf} [ns]",nbins,-1,-1,nbins,-1,-1);
  TH2D *hSBS_TRIGvRF = new TH2D("hSBS_TRIGvRF","Trigger and RF in F1; t_{trig} [ns]; t_{rf} [ns]",nbins,-1,-1,nbins,-1,-1);

  TH1D *hpred_mom = new TH1D("hpred_mom","; Predicted p_{N} [GeV]",100,0,3);
  TH1D *hpred_beta = new TH1D("hpred_beta","; Predicted #Beta",100,0.8,1);
  TH1D *hntof = new TH1D("hntof",";Predicted N_{TOF} [ns]",100,-1,1);
  TH1D *hW2_elastic = new TH1D("hW2_elastic",";W^{2} [GeV^{2}]",100,-1,-1);
  
  
  TH2D *hdxdy = new TH2D("hdxdy",";#Deltay [m]; #Deltax [m]",250,-4,4,250,-6,3);
  
  TH1D *htcoin = new TH1D("htcoin","",250,-50,50);
  TH1D *htcoin_corr = new TH1D("htcoin_corr","",250,-50,50);

  TH1D *htcoinRF = new TH1D("htcoinRF","",250,-50,50);
  TH1D *htcoinRF_corr = new TH1D("htcoinRF_corr","",250,-50,50);
  
  
  
  
  //event loop stuff
  skim_tree *T = new skim_tree(C);
  int treenum=-1,oldtreenum=-1;
  long nev = C->GetEntries();
  for( long ev=0; ev<nev; ev++ ){
    
    T->GetEntry(ev);
    if( ev % 100000 == 0 ) cout << "nevent = " << ev << " / " << nev << endl;
    
    //define some simple cuts here
    //skim file takes care of exclusive
    //condition of atleast 1 track and 1 hcal cluster
    if (T->bb_ntr<1) continue;
    //if (T->sbs_hcal_nclus<1) continue;
    if (T->hcal_blk_id->size()==0) continue;
    if(T->hcal_blk_id->size()==0 || T->hcal_tdctime == 0 || T->hodo_bar->size() == 0) continue;
    
    //3 hits is the minimum anyway.
    //this can be made tighter, 4 or 5 hits, if we want
    //if (T->bb_tr_nhits < 3) continue;
    
    //in principle we can clean up electron selection
    //jack has calibrated grinch timing for all kinematics
    //if (T->gr_trindex ==-1) continue;
    //if (T->gr_size < 3) continue;

    //quasi invariant mass squared cut
    //for GEN2 can leave fairly wide. even 0<W2<2 would probably be fine
    if (T->W2 < 0 || T->W2 > 2.0) continue;
    if(T->ps_e < 0.2) continue;
    if (fabs(T->dy - 0.63)>2*0.41) continue;
    //grab needed track parameters:
    double vz = T->bb_vz;
    double pathl = T->bb_pathl;
    double etof = pathl/0.299792458; //electron TOF from vertex to hodo.

    //for now predicted N pathl and tof from H2 elastics
    double npathl = T->sbs_pathl;
    double pred_beta = GetBeta(Mp,T->pred_mom);
    double ntof = npathl / (pred_beta * 0.2997);
    double ntof0 = dHCAL / (pred_beta * 0.2997);
    hpred_mom->Fill(T->pred_mom);
    hpred_beta->Fill(pred_beta);
    hntof->Fill(ntof);
    
    //kine stuff from skim tree
    double Q2 = T->Q2;
    double W2 = T->W2;
    double dx = T->dx;// + T->proton_deflection_x;
    double dy = T->dy;
    hW2_elastic->Fill(W2);
    
    //new trigger branches
    double bb_trigtime = T->bb_trigtime;
    double bb_rftime = T->bb_rftime;
    double sbs_trigtime = T->sbs_trigtime;
    double sbs_rftime = T->sbs_rftime;

    double T_rf = 2.008;
      
    //hodo fit stuff
    int ID;
    double tmean_old, tdiff_old;
    double tmean_CORR, tdiff_CORR;
    //if (T->hodo_bar->size()>0){
      //relevant hodo cluster parameters
      double tleft = T->hodo_tleft->at(0);
      double tright = T->hodo_tright->at(0);
      double totleft = T->hodo_totleft->at(0);
      double totright = T->hodo_totright->at(0);
      
      ID = int(T->hodo_bar->at(0));
      
      double yhodo = T->bb_y+zhodo*T->bb_ph_fp;
      
      // how to define our chi2 statistic? We want to correct all hodo PMT times to tvertex = 0:
      // tPMT = etof + t0 - walk * TOT + d/vscint
      // chi2 = sum_{i=1}^N \sum_{j=1}^Nhit (tPMT-etof-t0+walk*TOT - d/vscint)^2
      // dchi2/dt0_k = 2*\sum_i,j (tPMT-etof-t0+walk*TOT-d/vscint)*-1
      // dchi2/dwalk_k = 2*sum_i,j (tPMT-etof-t0+walk*TOT-d/vscint)*TOT
      // dchi2/d(1/v)_k = 2*sum_i,j (tPMT-etof-t0+walk*TOT-d/vscint)*-d

      double dLEFT = std::min(Lbar_hodo,std::max(0.0,Lbar_hodo/2.0 - yhodo));
      double dRIGHT = std::min(Lbar_hodo,std::max(0.0,Lbar_hodo/2.0 + yhodo));
      
      // we want to minimize
      // chi2 = sum_{i=1}^Nevent (tleft_corr^2 + tright_corr^2)
      // dchi2/dt0_j = sum_i (2*tleft_corr_i*dtleft_corr_i/dt0_j + 2*tright_corr_i*dtright_corr_i/dt0_j) = 0
      int ipar_t0 = ID;
      int ipar_vinv = ID+90;
      int ipar_wL = ID+180;
      int ipar_wR = ID+270;
    
      //sum_i (tL - etof - t0 + w*TOTL - dL/v) = 0  
      // tL-etof = t0L -w*TOTL + dL/v
      //there is a more elegant way to do this:
      //Left PMT:
      Mhodo(ipar_t0,ipar_t0) += 1.0;
      Mhodo(ipar_t0,ipar_vinv) += dLEFT;
      Mhodo(ipar_t0,ipar_wL) += -totleft;
      
      Mhodo(ipar_vinv, ipar_t0) += dLEFT;
      Mhodo(ipar_vinv, ipar_vinv) += pow(dLEFT,2);
      Mhodo(ipar_vinv, ipar_wL) += -totleft*dLEFT;
      
      Mhodo(ipar_wL, ipar_t0) += -totleft;
      Mhodo(ipar_wL, ipar_vinv) += -totleft*dLEFT;
      Mhodo(ipar_wL, ipar_wL) += pow(totleft,2);
      
      bhodo(ipar_t0) += (tleft - (etof-etof0))*1.0;
      bhodo(ipar_vinv) += (tleft -(etof-etof0))*dLEFT;
      bhodo(ipar_wL) += (tleft - (etof-etof0))*(-totleft);
      

      //Right PMT:
      Mhodo(ipar_t0,ipar_t0) += 1.0;
      Mhodo(ipar_t0,ipar_vinv) += dRIGHT;
      Mhodo(ipar_t0,ipar_wR) += -totright;
      
      Mhodo(ipar_vinv, ipar_t0) += dRIGHT;
      Mhodo(ipar_vinv, ipar_vinv) += pow(dRIGHT,2);
      Mhodo(ipar_vinv, ipar_wR) += -totright*dRIGHT;
      
      Mhodo(ipar_wR, ipar_t0) += -totright;
      Mhodo(ipar_wR, ipar_vinv) += -totright*dRIGHT;
      Mhodo(ipar_wR, ipar_wR) += pow(totright,2);
      
      bhodo(ipar_t0) += (tright - (etof-etof0))*1.0;
      bhodo(ipar_vinv) += (tright -(etof-etof0))*dRIGHT;
      bhodo(ipar_wR) += (tright - (etof-etof0))*(-totright);
    
      //if this is the second time running, calculate corrected variables, otherwise
      //these will just come out the same as the "old" variables since all new factors
      //will be zero upon first iteration.
      double tleft_CORR = tleft - (etof-etof0) - HODOt0[ID] + HODOwL[ID]*totleft - dLEFT/vscint[ID];
      double tright_CORR = tright - (etof-etof0) - HODOt0[ID] + HODOwR[ID]*totright - dRIGHT/vscint[ID];

      //calculate a corrected time difference without the propagation correction:
      tdiff_CORR = tleft_CORR - tright_CORR;
      tmean_CORR = 0.5 * (tleft_CORR + tright_CORR);
      double totmean = 0.5 * (totleft + totright);
      tmean_old = T->hodo_tmean->at(0);
      tdiff_old = T->hodo_tdiff->at(0);
    
      double tmean_raw = (tleft + tright)/2;
      double tdiff_raw = tleft - tright;
    
      //hall B vertex time calcs
      double rf_offset = -1.0;//1.0;
      double t_vz = tmean_CORR + bb_trigtime;
      double t_start = bb_rftime + vz/0.299792458;
      double dt_vz = t_vz - t_start;
      double dt_vz_corr = dt_vz - tpad_rf[ID];
      double tR = fmod(dt_vz,T_rf) + rf_offset;
      double tR_corr = fmod(dt_vz_corr,T_rf) + rf_offset;
      

      //Fill histograms
      htmean_hodoID_old->Fill( ID, tmean_old );
      htmean_hodoID_raw->Fill( ID, tmean_raw );
      htmean_hodoID_new->Fill( ID, tmean_CORR );

      htdiff_yhodo_old->Fill( yhodo, tdiff_old );
      htdiff_yhodo_raw->Fill( yhodo, tdiff_raw );
      htdiff_yhodo_new->Fill( yhodo, tdiff_CORR );
    
      //hodoscope resids (fmod times)
      hRes_all->Fill(tR);
      hRes_all_corr->Fill(tR_corr);
      h2d_ResTW_all->Fill(totmean, tR);
      h2d_ResTW_all_corr->Fill(totmean, tR_corr);
      hRes[ID]->Fill(tR);
      hRes_corr[ID]->Fill(tR_corr);

      //t_vz - t_hodo_start
      hdt_vz[ID]->Fill(dt_vz);
      hdt_vz_corr[ID]->Fill(dt_vz_corr);
      hdt_vz_all->Fill(dt_vz);
      hdt_vz_all_corr->Fill(dt_vz_corr);
      h2d_dt_vz->Fill(ID, dt_vz);
      h2d_dt_vz_corr->Fill(ID,dt_vz_corr);
      h2d_dt_vzTW->Fill(totmean, dt_vz);
      h2d_dt_vzTW_corr->Fill(totmean, dt_vz_corr);
      
      //}

    //hcal fit stuff
    //since we are doing timewalk, need to do it on the blks, then re-calculate the
    //new energy weighted cluster time average after individual blocks have TW corr
    double num_sum = 0;
    double denom_sum = 0;
    for(uint i=0; i<T->hcal_blk_id->size(); i++){
      int idblkHCAL = (int) T->hcal_blk_id->at(i) - 1;
      double tHCAL = T->hcal_blk_tdctime->at(i);
      double eblkHCAL = T->hcal_blk_e->at(i);

      //need these cuts to keep the minisation well constrained as a result of
      //unphysical noise outside these regions
      if(eblkHCAL>0.01 && tHCAL > 0 && tHCAL < 20){

	//note i have started with "all positive" definitions
	//tPMT-ntof = t0 + w0*E + w1/sqrt(E)
	// chi2 = sum_{i=1}^N \sum_{j=1}^Nhit (tPMT-ntof-t0 - w0*E - w1/sqrt(E))^2
	// dchi2/dt0_k = 2*\sum_i,j (tPMT-etof - t0 - w0*E - w1/sqrt(E))*-1
	// dchi2/dw0_k = 2*\sum_i,j (tPMT-etof - t0 - w0*E - w1/sqrt(E))*-E
	// dchi2/dw1_k = 2*\sum_i,j (tPMT-etof - t0 - w0*E - w1/sqrt(E))*-1/sqrt(E)
      
	int ipar_t0 = idblkHCAL;
	int ipar_w0 = idblkHCAL + 288; 
	int ipar_w1 = idblkHCAL + 2*288; 
      
	Mhcal(ipar_t0,ipar_t0) += 1.0;
	Mhcal(ipar_t0,ipar_w0) += eblkHCAL;
	Mhcal(ipar_t0,ipar_w1) += 1/sqrt(eblkHCAL);
      
	Mhcal(ipar_w0,ipar_t0) += eblkHCAL;
	Mhcal(ipar_w0,ipar_w0) += pow(eblkHCAL,2);
	Mhcal(ipar_w0,ipar_w1) += sqrt(eblkHCAL);
      
	Mhcal(ipar_w1,ipar_t0) += 1/sqrt(eblkHCAL);
	Mhcal(ipar_w1,ipar_w0) += sqrt(eblkHCAL);
	Mhcal(ipar_w1,ipar_w1) += 1/eblkHCAL;
      
      
	//(ntof - ntof0);
	//bhcal(ipar_t0) += (tHCAL) * 1.0;
	//bhcal(ipar_w0) += (tHCAL) * (eblkHCAL);
	//bhcal(ipar_w1) += (tHCAL) * (1/sqrt(eblkHCAL));
	bhcal(ipar_t0) += (tHCAL - (ntof - ntof0)) * 1.0;
	bhcal(ipar_w0) += (tHCAL - (ntof - ntof0)) * (eblkHCAL);
	bhcal(ipar_w1) += (tHCAL - (ntof - ntof0)) * (1/sqrt(eblkHCAL));
      }
      
      //double tHCAL_CORR = tHCAL - HCALt0[idblkHCAL] - HCALw0[idblkHCAL]*eblkHCAL - HCALw1[idblkHCAL]/sqrt(eblkHCAL);
      double tHCAL_CORR = tHCAL - (ntof - ntof0) - HCALt0[idblkHCAL] - HCALw0[idblkHCAL]*eblkHCAL - HCALw1[idblkHCAL]/sqrt(eblkHCAL);
	
      num_sum += tHCAL_CORR*eblkHCAL;
      denom_sum += eblkHCAL;
      
      hTWHCAL[idblkHCAL]->Fill(eblkHCAL,tHCAL);
      hTWHCAL_CORR[idblkHCAL]->Fill(eblkHCAL,tHCAL_CORR);
      
      double rf_sbsoffset = 0.0;
      double tHCAL_vz = tHCAL_CORR + sbs_trigtime;
      double tHCAL_start = sbs_rftime + vz/0.299792458; 
      double dtHCAL_vz = tHCAL_vz - tHCAL_start;
      double dtHCAL_vz_corr = dtHCAL_vz - tpad_rf_hcal[idblkHCAL];
      double tRHCAL = fmod(dtHCAL_vz, T_rf) + rf_sbsoffset;
      double tRHCAL_corr = fmod(dtHCAL_vz_corr, T_rf) + rf_sbsoffset;
      
      
      //hcal resids (fmods)? later
      hResHCAL_all->Fill(tRHCAL);
      hResHCAL_all_corr->Fill(tRHCAL_corr);
      hResHCAL[idblkHCAL]->Fill(tRHCAL);
      hResHCAL_corr[idblkHCAL]->Fill(tRHCAL_corr);
      
      //t_vz_hcal - t_hcal_start
      hdtHCAL_vz_all->Fill(dtHCAL_vz);
      hdtHCAL_vz_all_corr->Fill(dtHCAL_vz_corr);
      h2d_dtHCAL_vz->Fill(idblkHCAL,dtHCAL_vz);
      h2d_dtHCAL_vz_corr->Fill(idblkHCAL, dtHCAL_vz_corr);
      hdtHCAL_vz[idblkHCAL]->Fill(dtHCAL_vz);
      hdtHCAL_vz_corr[idblkHCAL]->Fill(dtHCAL_vz_corr);
      
      htmean_hcalID_old->Fill( idblkHCAL, tHCAL );
      htmean_hcalID_new->Fill( idblkHCAL, tHCAL_CORR );
      
      htmean_hcale_old->Fill( eblkHCAL, tHCAL );
      htmean_hcale_new->Fill( eblkHCAL, tHCAL_CORR );
    }
    //now calculate new cluster E weighted time average
    double hcal_tmean_corr;
    if(num_sum==0 || denom_sum ==0)
      hcal_tmean_corr = 0;
    else
      hcal_tmean_corr = num_sum / denom_sum;
    
    //draw v1190 vs F1 trigtime correlations
    htrig_BBvSBS->Fill(bb_trigtime, sbs_trigtime);
    hRF_BBvSBS->Fill(bb_rftime, sbs_rftime);
    hBB_TRIGvRF->Fill(bb_rftime, bb_trigtime);
    hSBS_TRIGvRF->Fill(sbs_rftime, sbs_trigtime);

    //draw dxdy corrected for proton deflection, should be a QE blob around (0,0)
    hdxdy->Fill( dy,dx );

    
    //if QE cuts passed, do coincidence time stuff (i assume?)
    //if( T->hcal_e>0.02 && W2>W2min && W2<W2max && sqrt(pow(dx,2)+pow(dy,2))<=0.24 ){
      // cout << "(Eprime,etheta,ephi)=("
      //      << ep.Mag() << ", " << etheta*57.3 << ", "
      //      << ephi*57.3 << ")" << endl;
      // cout << "(ptheta_etheta,pp_etheta)=("
      //      << ptheta_etheta*57.3 << ", " << pp_etheta << ")" << endl;
	
      // cout << "(ephi,pphi)=(" << ephi*57.3 << ", " << pphi_ephi*57.3 << ")"
      //      << endl;
      
      int idblkHCAL = (int) T->hcal_id - 1; 
      double eblkHCAL = T->hcal_blk_e->at(0);
      double tHCAL = T->hcal_tdctime;
      double tHCAL_CORR = tHCAL - (HCALt0[idblkHCAL] + HCALw0[idblkHCAL]*eblkHCAL + HCALw1[idblkHCAL]/sqrt(eblkHCAL));
      double rf_sbsoffset = 0.0;
      double tHCAL_vz = tHCAL_CORR + sbs_trigtime;
      double tHCAL_start = sbs_rftime + vz/0.299792458; 
      double dtHCAL_vz = tHCAL_vz - tHCAL_start;
      double dtHCAL_vz_corr = dtHCAL_vz - tpad_rf_hcal[idblkHCAL];
      double tRHCAL = fmod(dtHCAL_vz, T_rf) + rf_sbsoffset;
      double tRHCAL_corr = fmod(dtHCAL_vz_corr, T_rf) + rf_sbsoffset;
      
      
      htcoin_vs_HCALID_oldold->Fill( idblkHCAL, T->hcal_tdctime-tmean_old);
      htcoin_vs_HCALID_oldnew->Fill( idblkHCAL, tHCAL-tmean_CORR);
      htcoin_vs_HCALID_newold->Fill( idblkHCAL, tHCAL_CORR-tmean_old);
      htcoin_vs_HCALID_newnew->Fill( idblkHCAL, tHCAL_CORR-tmean_CORR);

      htcoin_vs_HODOID_oldold->Fill( ID, tHCAL-tmean_old);
      htcoin_vs_HODOID_oldnew->Fill( ID, tHCAL-tmean_CORR);
      htcoin_vs_HODOID_newold->Fill( ID, tHCAL_CORR-tmean_old);
      htcoin_vs_HODOID_newnew->Fill( ID, tHCAL_CORR-tmean_CORR);

      htcoin->Fill(T->hcal_tdctime-tmean_old);
      htcoin_corr->Fill(tHCAL_CORR-tmean_CORR);

      htcoinRF->Fill(dtHCAL_vz - dt_vz);
      htcoinRF_corr->Fill(dtHCAL_vz - dt_vz_corr);
      
      //}
  }//end of event loop


  TDecompSVD Ahodo(Mhodo);
  Ahodo.Solve(bhodo);

  TDecompSVD Ahcal(Mhcal);
  Ahcal.Solve(bhcal);

  TF1 *fres[90]; 
  for(int i=0; i<90; i++){
    fres[i] = new TF1(Form("fres_%i",i),"gaus",-2,2);
    int maxbin = hRes_corr[i]->GetMaximumBin();
    double max = hRes_corr[i]->GetXaxis()->GetBinCenter(maxbin);
    hRes_corr[i]->Fit(fres[i],"Q","",max-0.5,max+0.5);
    bhodo[i+360] = tpad_rf[i] + fres[i]->GetParameter(1);
  }

  TF1 *fTW[288];
  for(int i=0; i<288; i++){
    fTW[i] = new TF1(Form("fTW{%i",i),"[0]+[1]*x+[2]/sqrt(x)",0,0.2);
    fTW[i]->FixParameter(0, bhcal[i]);
    fTW[i]->FixParameter(1, bhcal[i+288]);
    fTW[i]->FixParameter(2, bhcal[i+2*288]);
    hTWHCAL[i]->Fit(fTW[i],"QR","");
  }
  //bhodo.Print();
  //bhcal.Print();
  
  ofstream hodoconfig;
  hodoconfig.open("hodoparams.txt");
  for ( int i=0; i< nparams_hodo; i++){
    hodoconfig << bhodo[i] << endl;
  }
  hodoconfig.close();
  
  ofstream hcalconfig("hcalparams.txt");
  for ( int i=0; i< nparams_HCAL; i++){
    hcalconfig << bhcal[i] << endl;
  }
  hcalconfig.close();

  fout->Write();
      
}

