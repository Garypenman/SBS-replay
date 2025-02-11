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

void TOFcal(Int_t kin_no=2, TString Target = "H2"){
  
  gStyle->SetPalette(kRainBow);
  gROOT->SetBatch(kTRUE);
  SetKinematics(kin_no,Target);
  
  TString currentline;
  
  TChain *C = new TChain("T");
  if(Target == "H2")
    C->Add(Form("/volatile/halla/sbs/gpenman/data/out_gen%i_H2.root",kin_no));
  else
    C->Add(Form("/volatile/halla/sbs/gpenman/data/out_gen%i_He3_qeneutron.root",kin_no));
  if(C)
    cout << "Chain setup" << endl;
  SkimTree *T = new SkimTree(C);
  long nev = C->GetEntries();
  if(T)
    cout << "SkimTree setup" << endl;

  
  //Default reference PMTs for aligning all other PMTs:
  int hodorefID = 44; // reference BAR (not PMT); we'll set the LEFT PMT to the reference by convention:
  int hcalrefID = 199; //row 16 column 6
  
  double zhodo = 1.854454; //meters
  double Lbar_hodo = 0.6; //meters
  double etof0 = (1.63+3.0)/0.299792458; // "central" TOF value:
  
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

  cout << "basic experimental variables initialised" << endl;
  
  //parameters per hcal-hodo channel combo
  //90 * (1 offset, 2 Hodo walk correction slopes, 1 hodo propagation speed, 1 paddle rf offset)
  //288 * (1 offset, 2 Hcal walk params for hall B like fit, 1 rf offsets)
  int nparams_hodo = 90 * (1 + 2 + 1 + 1);
  int nparams_hcal = 288 * (1 + 2 + 1);
  int nparams = nparams_hodo + nparams_hcal;

  //read in previously generated config for iterative calibration
  vector<double> params;
  vector<double> t0, HODOwL, HODOwR, vscint, tpad_rf; 
  vector<double> t1, HCALw0, HCALw1, tpad_rf_hcal;
  TString configfile = Form("params_GEN%i_%s.txt",kin_no,Target.Data());
  
  //get hodo params if they exist otherwise fill vector of zeros
  ifstream configin(configfile);
  if(!configin){
    cout << "Global config file not found, assuming first run of calibration." << endl;
    for (int i=0; i<nparams; i++)
      params.push_back(0.0);
  }
  else{
    while(currentline.ReadLine(configin)){
      params.push_back(stod(currentline.Data()));
    }
  }
  configin.close();

  cout << "parameter vector initialised" << endl;
  
  //now fill the vectors associated with each set of parameters
  for( int i=0; i<90; i++){
    t0.push_back (params[i] );
    vscint.push_back( 1.0/params[i+90] );
    HODOwL.push_back( params[i+180] );
    HODOwR.push_back( params[i+270] );
    tpad_rf.push_back( params[i+360] );
  }
  for ( int i=0; i<288; i++){
    int j=nparams_hodo+i;
    t1.push_back( params[j] );
    HCALw0.push_back( params[j+288] );
    HCALw1.push_back( params[j+2*288] );
    tpad_rf_hcal.push_back( params[j+3*288] );
  }
  cout << "parameter vector filled" << endl;
  
  //Ordering of params:
  //t0
  // t0i, walkileft, walkiright, 1/vi
  // later, rfi
  TMatrixD Mh(nparams,nparams);
  TVectorD bh(nparams); 

  for( int ipar=0; ipar<nparams; ipar++ ){
    for( int jpar=0; jpar<nparams; jpar++ ){
      Mh(ipar,jpar) = 0.0;
    }
    bh(ipar) = 0.0;
  }
  cout << "Matrix setup" << endl;
  
  
  //declare outfile here incase we want to choose to abort execution later on due to missing input files
  //but not overwrite previously generated outputfiles etc
  TString outputfilename=Form("TOFcal_GEN%i_%s.root",kin_no,Target.Data());
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
  TH2D *htcoin_vs_HCALID_oldold = new TH2D("htcoin_vs_HCALID_oldold", "OLD hodo, OLD HCAL; HCAL ID; t_{HCAL}-t_{HODO,corr} [ns]",288,-0.5,287.5,250,-10,10); 
  TH2D *htcoin_vs_HCALID_oldnew = new TH2D("htcoin_vs_HCALID_oldnew", "NEW hodo, OLD HCAL; HCAL ID; t_{HCAL}-t_{HODO,corr} [ns]",288,-0.5,287.5,250,-10,10); 
  TH2D *htcoin_vs_HCALID_newold = new TH2D("htcoin_vs_HCALID_newold", "OLD hodo, NEW HCAL; HCAL ID; t_{HCAL}-t_{HODO,corr} [ns]",288,-0.5,287.5,250,-10,10); 
  TH2D *htcoin_vs_HCALID_newnew = new TH2D("htcoin_vs_HCALID_newnew", "NEW hodo, NEW HCAL; HCAL ID; t_{HCAL}-t_{HODO,corr} [ns]",288,-0.5,287.5,250,-10,10); 

  TH2D *htcoin_vs_HODOID_oldold = new TH2D("htcoin_vs_HODOID_oldold", "OLD hodo, OLD HCAL; Hodo Bar; t_{HCAL}-t_{HODO,corr} [ns]",90,-0.5,89.5,250,-10,10); 
  TH2D *htcoin_vs_HODOID_oldnew = new TH2D("htcoin_vs_HODOID_oldnew", "NEW hodo, OLD HCAL; Hodo Bar; t_{HCAL}-t_{HODO,corr} [ns]",90,-0.5,89.5,250,-10,10); 
  TH2D *htcoin_vs_HODOID_newold = new TH2D("htcoin_vs_HODOID_newold", "OLD hodo, NEW HCAL; Hodo Bar; t_{HCAL}-t_{HODO,corr} [ns]",90,-0.5,89.5,250,-10,10); 
  TH2D *htcoin_vs_HODOID_newnew = new TH2D("htcoin_vs_HODOID_newnew", "NEW hodo, NEW HCAL; Hodo Bar; t_{HCAL}-t_{HODO,corr} [ns]",90,-0.5,89.5,250,-10,10); 

  fout->mkdir("hTWHODO")->cd();
  TH2D *hTWHODOL[90];
  TH2D *hTWHODOL_corr[90];
  TH2D *hTWHODOR[90];
  TH2D *hTWHODOR_corr[90];
  for( int i=0; i<90; i++){
    hTWHODOL[i] = new TH2D(Form("hTWHODOL_%i",i),"",100,0,30,100,-25,25);
    hTWHODOL_corr[i] = new TH2D(Form("hTWHODOL_corr_%i",i),"",100,0,30,100,-25,25);
    hTWHODOR[i] = new TH2D(Form("hTWHODOR_%i",i),"",100,0,30,100,-25,25);
    hTWHODOR_corr[i] = new TH2D(Form("hTWHODOR_corr_%i",i),"",100,0,30,100,-25,25);
  }
  fout->cd();
  
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
  TH2D *hTWHCAL_corr[288];
  for( int i=0; i<288; i++){
    hTWHCAL[i] = new TH2D(Form("hTWHCAL_%i",i),"",100,0,0.2,100,-25,25);
    hTWHCAL_corr[i] = new TH2D(Form("hTWHCAL_corr_%i",i),"",100,0,0.2,100,-25,25);
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
  TH1D *hntof = new TH1D("hntof",";Predicted N_{TOF} [ns]",100,-1,-1);
  TH1D *hW2_elastic = new TH1D("hW2_elastic",";W^{2} [GeV^{2}]",100,-1,-1);
  
  
  TH2D *hdxdy = new TH2D("hdxdy",";#Deltay [m]; #Deltax [m]",250,-4,4,250,-6,3);
  
  TH1D *htcoin = new TH1D("htcoin","",250,-50,50);
  TH1D *htcoin_corr = new TH1D("htcoin_corr","",250,-50,50);

  TH1D *htcoinRF = new TH1D("htcoinRF","",250,-50,50);
  TH1D *htcoinRF_corr = new TH1D("htcoinRF_corr","",250,-50,50);
  
  
  
  //event loopx
  for( long ev=0; ev<nev; ev++ ){
    
    T->GetEntry(ev);
    if( ev % 100000 == 0 ) cout << "nevent = " << ev << " / " << nev << endl;
    
    //define some simple cuts here
    //skim file takes care of exclusive
    //condition of atleast 1 track and 1 hcal cluster
    if (T->bb_ntr<1) continue;
    //if (T->sbs_hcal_nclus<1) continue;
    if (T->hcal_id==0)continue;
    if(T->hcal_blk_id->size()==0 || T->hcal_tdctime == 0 || T->hodo_bar->size() == 0 || T->hodo_refright == -999) continue;
    
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
    if (fabs(T->dy)>0.5) continue;

    
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
    if (bb_trigtime == -999 || bb_rftime == -999 || sbs_trigtime == -999 || sbs_rftime == -999) continue;
    double T_rf = 2.008;
    
    //hodo fit stuff
    int hodo_id = int(T->hodo_bar->at(0));

    double tmean_old = T->hodo_tmean->at(0);
    double tdiff_old = T->hodo_tdiff->at(0);
    
    double tleft = T->hodo_tleft->at(0);
    double tright = T->hodo_tright->at(0);
    double totleft = T->hodo_totleft->at(0);
    double totright = T->hodo_totright->at(0);
    
    
    double yhodo = T->bb_y+zhodo*T->bb_ph_fp;
    
    double dLEFT = std::min(Lbar_hodo,std::max(0.0,Lbar_hodo/2.0 - yhodo));
    double dRIGHT = std::min(Lbar_hodo,std::max(0.0,Lbar_hodo/2.0 + yhodo));
    
    
    // we want to minimize
    //note i have started with "all positive" definitions
    // hodo
    // tPMT = etof + t0 - walk * TOT + d/vscint
    // chi2 = sum_{i=1}^N \sum_{j=1}^Nhit (tPMT-etof-t0+walk*TOT - d/vscint)^2
    // dchi2/dt0_k = 2*\sum_i,j (tPMT-etof-t0+walk*TOT-d/vscint)*-1
    // dchi2/dwalk_k = 2*sum_i,j (tPMT-etof-t0+walk*TOT-d/vscint)*TOT
    // dchi2/d(1/v)_k = 2*sum_i,j (tPMT-etof-t0+walk*TOT-d/vscint)*-d
    // hcal
    //tPMT-ntof = t0 + w0*E + w1/sqrt(E)
    // chi2 = sum_{i=1}^N \sum_{j=1}^Nhit (tPMT-ntof-t0 - w0*E - w1/sqrt(E))^2
    // dchi2/dt0_k = 2*\sum_i,j (tPMT-etof - t0 - w0*E - w1/sqrt(E))*-1
    // dchi2/dw0_k = 2*\sum_i,j (tPMT-etof - t0 - w0*E - w1/sqrt(E))*-E
    // dchi2/dw1_k = 2*\sum_i,j (tPMT-etof - t0 - w0*E - w1/sqrt(E))*-1/sqrt(E)
    // chi2 = sum_{i=1}^Nevent (tleft_corr^2 + tright_corr^2)
    // dchi2/dt0_j = sum_i (2*tleft_corr_i*dtleft_corr_i/dt0_j + 2*tright_corr_i*dtright_corr_i/dt0_j) = 0
      
    
      
    //if this is the second time running, calculate corrected variables, otherwise
    //these will just come out the same as the "old" variables since all new factors
    //will be zero upon first iteration.
    double tleft_corr = tleft - (etof-etof0) - t0[hodo_id] - HODOwL[hodo_id]*totleft - dLEFT/vscint[hodo_id];
    double tright_corr = tright - (etof-etof0) - t0[hodo_id] - HODOwR[hodo_id]*totright - dRIGHT/vscint[hodo_id];
    //double tleft_corr = tleft - (etof-etof0);
    //double tright_corr = tright - (etof-etof0);
    
    //calculate a corrected time difference without the propagation correction:
    double tdiff_corr = tleft_corr - tright_corr;
    double tmean_corr = 0.5 * (tleft_corr + tright_corr);
    double totmean = 0.5 * (totleft + totright);
    
    double tmean_raw = 0.5 * (tleft + tright) - (etof-etof0);
    double tdiff_raw = tleft - tright;
    
    //hall B vertex time calcs
    double rf_offset = -1.0;//1.0;
    double t_vz = tmean_corr + T->hodo_refleft;
    double t_start = bb_rftime + vz/0.299792458;
    double dt_vz = t_vz - t_start;
    double dt_vz_corr = dt_vz - tpad_rf[hodo_id];
    double tR = fmod(dt_vz,T_rf) + rf_offset;
    double tR_corr = fmod(dt_vz_corr,T_rf) + rf_offset;
    
    //Fill histograms
    htmean_hodoID_old->Fill( hodo_id, tmean_old );
    htmean_hodoID_raw->Fill( hodo_id, tmean_raw );
    htmean_hodoID_new->Fill( hodo_id, tmean_corr );

    htdiff_yhodo_old->Fill( yhodo, tdiff_old );
    htdiff_yhodo_raw->Fill( yhodo, tdiff_raw );
    htdiff_yhodo_new->Fill( yhodo, tdiff_corr );

    //hodoscope tw
    hTWHODOL[hodo_id]->Fill( totleft, tleft );
    hTWHODOL_corr[hodo_id]->Fill( totleft, tleft_corr );
    hTWHODOR[hodo_id]->Fill( totright, tright );
    hTWHODOR_corr[hodo_id]->Fill( totright, tright_corr );

    //hodoscope resids (fmod times)
    hRes_all->Fill(tR);
    hRes_all_corr->Fill(tR_corr);
    h2d_ResTW_all->Fill(totmean, tR);
    h2d_ResTW_all_corr->Fill(totmean, tR_corr);
    hRes[hodo_id]->Fill(tR);
    hRes_corr[hodo_id]->Fill(tR_corr);

    //t_vz - t_hodo_start
    hdt_vz[hodo_id]->Fill(dt_vz);
    hdt_vz_corr[hodo_id]->Fill(dt_vz_corr);
    hdt_vz_all->Fill(dt_vz);
    hdt_vz_all_corr->Fill(dt_vz_corr);
    h2d_dt_vz->Fill(hodo_id, dt_vz);
    h2d_dt_vz_corr->Fill(hodo_id,dt_vz_corr);
    h2d_dt_vzTW->Fill(totmean, dt_vz);
    h2d_dt_vzTW_corr->Fill(totmean, dt_vz_corr);
      
      
    //hcal fit stuff
    //cout << T->hcal_blk_id->size() << endl;
    for(uint i=0; i<T->hcal_blk_id->size(); i++){
      int hcal_id = (int) T->hcal_blk_id->at(i) - 1;
      if (hcal_id==-1) continue;
      double hcal_e = T->hcal_blk_e->at(i);
      double tHCAL = T->hcal_blk_tdctime->at(i);
      
      //need these cuts to keep the minisation well constrained as a result of
      //unphysical noise outside these regions
      if(!(hcal_e>0.01 && tHCAL > 0 && tHCAL < 20)) continue;
      
      int ipar_t0 = hodo_id;
      int ipar_vinv = hodo_id + 90;
      int ipar_wL = hodo_id + 180;
      int ipar_wR = hodo_id + 270;
      int ipar_rf_hodo = hodo_id + 360;
      int ipar_t1 = nparams_hodo + hcal_id;
      int ipar_w0 = nparams_hodo + hcal_id + 288; 
      int ipar_w1 = nparams_hodo  + hcal_id + 2*288; 
      int ipar_rf_hcal = nparams_hodo + hcal_id + 3*288;
      
      Mh(ipar_t0,ipar_t0) += 1.0;
      Mh(ipar_t1,ipar_t1) += 1.0;
      //Mh(ipar_t0,ipar_t1) += -1.0;
      //Mh(ipar_t1,ipar_t0) += -1.0;
      
      Mh(ipar_t1,ipar_w0) += hcal_e;
      Mh(ipar_t1,ipar_w1) += 1/sqrt(hcal_e);
      
      Mh(ipar_w0,ipar_t1) += hcal_e;
      Mh(ipar_w0,ipar_w0) += pow(hcal_e,2);
      Mh(ipar_w0,ipar_w1) += sqrt(hcal_e);
      
      Mh(ipar_w1,ipar_t1) += 1/sqrt(hcal_e);
      Mh(ipar_w1,ipar_w0) += sqrt(hcal_e);
      Mh(ipar_w1,ipar_w1) += 1/hcal_e;
      
      
      //tleft
      Mh(ipar_t0,ipar_vinv) += -dLEFT;
      Mh(ipar_t0,ipar_wL) += -totleft;
      
      Mh(ipar_vinv, ipar_t0) += -dLEFT;
      Mh(ipar_vinv, ipar_vinv) += pow(dLEFT,2);
      Mh(ipar_vinv, ipar_wL) += totleft*dLEFT;
      
      Mh(ipar_wL, ipar_t0) += -totleft;
      Mh(ipar_wL, ipar_vinv) += totleft*dLEFT;
      Mh(ipar_wL, ipar_wL) += pow(totleft,2);
      
      
      //tright
      Mh(ipar_t0,ipar_vinv) += -dRIGHT;
      Mh(ipar_t0,ipar_wR) += -totright;
      
      Mh(ipar_vinv, ipar_t0) += -dRIGHT;
      Mh(ipar_vinv, ipar_vinv) += pow(dRIGHT,2);
      Mh(ipar_vinv, ipar_wR) += totright*dRIGHT;
      
      Mh(ipar_wR, ipar_t0) += -totright;
      Mh(ipar_wR, ipar_vinv) += totright*dRIGHT;
      Mh(ipar_wR, ipar_wR) += pow(totright,2);
      
      bh(ipar_t1) += (tHCAL-tmean_raw) * 1.0;
      bh(ipar_t0) += -(tHCAL-tmean_raw) * 1.0;
      bh(ipar_w0) += (tHCAL-tmean_raw) * (hcal_e);
      bh(ipar_w1) += (tHCAL-tmean_raw) * (1/sqrt(hcal_e));
      //bh(ipar_w0) += (tHCAL - (ntof - ntof0)) * (hcal_e);
      //bh(ipar_w1) += (tHCAL - (ntof - ntof0)) * (1/sqrt(hcal_e));
      bh(ipar_vinv) += -(tHCAL-tmean_raw)*(dRIGHT+dLEFT);
      bh(ipar_wR) += -(tright - (etof-etof0))*totright;
bh(ipar_wL) += -(tleft - (etof-etof0))*totleft;


      //if this is the second time running, calculate corrected variables, otherwise
      //these will just come out the same as the "old" variables since all new factors
      //will be zero upon first iteration.
      //double tHCAL_corr = tHCAL - (ntof - ntof0) - HCALw0[hcal_id]*hcal_e - HCALw1[hcal_id]/sqrt(hcal_e);
      double tHCAL_corr = tHCAL - t1[hcal_id] - HCALw0[hcal_id]*hcal_e - HCALw1[hcal_id]/sqrt(hcal_e);
      //double tHCAL_corr = tHCAL;
    
      hTWHCAL[hcal_id]->Fill(hcal_e,tHCAL);
      hTWHCAL_corr[hcal_id]->Fill(hcal_e,tHCAL_corr);
    
      double rf_sbsoffset = 0.0;
      double tHCAL_vz = tHCAL_corr + sbs_trigtime;
      double tHCAL_start = sbs_rftime + vz/0.299792458; 
      double dtHCAL_vz = tHCAL_vz - tHCAL_start;
      double dtHCAL_vz_corr = dtHCAL_vz - tpad_rf_hcal[hcal_id];
      double tRHCAL = fmod(dtHCAL_vz, T_rf) + rf_sbsoffset;
      double tRHCAL_corr = fmod(dtHCAL_vz_corr, T_rf) + rf_sbsoffset;
    
      //fill hcal histos
      htmean_hcalID_old->Fill( hcal_id, tHCAL );
      htmean_hcalID_new->Fill( hcal_id, tHCAL_corr );
    
      htmean_hcale_old->Fill( hcal_e, tHCAL );
      htmean_hcale_new->Fill( hcal_e, tHCAL_corr );
    
      //hcal resids (fmods)? later
      hResHCAL_all->Fill(tRHCAL);
      hResHCAL_all_corr->Fill(tRHCAL_corr);
      hResHCAL[hcal_id]->Fill(tRHCAL);
      hResHCAL_corr[hcal_id]->Fill(tRHCAL_corr);
    
      //t_vz_hcal - t_hcal_start
      hdtHCAL_vz_all->Fill(dtHCAL_vz);
      hdtHCAL_vz_all_corr->Fill(dtHCAL_vz_corr);
      h2d_dtHCAL_vz->Fill(hcal_id,dtHCAL_vz);
      h2d_dtHCAL_vz_corr->Fill(hcal_id, dtHCAL_vz_corr);
      hdtHCAL_vz[hcal_id]->Fill(dtHCAL_vz);
      hdtHCAL_vz_corr[hcal_id]->Fill(dtHCAL_vz_corr);
    

      //draw v1190 vs F1 trigtime correlations
      htrig_BBvSBS->Fill(bb_trigtime, sbs_trigtime);
      hRF_BBvSBS->Fill(bb_rftime, sbs_rftime);
      hBB_TRIGvRF->Fill(bb_rftime, bb_trigtime);
      hSBS_TRIGvRF->Fill(sbs_rftime, sbs_trigtime);

      //draw dxdy corrected for proton deflection, should be a QE blob around (0,0)
      hdxdy->Fill( dy,dx );

      //coincidence time plots
      htcoin_vs_HCALID_oldold->Fill( hcal_id, T->hcal_tdctime-tmean_old);
      htcoin_vs_HCALID_oldnew->Fill( hcal_id, tHCAL-tmean_corr);
      htcoin_vs_HCALID_newold->Fill( hcal_id, tHCAL_corr-tmean_old);
      htcoin_vs_HCALID_newnew->Fill( hcal_id, tHCAL_corr-tmean_corr);

      htcoin_vs_HODOID_oldold->Fill( hodo_id, tHCAL-tmean_old);
      htcoin_vs_HODOID_oldnew->Fill( hodo_id, tHCAL-tmean_corr);
      htcoin_vs_HODOID_newold->Fill( hodo_id, tHCAL_corr-tmean_old);
      htcoin_vs_HODOID_newnew->Fill( hodo_id, tHCAL_corr - tmean_corr);

      htcoin->Fill(tHCAL-tmean_raw);
      htcoin_corr->Fill(tHCAL_corr - tmean_corr);

      htcoinRF->Fill(dtHCAL_vz - dt_vz);
      htcoinRF_corr->Fill(dtHCAL_vz - dt_vz_corr);
    }
  }//end of event loop


  TDecompSVD Ah(Mh);
  Ah.Solve(bh);

  
  TF1 *fres[90]; 
  for(int i=0; i<90; i++){
    fres[i] = new TF1(Form("fres_%i",i),"gaus",-2,2);
    int maxbin = hRes_corr[i]->GetMaximumBin();
    double max = hRes_corr[i]->GetXaxis()->GetBinCenter(maxbin);
    hRes_corr[i]->Fit(fres[i],"Q","",max-0.5,max+0.5);
    //bh[360+i] = tpad_rf[i] + fres[i]->GetParameter(1);
    bh[360+i] = 0;
  }

  TF1 *fTW[288];
  for(int i=0; i<288; i++){
    int j=nparams_hodo+i;
    fTW[i] = new TF1(Form("fTW{%i",i),"[0]+[1]*x+[2]/sqrt(x)",0,0.2);
    fTW[i]->SetParameter(0, bh[j]);
    fTW[i]->FixParameter(1, bh[j+288]);
    fTW[i]->FixParameter(2, bh[j+2*288]);
    hTWHCAL[i]->Fit(fTW[i],"QR","");
  }
  //bh.Print();
  
  
  ofstream config;
  config.open(configfile);
  for ( int i=0; i< nparams; i++){
    config << bh[i] << endl;
  }
  config.close();
  
  
  fout->Write();
      
}

