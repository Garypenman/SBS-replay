#include "/home/gpenman/Analysis/SBS-Analysis/GEn/include/gen_ana.h"

void tof_skim(const Int_t kin_no = 2, const TString Target = "H2", TString rootfile_dir = ""){
  ROOT::EnableImplicitMT(32);

  if (Target == "H2"){
    cout << "Target is H2 Reference Cell" << endl;
  }
  else if (Target == "He3"){
    cout << "Target is He3 Production Cell" << endl;
  }
  else{
    cout << "Target doesn't make sense. Try harder." << endl;
    exit(1);
  }
  
  TChain* C = new TChain("T");
  
  //cache location of pass1
  rootfile_dir = Form("/cache/halla/sbs/prod/GEnII/pass1/GEN%i/%s/rootfiles",kin_no,Target.Data());
  if(kin_no == 2  && Target == "H2")
    //rootfile_dir = "/cache/halla/sbs/prod/GEnII/pass1/GEN2/H2/SBS100/rootfiles/";
    rootfile_dir = "/volatile/halla/sbs/gpenman/GEN_REPLAY/rootfiles/pass1/GEN2/H2/SBS100/";
  
  if(kin_no==4) {
    char input;
    while (true) {
      cout << "Analyse GEN4b? [y/n]" << endl;
      cin >> input;
      
      if ((input == 'y') || (input == 'n')) {
        break;
      }
    }
    if (input == 'y'){
      cout << "Adding GEN4a files, and changing rootfile_dir to GEN4b location" << endl;
      gSystem->Setenv("ROOTFILE_DIR_A",rootfile_dir); 
      C->Add("$ROOTFILE_DIR_A/e1209016_fullreplay_*.root");
      rootfile_dir = Form("/cache/halla/sbs/prod/GEnII/pass1/GEN4b/%s/rootfiles",Target.Data());
    }
  }

  gSystem->Setenv("ROOTFILE_DIR",rootfile_dir); 
  
  C->Add("$ROOTFILE_DIR/e1209016_fullreplay_*.root");
  
  SetKinematics(kin_no,Target);
  
  double hcal_voffset = 0.1;
  double hcal_hoffset = 0.25;
  
  //setting the correct coordinates for hcal axis in the lab frame
  //sean method
  TVector3 hcal_zaxis(sin(-th_sbs),0,cos(-th_sbs)); // Clock-wise rotation about Y axis
  TVector3 hcal_xaxis(0,-1,0);                                  
  TVector3 hcal_yaxis = hcal_zaxis.Cross(hcal_xaxis).Unit();
  std::vector<TVector3> hcal_axis = {hcal_xaxis, hcal_yaxis, hcal_zaxis};
  TVector3 hcal_origin = hcal_dist*hcal_axis[2] + hcal_voffset*hcal_axis[0] + hcal_hoffset*hcal_axis[1];
  
  //my method
  //TVector3 hcal_origin( -hcal_dist*sin(th_sbs), 0, hcal_dist*cos(th_sbs) );
  //TVector3 hcal_zaxis = hcal_origin.Unit();
  //TVector3 hcal_xaxis(0,-1,0);
  //TVector3 hcal_yaxis = hcal_zaxis.Cross( hcal_xaxis ).Unit();
  
  //can reconstruct photon 4-vec in the event loop
  TLorentzVector Tp4(0,0,0,Mp); //target 4vec
  TLorentzVector kp4(0,0,Eb,Eb); //beam 4vec
  TLorentzVector Qp4, kpp4, Rp4; //q, recoil electron, recoil nucleon
  //End of setting kinematics
  
  
  //stopwatch for benchmarking analysis time on large dataset
  TStopwatch t;
  t.Start();
  
  //outfile and chain check
  cout << "Checking tree, this may take some time dont worry." << endl;
  Long64_t nev = C->GetEntries();
  if(nev == 0){
    cout << gSystem->Getenv("ROOTFILE_DIR") << " appears to be empty. Check files and/or paths." << endl;
    exit(1);
  }
  
  //setup ttreereader
  TTreeReader rd (C);

  //kinematic vars
  TTreeReaderValue<Double_t> e_kine_Q2 = {rd,"e.kine.Q2"}; 
  TTreeReaderValue<Double_t> e_kine_W2 = {rd,"e.kine.W2"}; 
  TTreeReaderValue<Double_t> e_kine_epsilon = {rd,"e.kine.epsilon"}; 
  TTreeReaderValue<Double_t> e_kine_nu = {rd,"e.kine.nu"}; 
  TTreeReaderValue<Double_t> e_kine_x_bj ={rd,"e.kine.x_bj"}; 

  //bb tracks
  TTreeReaderArray<Double_t> bb_tr_x = {rd,"bb.tr.x"};
  TTreeReaderArray<Double_t> bb_tr_y = {rd,"bb.tr.y"};
  TTreeReaderArray<Double_t> bb_tr_px = {rd,"bb.tr.px"};
  TTreeReaderArray<Double_t> bb_tr_py = {rd,"bb.tr.py"};
  TTreeReaderArray<Double_t> bb_tr_pz = {rd,"bb.tr.pz"};
  TTreeReaderArray<Double_t> bb_tr_p = {rd,"bb.tr.p"};
  TTreeReaderArray<Double_t> bb_tr_th = {rd,"bb.tr.th"};
  TTreeReaderArray<Double_t> bb_tr_ph = {rd,"bb.tr.ph"};
  TTreeReaderArray<Double_t> bb_tr_vx = {rd,"bb.tr.vx"};
  TTreeReaderArray<Double_t> bb_tr_vy = {rd,"bb.tr.vy"};
  TTreeReaderArray<Double_t> bb_tr_vz = {rd,"bb.tr.vz"};
  TTreeReaderArray<Double_t> bb_tr_pathl = {rd, "bb.tr.pathl"};
   
  TTreeReaderValue<Double_t> bb_tr_n = {rd,"bb.tr.n"};
 
  TTreeReaderArray<Double_t> bb_tr_tg_th = {rd,"bb.tr.tg_th"};
  TTreeReaderArray<Double_t> bb_tr_tg_ph = {rd,"bb.tr.tg_ph"};
  TTreeReaderArray<Double_t> bb_tr_tg_y = {rd,"bb.tr.tg_y"};
  
  //TTreeReaderValue<Double_t> bb_gem_trigtime = {rd, "bb.gem.trigtime"};
  TTreeReaderArray<Double_t> bb_gem_track_nhits = {rd, "bb.gem.track.nhits"};

  //BBCal main cluster
  TTreeReaderValue<Double_t> bb_ps_x = {rd,"bb.ps.x"};
  TTreeReaderValue<Double_t> bb_sh_x = {rd,"bb.sh.x"};
  TTreeReaderValue<Double_t> bb_ps_y = {rd,"bb.ps.y"};
  TTreeReaderValue<Double_t> bb_sh_y = {rd,"bb.sh.y"};
  TTreeReaderValue<Double_t> bb_ps_e = {rd,"bb.ps.e"};
  TTreeReaderValue<Double_t> bb_sh_e = {rd,"bb.sh.e"};
  TTreeReaderValue<Double_t> bb_ps_atimeblk = {rd, "bb.ps.atimeblk"};
  TTreeReaderValue<Double_t> bb_sh_atimeblk = {rd, "bb.sh.atimeblk"};


  //hcal clusters
  TTreeReaderValue<Double_t> sbs_hcal_nclus = {rd,"sbs.hcal.nclus"};
  //times of highest E blk in cluster
  TTreeReaderArray<Double_t> sbs_hcal_clus_atimeblk = {rd, "sbs.hcal.clus.atimeblk"};
  TTreeReaderArray<Double_t> sbs_hcal_clus_tdctimeblk = {rd, "sbs.hcal.clus.tdctimeblk"};
  //energy weighted times of cluster
  TTreeReaderArray<Double_t> sbs_hcal_clus_adctime = {rd, "sbs.hcal.clus.adctime"};
  TTreeReaderArray<Double_t> sbs_hcal_clus_tdctime = {rd, "sbs.hcal.clus.tdctime"};
  TTreeReaderArray<Double_t> sbs_hcal_clus_x = {rd,"sbs.hcal.clus.x"};
  TTreeReaderArray<Double_t> sbs_hcal_clus_y = {rd,"sbs.hcal.clus.y"};
  TTreeReaderArray<Double_t> sbs_hcal_clus_e = {rd,"sbs.hcal.clus.e"};
  TTreeReaderArray<Double_t> sbs_hcal_clus_eblk = {rd,"sbs.hcal.clus.eblk"};
  TTreeReaderArray<Double_t> sbs_hcal_clus_row = {rd, "sbs.hcal.clus.row"};
  TTreeReaderArray<Double_t> sbs_hcal_clus_col = {rd, "sbs.hcal.clus.col"};
  TTreeReaderArray<Double_t> sbs_hcal_clus_id = {rd, "sbs.hcal.clus.id"};

  //blks in cluster
  TTreeReaderValue<Int_t> Ndata_sbs_hcal_clus_blk_id = {rd, "Ndata.sbs.hcal.clus_blk.id"};
  TTreeReaderArray<Double_t> sbs_hcal_clus_blk_x = {rd,"sbs.hcal.clus_blk.x"};
  TTreeReaderArray<Double_t> sbs_hcal_clus_blk_y = {rd,"sbs.hcal.clus_blk.y"};
  TTreeReaderArray<Double_t> sbs_hcal_clus_blk_e = {rd,"sbs.hcal.clus_blk.e"};
  TTreeReaderArray<Double_t> sbs_hcal_clus_blk_row = {rd, "sbs.hcal.clus_blk.row"};
  TTreeReaderArray<Double_t> sbs_hcal_clus_blk_col = {rd, "sbs.hcal.clus_blk.col"};
  TTreeReaderArray<Double_t> sbs_hcal_clus_blk_id = {rd, "sbs.hcal.clus_blk.id"};
  TTreeReaderArray<Double_t> sbs_hcal_clus_blk_atime = {rd, "sbs.hcal.clus_blk.atime"};
  TTreeReaderArray<Double_t> sbs_hcal_clus_blk_tdctime = {rd, "sbs.hcal.clus_blk.tdctime"};
  
  //hcal trigger stuff
  TTreeReaderValue<Int_t> Ndata_sbs_hcal_Ref_tdcelemID = {rd, "Ndata.sbs.hcal.Ref.tdcelemID"};
  TTreeReaderArray<Double_t> sbs_hcal_Ref_tdcelemID = {rd, "sbs.hcal.Ref.tdcelemID"};
  TTreeReaderArray<Double_t> sbs_hcal_Ref_tdc = {rd, "sbs.hcal.Ref.tdc"};
  
  //Event info
  TTreeReaderValue<ULong64_t> fEvtHdr_fEvtTime = {rd, "fEvtHdr.fEvtTime"};
  TTreeReaderValue<unsigned int> fEvtHdr_fEvtNum = {rd, "fEvtHdr.fEvtNum"};
  TTreeReaderValue<unsigned int> fEvtHdr_fEvtType = {rd, "fEvtHdr.fEvtType"};
  TTreeReaderValue<unsigned int> fEvtHdr_fEvtLen = {rd, "fEvtHdr.fEvtLen"};
  TTreeReaderValue<Int_t> fEvtHdr_fHelicity = {rd, "fEvtHdr.fHelicity"};
  TTreeReaderValue<unsigned int> fEvtHdr_fTrigBits = {rd, "fEvtHdr.fTrigBits"};
  TTreeReaderValue<unsigned int> fEvtHdr_fRun = {rd, "fEvtHdr.fRun"};

  //Hodoscope main cluster bars
  TTreeReaderValue<Int_t> Ndata_bb_hodotdc_clus_bar_tdc_id = {rd, "Ndata.bb.hodotdc.clus.bar.tdc.id"};
  TTreeReaderArray<Double_t> bb_hodotdc_clus_bar_tdc_id = {rd, "bb.hodotdc.clus.bar.tdc.id"};
  TTreeReaderArray<Double_t> bb_hodotdc_clus_bar_tdc_itrack = {rd, "bb.hodotdc.clus.bar.tdc.itrack"};
  TTreeReaderArray<Double_t> bb_hodotdc_clus_bar_tdc_meantime = {rd, "bb.hodotdc.clus.bar.tdc.meantime"};
  TTreeReaderArray<Double_t> bb_hodotdc_clus_bar_tdc_meantot = {rd, "bb.hodotdc.clus.bar.tdc.meantot"};
  TTreeReaderArray<Double_t> bb_hodotdc_clus_bar_tdc_tleft = {rd, "bb.hodotdc.clus.bar.tdc.tleft"};
  TTreeReaderArray<Double_t> bb_hodotdc_clus_bar_tdc_tright = {rd, "bb.hodotdc.clus.bar.tdc.tright"};
  TTreeReaderArray<Double_t> bb_hodotdc_clus_bar_tdc_totleft = {rd, "bb.hodotdc.clus.bar.tdc.totleft"};
  TTreeReaderArray<Double_t> bb_hodotdc_clus_bar_tdc_totright = {rd, "bb.hodotdc.clus.bar.tdc.totright"};
  TTreeReaderArray<Double_t> bb_hodotdc_clus_bar_tdc_timediff = {rd, "bb.hodotdc.clus.bar.tdc.timediff"};
  TTreeReaderArray<Double_t> bb_hodotdc_clus_bar_tdc_vpos = {rd, "bb.hodotdc.clus.bar.tdc.vpos"};
  TTreeReaderArray<Double_t> bb_hodotdc_clus_bar_tdc_timehitpos = {rd, "bb.hodotdc.clus.bar.tdc.timehitpos"};
  
  //grinch
  //TTreeReaderValue<Double_t> bb_grinch_tdc_bestcluster = {rd, "bb.grinch_tdc.bestcluster"};
  TTreeReaderValue<Double_t> bb_grinch_tdc_clus_adc = {rd, "bb.grinch_tdc.clus.adc"};
  TTreeReaderValue<Double_t> bb_grinch_tdc_clus_size = {rd, "bb.grinch_tdc.clus.size"};
  TTreeReaderValue<Double_t> bb_grinch_tdc_clus_t_mean = {rd, "bb.grinch_tdc.clus.t_mean"};
  TTreeReaderValue<Double_t> bb_grinch_tdc_clus_t_rms = {rd, "bb.grinch_tdc.clus.t_rms"};
  TTreeReaderValue<Double_t> bb_grinch_tdc_clus_tot_mean = {rd, "bb.grinch_tdc.clus.tot_mean"};
  TTreeReaderValue<Double_t> bb_grinch_tdc_clus_trackindex = {rd, "bb.grinch_tdc.clus.trackindex"};
  TTreeReaderValue<Double_t> bb_grinch_tdc_clus_x_mean = {rd, "bb.grinch_tdc.clus.x_mean"};
  TTreeReaderValue<Double_t> bb_grinch_tdc_clus_y_mean = {rd, "bb.grinch_tdc.clus.y_mean"};
  //TTreeReaderValue<Double_t> bb_grinch_tdc_nclus = {rd, "bb.grinch_tdc.nclus"};
  //TTreeReaderValue<Double_t> bb_grinch_tdc_ngoodhits = {rd, "bb.grinch_tdc.ngoodhits"};
  //TTreeReaderValue<Double_t> bb_grinch_tdc_ntrackmatch = {rd, "bb.grinch_tdc.ntrackmatch"};
  
  TTreeReaderValue<Int_t> Ndata_bb_tdctrig_tdc = {rd, "Ndata.bb.tdctrig.tdc"};
  TTreeReaderArray<Double_t> bb_tdctrig_tdc = {rd, "bb.tdctrig.tdc"};
  TTreeReaderValue<Int_t> Ndata_bb_tdctrig_tdcelemID = {rd, "Ndata.bb.tdctrig.tdcelemID"};
  TTreeReaderArray<Double_t> bb_tdctrig_tdcelemID = {rd, "bb.tdctrig.tdcelemID"};
  
  TTreeReaderValue<Double_t> IGL1I00OD16_16 = {rd, "IGL1I00OD16_16"};
  TTreeReaderValue<Double_t> scalhel_hel = {rd, "scalhel.hel"};
  
  //setup output file
  TString outfilename = Form("~/vol/data/out_tof_gen%i_%s.root",kin_no,Target.Data());
  TFile *outfile = new TFile(outfilename,"RECREATE");
  
  TTree *tree = new TTree("T","Output of GEN tof Skim");
  double Q2, W2, eps, nu, x_bj;
  double bb_x,bb_y,bb_px,bb_py,bb_pz,bb_p,bb_vx,bb_vy,bb_vz,bb_th_fp,bb_ph_fp,bb_th_tg,bb_ph_tg,bb_y_tg;
  double bb_pathl,sbs_pathl;
  double sh_e,ps_e,sh_x,sh_y,ps_x,ps_y;
  double hcal_e,hcal_eblk,hcal_x,hcal_y, hcal_row, hcal_col, hcal_id;
  tree->Branch("Q2",&Q2);
  tree->Branch("W2",&W2);
  tree->Branch("eps",&eps);
  tree->Branch("nu",&nu);
  tree->Branch("x_bj",&x_bj);

  tree->Branch("bb_x",&bb_x);
  tree->Branch("bb_y",&bb_y);
  tree->Branch("bb_px",&bb_px);
  tree->Branch("bb_py",&bb_py);
  tree->Branch("bb_pz",&bb_pz);
  tree->Branch("bb_p",&bb_p);
  tree->Branch("bb_vx",&bb_vx);
  tree->Branch("bb_vy",&bb_vy);
  tree->Branch("bb_vz",&bb_vz);
  tree->Branch("bb_th_fp",&bb_th_fp);
  tree->Branch("bb_th_tg",&bb_th_tg);
  tree->Branch("bb_ph_fp",&bb_ph_fp);
  tree->Branch("bb_ph_tg",&bb_ph_tg);
  tree->Branch("bb_y_tg",&bb_y_tg);

  
  tree->Branch("bb_pathl",&bb_pathl);
  tree->Branch("sbs_pathl",&sbs_pathl);

  double proton_defl_x;
  tree->Branch("proton_defl_x",&proton_defl_x);
  
  double bb_tr_nhits, bb_ntr;
  tree->Branch("bb_tr_nhits",&bb_tr_nhits);
  tree->Branch("bb_ntr",&bb_ntr);
  
  tree->Branch("sh_e",&sh_e);
  tree->Branch("sh_x",&sh_x);
  tree->Branch("sh_y",&sh_y);
  
  tree->Branch("ps_e",&ps_e);
  tree->Branch("ps_x",&ps_x);
  tree->Branch("ps_y",&ps_y);
  
  tree->Branch("hcal_e",&hcal_e);
  tree->Branch("hcal_eblk",&hcal_eblk);
  tree->Branch("hcal_x",&hcal_x);
  tree->Branch("hcal_y",&hcal_y);
  tree->Branch("hcal_row",&hcal_row);
  tree->Branch("hcal_col",&hcal_col);
  tree->Branch("hcal_id",&hcal_id);

  vector<double> hcal_blk_e;
  vector<double> hcal_blk_x;
  vector<double> hcal_blk_y;
  vector<double> hcal_blk_row;
  vector<double> hcal_blk_col;
  vector<double> hcal_blk_id;
  vector<double> hcal_blk_adctime;
  vector<double> hcal_blk_tdctime;
  tree->Branch("hcal_blk_e",&hcal_blk_e);
  tree->Branch("hcal_blk_x",&hcal_blk_x);
  tree->Branch("hcal_blk_y",&hcal_blk_y);
  tree->Branch("hcal_blk_row",&hcal_blk_row);
  tree->Branch("hcal_blk_col",&hcal_blk_col);
  tree->Branch("hcal_blk_id",&hcal_blk_id);
  tree->Branch("hcal_blk_tdctime",&hcal_blk_tdctime);
  tree->Branch("hcal_blk_adctime",&hcal_blk_adctime);
  
  double pred_ang_horiz, pred_ang_vert, pred_y, pred_x, pred_mom, dx, dy, dev;
  tree->Branch("pred_ang_horiz",&pred_ang_horiz);
  tree->Branch("pred_ang_vert",&pred_ang_vert);
  tree->Branch("pred_y",&pred_y);
  tree->Branch("pred_x",&pred_x);
  tree->Branch("pred_mom",&pred_mom);
  tree->Branch("dx",&dx);
  tree->Branch("dy",&dy);
  tree->Branch("dev",&dev);
  
  double ps_atime, sh_atime, hcal_adctime, hcal_atimeblk, coin_atime, coin_atime_corr;
  tree->Branch("ps_atime",&ps_atime);
  tree->Branch("sh_atime",&sh_atime);
  tree->Branch("hcal_adctime",&hcal_adctime); 
  tree->Branch("hcal_atimeblk",&hcal_atimeblk); 
  tree->Branch("coin_atime",&coin_atime);
  //tree->Branch("coin_atime_corr",&coin_atime_corr);
  
  double hcal_tdctime, hcal_tdctimeblk, coin_time, coin_time_corr;
  tree->Branch("hcal_tdctime",&hcal_tdctime);
  tree->Branch("hcal_tdctimeblk",&hcal_tdctimeblk);
  tree->Branch("coin_time",&coin_time);
  //tree->Branch("coin_time_corr",&coin_time_corr);

  vector<double> hcal_refid;
  vector<double> hcal_ref;
  tree->Branch("hcal_refid",&hcal_refid);
  tree->Branch("hcal_ref",&hcal_ref);
  
  
  vector<double> hodo_tmean;
  vector<double> hodo_tdiff;
  vector<double> hodo_totmean;
  vector<double> hodo_tleft;
  vector<double> hodo_tright;
  vector<double> hodo_totleft;
  vector<double> hodo_totright;
  vector<double> hodo_timehitpos;
  vector<double> hodo_bar;
  double hodo_nclus, hodo_nbar;
  tree->Branch("hodo_tmean",&hodo_tmean);
  tree->Branch("hodo_totmean",&hodo_totmean);
  tree->Branch("hodo_tdiff",&hodo_tdiff);
  tree->Branch("hodo_tleft",&hodo_tleft);
  tree->Branch("hodo_tright",&hodo_tright);
  tree->Branch("hodo_totleft",&hodo_totleft);
  tree->Branch("hodo_totright",&hodo_totright);
  tree->Branch("hodo_timehitpos",&hodo_timehitpos);
  tree->Branch("hodo_bar",&hodo_bar);
  tree->Branch("hodo_nclus",&hodo_nclus);
  tree->Branch("hodo_nbar",&hodo_nbar);
  
  double gr_adc,gr_size,gr_tmean,gr_totmean,gr_trindex,gr_xmean,gr_ymean,gr_nclus,gr_ngoodhits,gr_ntrackmatch;
  tree->Branch("gr_adc",&gr_adc);
  tree->Branch("gr_size",&gr_size);
  tree->Branch("gr_tmean",&gr_tmean);
  tree->Branch("gr_totmean",&gr_totmean);
  tree->Branch("gr_trindex",&gr_trindex);
  tree->Branch("gr_xmean",&gr_xmean);
  tree->Branch("gr_ymean",&gr_ymean);
  //tree->Branch("gr_nclus",&gr_nclus);
  //tree->Branch("gr_ngoodhits",&gr_ngoodhits);
  //tree->Branch("gr_ntrackmatch",&gr_ntrackmatch);
  
  //beam/target helicity information
  double helicity;
  int helicity_status;
  tree->Branch("helicity",&helicity);
  tree->Branch("helicity_status",&helicity_status);
  double He3Pol;
  tree->Branch("He3Pol",&He3Pol);
  double BeamPol;
  tree->Branch("BeamPol",&BeamPol);
  
  ULong64_t evt_time;
  tree->Branch("evt_time",&evt_time);
  uint evt_num, evt_type, evt_len;
  tree->Branch("evt_num",&evt_num);
  tree->Branch("evt_type",&evt_type);
  tree->Branch("evt_len",&evt_len);
  int evt_hel;
  tree->Branch("evt_hel",&evt_hel);
  uint run_num, trigbits;
  tree->Branch("run_num",&run_num);
  tree->Branch("trigbits",&trigbits);
  
  vector<int> tdctrig_id;
  vector<double> tdctrig;
  tree->Branch("tdctrig_id",&tdctrig_id);
  tree->Branch("tdctrig",&tdctrig);

  double bb_trigtime,bb_rftime,sbs_trigtime,sbs_rftime;
  tree->Branch("bb_trigtime",&bb_trigtime);
  tree->Branch("bb_rftime",&bb_rftime);
  tree->Branch("sbs_trigtime",&sbs_trigtime);
  tree->Branch("sbs_rftime",&sbs_rftime);
  
  //unique event id for using brufit splot fitting
  Long64_t UID=-1;
  tree->Branch("UID",&UID);
  
  //Data loop
  Long64_t ev=-1;
  //nev = 1000000;
  int treenum = 0, currenttreenum = 0, currentrunnum = 0;
  time_t run_time_unix;
  cout << "Processing " << nev << " events." << endl;

  int nev_good=0;
  int no_bbrf=0;
  int no_bbtrig=0;
  int no_sbsrf=0;
  int no_sbstrig=0;
  
  TString lastfile = "";
  while(rd.Next() && ev < nev){
    
    //std::cout << "C->GetName()    " << C->GetName() << std::endl;
    //std::cout << "C->GetCurrentFile()    " << C->GetCurrentFile() << std::endl;
    //std::cout << "C->GetCurrentFile()->GetName()    " << C->GetCurrentFile()->GetName() << std::endl;
    TString thisfile = C->GetCurrentFile()->GetName();
    if(thisfile != lastfile){
      cout << thisfile << endl;
      lastfile = thisfile;
    }
    
    ev = rd.GetCurrentEntry();
    UID += 1;
    
    if( ev % 100000 == 0 )
      cout << ev << " / " << nev << ": " << *fEvtHdr_fRun << endl;
    // Pre-cuts -- these should be in the replay cdef
    //if( *bb_tr_n <= 0 ) continue; 
    
    //hcal 2 arm cuts
    //if ( *sbs_hcal_nclus <= 0) continue;
    //if ( sbs_hcal_clus_tdctimeblk[0] == -1000  ||  sbs_hcal_clus_tdctime[0] == 0.00) continue;

    nev_good++;
    
    //target z vertex cut (get rid of end caps)
    //if (Target == "He3")
    //if(fabs(bb_tr_vz[0]) > 0.27) continue;
    
    //if(*fEvtHdr_fTrigBits != 4) continue;
    
    //bbgem track hits cut
    //if( bb_gem_track_nhits[0] <= 3 ) continue;
    
    //no hodo time cut?
    //if (*Ndata_bb_hodotdc_clus_bar_tdc_id <= 0) continue;
    //bad hcal tdc time cut
    //if (*sbs_hcal_tdctimeblk <= -999) continue;
    
    evt_time = *fEvtHdr_fEvtTime;
    evt_num = *fEvtHdr_fEvtNum;
    evt_type = *fEvtHdr_fEvtType;
    evt_len = *fEvtHdr_fEvtLen;
    evt_hel = *fEvtHdr_fHelicity;
    trigbits = *fEvtHdr_fTrigBits;
    run_num = *fEvtHdr_fRun;
    
    // _tg = target plane, _fp = focal plane or in-plane, tr_ = local transport coord system 
    //set variables in bigbite
    bb_x = bb_tr_x[0];
    bb_y = bb_tr_y[0];
    
    bb_px = bb_tr_px[0];
    bb_py = bb_tr_py[0];
    bb_pz = bb_tr_pz[0];
    bb_p = bb_tr_p[0];
    
    bb_vx = bb_tr_vx[0];
    bb_vy = bb_tr_vy[0];
    bb_vz = bb_tr_vz[0];
    
    bb_th_fp = bb_tr_th[0];
    bb_ph_fp = bb_tr_ph[0];
    bb_th_tg = bb_tr_tg_th[0];
    bb_ph_tg = bb_tr_tg_ph[0];
    bb_y_tg = bb_tr_tg_y[0];
    
    bb_pathl = bb_tr_pathl[0];

    //trigtime = *bb_gem_trigtime;
    bb_tr_nhits = bb_gem_track_nhits[0];
    bb_ntr = *bb_tr_n;
    
    sh_e = *bb_sh_e;  
    ps_e = *bb_ps_e;
    sh_x = *bb_sh_x;
    sh_y = *bb_sh_y;
    ps_x = *bb_ps_x;
    ps_y = *bb_ps_y;
    
    //single arm calculations
    kpp4.SetPxPyPzE(bb_px,bb_py,bb_pz,bb_p);
    Qp4 = kp4 - kpp4;
    TVector3 qdir = Qp4.Vect().Unit();
    Rp4 = Tp4 + Qp4;
    
    Q2 = *e_kine_Q2;
    W2 = *e_kine_W2;
    eps = *e_kine_epsilon;
    nu = *e_kine_nu;
    x_bj = *e_kine_x_bj;
    //Q2 = -Qp4.Mag2();
    //W2 = Rp4.Mag2(); 
    //eps = 0;
    //nu = kp4.E() - kpp4.E();
    //x_bj = Q2/ (2*Tp4*Qp4);
    
    //if (W2 > 2.0) continue;
    
    Double_t trx_sh  = (bb_x + (show_dist) * bb_th_fp );
    Double_t try_sh  = (bb_y + (show_dist) * bb_ph_fp );
    
    Rp4.RotateY(th_sbs);
    pred_ang_horiz = TMath::ATan(Rp4.Px() / Rp4.Pz());
    pred_ang_vert = TMath::ATan(Rp4.Py() / sqrt(Rp4.Px()*Rp4.Px()  +  Rp4.Pz()*Rp4.Pz()));

    //Calculate expected proton deflection using crude model:
    double GEMpitch = 10.0*TMath::DegToRad();
    double sbsfield = 1.0; //only 0.3 for GEN2 H2 30% data
    double sbsmaxfield = 1.3; //Tesla
    double sbsdist = 2.8;
    double Dgap = 48.0*2.54/100.0; //about 1.22 m
    double BdL = sbsfield * sbsmaxfield * Dgap;
    double thetabend_prot = 0.3 * BdL / Rp4.Rho(); 
    proton_defl_x = tan(thetabend_prot) * (hcal_dist - (sbsdist + Dgap/2.0) );
    //double proton_deflection_pathl = (hcal_dist - (sbsdist + Dgap/2.0)) / cos(thetabend_prot);
    
    //Now we need to calculate the "true" trajectory bend angle for the electron from the reconstructed angles:
    TVector3 enhat_tgt( bb_th_tg, bb_ph_tg, 1.0 );
    enhat_tgt = enhat_tgt.Unit();
	
    TVector3 enhat_fp( bb_th_fp, bb_ph_fp, 1.0 );
    enhat_fp = enhat_fp.Unit();

    TVector3 GEMzaxis(-sin(GEMpitch),0,cos(GEMpitch));
    TVector3 GEMyaxis(0,1,0);
    TVector3 GEMxaxis = (GEMyaxis.Cross(GEMzaxis)).Unit();
	
    TVector3 enhat_fp_rot = enhat_fp.X() * GEMxaxis + enhat_fp.Y() * GEMyaxis + enhat_fp.Z() * GEMzaxis;

    double thetabend = acos( enhat_fp_rot.Dot( enhat_tgt ) );
    
    //set variables in sbs
    hcal_e = sbs_hcal_clus_e[0];
    hcal_eblk = sbs_hcal_clus_eblk[0];
    hcal_x = sbs_hcal_clus_x[0];
    hcal_y = sbs_hcal_clus_y[0];
    hcal_row = sbs_hcal_clus_row[0];
    hcal_col = sbs_hcal_clus_col[0];
    hcal_id = sbs_hcal_clus_id[0];
    //energy weighted adc,tdc times of cluster
    hcal_adctime = sbs_hcal_clus_adctime[0];
    hcal_tdctime = sbs_hcal_clus_tdctime[0];
    //adc,tdc time of highest energy block in main cluster
    hcal_atimeblk = sbs_hcal_clus_atimeblk[0];
    hcal_tdctimeblk = sbs_hcal_clus_tdctimeblk[0];
    
    
    //hcal blks in main cluster
    int nhcal = *Ndata_sbs_hcal_clus_blk_id;
    hcal_blk_e.clear();
    hcal_blk_x.clear();
    hcal_blk_y.clear();
    hcal_blk_row.clear();
    hcal_blk_col.clear();
    hcal_blk_id.clear();
    hcal_blk_adctime.clear();
    hcal_blk_tdctime.clear();
    for(int i=0; i<nhcal; i++){
      hcal_blk_e.push_back(sbs_hcal_clus_blk_e[i]);
      hcal_blk_x.push_back(sbs_hcal_clus_blk_x[i]);
      hcal_blk_y.push_back(sbs_hcal_clus_blk_y[i]);
      hcal_blk_row.push_back(sbs_hcal_clus_blk_row[i]);
      hcal_blk_col.push_back(sbs_hcal_clus_blk_col[i]);
      hcal_blk_id.push_back(sbs_hcal_clus_blk_id[i]);
      hcal_blk_adctime.push_back(sbs_hcal_clus_blk_atime[i]);
      hcal_blk_tdctime.push_back(sbs_hcal_clus_blk_tdctime[i]);
    }
    
    //timing variables
    ps_atime = *bb_ps_atimeblk;
    sh_atime = *bb_sh_atimeblk;
    coin_atime = hcal_adctime - sh_atime;
    //int hcal_rownum = (int) sbs_hcal_clus_rowblk;
    //coin_atime_corr = coin_atime - hcal_atcoin_means[hcal_rownum];
    coin_time = hcal_tdctime - bb_hodotdc_clus_bar_tdc_meantime[0];
    //int hodo_barid = (int) bb_hodotdc_clus_bar_tdc_id[0];
    //coin_time_corr = coin_time - hodo_tcoin_means[hodo_barid]; 
    hodo_nclus = *Ndata_bb_hodotdc_clus_bar_tdc_id;
    int nhodo = *Ndata_bb_hodotdc_clus_bar_tdc_id;
    hodo_nbar = nhodo;
    hodo_tmean.clear();
    hodo_tdiff.clear();
    hodo_totmean.clear();
    hodo_tleft.clear();
    hodo_tright.clear();
    hodo_totleft.clear();
    hodo_totright.clear();
    hodo_timehitpos.clear();
    hodo_bar.clear();
    
    for (int i=0; i<nhodo; i++){
      hodo_tmean.push_back(bb_hodotdc_clus_bar_tdc_meantime[i]);
      hodo_tdiff.push_back(bb_hodotdc_clus_bar_tdc_timediff[i]);
      hodo_totmean.push_back(bb_hodotdc_clus_bar_tdc_meantot[i]);
      hodo_tleft.push_back(bb_hodotdc_clus_bar_tdc_tleft[i]);
      hodo_tright.push_back(bb_hodotdc_clus_bar_tdc_tright[i]);
      hodo_totleft.push_back(bb_hodotdc_clus_bar_tdc_totleft[i]);
      hodo_totright.push_back(bb_hodotdc_clus_bar_tdc_totright[i]);
      hodo_timehitpos.push_back(bb_hodotdc_clus_bar_tdc_timehitpos[i]);
      hodo_bar.push_back(bb_hodotdc_clus_bar_tdc_id[i]);
    }
    
    tdctrig_id.clear();
    tdctrig.clear();
    int ntdctrig = *Ndata_bb_tdctrig_tdcelemID;
    bb_trigtime=-999;
    bb_rftime=-999;
    for (int i=0; i<ntdctrig; i++){
      tdctrig_id.push_back(bb_tdctrig_tdcelemID[i]);
      tdctrig.push_back(bb_tdctrig_tdc[i]);
      if(bb_tdctrig_tdcelemID[i]==4) bb_rftime=bb_tdctrig_tdc[i];
      if(bb_tdctrig_tdcelemID[i]==5) bb_trigtime=bb_tdctrig_tdc[i];
    }
    if(bb_rftime == -999) no_bbrf++;
    if(bb_trigtime == -999) no_bbtrig++;
    
    hcal_refid.clear();
    hcal_ref.clear();
    int nhcalref = *Ndata_sbs_hcal_Ref_tdcelemID;
    sbs_trigtime=-999;
    sbs_rftime=-999;
    for (int i=0; i<nhcalref; i++){
      hcal_refid.push_back(sbs_hcal_Ref_tdcelemID[i]);
      hcal_ref.push_back(sbs_hcal_Ref_tdc[i]);
      if(sbs_hcal_Ref_tdcelemID[i]==2) sbs_trigtime=sbs_hcal_Ref_tdc[i];
      if(sbs_hcal_Ref_tdcelemID[i]==3) sbs_rftime=sbs_hcal_Ref_tdc[i];
    }
    if(sbs_rftime == -999) no_sbsrf++;
    if(sbs_trigtime == -999) no_sbstrig++;

    //if(sbs_rftime == -999 || sbs_trigtime == -999 || bb_rftime == -999 || bb_trigtime == -999) continue;
    
    gr_adc = *bb_grinch_tdc_clus_adc;
    gr_size = *bb_grinch_tdc_clus_size;
    gr_tmean = *bb_grinch_tdc_clus_t_mean;
    gr_totmean = *bb_grinch_tdc_clus_tot_mean;
    gr_trindex = *bb_grinch_tdc_clus_trackindex;
    gr_xmean = *bb_grinch_tdc_clus_x_mean;
    gr_ymean = *bb_grinch_tdc_clus_y_mean;
    //gr_nclus = *bb_grinch_tdc_nclus;
    //gr_ngoodhits = *bb_grinch_tdc_ngoodhits;
    //gr_ntrackmatch = *bb_grinch_tdc_ntrackmatch;
  
    //global hall coords
    //x is -ve vertical
    //y positive left
    //z in dir beam
    //theta polar therefore dx/dz
    //phi dy/dz - in same plane as angle of spectrometers ("in-plane")
    //vertex corrections
    TVector3 vertex(bb_vx,bb_vy,bb_vz);
    double sintersect = (hcal_origin-vertex).Dot( hcal_zaxis )/qdir.Dot( hcal_zaxis );
    TVector3 hcal_intersect = vertex + sintersect * qdir;
    
    //pred_x = hcal_intersect.Dot( hcal_xaxis );
    //pred_y = hcal_intersect.Dot( hcal_yaxis );
    //correction from sean:
    pred_x = (hcal_intersect-hcal_origin).Dot( hcal_xaxis );
    pred_y = (hcal_intersect-hcal_origin).Dot( hcal_yaxis );

    sbs_pathl = hcal_intersect.Mag();
    //sbs_pathl_defl =(dSBS-(Dgap/2.0)) + proton_deflection_pathl;

    dx = hcal_x - pred_x;
    dy = hcal_y - pred_y;
    
    pred_mom = Rp4.Rho();
    dev = (sqrt(dx*dx + dy*dy) / hcal_dist) * pred_mom; //units gev
    
    Double_t IHWP = *IGL1I00OD16_16;
    //moller asymmetry flipped for kin3,4 compared to kin2
    if(fabs(IHWP) != 1) IHWP = 0;
    if (kin_no == 3 || kin_no == 4 ) IHWP *= -1;
    helicity = IHWP * (*scalhel_hel); //helicity from tree
    
    tree->Fill();
  }//End of event loop
  cout << "End of event loop" << endl;
  cout << "Total Good Events: " << nev_good << endl;
  cout << "Events with no BB Trig :" << no_bbtrig << endl;
  cout << "Events with no BB RF :" << no_bbrf << endl;
  cout << "Events with no SBS Trig :" << no_sbstrig << endl;
  cout << "Events with no SBS RF :" << no_sbsrf << endl;
  cout <<"Saving rootfiles" << endl;
  outfile->cd();
  outfile->Write();
  outfile->Close();
  t.Stop();
  cout << "Total Time: ";
  t.Print();
  cout << " s" << endl;
  
}//End of macro
  
