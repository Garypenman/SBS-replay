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
    rootfile_dir = "/cache/halla/sbs/prod/GEnII/pass1/GEN2/H2/SBS100/rootfiles/";
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

  //my replays
  rootfile_dir = Form("/volatile/halla/sbs/gpenman/GEN_REPLAY/rootfiles/");
  
  gSystem->Setenv("ROOTFILE_DIR",rootfile_dir); 
  
  C->Add("$ROOTFILE_DIR/e1209016_fullreplay_*.root");
  
  SetKinematics(kin_no);
  
  double hcal_voffset = 0.0;
  double hcal_hoffset = 0.0;
  
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
  
  TTreeReaderValue<Double_t> e_kine_Q2 = {rd,"e.kine.Q2"}; 
  TTreeReaderValue<Double_t> e_kine_W2 = {rd,"e.kine.W2"}; 
  TTreeReaderValue<Double_t> e_kine_epsilon = {rd,"e.kine.epsilon"}; 
  TTreeReaderValue<Double_t> e_kine_nu = {rd,"e.kine.nu"}; 
  TTreeReaderValue<Double_t> e_kine_x_bj ={rd,"e.kine.x_bj"}; 
  
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
  
  TTreeReaderValue<Double_t> bb_gem_trigtime = {rd, "bb.gem.trigtime"};
  TTreeReaderArray<Double_t> bb_gem_track_nhits = {rd, "bb.gem.track.nhits"};
  
  /*
  TTreeReaderArray<Double_t> sbs_tr_x = {rd,"sbs.tr.x"};
  TTreeReaderArray<Double_t> sbs_tr_y = {rd,"sbs.tr.y"};
  TTreeReaderArray<Double_t> sbs_tr_px = {rd,"sbs.tr.px"};
  TTreeReaderArray<Double_t> sbs_tr_py = {rd,"sbs.tr.py"};
  TTreeReaderArray<Double_t> sbs_tr_pz = {rd,"sbs.tr.pz"};
  TTreeReaderArray<Double_t> sbs_tr_p = {rd,"sbs.tr.p"};
  TTreeReaderArray<Double_t> sbs_tr_th = {rd,"sbs.tr.th"};
  TTreeReaderArray<Double_t> sbs_tr_ph = {rd,"sbs.tr.ph"};
  TTreeReaderArray<Double_t> sbs_tr_vx = {rd,"sbs.tr.vx"};
  TTreeReaderArray<Double_t> sbs_tr_vy = {rd,"sbs.tr.vy"};
  TTreeReaderArray<Double_t> sbs_tr_vz = {rd,"sbs.tr.vz"};

  TTreeReaderValue<Double_t> sbs_tr_n = {rd,"sbs.tr.n"};
 
  TTreeReaderArray<Double_t> sbs_tr_tg_th = {rd,"sbs.tr.tg_th"};
  TTreeReaderArray<Double_t> sbs_tr_tg_ph = {rd,"sbs.tr.tg_ph"};
  TTreeReaderArray<Double_t> sbs_tr_tg_y = {rd,"sbs.tr.tg_y"};
  */
  
   
  TTreeReaderValue<Double_t> bb_ps_x = {rd,"bb.ps.x"};
  TTreeReaderValue<Double_t> bb_sh_x = {rd,"bb.sh.x"};
  TTreeReaderValue<Double_t> bb_ps_y = {rd,"bb.ps.y"};
  TTreeReaderValue<Double_t> bb_sh_y = {rd,"bb.sh.y"};
  TTreeReaderValue<Double_t> bb_ps_e = {rd,"bb.ps.e"};
  TTreeReaderValue<Double_t> bb_sh_e = {rd,"bb.sh.e"};
  TTreeReaderValue<Double_t> bb_ps_atimeblk = {rd, "bb.ps.atimeblk"};
  TTreeReaderValue<Double_t> bb_sh_atimeblk = {rd, "bb.sh.atimeblk"};
  
  TTreeReaderValue<Double_t> sbs_hcal_atimeblk = {rd, "sbs.hcal.atimeblk"};
  TTreeReaderValue<Double_t> sbs_hcal_tdctimeblk = {rd, "sbs.hcal.tdctimeblk"};
  TTreeReaderValue<Double_t> sbs_hcal_x = {rd,"sbs.hcal.x"};
  TTreeReaderValue<Double_t> sbs_hcal_y = {rd,"sbs.hcal.y"};
  TTreeReaderValue<Double_t> sbs_hcal_e = {rd,"sbs.hcal.e"};
  TTreeReaderValue<Double_t> sbs_hcal_eblk = {rd,"sbs.hcal.eblk"};
  TTreeReaderValue<Double_t> sbs_hcal_nclus = {rd,"sbs.hcal.nclus"};
  TTreeReaderValue<Double_t> sbs_hcal_rowblk = {rd, "sbs.hcal.rowblk"};
  TTreeReaderValue<Double_t> sbs_hcal_colblk = {rd, "sbs.hcal.colblk"};
  TTreeReaderValue<Double_t> sbs_hcal_idblk = {rd, "sbs.hcal.idblk"};
  
  TTreeReaderValue<ULong64_t> fEvtHdr_fEvtTime = {rd, "fEvtHdr.fEvtTime"};
  TTreeReaderValue<unsigned int> fEvtHdr_fEvtNum = {rd, "fEvtHdr.fEvtNum"};
  TTreeReaderValue<unsigned int> fEvtHdr_fEvtType = {rd, "fEvtHdr.fEvtType"};
  TTreeReaderValue<unsigned int> fEvtHdr_fEvtLen = {rd, "fEvtHdr.fEvtLen"};
  TTreeReaderValue<Int_t> fEvtHdr_fHelicity = {rd, "fEvtHdr.fHelicity"};
  TTreeReaderValue<unsigned int> fEvtHdr_fTrigBits = {rd, "fEvtHdr.fTrigBits"};
  TTreeReaderValue<unsigned int> fEvtHdr_fRun = {rd, "fEvtHdr.fRun"};
  
  TTreeReaderValue<Int_t> Ndata_bb_hodotdc_clus_bar_tdc_id = {rd, "Ndata.bb.hodotdc.clus.bar.tdc.id"};
  TTreeReaderArray<Double_t> bb_hodotdc_clus_bar_tdc_id = {rd, "bb.hodotdc.clus.bar.tdc.id"};
  TTreeReaderArray<Double_t> bb_hodotdc_clus_bar_tdc_itrack = {rd, "bb.hodotdc.clus.bar.tdc.itrack"};
  TTreeReaderArray<Double_t> bb_hodotdc_clus_bar_tdc_meantime = {rd, "bb.hodotdc.clus.bar.tdc.meantime"};
  TTreeReaderArray<Double_t> bb_hodotdc_clus_bar_tdc_tleft = {rd, "bb.hodotdc.clus.bar.tdc.tleft"};
  TTreeReaderArray<Double_t> bb_hodotdc_clus_bar_tdc_tright = {rd, "bb.hodotdc.clus.bar.tdc.tright"};
  TTreeReaderArray<Double_t> bb_hodotdc_clus_bar_tdc_totleft = {rd, "bb.hodotdc.clus.bar.tdc.totleft"};
  TTreeReaderArray<Double_t> bb_hodotdc_clus_bar_tdc_totright = {rd, "bb.hodotdc.clus.bar.tdc.totright"};
  TTreeReaderArray<Double_t> bb_hodotdc_clus_bar_tdc_meantot = {rd, "bb.hodotdc.clus.bar.tdc.meantot"};
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
  
  
  //paste any variables from below that you need otherwise keep them commented for speed
  //to quickly jump down search for string "//end of reader"
  /*
  TTreeReaderValue<Double_t> rast12x = {rd, "rast12x"};
   TTreeReaderValue<Double_t> rast12y = {rd, "rast12y"};
   TTreeReaderValue<Double_t> rast12xmm = {rd, "rast12xmm"};
   TTreeReaderValue<Double_t> rast12ymm = {rd, "rast12ymm"};
   TTreeReaderValue<Double_t> ax_new = {rd, "ax_new"};
   TTreeReaderValue<Double_t> bx_new = {rd, "bx_new"};
   TTreeReaderValue<Double_t> ay_new = {rd, "ay_new"};
   TTreeReaderValue<Double_t> by_new = {rd, "by_new"};
   TTreeReaderValue<Double_t> targx = {rd, "targx"};
   TTreeReaderValue<Double_t> targy = {rd, "targy"};
   TTreeReaderValue<Double_t> scalhel_true_hel = {rd, "scalhel.true_hel"};
   TTreeReaderValue<Double_t> scalhel_hel_stable = {rd, "scalhel.hel_stable"};
   
   TTreeReaderValue<Int_t> Ndata_bb_bbtrig_a_amp_p = {rd, "Ndata.bb.bbtrig.a_amp_p"};
   TTreeReaderArray<Double_t> bb_bbtrig_a_amp_p = {rd, "bb.bbtrig.a_amp_p"};
   TTreeReaderValue<Int_t> Ndata_bb_bbtrig_a_time = {rd, "Ndata.bb.bbtrig.a_time"};
   TTreeReaderArray<Double_t> bb_bbtrig_a_time = {rd, "bb.bbtrig.a_time"};
   TTreeReaderValue<Int_t> Ndata_bb_bbtrig_adcelemID = {rd, "Ndata.bb.bbtrig.adcelemID"};
   TTreeReaderArray<Double_t> bb_bbtrig_adcelemID = {rd, "bb.bbtrig.adcelemID"};
   TTreeReaderValue<Int_t> Ndata_bb_bbtrig_tdc = {rd, "Ndata.bb.bbtrig.tdc"};
   TTreeReaderArray<Double_t> bb_bbtrig_tdc = {rd, "bb.bbtrig.tdc"};
   TTreeReaderValue<Int_t> Ndata_bb_bbtrig_tdcelemID = {rd, "Ndata.bb.bbtrig.tdcelemID"};
   TTreeReaderArray<Double_t> bb_bbtrig_tdcelemID = {rd, "bb.bbtrig.tdcelemID"};
   TTreeReaderValue<Int_t> Ndata_bb_eps_over_etot = {rd, "Ndata.bb.eps_over_etot"};
   TTreeReaderArray<Double_t> bb_eps_over_etot = {rd, "bb.eps_over_etot"};
   TTreeReaderValue<Int_t> Ndata_bb_etot_over_p = {rd, "Ndata.bb.etot_over_p"};
   TTreeReaderArray<Double_t> bb_etot_over_p = {rd, "bb.etot_over_p"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_ADCU = {rd, "Ndata.bb.gem.hit.ADCU"};
   TTreeReaderArray<Double_t> bb_gem_hit_ADCU = {rd, "bb.gem.hit.ADCU"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_ADCU_deconv = {rd, "Ndata.bb.gem.hit.ADCU_deconv"};
   TTreeReaderArray<Double_t> bb_gem_hit_ADCU_deconv = {rd, "bb.gem.hit.ADCU_deconv"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_ADCV = {rd, "Ndata.bb.gem.hit.ADCV"};
   TTreeReaderArray<Double_t> bb_gem_hit_ADCV = {rd, "bb.gem.hit.ADCV"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_ADCV_deconv = {rd, "Ndata.bb.gem.hit.ADCV_deconv"};
   TTreeReaderArray<Double_t> bb_gem_hit_ADCV_deconv = {rd, "bb.gem.hit.ADCV_deconv"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_ADCasym = {rd, "Ndata.bb.gem.hit.ADCasym"};
   TTreeReaderArray<Double_t> bb_gem_hit_ADCasym = {rd, "bb.gem.hit.ADCasym"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_ADCasym_deconv = {rd, "Ndata.bb.gem.hit.ADCasym_deconv"};
   TTreeReaderArray<Double_t> bb_gem_hit_ADCasym_deconv = {rd, "bb.gem.hit.ADCasym_deconv"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_ADCavg = {rd, "Ndata.bb.gem.hit.ADCavg"};
   TTreeReaderArray<Double_t> bb_gem_hit_ADCavg = {rd, "bb.gem.hit.ADCavg"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_ADCavg_deconv = {rd, "Ndata.bb.gem.hit.ADCavg_deconv"};
   TTreeReaderArray<Double_t> bb_gem_hit_ADCavg_deconv = {rd, "bb.gem.hit.ADCavg_deconv"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_ADCfrac0_Umax = {rd, "Ndata.bb.gem.hit.ADCfrac0_Umax"};
   TTreeReaderArray<Double_t> bb_gem_hit_ADCfrac0_Umax = {rd, "bb.gem.hit.ADCfrac0_Umax"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_ADCfrac0_Vmax = {rd, "Ndata.bb.gem.hit.ADCfrac0_Vmax"};
   TTreeReaderArray<Double_t> bb_gem_hit_ADCfrac0_Vmax = {rd, "bb.gem.hit.ADCfrac0_Vmax"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_ADCfrac1_Umax = {rd, "Ndata.bb.gem.hit.ADCfrac1_Umax"};
   TTreeReaderArray<Double_t> bb_gem_hit_ADCfrac1_Umax = {rd, "bb.gem.hit.ADCfrac1_Umax"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_ADCfrac1_Vmax = {rd, "Ndata.bb.gem.hit.ADCfrac1_Vmax"};
   TTreeReaderArray<Double_t> bb_gem_hit_ADCfrac1_Vmax = {rd, "bb.gem.hit.ADCfrac1_Vmax"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_ADCfrac2_Umax = {rd, "Ndata.bb.gem.hit.ADCfrac2_Umax"};
   TTreeReaderArray<Double_t> bb_gem_hit_ADCfrac2_Umax = {rd, "bb.gem.hit.ADCfrac2_Umax"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_ADCfrac2_Vmax = {rd, "Ndata.bb.gem.hit.ADCfrac2_Vmax"};
   TTreeReaderArray<Double_t> bb_gem_hit_ADCfrac2_Vmax = {rd, "bb.gem.hit.ADCfrac2_Vmax"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_ADCfrac3_Umax = {rd, "Ndata.bb.gem.hit.ADCfrac3_Umax"};
   TTreeReaderArray<Double_t> bb_gem_hit_ADCfrac3_Umax = {rd, "bb.gem.hit.ADCfrac3_Umax"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_ADCfrac3_Vmax = {rd, "Ndata.bb.gem.hit.ADCfrac3_Vmax"};
   TTreeReaderArray<Double_t> bb_gem_hit_ADCfrac3_Vmax = {rd, "bb.gem.hit.ADCfrac3_Vmax"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_ADCfrac4_Umax = {rd, "Ndata.bb.gem.hit.ADCfrac4_Umax"};
   TTreeReaderArray<Double_t> bb_gem_hit_ADCfrac4_Umax = {rd, "bb.gem.hit.ADCfrac4_Umax"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_ADCfrac4_Vmax = {rd, "Ndata.bb.gem.hit.ADCfrac4_Vmax"};
   TTreeReaderArray<Double_t> bb_gem_hit_ADCfrac4_Vmax = {rd, "bb.gem.hit.ADCfrac4_Vmax"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_ADCfrac5_Umax = {rd, "Ndata.bb.gem.hit.ADCfrac5_Umax"};
   TTreeReaderArray<Double_t> bb_gem_hit_ADCfrac5_Umax = {rd, "bb.gem.hit.ADCfrac5_Umax"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_ADCfrac5_Vmax = {rd, "Ndata.bb.gem.hit.ADCfrac5_Vmax"};
   TTreeReaderArray<Double_t> bb_gem_hit_ADCfrac5_Vmax = {rd, "bb.gem.hit.ADCfrac5_Vmax"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_ADCmaxcomboUclust_deconv = {rd, "Ndata.bb.gem.hit.ADCmaxcomboUclust_deconv"};
   TTreeReaderArray<Double_t> bb_gem_hit_ADCmaxcomboUclust_deconv = {rd, "bb.gem.hit.ADCmaxcomboUclust_deconv"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_ADCmaxcomboVclust_deconv = {rd, "Ndata.bb.gem.hit.ADCmaxcomboVclust_deconv"};
   TTreeReaderArray<Double_t> bb_gem_hit_ADCmaxcomboVclust_deconv = {rd, "bb.gem.hit.ADCmaxcomboVclust_deconv"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_ADCmaxsampU = {rd, "Ndata.bb.gem.hit.ADCmaxsampU"};
   TTreeReaderArray<Double_t> bb_gem_hit_ADCmaxsampU = {rd, "bb.gem.hit.ADCmaxsampU"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_ADCmaxsampUclust = {rd, "Ndata.bb.gem.hit.ADCmaxsampUclust"};
   TTreeReaderArray<Double_t> bb_gem_hit_ADCmaxsampUclust = {rd, "bb.gem.hit.ADCmaxsampUclust"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_ADCmaxsampUclust_deconv = {rd, "Ndata.bb.gem.hit.ADCmaxsampUclust_deconv"};
   TTreeReaderArray<Double_t> bb_gem_hit_ADCmaxsampUclust_deconv = {rd, "bb.gem.hit.ADCmaxsampUclust_deconv"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_ADCmaxsampV = {rd, "Ndata.bb.gem.hit.ADCmaxsampV"};
   TTreeReaderArray<Double_t> bb_gem_hit_ADCmaxsampV = {rd, "bb.gem.hit.ADCmaxsampV"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_ADCmaxsampVclust = {rd, "Ndata.bb.gem.hit.ADCmaxsampVclust"};
   TTreeReaderArray<Double_t> bb_gem_hit_ADCmaxsampVclust = {rd, "bb.gem.hit.ADCmaxsampVclust"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_ADCmaxsampVclust_deconv = {rd, "Ndata.bb.gem.hit.ADCmaxsampVclust_deconv"};
   TTreeReaderArray<Double_t> bb_gem_hit_ADCmaxsampVclust_deconv = {rd, "bb.gem.hit.ADCmaxsampVclust_deconv"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_ADCmaxstripU = {rd, "Ndata.bb.gem.hit.ADCmaxstripU"};
   TTreeReaderArray<Double_t> bb_gem_hit_ADCmaxstripU = {rd, "bb.gem.hit.ADCmaxstripU"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_ADCmaxstripV = {rd, "Ndata.bb.gem.hit.ADCmaxstripV"};
   TTreeReaderArray<Double_t> bb_gem_hit_ADCmaxstripV = {rd, "bb.gem.hit.ADCmaxstripV"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_BUILD_ALL_SAMPLES_U = {rd, "Ndata.bb.gem.hit.BUILD_ALL_SAMPLES_U"};
   TTreeReaderArray<Double_t> bb_gem_hit_BUILD_ALL_SAMPLES_U = {rd, "bb.gem.hit.BUILD_ALL_SAMPLES_U"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_BUILD_ALL_SAMPLES_V = {rd, "Ndata.bb.gem.hit.BUILD_ALL_SAMPLES_V"};
   TTreeReaderArray<Double_t> bb_gem_hit_BUILD_ALL_SAMPLES_V = {rd, "bb.gem.hit.BUILD_ALL_SAMPLES_V"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_CM_GOOD_U = {rd, "Ndata.bb.gem.hit.CM_GOOD_U"};
   TTreeReaderArray<Double_t> bb_gem_hit_CM_GOOD_U = {rd, "bb.gem.hit.CM_GOOD_U"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_CM_GOOD_V = {rd, "Ndata.bb.gem.hit.CM_GOOD_V"};
   TTreeReaderArray<Double_t> bb_gem_hit_CM_GOOD_V = {rd, "bb.gem.hit.CM_GOOD_V"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_DeconvADC0_Umax = {rd, "Ndata.bb.gem.hit.DeconvADC0_Umax"};
   TTreeReaderArray<Double_t> bb_gem_hit_DeconvADC0_Umax = {rd, "bb.gem.hit.DeconvADC0_Umax"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_DeconvADC0_Vmax = {rd, "Ndata.bb.gem.hit.DeconvADC0_Vmax"};
   TTreeReaderArray<Double_t> bb_gem_hit_DeconvADC0_Vmax = {rd, "bb.gem.hit.DeconvADC0_Vmax"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_DeconvADC1_Umax = {rd, "Ndata.bb.gem.hit.DeconvADC1_Umax"};
   TTreeReaderArray<Double_t> bb_gem_hit_DeconvADC1_Umax = {rd, "bb.gem.hit.DeconvADC1_Umax"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_DeconvADC1_Vmax = {rd, "Ndata.bb.gem.hit.DeconvADC1_Vmax"};
   TTreeReaderArray<Double_t> bb_gem_hit_DeconvADC1_Vmax = {rd, "bb.gem.hit.DeconvADC1_Vmax"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_DeconvADC2_Umax = {rd, "Ndata.bb.gem.hit.DeconvADC2_Umax"};
   TTreeReaderArray<Double_t> bb_gem_hit_DeconvADC2_Umax = {rd, "bb.gem.hit.DeconvADC2_Umax"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_DeconvADC2_Vmax = {rd, "Ndata.bb.gem.hit.DeconvADC2_Vmax"};
   TTreeReaderArray<Double_t> bb_gem_hit_DeconvADC2_Vmax = {rd, "bb.gem.hit.DeconvADC2_Vmax"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_DeconvADC3_Umax = {rd, "Ndata.bb.gem.hit.DeconvADC3_Umax"};
   TTreeReaderArray<Double_t> bb_gem_hit_DeconvADC3_Umax = {rd, "bb.gem.hit.DeconvADC3_Umax"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_DeconvADC3_Vmax = {rd, "Ndata.bb.gem.hit.DeconvADC3_Vmax"};
   TTreeReaderArray<Double_t> bb_gem_hit_DeconvADC3_Vmax = {rd, "bb.gem.hit.DeconvADC3_Vmax"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_DeconvADC4_Umax = {rd, "Ndata.bb.gem.hit.DeconvADC4_Umax"};
   TTreeReaderArray<Double_t> bb_gem_hit_DeconvADC4_Umax = {rd, "bb.gem.hit.DeconvADC4_Umax"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_DeconvADC4_Vmax = {rd, "Ndata.bb.gem.hit.DeconvADC4_Vmax"};
   TTreeReaderArray<Double_t> bb_gem_hit_DeconvADC4_Vmax = {rd, "bb.gem.hit.DeconvADC4_Vmax"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_DeconvADC5_Umax = {rd, "Ndata.bb.gem.hit.DeconvADC5_Umax"};
   TTreeReaderArray<Double_t> bb_gem_hit_DeconvADC5_Umax = {rd, "bb.gem.hit.DeconvADC5_Umax"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_DeconvADC5_Vmax = {rd, "Ndata.bb.gem.hit.DeconvADC5_Vmax"};
   TTreeReaderArray<Double_t> bb_gem_hit_DeconvADC5_Vmax = {rd, "bb.gem.hit.DeconvADC5_Vmax"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_DeconvADCmaxcomboU = {rd, "Ndata.bb.gem.hit.DeconvADCmaxcomboU"};
   TTreeReaderArray<Double_t> bb_gem_hit_DeconvADCmaxcomboU = {rd, "bb.gem.hit.DeconvADCmaxcomboU"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_DeconvADCmaxcomboV = {rd, "Ndata.bb.gem.hit.DeconvADCmaxcomboV"};
   TTreeReaderArray<Double_t> bb_gem_hit_DeconvADCmaxcomboV = {rd, "bb.gem.hit.DeconvADCmaxcomboV"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_DeconvADCmaxsampU = {rd, "Ndata.bb.gem.hit.DeconvADCmaxsampU"};
   TTreeReaderArray<Double_t> bb_gem_hit_DeconvADCmaxsampU = {rd, "bb.gem.hit.DeconvADCmaxsampU"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_DeconvADCmaxsampV = {rd, "Ndata.bb.gem.hit.DeconvADCmaxsampV"};
   TTreeReaderArray<Double_t> bb_gem_hit_DeconvADCmaxsampV = {rd, "bb.gem.hit.DeconvADCmaxsampV"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_DeconvADCmaxstripU = {rd, "Ndata.bb.gem.hit.DeconvADCmaxstripU"};
   TTreeReaderArray<Double_t> bb_gem_hit_DeconvADCmaxstripU = {rd, "bb.gem.hit.DeconvADCmaxstripU"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_DeconvADCmaxstripV = {rd, "Ndata.bb.gem.hit.DeconvADCmaxstripV"};
   TTreeReaderArray<Double_t> bb_gem_hit_DeconvADCmaxstripV = {rd, "bb.gem.hit.DeconvADCmaxstripV"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_ENABLE_CM_U = {rd, "Ndata.bb.gem.hit.ENABLE_CM_U"};
   TTreeReaderArray<Double_t> bb_gem_hit_ENABLE_CM_U = {rd, "bb.gem.hit.ENABLE_CM_U"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_ENABLE_CM_V = {rd, "Ndata.bb.gem.hit.ENABLE_CM_V"};
   TTreeReaderArray<Double_t> bb_gem_hit_ENABLE_CM_V = {rd, "bb.gem.hit.ENABLE_CM_V"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_TSchi2_Umax = {rd, "Ndata.bb.gem.hit.TSchi2_Umax"};
   TTreeReaderArray<Double_t> bb_gem_hit_TSchi2_Umax = {rd, "bb.gem.hit.TSchi2_Umax"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_TSchi2_Vmax = {rd, "Ndata.bb.gem.hit.TSchi2_Vmax"};
   TTreeReaderArray<Double_t> bb_gem_hit_TSchi2_Vmax = {rd, "bb.gem.hit.TSchi2_Vmax"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_TSprob_Umax = {rd, "Ndata.bb.gem.hit.TSprob_Umax"};
   TTreeReaderArray<Double_t> bb_gem_hit_TSprob_Umax = {rd, "bb.gem.hit.TSprob_Umax"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_TSprob_Vmax = {rd, "Ndata.bb.gem.hit.TSprob_Vmax"};
   TTreeReaderArray<Double_t> bb_gem_hit_TSprob_Vmax = {rd, "bb.gem.hit.TSprob_Vmax"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_Tavg = {rd, "Ndata.bb.gem.hit.Tavg"};
   TTreeReaderArray<Double_t> bb_gem_hit_Tavg = {rd, "bb.gem.hit.Tavg"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_Tavg_corr = {rd, "Ndata.bb.gem.hit.Tavg_corr"};
   TTreeReaderArray<Double_t> bb_gem_hit_Tavg_corr = {rd, "bb.gem.hit.Tavg_corr"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_Tavg_deconv = {rd, "Ndata.bb.gem.hit.Tavg_deconv"};
   TTreeReaderArray<Double_t> bb_gem_hit_Tavg_deconv = {rd, "bb.gem.hit.Tavg_deconv"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_Tavg_fit = {rd, "Ndata.bb.gem.hit.Tavg_fit"};
   TTreeReaderArray<Double_t> bb_gem_hit_Tavg_fit = {rd, "bb.gem.hit.Tavg_fit"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_Utime = {rd, "Ndata.bb.gem.hit.Utime"};
   TTreeReaderArray<Double_t> bb_gem_hit_Utime = {rd, "bb.gem.hit.Utime"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_UtimeDeconv = {rd, "Ndata.bb.gem.hit.UtimeDeconv"};
   TTreeReaderArray<Double_t> bb_gem_hit_UtimeDeconv = {rd, "bb.gem.hit.UtimeDeconv"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_UtimeFit = {rd, "Ndata.bb.gem.hit.UtimeFit"};
   TTreeReaderArray<Double_t> bb_gem_hit_UtimeFit = {rd, "bb.gem.hit.UtimeFit"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_UtimeMaxStrip = {rd, "Ndata.bb.gem.hit.UtimeMaxStrip"};
   TTreeReaderArray<Double_t> bb_gem_hit_UtimeMaxStrip = {rd, "bb.gem.hit.UtimeMaxStrip"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_UtimeMaxStripDeconv = {rd, "Ndata.bb.gem.hit.UtimeMaxStripDeconv"};
   TTreeReaderArray<Double_t> bb_gem_hit_UtimeMaxStripDeconv = {rd, "bb.gem.hit.UtimeMaxStripDeconv"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_UtimeMaxStripFit = {rd, "Ndata.bb.gem.hit.UtimeMaxStripFit"};
   TTreeReaderArray<Double_t> bb_gem_hit_UtimeMaxStripFit = {rd, "bb.gem.hit.UtimeMaxStripFit"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_Vtime = {rd, "Ndata.bb.gem.hit.Vtime"};
   TTreeReaderArray<Double_t> bb_gem_hit_Vtime = {rd, "bb.gem.hit.Vtime"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_VtimeDeconv = {rd, "Ndata.bb.gem.hit.VtimeDeconv"};
   TTreeReaderArray<Double_t> bb_gem_hit_VtimeDeconv = {rd, "bb.gem.hit.VtimeDeconv"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_VtimeFit = {rd, "Ndata.bb.gem.hit.VtimeFit"};
   TTreeReaderArray<Double_t> bb_gem_hit_VtimeFit = {rd, "bb.gem.hit.VtimeFit"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_VtimeMaxStrip = {rd, "Ndata.bb.gem.hit.VtimeMaxStrip"};
   TTreeReaderArray<Double_t> bb_gem_hit_VtimeMaxStrip = {rd, "bb.gem.hit.VtimeMaxStrip"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_VtimeMaxStripDeconv = {rd, "Ndata.bb.gem.hit.VtimeMaxStripDeconv"};
   TTreeReaderArray<Double_t> bb_gem_hit_VtimeMaxStripDeconv = {rd, "bb.gem.hit.VtimeMaxStripDeconv"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_VtimeMaxStripFit = {rd, "Ndata.bb.gem.hit.VtimeMaxStripFit"};
   TTreeReaderArray<Double_t> bb_gem_hit_VtimeMaxStripFit = {rd, "bb.gem.hit.VtimeMaxStripFit"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_ccor_clust = {rd, "Ndata.bb.gem.hit.ccor_clust"};
   TTreeReaderArray<Double_t> bb_gem_hit_ccor_clust = {rd, "bb.gem.hit.ccor_clust"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_ccor_clust_deconv = {rd, "Ndata.bb.gem.hit.ccor_clust_deconv"};
   TTreeReaderArray<Double_t> bb_gem_hit_ccor_clust_deconv = {rd, "bb.gem.hit.ccor_clust_deconv"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_ccor_strip = {rd, "Ndata.bb.gem.hit.ccor_strip"};
   TTreeReaderArray<Double_t> bb_gem_hit_ccor_strip = {rd, "bb.gem.hit.ccor_strip"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_ccor_strip_deconv = {rd, "Ndata.bb.gem.hit.ccor_strip_deconv"};
   TTreeReaderArray<Double_t> bb_gem_hit_ccor_strip_deconv = {rd, "bb.gem.hit.ccor_strip_deconv"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_deltat = {rd, "Ndata.bb.gem.hit.deltat"};
   TTreeReaderArray<Double_t> bb_gem_hit_deltat = {rd, "bb.gem.hit.deltat"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_deltat_deconv = {rd, "Ndata.bb.gem.hit.deltat_deconv"};
   TTreeReaderArray<Double_t> bb_gem_hit_deltat_deconv = {rd, "bb.gem.hit.deltat_deconv"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_deltat_fit = {rd, "Ndata.bb.gem.hit.deltat_fit"};
   TTreeReaderArray<Double_t> bb_gem_hit_deltat_fit = {rd, "bb.gem.hit.deltat_fit"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_eresidu = {rd, "Ndata.bb.gem.hit.eresidu"};
   TTreeReaderArray<Double_t> bb_gem_hit_eresidu = {rd, "bb.gem.hit.eresidu"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_eresidv = {rd, "Ndata.bb.gem.hit.eresidv"};
   TTreeReaderArray<Double_t> bb_gem_hit_eresidv = {rd, "bb.gem.hit.eresidv"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_icombomaxUclustDeconv = {rd, "Ndata.bb.gem.hit.icombomaxUclustDeconv"};
   TTreeReaderArray<Double_t> bb_gem_hit_icombomaxUclustDeconv = {rd, "bb.gem.hit.icombomaxUclustDeconv"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_icombomaxUstripDeconv = {rd, "Ndata.bb.gem.hit.icombomaxUstripDeconv"};
   TTreeReaderArray<Double_t> bb_gem_hit_icombomaxUstripDeconv = {rd, "bb.gem.hit.icombomaxUstripDeconv"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_icombomaxVclustDeconv = {rd, "Ndata.bb.gem.hit.icombomaxVclustDeconv"};
   TTreeReaderArray<Double_t> bb_gem_hit_icombomaxVclustDeconv = {rd, "bb.gem.hit.icombomaxVclustDeconv"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_icombomaxVstripDeconv = {rd, "Ndata.bb.gem.hit.icombomaxVstripDeconv"};
   TTreeReaderArray<Double_t> bb_gem_hit_icombomaxVstripDeconv = {rd, "bb.gem.hit.icombomaxVstripDeconv"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_isampmaxUclust = {rd, "Ndata.bb.gem.hit.isampmaxUclust"};
   TTreeReaderArray<Double_t> bb_gem_hit_isampmaxUclust = {rd, "bb.gem.hit.isampmaxUclust"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_isampmaxUclustDeconv = {rd, "Ndata.bb.gem.hit.isampmaxUclustDeconv"};
   TTreeReaderArray<Double_t> bb_gem_hit_isampmaxUclustDeconv = {rd, "bb.gem.hit.isampmaxUclustDeconv"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_isampmaxUstrip = {rd, "Ndata.bb.gem.hit.isampmaxUstrip"};
   TTreeReaderArray<Double_t> bb_gem_hit_isampmaxUstrip = {rd, "bb.gem.hit.isampmaxUstrip"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_isampmaxUstripDeconv = {rd, "Ndata.bb.gem.hit.isampmaxUstripDeconv"};
   TTreeReaderArray<Double_t> bb_gem_hit_isampmaxUstripDeconv = {rd, "bb.gem.hit.isampmaxUstripDeconv"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_isampmaxVclust = {rd, "Ndata.bb.gem.hit.isampmaxVclust"};
   TTreeReaderArray<Double_t> bb_gem_hit_isampmaxVclust = {rd, "bb.gem.hit.isampmaxVclust"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_isampmaxVclustDeconv = {rd, "Ndata.bb.gem.hit.isampmaxVclustDeconv"};
   TTreeReaderArray<Double_t> bb_gem_hit_isampmaxVclustDeconv = {rd, "bb.gem.hit.isampmaxVclustDeconv"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_isampmaxVstrip = {rd, "Ndata.bb.gem.hit.isampmaxVstrip"};
   TTreeReaderArray<Double_t> bb_gem_hit_isampmaxVstrip = {rd, "bb.gem.hit.isampmaxVstrip"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_isampmaxVstripDeconv = {rd, "Ndata.bb.gem.hit.isampmaxVstripDeconv"};
   TTreeReaderArray<Double_t> bb_gem_hit_isampmaxVstripDeconv = {rd, "bb.gem.hit.isampmaxVstripDeconv"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_layer = {rd, "Ndata.bb.gem.hit.layer"};
   TTreeReaderArray<Double_t> bb_gem_hit_layer = {rd, "bb.gem.hit.layer"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_module = {rd, "Ndata.bb.gem.hit.module"};
   TTreeReaderArray<Double_t> bb_gem_hit_module = {rd, "bb.gem.hit.module"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_nstripu = {rd, "Ndata.bb.gem.hit.nstripu"};
   TTreeReaderArray<Double_t> bb_gem_hit_nstripu = {rd, "bb.gem.hit.nstripu"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_nstripv = {rd, "Ndata.bb.gem.hit.nstripv"};
   TTreeReaderArray<Double_t> bb_gem_hit_nstripv = {rd, "bb.gem.hit.nstripv"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_residu = {rd, "Ndata.bb.gem.hit.residu"};
   TTreeReaderArray<Double_t> bb_gem_hit_residu = {rd, "bb.gem.hit.residu"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_residv = {rd, "Ndata.bb.gem.hit.residv"};
   TTreeReaderArray<Double_t> bb_gem_hit_residv = {rd, "bb.gem.hit.residv"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_trackindex = {rd, "Ndata.bb.gem.hit.trackindex"};
   TTreeReaderArray<Double_t> bb_gem_hit_trackindex = {rd, "bb.gem.hit.trackindex"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_u = {rd, "Ndata.bb.gem.hit.u"};
   TTreeReaderArray<Double_t> bb_gem_hit_u = {rd, "bb.gem.hit.u"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_umoment = {rd, "Ndata.bb.gem.hit.umoment"};
   TTreeReaderArray<Double_t> bb_gem_hit_umoment = {rd, "bb.gem.hit.umoment"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_usigma = {rd, "Ndata.bb.gem.hit.usigma"};
   TTreeReaderArray<Double_t> bb_gem_hit_usigma = {rd, "bb.gem.hit.usigma"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_ustriphi = {rd, "Ndata.bb.gem.hit.ustriphi"};
   TTreeReaderArray<Double_t> bb_gem_hit_ustriphi = {rd, "bb.gem.hit.ustriphi"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_ustriplo = {rd, "Ndata.bb.gem.hit.ustriplo"};
   TTreeReaderArray<Double_t> bb_gem_hit_ustriplo = {rd, "bb.gem.hit.ustriplo"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_ustripmax = {rd, "Ndata.bb.gem.hit.ustripmax"};
   TTreeReaderArray<Double_t> bb_gem_hit_ustripmax = {rd, "bb.gem.hit.ustripmax"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_v = {rd, "Ndata.bb.gem.hit.v"};
   TTreeReaderArray<Double_t> bb_gem_hit_v = {rd, "bb.gem.hit.v"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_vmoment = {rd, "Ndata.bb.gem.hit.vmoment"};
   TTreeReaderArray<Double_t> bb_gem_hit_vmoment = {rd, "bb.gem.hit.vmoment"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_vsigma = {rd, "Ndata.bb.gem.hit.vsigma"};
   TTreeReaderArray<Double_t> bb_gem_hit_vsigma = {rd, "bb.gem.hit.vsigma"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_vstriphi = {rd, "Ndata.bb.gem.hit.vstriphi"};
   TTreeReaderArray<Double_t> bb_gem_hit_vstriphi = {rd, "bb.gem.hit.vstriphi"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_vstriplo = {rd, "Ndata.bb.gem.hit.vstriplo"};
   TTreeReaderArray<Double_t> bb_gem_hit_vstriplo = {rd, "bb.gem.hit.vstriplo"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_vstripmax = {rd, "Ndata.bb.gem.hit.vstripmax"};
   TTreeReaderArray<Double_t> bb_gem_hit_vstripmax = {rd, "bb.gem.hit.vstripmax"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_xglobal = {rd, "Ndata.bb.gem.hit.xglobal"};
   TTreeReaderArray<Double_t> bb_gem_hit_xglobal = {rd, "bb.gem.hit.xglobal"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_xlocal = {rd, "Ndata.bb.gem.hit.xlocal"};
   TTreeReaderArray<Double_t> bb_gem_hit_xlocal = {rd, "bb.gem.hit.xlocal"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_yglobal = {rd, "Ndata.bb.gem.hit.yglobal"};
   TTreeReaderArray<Double_t> bb_gem_hit_yglobal = {rd, "bb.gem.hit.yglobal"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_ylocal = {rd, "Ndata.bb.gem.hit.ylocal"};
   TTreeReaderArray<Double_t> bb_gem_hit_ylocal = {rd, "bb.gem.hit.ylocal"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_hit_zglobal = {rd, "Ndata.bb.gem.hit.zglobal"};
   TTreeReaderArray<Double_t> bb_gem_hit_zglobal = {rd, "bb.gem.hit.zglobal"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_n2Dhit_layer = {rd, "Ndata.bb.gem.n2Dhit_layer"};
   TTreeReaderArray<Double_t> bb_gem_n2Dhit_layer = {rd, "bb.gem.n2Dhit_layer"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_nclustu_layer = {rd, "Ndata.bb.gem.nclustu_layer"};
   TTreeReaderArray<Double_t> bb_gem_nclustu_layer = {rd, "bb.gem.nclustu_layer"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_nclustv_layer = {rd, "Ndata.bb.gem.nclustv_layer"};
   TTreeReaderArray<Double_t> bb_gem_nclustv_layer = {rd, "bb.gem.nclustv_layer"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_nstripsu_layer = {rd, "Ndata.bb.gem.nstripsu_layer"};
   TTreeReaderArray<Double_t> bb_gem_nstripsu_layer = {rd, "bb.gem.nstripsu_layer"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_nstripsv_layer = {rd, "Ndata.bb.gem.nstripsv_layer"};
   TTreeReaderArray<Double_t> bb_gem_nstripsv_layer = {rd, "bb.gem.nstripsv_layer"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_track_chi2ndf = {rd, "Ndata.bb.gem.track.chi2ndf"};
   TTreeReaderArray<Double_t> bb_gem_track_chi2ndf = {rd, "bb.gem.track.chi2ndf"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_track_chi2ndf_hitquality = {rd, "Ndata.bb.gem.track.chi2ndf_hitquality"};
   TTreeReaderArray<Double_t> bb_gem_track_chi2ndf_hitquality = {rd, "bb.gem.track.chi2ndf_hitquality"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_track_ngoodhits = {rd, "Ndata.bb.gem.track.ngoodhits"};
   TTreeReaderArray<Double_t> bb_gem_track_ngoodhits = {rd, "bb.gem.track.ngoodhits"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_track_nhits = {rd, "Ndata.bb.gem.track.nhits"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_track_t0 = {rd, "Ndata.bb.gem.track.t0"};
   TTreeReaderArray<Double_t> bb_gem_track_t0 = {rd, "bb.gem.track.t0"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_track_x = {rd, "Ndata.bb.gem.track.x"};
   TTreeReaderArray<Double_t> bb_gem_track_x = {rd, "bb.gem.track.x"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_track_xp = {rd, "Ndata.bb.gem.track.xp"};
   TTreeReaderArray<Double_t> bb_gem_track_xp = {rd, "bb.gem.track.xp"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_track_y = {rd, "Ndata.bb.gem.track.y"};
   TTreeReaderArray<Double_t> bb_gem_track_y = {rd, "bb.gem.track.y"};
   TTreeReaderValue<Int_t> Ndata_bb_gem_track_yp = {rd, "Ndata.bb.gem.track.yp"};
   TTreeReaderArray<Double_t> bb_gem_track_yp = {rd, "bb.gem.track.yp"};
   TTreeReaderValue<Int_t> Ndata_bb_grinch_tdc_allclus_adc = {rd, "Ndata.bb.grinch_tdc.allclus.adc"};
   TTreeReaderArray<Double_t> bb_grinch_tdc_allclus_adc = {rd, "bb.grinch_tdc.allclus.adc"};
   TTreeReaderValue<Int_t> Ndata_bb_grinch_tdc_allclus_size = {rd, "Ndata.bb.grinch_tdc.allclus.size"};
   TTreeReaderArray<Double_t> bb_grinch_tdc_allclus_size = {rd, "bb.grinch_tdc.allclus.size"};
   TTreeReaderValue<Int_t> Ndata_bb_grinch_tdc_allclus_t_mean = {rd, "Ndata.bb.grinch_tdc.allclus.t_mean"};
   TTreeReaderArray<Double_t> bb_grinch_tdc_allclus_t_mean = {rd, "bb.grinch_tdc.allclus.t_mean"};
   TTreeReaderValue<Int_t> Ndata_bb_grinch_tdc_allclus_t_rms = {rd, "Ndata.bb.grinch_tdc.allclus.t_rms"};
   TTreeReaderArray<Double_t> bb_grinch_tdc_allclus_t_rms = {rd, "bb.grinch_tdc.allclus.t_rms"};
   TTreeReaderValue<Int_t> Ndata_bb_grinch_tdc_allclus_tot_mean = {rd, "Ndata.bb.grinch_tdc.allclus.tot_mean"};
   TTreeReaderArray<Double_t> bb_grinch_tdc_allclus_tot_mean = {rd, "bb.grinch_tdc.allclus.tot_mean"};
   TTreeReaderValue<Int_t> Ndata_bb_grinch_tdc_allclus_trackindex = {rd, "Ndata.bb.grinch_tdc.allclus.trackindex"};
   TTreeReaderArray<Double_t> bb_grinch_tdc_allclus_trackindex = {rd, "bb.grinch_tdc.allclus.trackindex"};
   TTreeReaderValue<Int_t> Ndata_bb_grinch_tdc_allclus_x_mean = {rd, "Ndata.bb.grinch_tdc.allclus.x_mean"};
   TTreeReaderArray<Double_t> bb_grinch_tdc_allclus_x_mean = {rd, "bb.grinch_tdc.allclus.x_mean"};
   TTreeReaderValue<Int_t> Ndata_bb_grinch_tdc_allclus_y_mean = {rd, "Ndata.bb.grinch_tdc.allclus.y_mean"};
   TTreeReaderArray<Double_t> bb_grinch_tdc_allclus_y_mean = {rd, "bb.grinch_tdc.allclus.y_mean"};
   TTreeReaderValue<Int_t> Ndata_bb_grinch_tdc_hit_amp = {rd, "Ndata.bb.grinch_tdc.hit.amp"};
   TTreeReaderArray<Double_t> bb_grinch_tdc_hit_amp = {rd, "bb.grinch_tdc.hit.amp"};
   TTreeReaderValue<Int_t> Ndata_bb_grinch_tdc_hit_clustindex = {rd, "Ndata.bb.grinch_tdc.hit.clustindex"};
   TTreeReaderArray<Double_t> bb_grinch_tdc_hit_clustindex = {rd, "bb.grinch_tdc.hit.clustindex"};
   TTreeReaderValue<Int_t> Ndata_bb_grinch_tdc_hit_col = {rd, "Ndata.bb.grinch_tdc.hit.col"};
   TTreeReaderArray<Double_t> bb_grinch_tdc_hit_col = {rd, "bb.grinch_tdc.hit.col"};
   TTreeReaderValue<Int_t> Ndata_bb_grinch_tdc_hit_pmtnum = {rd, "Ndata.bb.grinch_tdc.hit.pmtnum"};
   TTreeReaderArray<Double_t> bb_grinch_tdc_hit_pmtnum = {rd, "bb.grinch_tdc.hit.pmtnum"};
   TTreeReaderValue<Int_t> Ndata_bb_grinch_tdc_hit_row = {rd, "Ndata.bb.grinch_tdc.hit.row"};
   TTreeReaderArray<Double_t> bb_grinch_tdc_hit_row = {rd, "bb.grinch_tdc.hit.row"};
   TTreeReaderValue<Int_t> Ndata_bb_grinch_tdc_hit_time = {rd, "Ndata.bb.grinch_tdc.hit.time"};
   TTreeReaderArray<Double_t> bb_grinch_tdc_hit_time = {rd, "bb.grinch_tdc.hit.time"};
   TTreeReaderValue<Int_t> Ndata_bb_grinch_tdc_hit_trackindex = {rd, "Ndata.bb.grinch_tdc.hit.trackindex"};
   TTreeReaderArray<Double_t> bb_grinch_tdc_hit_trackindex = {rd, "bb.grinch_tdc.hit.trackindex"};
   TTreeReaderValue<Int_t> Ndata_bb_grinch_tdc_hit_xhit = {rd, "Ndata.bb.grinch_tdc.hit.xhit"};
   TTreeReaderArray<Double_t> bb_grinch_tdc_hit_xhit = {rd, "bb.grinch_tdc.hit.xhit"};
   TTreeReaderValue<Int_t> Ndata_bb_grinch_tdc_hit_yhit = {rd, "Ndata.bb.grinch_tdc.hit.yhit"};
   TTreeReaderArray<Double_t> bb_grinch_tdc_hit_yhit = {rd, "bb.grinch_tdc.hit.yhit"};
   TTreeReaderValue<Int_t> Ndata_bb_hodotdc_clus_bar_tdc_id = {rd, "Ndata.bb.hodotdc.clus.bar.tdc.id"};
   TTreeReaderArray<Double_t> bb_hodotdc_clus_bar_tdc_id = {rd, "bb.hodotdc.clus.bar.tdc.id"};
   TTreeReaderValue<Int_t> Ndata_bb_hodotdc_clus_bar_tdc_itrack = {rd, "Ndata.bb.hodotdc.clus.bar.tdc.itrack"};
   TTreeReaderArray<Double_t> bb_hodotdc_clus_bar_tdc_itrack = {rd, "bb.hodotdc.clus.bar.tdc.itrack"};
   TTreeReaderValue<Int_t> Ndata_bb_hodotdc_clus_bar_tdc_meantime = {rd, "Ndata.bb.hodotdc.clus.bar.tdc.meantime"};
   TTreeReaderArray<Double_t> bb_hodotdc_clus_bar_tdc_meantime = {rd, "bb.hodotdc.clus.bar.tdc.meantime"};
   TTreeReaderValue<Int_t> Ndata_bb_hodotdc_clus_bar_tdc_meantot = {rd, "Ndata.bb.hodotdc.clus.bar.tdc.meantot"};
   TTreeReaderArray<Double_t> bb_hodotdc_clus_bar_tdc_meantot = {rd, "bb.hodotdc.clus.bar.tdc.meantot"};
   TTreeReaderValue<Int_t> Ndata_bb_hodotdc_clus_bar_tdc_timediff = {rd, "Ndata.bb.hodotdc.clus.bar.tdc.timediff"};
   TTreeReaderArray<Double_t> bb_hodotdc_clus_bar_tdc_timediff = {rd, "bb.hodotdc.clus.bar.tdc.timediff"};
   TTreeReaderValue<Int_t> Ndata_bb_hodotdc_clus_bar_tdc_timehitpos = {rd, "Ndata.bb.hodotdc.clus.bar.tdc.timehitpos"};
   TTreeReaderValue<Int_t> Ndata_bb_hodotdc_clus_bar_tdc_vpos = {rd, "Ndata.bb.hodotdc.clus.bar.tdc.vpos"};
   TTreeReaderArray<Double_t> bb_hodotdc_clus_bar_tdc_vpos = {rd, "bb.hodotdc.clus.bar.tdc.vpos"};
   TTreeReaderValue<Int_t> Ndata_bb_hodotdc_clus_id = {rd, "Ndata.bb.hodotdc.clus.id"};
   TTreeReaderValue<Int_t> Ndata_bb_hodotdc_clus_size = {rd, "Ndata.bb.hodotdc.clus.size"};
   TTreeReaderArray<Double_t> bb_hodotdc_clus_size = {rd, "bb.hodotdc.clus.size"};
   TTreeReaderValue<Int_t> Ndata_bb_hodotdc_clus_tdiff = {rd, "Ndata.bb.hodotdc.clus.tdiff"};
   TTreeReaderArray<Double_t> bb_hodotdc_clus_tdiff = {rd, "bb.hodotdc.clus.tdiff"};
   TTreeReaderValue<Int_t> Ndata_bb_hodotdc_clus_tmean = {rd, "Ndata.bb.hodotdc.clus.tmean"};
   TTreeReaderValue<Int_t> Ndata_bb_hodotdc_clus_totmean = {rd, "Ndata.bb.hodotdc.clus.totmean"};
   TTreeReaderArray<Double_t> bb_hodotdc_clus_totmean = {rd, "bb.hodotdc.clus.totmean"};
   TTreeReaderValue<Int_t> Ndata_bb_hodotdc_clus_trackindex = {rd, "Ndata.bb.hodotdc.clus.trackindex"};
   TTreeReaderArray<Double_t> bb_hodotdc_clus_trackindex = {rd, "bb.hodotdc.clus.trackindex"};
   TTreeReaderValue<Int_t> Ndata_bb_hodotdc_clus_xmean = {rd, "Ndata.bb.hodotdc.clus.xmean"};
   TTreeReaderValue<Int_t> Ndata_bb_hodotdc_clus_ymean = {rd, "Ndata.bb.hodotdc.clus.ymean"};
   TTreeReaderValue<Int_t> Ndata_bb_ps_clus_again = {rd, "Ndata.bb.ps.clus.again"};
   TTreeReaderArray<Double_t> bb_ps_clus_again = {rd, "bb.ps.clus.again"};
   TTreeReaderValue<Int_t> Ndata_bb_ps_clus_atime = {rd, "Ndata.bb.ps.clus.atime"};
   TTreeReaderArray<Double_t> bb_ps_clus_atime = {rd, "bb.ps.clus.atime"};
   TTreeReaderValue<Int_t> Ndata_bb_ps_clus_col = {rd, "Ndata.bb.ps.clus.col"};
   TTreeReaderArray<Double_t> bb_ps_clus_col = {rd, "bb.ps.clus.col"};
   TTreeReaderValue<Int_t> Ndata_bb_ps_clus_e = {rd, "Ndata.bb.ps.clus.e"};
   TTreeReaderArray<Double_t> bb_ps_clus_e = {rd, "bb.ps.clus.e"};
   TTreeReaderValue<Int_t> Ndata_bb_ps_clus_eblk = {rd, "Ndata.bb.ps.clus.eblk"};
   TTreeReaderArray<Double_t> bb_ps_clus_eblk = {rd, "bb.ps.clus.eblk"};
   TTreeReaderValue<Int_t> Ndata_bb_ps_clus_id = {rd, "Ndata.bb.ps.clus.id"};
   TTreeReaderArray<Double_t> bb_ps_clus_id = {rd, "bb.ps.clus.id"};
   TTreeReaderValue<Int_t> Ndata_bb_ps_clus_nblk = {rd, "Ndata.bb.ps.clus.nblk"};
   TTreeReaderArray<Double_t> bb_ps_clus_nblk = {rd, "bb.ps.clus.nblk"};
   TTreeReaderValue<Int_t> Ndata_bb_ps_clus_row = {rd, "Ndata.bb.ps.clus.row"};
   TTreeReaderArray<Double_t> bb_ps_clus_row = {rd, "bb.ps.clus.row"};
   TTreeReaderValue<Int_t> Ndata_bb_ps_clus_tdctime = {rd, "Ndata.bb.ps.clus.tdctime"};
   TTreeReaderArray<Double_t> bb_ps_clus_tdctime = {rd, "bb.ps.clus.tdctime"};
   TTreeReaderValue<Int_t> Ndata_bb_ps_clus_x = {rd, "Ndata.bb.ps.clus.x"};
   TTreeReaderArray<Double_t> bb_ps_clus_x = {rd, "bb.ps.clus.x"};
   TTreeReaderValue<Int_t> Ndata_bb_ps_clus_y = {rd, "Ndata.bb.ps.clus.y"};
   TTreeReaderArray<Double_t> bb_ps_clus_y = {rd, "bb.ps.clus.y"};
   TTreeReaderValue<Int_t> Ndata_bb_ps_clus_blk_again = {rd, "Ndata.bb.ps.clus_blk.again"};
   TTreeReaderArray<Double_t> bb_ps_clus_blk_again = {rd, "bb.ps.clus_blk.again"};
   TTreeReaderValue<Int_t> Ndata_bb_ps_clus_blk_atime = {rd, "Ndata.bb.ps.clus_blk.atime"};
   TTreeReaderArray<Double_t> bb_ps_clus_blk_atime = {rd, "bb.ps.clus_blk.atime"};
   TTreeReaderValue<Int_t> Ndata_bb_ps_clus_blk_col = {rd, "Ndata.bb.ps.clus_blk.col"};
   TTreeReaderArray<Double_t> bb_ps_clus_blk_col = {rd, "bb.ps.clus_blk.col"};
   TTreeReaderValue<Int_t> Ndata_bb_ps_clus_blk_e = {rd, "Ndata.bb.ps.clus_blk.e"};
   TTreeReaderArray<Double_t> bb_ps_clus_blk_e = {rd, "bb.ps.clus_blk.e"};
   TTreeReaderValue<Int_t> Ndata_bb_ps_clus_blk_id = {rd, "Ndata.bb.ps.clus_blk.id"};
   TTreeReaderArray<Double_t> bb_ps_clus_blk_id = {rd, "bb.ps.clus_blk.id"};
   TTreeReaderValue<Int_t> Ndata_bb_ps_clus_blk_row = {rd, "Ndata.bb.ps.clus_blk.row"};
   TTreeReaderArray<Double_t> bb_ps_clus_blk_row = {rd, "bb.ps.clus_blk.row"};
   TTreeReaderValue<Int_t> Ndata_bb_ps_clus_blk_tdctime = {rd, "Ndata.bb.ps.clus_blk.tdctime"};
   TTreeReaderArray<Double_t> bb_ps_clus_blk_tdctime = {rd, "bb.ps.clus_blk.tdctime"};
   TTreeReaderValue<Int_t> Ndata_bb_ps_clus_blk_x = {rd, "Ndata.bb.ps.clus_blk.x"};
   TTreeReaderArray<Double_t> bb_ps_clus_blk_x = {rd, "bb.ps.clus_blk.x"};
   TTreeReaderValue<Int_t> Ndata_bb_ps_clus_blk_y = {rd, "Ndata.bb.ps.clus_blk.y"};
   TTreeReaderArray<Double_t> bb_ps_clus_blk_y = {rd, "bb.ps.clus_blk.y"};
   TTreeReaderValue<Int_t> Ndata_bb_sh_clus_again = {rd, "Ndata.bb.sh.clus.again"};
   TTreeReaderArray<Double_t> bb_sh_clus_again = {rd, "bb.sh.clus.again"};
   TTreeReaderValue<Int_t> Ndata_bb_sh_clus_atime = {rd, "Ndata.bb.sh.clus.atime"};
   TTreeReaderArray<Double_t> bb_sh_clus_atime = {rd, "bb.sh.clus.atime"};
   TTreeReaderValue<Int_t> Ndata_bb_sh_clus_col = {rd, "Ndata.bb.sh.clus.col"};
   TTreeReaderArray<Double_t> bb_sh_clus_col = {rd, "bb.sh.clus.col"};
   TTreeReaderValue<Int_t> Ndata_bb_sh_clus_e = {rd, "Ndata.bb.sh.clus.e"};
   TTreeReaderArray<Double_t> bb_sh_clus_e = {rd, "bb.sh.clus.e"};
   TTreeReaderValue<Int_t> Ndata_bb_sh_clus_eblk = {rd, "Ndata.bb.sh.clus.eblk"};
   TTreeReaderArray<Double_t> bb_sh_clus_eblk = {rd, "bb.sh.clus.eblk"};
   TTreeReaderValue<Int_t> Ndata_bb_sh_clus_id = {rd, "Ndata.bb.sh.clus.id"};
   TTreeReaderArray<Double_t> bb_sh_clus_id = {rd, "bb.sh.clus.id"};
   TTreeReaderValue<Int_t> Ndata_bb_sh_clus_nblk = {rd, "Ndata.bb.sh.clus.nblk"};
   TTreeReaderArray<Double_t> bb_sh_clus_nblk = {rd, "bb.sh.clus.nblk"};
   TTreeReaderValue<Int_t> Ndata_bb_sh_clus_row = {rd, "Ndata.bb.sh.clus.row"};
   TTreeReaderArray<Double_t> bb_sh_clus_row = {rd, "bb.sh.clus.row"};
   TTreeReaderValue<Int_t> Ndata_bb_sh_clus_tdctime = {rd, "Ndata.bb.sh.clus.tdctime"};
   TTreeReaderArray<Double_t> bb_sh_clus_tdctime = {rd, "bb.sh.clus.tdctime"};
   TTreeReaderValue<Int_t> Ndata_bb_sh_clus_x = {rd, "Ndata.bb.sh.clus.x"};
   TTreeReaderArray<Double_t> bb_sh_clus_x = {rd, "bb.sh.clus.x"};
   TTreeReaderValue<Int_t> Ndata_bb_sh_clus_y = {rd, "Ndata.bb.sh.clus.y"};
   TTreeReaderArray<Double_t> bb_sh_clus_y = {rd, "bb.sh.clus.y"};
   TTreeReaderValue<Int_t> Ndata_bb_sh_clus_blk_again = {rd, "Ndata.bb.sh.clus_blk.again"};
   TTreeReaderArray<Double_t> bb_sh_clus_blk_again = {rd, "bb.sh.clus_blk.again"};
   TTreeReaderValue<Int_t> Ndata_bb_sh_clus_blk_atime = {rd, "Ndata.bb.sh.clus_blk.atime"};
   TTreeReaderArray<Double_t> bb_sh_clus_blk_atime = {rd, "bb.sh.clus_blk.atime"};
   TTreeReaderValue<Int_t> Ndata_bb_sh_clus_blk_col = {rd, "Ndata.bb.sh.clus_blk.col"};
   TTreeReaderArray<Double_t> bb_sh_clus_blk_col = {rd, "bb.sh.clus_blk.col"};
   TTreeReaderValue<Int_t> Ndata_bb_sh_clus_blk_e = {rd, "Ndata.bb.sh.clus_blk.e"};
   TTreeReaderArray<Double_t> bb_sh_clus_blk_e = {rd, "bb.sh.clus_blk.e"};
   TTreeReaderValue<Int_t> Ndata_bb_sh_clus_blk_id = {rd, "Ndata.bb.sh.clus_blk.id"};
   TTreeReaderArray<Double_t> bb_sh_clus_blk_id = {rd, "bb.sh.clus_blk.id"};
   TTreeReaderValue<Int_t> Ndata_bb_sh_clus_blk_row = {rd, "Ndata.bb.sh.clus_blk.row"};
   TTreeReaderArray<Double_t> bb_sh_clus_blk_row = {rd, "bb.sh.clus_blk.row"};
   TTreeReaderValue<Int_t> Ndata_bb_sh_clus_blk_tdctime = {rd, "Ndata.bb.sh.clus_blk.tdctime"};
   TTreeReaderArray<Double_t> bb_sh_clus_blk_tdctime = {rd, "bb.sh.clus_blk.tdctime"};
   TTreeReaderValue<Int_t> Ndata_bb_sh_clus_blk_x = {rd, "Ndata.bb.sh.clus_blk.x"};
   TTreeReaderArray<Double_t> bb_sh_clus_blk_x = {rd, "bb.sh.clus_blk.x"};
   TTreeReaderValue<Int_t> Ndata_bb_sh_clus_blk_y = {rd, "Ndata.bb.sh.clus_blk.y"};
   TTreeReaderArray<Double_t> bb_sh_clus_blk_y = {rd, "bb.sh.clus_blk.y"};
   TTreeReaderValue<Int_t> Ndata_bb_tr_beta = {rd, "Ndata.bb.tr.beta"};
   TTreeReaderArray<Double_t> bb_tr_beta = {rd, "bb.tr.beta"};
   TTreeReaderValue<Int_t> Ndata_bb_tr_chi2 = {rd, "Ndata.bb.tr.chi2"};
   TTreeReaderArray<Double_t> bb_tr_chi2 = {rd, "bb.tr.chi2"};
   TTreeReaderValue<Int_t> Ndata_bb_tr_d_ph = {rd, "Ndata.bb.tr.d_ph"};
   TTreeReaderArray<Double_t> bb_tr_d_ph = {rd, "bb.tr.d_ph"};
   TTreeReaderValue<Int_t> Ndata_bb_tr_d_th = {rd, "Ndata.bb.tr.d_th"};
   TTreeReaderArray<Double_t> bb_tr_d_th = {rd, "bb.tr.d_th"};
   TTreeReaderValue<Int_t> Ndata_bb_tr_d_x = {rd, "Ndata.bb.tr.d_x"};
   TTreeReaderArray<Double_t> bb_tr_d_x = {rd, "bb.tr.d_x"};
   TTreeReaderValue<Int_t> Ndata_bb_tr_d_y = {rd, "Ndata.bb.tr.d_y"};
   TTreeReaderArray<Double_t> bb_tr_d_y = {rd, "bb.tr.d_y"};
   TTreeReaderValue<Int_t> Ndata_bb_tr_dbeta = {rd, "Ndata.bb.tr.dbeta"};
   TTreeReaderArray<Double_t> bb_tr_dbeta = {rd, "bb.tr.dbeta"};
   TTreeReaderValue<Int_t> Ndata_bb_tr_dtime = {rd, "Ndata.bb.tr.dtime"};
   TTreeReaderArray<Double_t> bb_tr_dtime = {rd, "bb.tr.dtime"};
   TTreeReaderValue<Int_t> Ndata_bb_tr_flag = {rd, "Ndata.bb.tr.flag"};
   TTreeReaderArray<Double_t> bb_tr_flag = {rd, "bb.tr.flag"};
   TTreeReaderValue<Int_t> Ndata_bb_tr_ndof = {rd, "Ndata.bb.tr.ndof"};
   TTreeReaderArray<Double_t> bb_tr_ndof = {rd, "bb.tr.ndof"};
   TTreeReaderValue<Int_t> Ndata_bb_tr_p = {rd, "Ndata.bb.tr.p"};
   TTreeReaderArray<Double_t> bb_tr_p = {rd, "bb.tr.p"};
   TTreeReaderValue<Int_t> Ndata_bb_tr_pathl = {rd, "Ndata.bb.tr.pathl"};
   TTreeReaderValue<Int_t> Ndata_bb_tr_ph = {rd, "Ndata.bb.tr.ph"};
   TTreeReaderArray<Double_t> bb_tr_ph = {rd, "bb.tr.ph"};
   TTreeReaderValue<Int_t> Ndata_bb_tr_px = {rd, "Ndata.bb.tr.px"};
   TTreeReaderArray<Double_t> bb_tr_px = {rd, "bb.tr.px"};
   TTreeReaderValue<Int_t> Ndata_bb_tr_py = {rd, "Ndata.bb.tr.py"};
   TTreeReaderArray<Double_t> bb_tr_py = {rd, "bb.tr.py"};
   TTreeReaderValue<Int_t> Ndata_bb_tr_pz = {rd, "Ndata.bb.tr.pz"};
   TTreeReaderArray<Double_t> bb_tr_pz = {rd, "bb.tr.pz"};
   TTreeReaderValue<Int_t> Ndata_bb_tr_r_ph = {rd, "Ndata.bb.tr.r_ph"};
   TTreeReaderArray<Double_t> bb_tr_r_ph = {rd, "bb.tr.r_ph"};
   TTreeReaderValue<Int_t> Ndata_bb_tr_r_th = {rd, "Ndata.bb.tr.r_th"};
   TTreeReaderArray<Double_t> bb_tr_r_th = {rd, "bb.tr.r_th"};
   TTreeReaderValue<Int_t> Ndata_bb_tr_r_x = {rd, "Ndata.bb.tr.r_x"};
   TTreeReaderArray<Double_t> bb_tr_r_x = {rd, "bb.tr.r_x"};
   TTreeReaderValue<Int_t> Ndata_bb_tr_r_y = {rd, "Ndata.bb.tr.r_y"};
   TTreeReaderArray<Double_t> bb_tr_r_y = {rd, "bb.tr.r_y"};
   TTreeReaderValue<Int_t> Ndata_bb_tr_tg_dp = {rd, "Ndata.bb.tr.tg_dp"};
   TTreeReaderArray<Double_t> bb_tr_tg_dp = {rd, "bb.tr.tg_dp"};
   TTreeReaderValue<Int_t> Ndata_bb_tr_tg_ph = {rd, "Ndata.bb.tr.tg_ph"};
   TTreeReaderArray<Double_t> bb_tr_tg_ph = {rd, "bb.tr.tg_ph"};
   TTreeReaderValue<Int_t> Ndata_bb_tr_tg_th = {rd, "Ndata.bb.tr.tg_th"};
   TTreeReaderArray<Double_t> bb_tr_tg_th = {rd, "bb.tr.tg_th"};
   TTreeReaderValue<Int_t> Ndata_bb_tr_tg_x = {rd, "Ndata.bb.tr.tg_x"};
   TTreeReaderArray<Double_t> bb_tr_tg_x = {rd, "bb.tr.tg_x"};
   TTreeReaderValue<Int_t> Ndata_bb_tr_tg_y = {rd, "Ndata.bb.tr.tg_y"};
   TTreeReaderArray<Double_t> bb_tr_tg_y = {rd, "bb.tr.tg_y"};
   TTreeReaderValue<Int_t> Ndata_bb_tr_th = {rd, "Ndata.bb.tr.th"};
   TTreeReaderArray<Double_t> bb_tr_th = {rd, "bb.tr.th"};
   TTreeReaderValue<Int_t> Ndata_bb_tr_time = {rd, "Ndata.bb.tr.time"};
   TTreeReaderArray<Double_t> bb_tr_time = {rd, "bb.tr.time"};
   TTreeReaderValue<Int_t> Ndata_bb_tr_vx = {rd, "Ndata.bb.tr.vx"};
   TTreeReaderArray<Double_t> bb_tr_vx = {rd, "bb.tr.vx"};
   TTreeReaderValue<Int_t> Ndata_bb_tr_vy = {rd, "Ndata.bb.tr.vy"};
   TTreeReaderArray<Double_t> bb_tr_vy = {rd, "bb.tr.vy"};
   TTreeReaderValue<Int_t> Ndata_bb_tr_vz = {rd, "Ndata.bb.tr.vz"};
   TTreeReaderArray<Double_t> bb_tr_vz = {rd, "bb.tr.vz"};
   TTreeReaderValue<Int_t> Ndata_bb_tr_x = {rd, "Ndata.bb.tr.x"};
   TTreeReaderArray<Double_t> bb_tr_x = {rd, "bb.tr.x"};
   TTreeReaderValue<Int_t> Ndata_bb_tr_y = {rd, "Ndata.bb.tr.y"};
   TTreeReaderArray<Double_t> bb_tr_y = {rd, "bb.tr.y"};
   TTreeReaderValue<Int_t> Ndata_bb_x_bcp = {rd, "Ndata.bb.x_bcp"};
   TTreeReaderArray<Double_t> bb_x_bcp = {rd, "bb.x_bcp"};
   TTreeReaderValue<Int_t> Ndata_bb_x_fcp = {rd, "Ndata.bb.x_fcp"};
   TTreeReaderArray<Double_t> bb_x_fcp = {rd, "bb.x_fcp"};
   TTreeReaderValue<Int_t> Ndata_bb_y_bcp = {rd, "Ndata.bb.y_bcp"};
   TTreeReaderArray<Double_t> bb_y_bcp = {rd, "bb.y_bcp"};
   TTreeReaderValue<Int_t> Ndata_bb_y_fcp = {rd, "Ndata.bb.y_fcp"};
   TTreeReaderArray<Double_t> bb_y_fcp = {rd, "bb.y_fcp"};
   TTreeReaderValue<Int_t> Ndata_bb_z_bcp = {rd, "Ndata.bb.z_bcp"};
   TTreeReaderArray<Double_t> bb_z_bcp = {rd, "bb.z_bcp"};
   TTreeReaderValue<Int_t> Ndata_bb_z_fcp = {rd, "Ndata.bb.z_fcp"};
   TTreeReaderArray<Double_t> bb_z_fcp = {rd, "bb.z_fcp"};
   TTreeReaderValue<Int_t> Ndata_sbs_hcal_clus_again = {rd, "Ndata.sbs.hcal.clus.again"};
   TTreeReaderArray<Double_t> sbs_hcal_clus_again = {rd, "sbs.hcal.clus.again"};
   TTreeReaderValue<Int_t> Ndata_sbs_hcal_clus_atime = {rd, "Ndata.sbs.hcal.clus.atime"};
   TTreeReaderArray<Double_t> sbs_hcal_clus_atime = {rd, "sbs.hcal.clus.atime"};
   TTreeReaderValue<Int_t> Ndata_sbs_hcal_clus_col = {rd, "Ndata.sbs.hcal.clus.col"};
   TTreeReaderArray<Double_t> sbs_hcal_clus_col = {rd, "sbs.hcal.clus.col"};
   TTreeReaderValue<Int_t> Ndata_sbs_hcal_clus_e = {rd, "Ndata.sbs.hcal.clus.e"};
   TTreeReaderArray<Double_t> sbs_hcal_clus_e = {rd, "sbs.hcal.clus.e"};
   TTreeReaderValue<Int_t> Ndata_sbs_hcal_clus_eblk = {rd, "Ndata.sbs.hcal.clus.eblk"};
   TTreeReaderArray<Double_t> sbs_hcal_clus_eblk = {rd, "sbs.hcal.clus.eblk"};
   TTreeReaderValue<Int_t> Ndata_sbs_hcal_clus_id = {rd, "Ndata.sbs.hcal.clus.id"};
   TTreeReaderArray<Double_t> sbs_hcal_clus_id = {rd, "sbs.hcal.clus.id"};
   TTreeReaderValue<Int_t> Ndata_sbs_hcal_clus_nblk = {rd, "Ndata.sbs.hcal.clus.nblk"};
   TTreeReaderArray<Double_t> sbs_hcal_clus_nblk = {rd, "sbs.hcal.clus.nblk"};
   TTreeReaderValue<Int_t> Ndata_sbs_hcal_clus_row = {rd, "Ndata.sbs.hcal.clus.row"};
   TTreeReaderArray<Double_t> sbs_hcal_clus_row = {rd, "sbs.hcal.clus.row"};
   TTreeReaderValue<Int_t> Ndata_sbs_hcal_clus_tdctime = {rd, "Ndata.sbs.hcal.clus.tdctime"};
   TTreeReaderArray<Double_t> sbs_hcal_clus_tdctime = {rd, "sbs.hcal.clus.tdctime"};
   TTreeReaderValue<Int_t> Ndata_sbs_hcal_clus_x = {rd, "Ndata.sbs.hcal.clus.x"};
   TTreeReaderArray<Double_t> sbs_hcal_clus_x = {rd, "sbs.hcal.clus.x"};
   TTreeReaderValue<Int_t> Ndata_sbs_hcal_clus_y = {rd, "Ndata.sbs.hcal.clus.y"};
   TTreeReaderArray<Double_t> sbs_hcal_clus_y = {rd, "sbs.hcal.clus.y"};
   TTreeReaderValue<Int_t> Ndata_sbs_hcal_clus_blk_again = {rd, "Ndata.sbs.hcal.clus_blk.again"};
   TTreeReaderArray<Double_t> sbs_hcal_clus_blk_again = {rd, "sbs.hcal.clus_blk.again"};
   TTreeReaderValue<Int_t> Ndata_sbs_hcal_clus_blk_atime = {rd, "Ndata.sbs.hcal.clus_blk.atime"};
   TTreeReaderArray<Double_t> sbs_hcal_clus_blk_atime = {rd, "sbs.hcal.clus_blk.atime"};
   TTreeReaderValue<Int_t> Ndata_sbs_hcal_clus_blk_col = {rd, "Ndata.sbs.hcal.clus_blk.col"};
   TTreeReaderArray<Double_t> sbs_hcal_clus_blk_col = {rd, "sbs.hcal.clus_blk.col"};
   TTreeReaderValue<Int_t> Ndata_sbs_hcal_clus_blk_e = {rd, "Ndata.sbs.hcal.clus_blk.e"};
   TTreeReaderArray<Double_t> sbs_hcal_clus_blk_e = {rd, "sbs.hcal.clus_blk.e"};
   TTreeReaderValue<Int_t> Ndata_sbs_hcal_clus_blk_id = {rd, "Ndata.sbs.hcal.clus_blk.id"};
   TTreeReaderArray<Double_t> sbs_hcal_clus_blk_id = {rd, "sbs.hcal.clus_blk.id"};
   TTreeReaderValue<Int_t> Ndata_sbs_hcal_clus_blk_row = {rd, "Ndata.sbs.hcal.clus_blk.row"};
   TTreeReaderArray<Double_t> sbs_hcal_clus_blk_row = {rd, "sbs.hcal.clus_blk.row"};
   TTreeReaderValue<Int_t> Ndata_sbs_hcal_clus_blk_tdctime = {rd, "Ndata.sbs.hcal.clus_blk.tdctime"};
   TTreeReaderArray<Double_t> sbs_hcal_clus_blk_tdctime = {rd, "sbs.hcal.clus_blk.tdctime"};
   TTreeReaderValue<Int_t> Ndata_sbs_hcal_clus_blk_x = {rd, "Ndata.sbs.hcal.clus_blk.x"};
   TTreeReaderArray<Double_t> sbs_hcal_clus_blk_x = {rd, "sbs.hcal.clus_blk.x"};
   TTreeReaderValue<Int_t> Ndata_sbs_hcal_clus_blk_y = {rd, "Ndata.sbs.hcal.clus_blk.y"};
   TTreeReaderArray<Double_t> sbs_hcal_clus_blk_y = {rd, "sbs.hcal.clus_blk.y"};
   TTreeReaderValue<Int_t> Ndata_sbs_tr_beta = {rd, "Ndata.sbs.tr.beta"};
   TTreeReaderArray<Double_t> sbs_tr_beta = {rd, "sbs.tr.beta"};
   TTreeReaderValue<Int_t> Ndata_sbs_tr_chi2 = {rd, "Ndata.sbs.tr.chi2"};
   TTreeReaderArray<Double_t> sbs_tr_chi2 = {rd, "sbs.tr.chi2"};
   TTreeReaderValue<Int_t> Ndata_sbs_tr_d_ph = {rd, "Ndata.sbs.tr.d_ph"};
   TTreeReaderArray<Double_t> sbs_tr_d_ph = {rd, "sbs.tr.d_ph"};
   TTreeReaderValue<Int_t> Ndata_sbs_tr_d_th = {rd, "Ndata.sbs.tr.d_th"};
   TTreeReaderArray<Double_t> sbs_tr_d_th = {rd, "sbs.tr.d_th"};
   TTreeReaderValue<Int_t> Ndata_sbs_tr_d_x = {rd, "Ndata.sbs.tr.d_x"};
   TTreeReaderArray<Double_t> sbs_tr_d_x = {rd, "sbs.tr.d_x"};
   TTreeReaderValue<Int_t> Ndata_sbs_tr_d_y = {rd, "Ndata.sbs.tr.d_y"};
   TTreeReaderArray<Double_t> sbs_tr_d_y = {rd, "sbs.tr.d_y"};
   TTreeReaderValue<Int_t> Ndata_sbs_tr_dbeta = {rd, "Ndata.sbs.tr.dbeta"};
   TTreeReaderArray<Double_t> sbs_tr_dbeta = {rd, "sbs.tr.dbeta"};
   TTreeReaderValue<Int_t> Ndata_sbs_tr_dtime = {rd, "Ndata.sbs.tr.dtime"};
   TTreeReaderArray<Double_t> sbs_tr_dtime = {rd, "sbs.tr.dtime"};
   TTreeReaderValue<Int_t> Ndata_sbs_tr_flag = {rd, "Ndata.sbs.tr.flag"};
   TTreeReaderArray<Double_t> sbs_tr_flag = {rd, "sbs.tr.flag"};
   TTreeReaderValue<Int_t> Ndata_sbs_tr_ndof = {rd, "Ndata.sbs.tr.ndof"};
   TTreeReaderArray<Double_t> sbs_tr_ndof = {rd, "sbs.tr.ndof"};
   TTreeReaderValue<Int_t> Ndata_sbs_tr_p = {rd, "Ndata.sbs.tr.p"};
   TTreeReaderArray<Double_t> sbs_tr_p = {rd, "sbs.tr.p"};
   TTreeReaderValue<Int_t> Ndata_sbs_tr_pathl = {rd, "Ndata.sbs.tr.pathl"};
   TTreeReaderArray<Double_t> sbs_tr_pathl = {rd, "sbs.tr.pathl"};
   TTreeReaderValue<Int_t> Ndata_sbs_tr_ph = {rd, "Ndata.sbs.tr.ph"};
   TTreeReaderArray<Double_t> sbs_tr_ph = {rd, "sbs.tr.ph"};
   TTreeReaderValue<Int_t> Ndata_sbs_tr_px = {rd, "Ndata.sbs.tr.px"};
   TTreeReaderArray<Double_t> sbs_tr_px = {rd, "sbs.tr.px"};
   TTreeReaderValue<Int_t> Ndata_sbs_tr_py = {rd, "Ndata.sbs.tr.py"};
   TTreeReaderArray<Double_t> sbs_tr_py = {rd, "sbs.tr.py"};
   TTreeReaderValue<Int_t> Ndata_sbs_tr_pz = {rd, "Ndata.sbs.tr.pz"};
   TTreeReaderArray<Double_t> sbs_tr_pz = {rd, "sbs.tr.pz"};
   TTreeReaderValue<Int_t> Ndata_sbs_tr_r_ph = {rd, "Ndata.sbs.tr.r_ph"};
   TTreeReaderArray<Double_t> sbs_tr_r_ph = {rd, "sbs.tr.r_ph"};
   TTreeReaderValue<Int_t> Ndata_sbs_tr_r_th = {rd, "Ndata.sbs.tr.r_th"};
   TTreeReaderArray<Double_t> sbs_tr_r_th = {rd, "sbs.tr.r_th"};
   TTreeReaderValue<Int_t> Ndata_sbs_tr_r_x = {rd, "Ndata.sbs.tr.r_x"};
   TTreeReaderArray<Double_t> sbs_tr_r_x = {rd, "sbs.tr.r_x"};
   TTreeReaderValue<Int_t> Ndata_sbs_tr_r_y = {rd, "Ndata.sbs.tr.r_y"};
   TTreeReaderArray<Double_t> sbs_tr_r_y = {rd, "sbs.tr.r_y"};
   TTreeReaderValue<Int_t> Ndata_sbs_tr_tg_dp = {rd, "Ndata.sbs.tr.tg_dp"};
   TTreeReaderArray<Double_t> sbs_tr_tg_dp = {rd, "sbs.tr.tg_dp"};
   TTreeReaderValue<Int_t> Ndata_sbs_tr_tg_ph = {rd, "Ndata.sbs.tr.tg_ph"};
   TTreeReaderArray<Double_t> sbs_tr_tg_ph = {rd, "sbs.tr.tg_ph"};
   TTreeReaderValue<Int_t> Ndata_sbs_tr_tg_th = {rd, "Ndata.sbs.tr.tg_th"};
   TTreeReaderArray<Double_t> sbs_tr_tg_th = {rd, "sbs.tr.tg_th"};
   TTreeReaderValue<Int_t> Ndata_sbs_tr_tg_y = {rd, "Ndata.sbs.tr.tg_y"};
   TTreeReaderArray<Double_t> sbs_tr_tg_y = {rd, "sbs.tr.tg_y"};
   TTreeReaderValue<Int_t> Ndata_sbs_tr_th = {rd, "Ndata.sbs.tr.th"};
   TTreeReaderArray<Double_t> sbs_tr_th = {rd, "sbs.tr.th"};
   TTreeReaderValue<Int_t> Ndata_sbs_tr_time = {rd, "Ndata.sbs.tr.time"};
   TTreeReaderArray<Double_t> sbs_tr_time = {rd, "sbs.tr.time"};
   TTreeReaderValue<Int_t> Ndata_sbs_tr_vx = {rd, "Ndata.sbs.tr.vx"};
   TTreeReaderArray<Double_t> sbs_tr_vx = {rd, "sbs.tr.vx"};
   TTreeReaderValue<Int_t> Ndata_sbs_tr_vy = {rd, "Ndata.sbs.tr.vy"};
   TTreeReaderArray<Double_t> sbs_tr_vy = {rd, "sbs.tr.vy"};
   TTreeReaderValue<Int_t> Ndata_sbs_tr_vz = {rd, "Ndata.sbs.tr.vz"};
   TTreeReaderArray<Double_t> sbs_tr_vz = {rd, "sbs.tr.vz"};
   TTreeReaderValue<Int_t> Ndata_sbs_tr_x = {rd, "Ndata.sbs.tr.x"};
   TTreeReaderArray<Double_t> sbs_tr_x = {rd, "sbs.tr.x"};
   TTreeReaderValue<Int_t> Ndata_sbs_tr_y = {rd, "Ndata.sbs.tr.y"};
   TTreeReaderArray<Double_t> sbs_tr_y = {rd, "sbs.tr.y"};
   TTreeReaderValue<Int_t> Ndata_sbs_trig_a_amp_p = {rd, "Ndata.sbs.trig.a_amp_p"};
   TTreeReaderArray<Double_t> sbs_trig_a_amp_p = {rd, "sbs.trig.a_amp_p"};
   TTreeReaderValue<Int_t> Ndata_sbs_trig_a_p = {rd, "Ndata.sbs.trig.a_p"};
   TTreeReaderArray<Double_t> sbs_trig_a_p = {rd, "sbs.trig.a_p"};
   TTreeReaderValue<Int_t> Ndata_sbs_trig_a_time = {rd, "Ndata.sbs.trig.a_time"};
   TTreeReaderArray<Double_t> sbs_trig_a_time = {rd, "sbs.trig.a_time"};
   TTreeReaderValue<Int_t> Ndata_sbs_trig_adcelemID = {rd, "Ndata.sbs.trig.adcelemID"};
   TTreeReaderArray<Double_t> sbs_trig_adcelemID = {rd, "sbs.trig.adcelemID"};
   TTreeReaderValue<Int_t> Ndata_sbs_x_bcp = {rd, "Ndata.sbs.x_bcp"};
   TTreeReaderArray<Double_t> sbs_x_bcp = {rd, "sbs.x_bcp"};
   TTreeReaderValue<Int_t> Ndata_sbs_x_fcp = {rd, "Ndata.sbs.x_fcp"};
   TTreeReaderArray<Double_t> sbs_x_fcp = {rd, "sbs.x_fcp"};
   TTreeReaderValue<Int_t> Ndata_sbs_y_bcp = {rd, "Ndata.sbs.y_bcp"};
   TTreeReaderArray<Double_t> sbs_y_bcp = {rd, "sbs.y_bcp"};
   TTreeReaderValue<Int_t> Ndata_sbs_y_fcp = {rd, "Ndata.sbs.y_fcp"};
   TTreeReaderArray<Double_t> sbs_y_fcp = {rd, "sbs.y_fcp"};
   TTreeReaderValue<Int_t> Ndata_sbs_z_bcp = {rd, "Ndata.sbs.z_bcp"};
   TTreeReaderArray<Double_t> sbs_z_bcp = {rd, "sbs.z_bcp"};
   TTreeReaderValue<Int_t> Ndata_sbs_z_fcp = {rd, "Ndata.sbs.z_fcp"};
   TTreeReaderArray<Double_t> sbs_z_fcp = {rd, "sbs.z_fcp"};
   TTreeReaderValue<Double_t> BB_gold_beta = {rd, "BB.gold.beta"};
   TTreeReaderValue<Double_t> BB_gold_dp = {rd, "BB.gold.dp"};
   TTreeReaderValue<Double_t> BB_gold_index = {rd, "BB.gold.index"};
   TTreeReaderValue<Double_t> BB_gold_ok = {rd, "BB.gold.ok"};
   TTreeReaderValue<Double_t> BB_gold_p = {rd, "BB.gold.p"};
   TTreeReaderValue<Double_t> BB_gold_ph = {rd, "BB.gold.ph"};
   TTreeReaderValue<Double_t> BB_gold_px = {rd, "BB.gold.px"};
   TTreeReaderValue<Double_t> BB_gold_py = {rd, "BB.gold.py"};
   TTreeReaderValue<Double_t> BB_gold_pz = {rd, "BB.gold.pz"};
   TTreeReaderValue<Double_t> BB_gold_th = {rd, "BB.gold.th"};
   TTreeReaderValue<Double_t> BB_gold_x = {rd, "BB.gold.x"};
   TTreeReaderValue<Double_t> BB_gold_y = {rd, "BB.gold.y"};
   TTreeReaderValue<Double_t> Lrb_BPMA_rawcur_1 = {rd, "Lrb.BPMA.rawcur.1"};
   TTreeReaderValue<Double_t> Lrb_BPMA_rawcur_2 = {rd, "Lrb.BPMA.rawcur.2"};
   TTreeReaderValue<Double_t> Lrb_BPMA_rawcur_3 = {rd, "Lrb.BPMA.rawcur.3"};
   TTreeReaderValue<Double_t> Lrb_BPMA_rawcur_4 = {rd, "Lrb.BPMA.rawcur.4"};
   TTreeReaderValue<Double_t> Lrb_BPMA_rotpos1 = {rd, "Lrb.BPMA.rotpos1"};
   TTreeReaderValue<Double_t> Lrb_BPMA_rotpos2 = {rd, "Lrb.BPMA.rotpos2"};
   TTreeReaderValue<Double_t> Lrb_BPMA_x = {rd, "Lrb.BPMA.x"};
   TTreeReaderValue<Double_t> Lrb_BPMA_y = {rd, "Lrb.BPMA.y"};
   TTreeReaderValue<Double_t> Lrb_BPMA_z = {rd, "Lrb.BPMA.z"};
   TTreeReaderValue<Double_t> Lrb_BPMB_rawcur_1 = {rd, "Lrb.BPMB.rawcur.1"};
   TTreeReaderValue<Double_t> Lrb_BPMB_rawcur_2 = {rd, "Lrb.BPMB.rawcur.2"};
   TTreeReaderValue<Double_t> Lrb_BPMB_rawcur_3 = {rd, "Lrb.BPMB.rawcur.3"};
   TTreeReaderValue<Double_t> Lrb_BPMB_rawcur_4 = {rd, "Lrb.BPMB.rawcur.4"};
   TTreeReaderValue<Double_t> Lrb_BPMB_rotpos1 = {rd, "Lrb.BPMB.rotpos1"};
   TTreeReaderValue<Double_t> Lrb_BPMB_rotpos2 = {rd, "Lrb.BPMB.rotpos2"};
   TTreeReaderValue<Double_t> Lrb_BPMB_x = {rd, "Lrb.BPMB.x"};
   TTreeReaderValue<Double_t> Lrb_BPMB_y = {rd, "Lrb.BPMB.y"};
   TTreeReaderValue<Double_t> Lrb_BPMB_z = {rd, "Lrb.BPMB.z"};
   TTreeReaderValue<Double_t> Lrb_Raster_bpma_x = {rd, "Lrb.Raster.bpma.x"};
   TTreeReaderValue<Double_t> Lrb_Raster_bpma_y = {rd, "Lrb.Raster.bpma.y"};
   TTreeReaderValue<Double_t> Lrb_Raster_bpma_z = {rd, "Lrb.Raster.bpma.z"};
   TTreeReaderValue<Double_t> Lrb_Raster_bpmb_x = {rd, "Lrb.Raster.bpmb.x"};
   TTreeReaderValue<Double_t> Lrb_Raster_bpmb_y = {rd, "Lrb.Raster.bpmb.y"};
   TTreeReaderValue<Double_t> Lrb_Raster_bpmb_z = {rd, "Lrb.Raster.bpmb.z"};
   TTreeReaderValue<Double_t> Lrb_Raster_rawcur_x = {rd, "Lrb.Raster.rawcur.x"};
   TTreeReaderValue<Double_t> Lrb_Raster_rawcur_y = {rd, "Lrb.Raster.rawcur.y"};
   TTreeReaderValue<Double_t> Lrb_Raster_rawslope_x = {rd, "Lrb.Raster.rawslope.x"};
   TTreeReaderValue<Double_t> Lrb_Raster_rawslope_y = {rd, "Lrb.Raster.rawslope.y"};
   TTreeReaderValue<Double_t> Lrb_Raster_target_dir_x = {rd, "Lrb.Raster.target.dir.x"};
   TTreeReaderValue<Double_t> Lrb_Raster_target_dir_y = {rd, "Lrb.Raster.target.dir.y"};
   TTreeReaderValue<Double_t> Lrb_Raster_target_dir_z = {rd, "Lrb.Raster.target.dir.z"};
   TTreeReaderValue<Double_t> Lrb_Raster_target_x = {rd, "Lrb.Raster.target.x"};
   TTreeReaderValue<Double_t> Lrb_Raster_target_y = {rd, "Lrb.Raster.target.y"};
   TTreeReaderValue<Double_t> Lrb_Raster_target_z = {rd, "Lrb.Raster.target.z"};
   TTreeReaderValue<Double_t> Lrb_Raster2_bpma_x = {rd, "Lrb.Raster2.bpma.x"};
   TTreeReaderValue<Double_t> Lrb_Raster2_bpma_y = {rd, "Lrb.Raster2.bpma.y"};
   TTreeReaderValue<Double_t> Lrb_Raster2_bpma_z = {rd, "Lrb.Raster2.bpma.z"};
   TTreeReaderValue<Double_t> Lrb_Raster2_bpmb_x = {rd, "Lrb.Raster2.bpmb.x"};
   TTreeReaderValue<Double_t> Lrb_Raster2_bpmb_y = {rd, "Lrb.Raster2.bpmb.y"};
   TTreeReaderValue<Double_t> Lrb_Raster2_bpmb_z = {rd, "Lrb.Raster2.bpmb.z"};
   TTreeReaderValue<Double_t> Lrb_Raster2_rawcur_x = {rd, "Lrb.Raster2.rawcur.x"};
   TTreeReaderValue<Double_t> Lrb_Raster2_rawcur_y = {rd, "Lrb.Raster2.rawcur.y"};
   TTreeReaderValue<Double_t> Lrb_Raster2_rawslope_x = {rd, "Lrb.Raster2.rawslope.x"};
   TTreeReaderValue<Double_t> Lrb_Raster2_rawslope_y = {rd, "Lrb.Raster2.rawslope.y"};
   TTreeReaderValue<Double_t> Lrb_Raster2_target_dir_x = {rd, "Lrb.Raster2.target.dir.x"};
   TTreeReaderValue<Double_t> Lrb_Raster2_target_dir_y = {rd, "Lrb.Raster2.target.dir.y"};
   TTreeReaderValue<Double_t> Lrb_Raster2_target_dir_z = {rd, "Lrb.Raster2.target.dir.z"};
   TTreeReaderValue<Double_t> Lrb_Raster2_target_x = {rd, "Lrb.Raster2.target.x"};
   TTreeReaderValue<Double_t> Lrb_Raster2_target_y = {rd, "Lrb.Raster2.target.y"};
   TTreeReaderValue<Double_t> Lrb_Raster2_target_z = {rd, "Lrb.Raster2.target.z"};
   TTreeReaderValue<Double_t> SBSrb_BPMA_rawcur_1 = {rd, "SBSrb.BPMA.rawcur.1"};
   TTreeReaderValue<Double_t> SBSrb_BPMA_rawcur_2 = {rd, "SBSrb.BPMA.rawcur.2"};
   TTreeReaderValue<Double_t> SBSrb_BPMA_rawcur_3 = {rd, "SBSrb.BPMA.rawcur.3"};
   TTreeReaderValue<Double_t> SBSrb_BPMA_rawcur_4 = {rd, "SBSrb.BPMA.rawcur.4"};
   TTreeReaderValue<Double_t> SBSrb_BPMA_rotpos1 = {rd, "SBSrb.BPMA.rotpos1"};
   TTreeReaderValue<Double_t> SBSrb_BPMA_rotpos2 = {rd, "SBSrb.BPMA.rotpos2"};
   TTreeReaderValue<Double_t> SBSrb_BPMA_x = {rd, "SBSrb.BPMA.x"};
   TTreeReaderValue<Double_t> SBSrb_BPMA_y = {rd, "SBSrb.BPMA.y"};
   TTreeReaderValue<Double_t> SBSrb_BPMA_z = {rd, "SBSrb.BPMA.z"};
   TTreeReaderValue<Double_t> SBSrb_BPMB_rawcur_1 = {rd, "SBSrb.BPMB.rawcur.1"};
   TTreeReaderValue<Double_t> SBSrb_BPMB_rawcur_2 = {rd, "SBSrb.BPMB.rawcur.2"};
   TTreeReaderValue<Double_t> SBSrb_BPMB_rawcur_3 = {rd, "SBSrb.BPMB.rawcur.3"};
   TTreeReaderValue<Double_t> SBSrb_BPMB_rawcur_4 = {rd, "SBSrb.BPMB.rawcur.4"};
   TTreeReaderValue<Double_t> SBSrb_BPMB_rotpos1 = {rd, "SBSrb.BPMB.rotpos1"};
   TTreeReaderValue<Double_t> SBSrb_BPMB_rotpos2 = {rd, "SBSrb.BPMB.rotpos2"};
   TTreeReaderValue<Double_t> SBSrb_BPMB_x = {rd, "SBSrb.BPMB.x"};
   TTreeReaderValue<Double_t> SBSrb_BPMB_y = {rd, "SBSrb.BPMB.y"};
   TTreeReaderValue<Double_t> SBSrb_BPMB_z = {rd, "SBSrb.BPMB.z"};
   TTreeReaderValue<Double_t> SBSrb_Raster_bpma_x = {rd, "SBSrb.Raster.bpma.x"};
   TTreeReaderValue<Double_t> SBSrb_Raster_bpma_y = {rd, "SBSrb.Raster.bpma.y"};
   TTreeReaderValue<Double_t> SBSrb_Raster_bpma_z = {rd, "SBSrb.Raster.bpma.z"};
   TTreeReaderValue<Double_t> SBSrb_Raster_bpmb_x = {rd, "SBSrb.Raster.bpmb.x"};
   TTreeReaderValue<Double_t> SBSrb_Raster_bpmb_y = {rd, "SBSrb.Raster.bpmb.y"};
   TTreeReaderValue<Double_t> SBSrb_Raster_bpmb_z = {rd, "SBSrb.Raster.bpmb.z"};
   TTreeReaderValue<Double_t> SBSrb_Raster_rawcur_x = {rd, "SBSrb.Raster.rawcur.x"};
   TTreeReaderValue<Double_t> SBSrb_Raster_rawcur_y = {rd, "SBSrb.Raster.rawcur.y"};
   TTreeReaderValue<Double_t> SBSrb_Raster_rawslope_x = {rd, "SBSrb.Raster.rawslope.x"};
   TTreeReaderValue<Double_t> SBSrb_Raster_rawslope_y = {rd, "SBSrb.Raster.rawslope.y"};
   TTreeReaderValue<Double_t> SBSrb_Raster_target_dir_x = {rd, "SBSrb.Raster.target.dir.x"};
   TTreeReaderValue<Double_t> SBSrb_Raster_target_dir_y = {rd, "SBSrb.Raster.target.dir.y"};
   TTreeReaderValue<Double_t> SBSrb_Raster_target_dir_z = {rd, "SBSrb.Raster.target.dir.z"};
   TTreeReaderValue<Double_t> SBSrb_Raster_target_x = {rd, "SBSrb.Raster.target.x"};
   TTreeReaderValue<Double_t> SBSrb_Raster_target_y = {rd, "SBSrb.Raster.target.y"};
   TTreeReaderValue<Double_t> SBSrb_Raster_target_z = {rd, "SBSrb.Raster.target.z"};
   TTreeReaderValue<Double_t> SBSrb_Raster2_bpma_x = {rd, "SBSrb.Raster2.bpma.x"};
   TTreeReaderValue<Double_t> SBSrb_Raster2_bpma_y = {rd, "SBSrb.Raster2.bpma.y"};
   TTreeReaderValue<Double_t> SBSrb_Raster2_bpma_z = {rd, "SBSrb.Raster2.bpma.z"};
   TTreeReaderValue<Double_t> SBSrb_Raster2_bpmb_x = {rd, "SBSrb.Raster2.bpmb.x"};
   TTreeReaderValue<Double_t> SBSrb_Raster2_bpmb_y = {rd, "SBSrb.Raster2.bpmb.y"};
   TTreeReaderValue<Double_t> SBSrb_Raster2_bpmb_z = {rd, "SBSrb.Raster2.bpmb.z"};
   TTreeReaderValue<Double_t> SBSrb_Raster2_rawcur_x = {rd, "SBSrb.Raster2.rawcur.x"};
   TTreeReaderValue<Double_t> SBSrb_Raster2_rawcur_y = {rd, "SBSrb.Raster2.rawcur.y"};
   TTreeReaderValue<Double_t> SBSrb_Raster2_rawslope_x = {rd, "SBSrb.Raster2.rawslope.x"};
   TTreeReaderValue<Double_t> SBSrb_Raster2_rawslope_y = {rd, "SBSrb.Raster2.rawslope.y"};
   TTreeReaderValue<Double_t> SBSrb_Raster2_target_dir_x = {rd, "SBSrb.Raster2.target.dir.x"};
   TTreeReaderValue<Double_t> SBSrb_Raster2_target_dir_y = {rd, "SBSrb.Raster2.target.dir.y"};
   TTreeReaderValue<Double_t> SBSrb_Raster2_target_dir_z = {rd, "SBSrb.Raster2.target.dir.z"};
   TTreeReaderValue<Double_t> SBSrb_Raster2_target_x = {rd, "SBSrb.Raster2.target.x"};
   TTreeReaderValue<Double_t> SBSrb_Raster2_target_y = {rd, "SBSrb.Raster2.target.y"};
   TTreeReaderValue<Double_t> SBSrb_Raster2_target_z = {rd, "SBSrb.Raster2.target.z"};
   TTreeReaderValue<Double_t> bb_gem_hit_ngoodhits = {rd, "bb.gem.hit.ngoodhits"};
   TTreeReaderValue<Double_t> bb_gem_m0_clust_nclustu = {rd, "bb.gem.m0.clust.nclustu"};
   TTreeReaderValue<Double_t> bb_gem_m0_clust_nclustu_tot = {rd, "bb.gem.m0.clust.nclustu_tot"};
   TTreeReaderValue<Double_t> bb_gem_m0_clust_nclustv = {rd, "bb.gem.m0.clust.nclustv"};
   TTreeReaderValue<Double_t> bb_gem_m0_clust_nclustv_tot = {rd, "bb.gem.m0.clust.nclustv_tot"};
   TTreeReaderValue<Double_t> bb_gem_m0_strip_nstrips_keep = {rd, "bb.gem.m0.strip.nstrips_keep"};
   TTreeReaderValue<Double_t> bb_gem_m0_strip_nstrips_keepU = {rd, "bb.gem.m0.strip.nstrips_keepU"};
   TTreeReaderValue<Double_t> bb_gem_m0_strip_nstrips_keepV = {rd, "bb.gem.m0.strip.nstrips_keepV"};
   TTreeReaderValue<Double_t> bb_gem_m0_strip_nstrips_keep_lmax = {rd, "bb.gem.m0.strip.nstrips_keep_lmax"};
   TTreeReaderValue<Double_t> bb_gem_m0_strip_nstrips_keep_lmaxU = {rd, "bb.gem.m0.strip.nstrips_keep_lmaxU"};
   TTreeReaderValue<Double_t> bb_gem_m0_strip_nstrips_keep_lmaxV = {rd, "bb.gem.m0.strip.nstrips_keep_lmaxV"};
   TTreeReaderValue<Double_t> bb_gem_m0_strip_nstripsfired = {rd, "bb.gem.m0.strip.nstripsfired"};
   TTreeReaderValue<Double_t> bb_gem_m1_clust_nclustu = {rd, "bb.gem.m1.clust.nclustu"};
   TTreeReaderValue<Double_t> bb_gem_m1_clust_nclustu_tot = {rd, "bb.gem.m1.clust.nclustu_tot"};
   TTreeReaderValue<Double_t> bb_gem_m1_clust_nclustv = {rd, "bb.gem.m1.clust.nclustv"};
   TTreeReaderValue<Double_t> bb_gem_m1_clust_nclustv_tot = {rd, "bb.gem.m1.clust.nclustv_tot"};
   TTreeReaderValue<Double_t> bb_gem_m1_strip_nstrips_keep = {rd, "bb.gem.m1.strip.nstrips_keep"};
   TTreeReaderValue<Double_t> bb_gem_m1_strip_nstrips_keepU = {rd, "bb.gem.m1.strip.nstrips_keepU"};
   TTreeReaderValue<Double_t> bb_gem_m1_strip_nstrips_keepV = {rd, "bb.gem.m1.strip.nstrips_keepV"};
   TTreeReaderValue<Double_t> bb_gem_m1_strip_nstrips_keep_lmax = {rd, "bb.gem.m1.strip.nstrips_keep_lmax"};
   TTreeReaderValue<Double_t> bb_gem_m1_strip_nstrips_keep_lmaxU = {rd, "bb.gem.m1.strip.nstrips_keep_lmaxU"};
   TTreeReaderValue<Double_t> bb_gem_m1_strip_nstrips_keep_lmaxV = {rd, "bb.gem.m1.strip.nstrips_keep_lmaxV"};
   TTreeReaderValue<Double_t> bb_gem_m1_strip_nstripsfired = {rd, "bb.gem.m1.strip.nstripsfired"};
   TTreeReaderValue<Double_t> bb_gem_m2_clust_nclustu = {rd, "bb.gem.m2.clust.nclustu"};
   TTreeReaderValue<Double_t> bb_gem_m2_clust_nclustu_tot = {rd, "bb.gem.m2.clust.nclustu_tot"};
   TTreeReaderValue<Double_t> bb_gem_m2_clust_nclustv = {rd, "bb.gem.m2.clust.nclustv"};
   TTreeReaderValue<Double_t> bb_gem_m2_clust_nclustv_tot = {rd, "bb.gem.m2.clust.nclustv_tot"};
   TTreeReaderValue<Double_t> bb_gem_m2_strip_nstrips_keep = {rd, "bb.gem.m2.strip.nstrips_keep"};
   TTreeReaderValue<Double_t> bb_gem_m2_strip_nstrips_keepU = {rd, "bb.gem.m2.strip.nstrips_keepU"};
   TTreeReaderValue<Double_t> bb_gem_m2_strip_nstrips_keepV = {rd, "bb.gem.m2.strip.nstrips_keepV"};
   TTreeReaderValue<Double_t> bb_gem_m2_strip_nstrips_keep_lmax = {rd, "bb.gem.m2.strip.nstrips_keep_lmax"};
   TTreeReaderValue<Double_t> bb_gem_m2_strip_nstrips_keep_lmaxU = {rd, "bb.gem.m2.strip.nstrips_keep_lmaxU"};
   TTreeReaderValue<Double_t> bb_gem_m2_strip_nstrips_keep_lmaxV = {rd, "bb.gem.m2.strip.nstrips_keep_lmaxV"};
   TTreeReaderValue<Double_t> bb_gem_m2_strip_nstripsfired = {rd, "bb.gem.m2.strip.nstripsfired"};
   TTreeReaderValue<Double_t> bb_gem_m3_clust_nclustu = {rd, "bb.gem.m3.clust.nclustu"};
   TTreeReaderValue<Double_t> bb_gem_m3_clust_nclustu_tot = {rd, "bb.gem.m3.clust.nclustu_tot"};
   TTreeReaderValue<Double_t> bb_gem_m3_clust_nclustv = {rd, "bb.gem.m3.clust.nclustv"};
   TTreeReaderValue<Double_t> bb_gem_m3_clust_nclustv_tot = {rd, "bb.gem.m3.clust.nclustv_tot"};
   TTreeReaderValue<Double_t> bb_gem_m3_strip_nstrips_keep = {rd, "bb.gem.m3.strip.nstrips_keep"};
   TTreeReaderValue<Double_t> bb_gem_m3_strip_nstrips_keepU = {rd, "bb.gem.m3.strip.nstrips_keepU"};
   TTreeReaderValue<Double_t> bb_gem_m3_strip_nstrips_keepV = {rd, "bb.gem.m3.strip.nstrips_keepV"};
   TTreeReaderValue<Double_t> bb_gem_m3_strip_nstrips_keep_lmax = {rd, "bb.gem.m3.strip.nstrips_keep_lmax"};
   TTreeReaderValue<Double_t> bb_gem_m3_strip_nstrips_keep_lmaxU = {rd, "bb.gem.m3.strip.nstrips_keep_lmaxU"};
   TTreeReaderValue<Double_t> bb_gem_m3_strip_nstrips_keep_lmaxV = {rd, "bb.gem.m3.strip.nstrips_keep_lmaxV"};
   TTreeReaderValue<Double_t> bb_gem_m3_strip_nstripsfired = {rd, "bb.gem.m3.strip.nstripsfired"};
   TTreeReaderValue<Double_t> bb_gem_m4_clust_nclustu = {rd, "bb.gem.m4.clust.nclustu"};
   TTreeReaderValue<Double_t> bb_gem_m4_clust_nclustu_tot = {rd, "bb.gem.m4.clust.nclustu_tot"};
   TTreeReaderValue<Double_t> bb_gem_m4_clust_nclustv = {rd, "bb.gem.m4.clust.nclustv"};
   TTreeReaderValue<Double_t> bb_gem_m4_clust_nclustv_tot = {rd, "bb.gem.m4.clust.nclustv_tot"};
   TTreeReaderValue<Double_t> bb_gem_m4_strip_nstrips_keep = {rd, "bb.gem.m4.strip.nstrips_keep"};
   TTreeReaderValue<Double_t> bb_gem_m4_strip_nstrips_keepU = {rd, "bb.gem.m4.strip.nstrips_keepU"};
   TTreeReaderValue<Double_t> bb_gem_m4_strip_nstrips_keepV = {rd, "bb.gem.m4.strip.nstrips_keepV"};
   TTreeReaderValue<Double_t> bb_gem_m4_strip_nstrips_keep_lmax = {rd, "bb.gem.m4.strip.nstrips_keep_lmax"};
   TTreeReaderValue<Double_t> bb_gem_m4_strip_nstrips_keep_lmaxU = {rd, "bb.gem.m4.strip.nstrips_keep_lmaxU"};
   TTreeReaderValue<Double_t> bb_gem_m4_strip_nstrips_keep_lmaxV = {rd, "bb.gem.m4.strip.nstrips_keep_lmaxV"};
   TTreeReaderValue<Double_t> bb_gem_m4_strip_nstripsfired = {rd, "bb.gem.m4.strip.nstripsfired"};
   TTreeReaderValue<Double_t> bb_gem_m5_clust_nclustu = {rd, "bb.gem.m5.clust.nclustu"};
   TTreeReaderValue<Double_t> bb_gem_m5_clust_nclustu_tot = {rd, "bb.gem.m5.clust.nclustu_tot"};
   TTreeReaderValue<Double_t> bb_gem_m5_clust_nclustv = {rd, "bb.gem.m5.clust.nclustv"};
   TTreeReaderValue<Double_t> bb_gem_m5_clust_nclustv_tot = {rd, "bb.gem.m5.clust.nclustv_tot"};
   TTreeReaderValue<Double_t> bb_gem_m5_strip_nstrips_keep = {rd, "bb.gem.m5.strip.nstrips_keep"};
   TTreeReaderValue<Double_t> bb_gem_m5_strip_nstrips_keepU = {rd, "bb.gem.m5.strip.nstrips_keepU"};
   TTreeReaderValue<Double_t> bb_gem_m5_strip_nstrips_keepV = {rd, "bb.gem.m5.strip.nstrips_keepV"};
   TTreeReaderValue<Double_t> bb_gem_m5_strip_nstrips_keep_lmax = {rd, "bb.gem.m5.strip.nstrips_keep_lmax"};
   TTreeReaderValue<Double_t> bb_gem_m5_strip_nstrips_keep_lmaxU = {rd, "bb.gem.m5.strip.nstrips_keep_lmaxU"};
   TTreeReaderValue<Double_t> bb_gem_m5_strip_nstrips_keep_lmaxV = {rd, "bb.gem.m5.strip.nstrips_keep_lmaxV"};
   TTreeReaderValue<Double_t> bb_gem_m5_strip_nstripsfired = {rd, "bb.gem.m5.strip.nstripsfired"};
   TTreeReaderValue<Double_t> bb_gem_m6_clust_nclustu = {rd, "bb.gem.m6.clust.nclustu"};
   TTreeReaderValue<Double_t> bb_gem_m6_clust_nclustu_tot = {rd, "bb.gem.m6.clust.nclustu_tot"};
   TTreeReaderValue<Double_t> bb_gem_m6_clust_nclustv = {rd, "bb.gem.m6.clust.nclustv"};
   TTreeReaderValue<Double_t> bb_gem_m6_clust_nclustv_tot = {rd, "bb.gem.m6.clust.nclustv_tot"};
   TTreeReaderValue<Double_t> bb_gem_m6_strip_nstrips_keep = {rd, "bb.gem.m6.strip.nstrips_keep"};
   TTreeReaderValue<Double_t> bb_gem_m6_strip_nstrips_keepU = {rd, "bb.gem.m6.strip.nstrips_keepU"};
   TTreeReaderValue<Double_t> bb_gem_m6_strip_nstrips_keepV = {rd, "bb.gem.m6.strip.nstrips_keepV"};
   TTreeReaderValue<Double_t> bb_gem_m6_strip_nstrips_keep_lmax = {rd, "bb.gem.m6.strip.nstrips_keep_lmax"};
   TTreeReaderValue<Double_t> bb_gem_m6_strip_nstrips_keep_lmaxU = {rd, "bb.gem.m6.strip.nstrips_keep_lmaxU"};
   TTreeReaderValue<Double_t> bb_gem_m6_strip_nstrips_keep_lmaxV = {rd, "bb.gem.m6.strip.nstrips_keep_lmaxV"};
   TTreeReaderValue<Double_t> bb_gem_m6_strip_nstripsfired = {rd, "bb.gem.m6.strip.nstripsfired"};
   TTreeReaderValue<Double_t> bb_gem_m7_clust_nclustu = {rd, "bb.gem.m7.clust.nclustu"};
   TTreeReaderValue<Double_t> bb_gem_m7_clust_nclustu_tot = {rd, "bb.gem.m7.clust.nclustu_tot"};
   TTreeReaderValue<Double_t> bb_gem_m7_clust_nclustv = {rd, "bb.gem.m7.clust.nclustv"};
   TTreeReaderValue<Double_t> bb_gem_m7_clust_nclustv_tot = {rd, "bb.gem.m7.clust.nclustv_tot"};
   TTreeReaderValue<Double_t> bb_gem_m7_strip_nstrips_keep = {rd, "bb.gem.m7.strip.nstrips_keep"};
   TTreeReaderValue<Double_t> bb_gem_m7_strip_nstrips_keepU = {rd, "bb.gem.m7.strip.nstrips_keepU"};
   TTreeReaderValue<Double_t> bb_gem_m7_strip_nstrips_keepV = {rd, "bb.gem.m7.strip.nstrips_keepV"};
   TTreeReaderValue<Double_t> bb_gem_m7_strip_nstrips_keep_lmax = {rd, "bb.gem.m7.strip.nstrips_keep_lmax"};
   TTreeReaderValue<Double_t> bb_gem_m7_strip_nstrips_keep_lmaxU = {rd, "bb.gem.m7.strip.nstrips_keep_lmaxU"};
   TTreeReaderValue<Double_t> bb_gem_m7_strip_nstrips_keep_lmaxV = {rd, "bb.gem.m7.strip.nstrips_keep_lmaxV"};
   TTreeReaderValue<Double_t> bb_gem_m7_strip_nstripsfired = {rd, "bb.gem.m7.strip.nstripsfired"};
   TTreeReaderValue<Double_t> bb_gem_nlayershit = {rd, "bb.gem.nlayershit"};
   TTreeReaderValue<Double_t> bb_gem_nlayershitu = {rd, "bb.gem.nlayershitu"};
   TTreeReaderValue<Double_t> bb_gem_nlayershituv = {rd, "bb.gem.nlayershituv"};
   TTreeReaderValue<Double_t> bb_gem_nlayershitv = {rd, "bb.gem.nlayershitv"};
   TTreeReaderValue<Double_t> bb_gem_track_besttrack = {rd, "bb.gem.track.besttrack"};
   TTreeReaderValue<Double_t> bb_gem_track_ntrack = {rd, "bb.gem.track.ntrack"};
   
   TTreeReaderValue<Double_t> bb_hodotdc_nclus = {rd, "bb.hodotdc.nclus"};
   TTreeReaderValue<Double_t> bb_ps_againblk = {rd, "bb.ps.againblk"};
   TTreeReaderValue<Double_t> bb_ps_colblk = {rd, "bb.ps.colblk"};
   TTreeReaderValue<Double_t> bb_ps_e = {rd, "bb.ps.e"};
   TTreeReaderValue<Double_t> bb_ps_eblk = {rd, "bb.ps.eblk"};
   TTreeReaderValue<Double_t> bb_ps_idblk = {rd, "bb.ps.idblk"};
   TTreeReaderValue<Double_t> bb_ps_index = {rd, "bb.ps.index"};
   TTreeReaderValue<Double_t> bb_ps_nblk = {rd, "bb.ps.nblk"};
   TTreeReaderValue<Double_t> bb_ps_nclus = {rd, "bb.ps.nclus"};
   TTreeReaderValue<Double_t> bb_ps_ngoodADChits = {rd, "bb.ps.ngoodADChits"};
   TTreeReaderValue<Double_t> bb_ps_rowblk = {rd, "bb.ps.rowblk"};
   TTreeReaderValue<Double_t> bb_ps_x = {rd, "bb.ps.x"};
   TTreeReaderValue<Double_t> bb_ps_y = {rd, "bb.ps.y"};
   TTreeReaderValue<Double_t> bb_sh_againblk = {rd, "bb.sh.againblk"};
   TTreeReaderValue<Double_t> bb_sh_colblk = {rd, "bb.sh.colblk"};
   TTreeReaderValue<Double_t> bb_sh_e = {rd, "bb.sh.e"};
   TTreeReaderValue<Double_t> bb_sh_eblk = {rd, "bb.sh.eblk"};
   TTreeReaderValue<Double_t> bb_sh_idblk = {rd, "bb.sh.idblk"};
   TTreeReaderValue<Double_t> bb_sh_index = {rd, "bb.sh.index"};
   TTreeReaderValue<Double_t> bb_sh_nblk = {rd, "bb.sh.nblk"};
   TTreeReaderValue<Double_t> bb_sh_nclus = {rd, "bb.sh.nclus"};
   TTreeReaderValue<Double_t> bb_sh_ngoodADChits = {rd, "bb.sh.ngoodADChits"};
   TTreeReaderValue<Double_t> bb_sh_rowblk = {rd, "bb.sh.rowblk"};
   TTreeReaderValue<Double_t> bb_sh_x = {rd, "bb.sh.x"};
   TTreeReaderValue<Double_t> bb_sh_y = {rd, "bb.sh.y"};
   TTreeReaderValue<Double_t> bb_tr_n = {rd, "bb.tr.n"};
   TTreeReaderValue<Double_t> bb_ts_over_threshold = {rd, "bb.ts.over_threshold"};
   TTreeReaderValue<Double_t> e_kine_omega = {rd, "e.kine.omega"};
   TTreeReaderValue<Double_t> e_kine_ph_q = {rd, "e.kine.ph_q"};
   TTreeReaderValue<Double_t> e_kine_q3m = {rd, "e.kine.q3m"};
   TTreeReaderValue<Double_t> e_kine_q_x = {rd, "e.kine.q_x"};
   TTreeReaderValue<Double_t> e_kine_q_y = {rd, "e.kine.q_y"};
   TTreeReaderValue<Double_t> e_kine_q_z = {rd, "e.kine.q_z"};
   TTreeReaderValue<Double_t> e_kine_th_q = {rd, "e.kine.th_q"};
   TTreeReaderValue<Double_t> e_kine_x_bj = {rd, "e.kine.x_bj"};
   TTreeReaderValue<Double_t> g_datatype = {rd, "g.datatype"};
   TTreeReaderValue<Double_t> g_evlen = {rd, "g.evlen"};
   TTreeReaderValue<Double_t> g_evnum = {rd, "g.evnum"};
   TTreeReaderValue<Double_t> g_evtime = {rd, "g.evtime"};
   TTreeReaderValue<Double_t> g_evtyp = {rd, "g.evtyp"};
   TTreeReaderValue<Double_t> g_runnum = {rd, "g.runnum"};
   TTreeReaderValue<Double_t> g_runtime = {rd, "g.runtime"};
   TTreeReaderValue<Double_t> g_runtype = {rd, "g.runtype"};
   TTreeReaderValue<Double_t> g_trigbits = {rd, "g.trigbits"};
   TTreeReaderValue<Double_t> sbs_hcal_againblk = {rd, "sbs.hcal.againblk"};
   TTreeReaderValue<Double_t> sbs_hcal_eblk = {rd, "sbs.hcal.eblk"};
   TTreeReaderValue<Double_t> sbs_hcal_idblk = {rd, "sbs.hcal.idblk"};
   TTreeReaderValue<Double_t> sbs_hcal_index = {rd, "sbs.hcal.index"};
   TTreeReaderValue<Double_t> sbs_hcal_nblk = {rd, "sbs.hcal.nblk"};
   TTreeReaderValue<Double_t> sbs_hcal_nclus = {rd, "sbs.hcal.nclus"};
   TTreeReaderValue<Double_t> sbs_hcal_ngoodADChits = {rd, "sbs.hcal.ngoodADChits"};
   TTreeReaderValue<Double_t> sbs_tr_n = {rd, "sbs.tr.n"};
   TTreeReaderValue<Double_t> scalhel_errcode = {rd, "scalhel.errcode"};
   TTreeReaderValue<Double_t> scalhel_evtcount = {rd, "scalhel.evtcount"};
   TTreeReaderValue<Double_t> scalhel_hel = {rd, "scalhel.hel"};
   TTreeReaderValue<Double_t> scalhel_lhrs_fadc_hel = {rd, "scalhel.lhrs.fadc.hel"};
   TTreeReaderValue<Double_t> scalhel_lhrs_fadc_pat = {rd, "scalhel.lhrs.fadc.pat"};
   TTreeReaderValue<Double_t> scalhel_lhrs_fadc_tsettle = {rd, "scalhel.lhrs.fadc.tsettle"};
   TTreeReaderValue<Double_t> scalhel_patcount = {rd, "scalhel.patcount"};
   TTreeReaderValue<Double_t> scalhel_patphase = {rd, "scalhel.patphase"};
   TTreeReaderValue<Double_t> scalhel_seed = {rd, "scalhel.seed"};
   TTreeReaderValue<Double_t> singletrack_bb = {rd, "singletrack_bb"};
   TTreeReaderValue<Double_t> anytrack_bb = {rd, "anytrack_bb"};
   TTreeReaderValue<Double_t> singletrack = {rd, "singletrack"};
   TTreeReaderValue<Double_t> anytrack = {rd, "anytrack"};
   TTreeReaderValue<Double_t> HALLA_p = {rd, "HALLA_p"};
   TTreeReaderValue<Double_t> hac_bcm_average = {rd, "hac_bcm_average"};
   TTreeReaderValue<Double_t> IPM1H04A_XPOS = {rd, "IPM1H04A.XPOS"};
   TTreeReaderValue<Double_t> IPM1H04A_YPOS = {rd, "IPM1H04A.YPOS"};
   TTreeReaderValue<Double_t> IPM1H04E_XPOS = {rd, "IPM1H04E.XPOS"};
   TTreeReaderValue<Double_t> IPM1H04E_YPOS = {rd, "IPM1H04E.YPOS"};
   TTreeReaderValue<Double_t> IGL1I00OD16_16 = {rd, "IGL1I00OD16_16"};
   */
  //end of reader
  
  
  //setup output file
  TString outfilename = Form("~/vol/data/out_tof_gen%i_%s.root",kin_no,Target.Data());
  TFile *outfile = new TFile(outfilename,"RECREATE");
  
  TTree *tree = new TTree("T","Output of GEN tof Skim");
  double Q2, W2, eps, nu, x_bj;
  double bb_x,bb_y,bb_px,bb_py,bb_pz,bb_p,bb_vx,bb_vy,bb_vz,bb_th_fp,bb_ph_fp,bb_th_tg,bb_ph_tg,bb_y_tg;
  double bb_pathl;
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
  
  double bb_tr_nhits, bb_ntr, trigtime;
  tree->Branch("bb_tr_nhits",&bb_tr_nhits);
  tree->Branch("bb_ntr",&bb_ntr);
  tree->Branch("trigtime",&trigtime);

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
  
  double pred_ang_horiz, pred_ang_vert, pred_y, pred_x, pred_mom, dx, dy, dev;
  tree->Branch("pred_ang_horiz",&pred_ang_horiz);
  tree->Branch("pred_ang_vert",&pred_ang_vert);
  tree->Branch("pred_y",&pred_y);
  tree->Branch("pred_x",&pred_x);
  tree->Branch("pred_mom",&pred_mom);
  tree->Branch("dx",&dx);
  tree->Branch("dy",&dy);
  tree->Branch("dev",&dev);
  
  double ps_atime, sh_atime, hcal_atime, coin_atime, coin_atime_corr;
  tree->Branch("ps_atime",&ps_atime);
  tree->Branch("sh_atime",&sh_atime);
  tree->Branch("hcal_atime",&hcal_atime); 
  tree->Branch("coin_atime",&coin_atime);
  //tree->Branch("coin_atime_corr",&coin_atime_corr);
  
  double hcal_time, coin_time, coin_time_corr;
  tree->Branch("hcal_time",&hcal_time);
  tree->Branch("coin_time",&coin_time);
  //tree->Branch("coin_time_corr",&coin_time_corr);
  
  vector<double> hodo_tmean;
  vector<double> hodo_tleft;
  vector<double> hodo_tright;
  vector<double> hodo_totleft;
  vector<double> hodo_totright;
  vector<double> hodo_timehitpos;
  vector<double> hodo_bar;
  tree->Branch("hodo_tmean",&hodo_tmean);
  tree->Branch("hodo_tleft",&hodo_tleft);
  tree->Branch("hodo_tright",&hodo_tright);
  tree->Branch("hodo_totleft",&hodo_totleft);
  tree->Branch("hodo_totright",&hodo_totright);
  tree->Branch("hodo_timehitpos",&hodo_timehitpos);
  tree->Branch("hodo_bar",&hodo_bar);
  
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

  //unique event id for using brufit splot fitting
  Long64_t UID=-1;
  tree->Branch("UID",&UID);
  
  //Data loop
  Long64_t ev=-1;
  //nev = 1000000;
  int treenum = 0, currenttreenum = 0, currentrunnum = 0;
  time_t run_time_unix;
  cout << "Processing " << nev << " events." << endl;
  
  while(rd.Next() && ev < nev){
    
    ev = rd.GetCurrentEntry();
    UID += 1;
    
    if( ev % 100000 == 0 )
      cout << ev << " / " << nev << endl;
    
    // Pre-cuts -- these should be in the replay cdef
    if( *bb_tr_n <= 0 ) continue; 
    
    //hcal 2 arm cuts
    if ( *sbs_hcal_nclus <= 0) continue;
    
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

    trigtime = *bb_gem_trigtime;
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
    
    //set variables in sbs
    hcal_e = *sbs_hcal_e;
    hcal_eblk = *sbs_hcal_eblk;
    hcal_x = *sbs_hcal_x;
    hcal_y = *sbs_hcal_y;
    hcal_row = *sbs_hcal_rowblk;
    hcal_col = *sbs_hcal_colblk;
    hcal_id = *sbs_hcal_idblk;

    //timing variables
    ps_atime = *bb_ps_atimeblk;
    sh_atime = *bb_sh_atimeblk;
    hcal_atime = *sbs_hcal_atimeblk;
    coin_atime = hcal_atime - sh_atime;
    //int hcal_rownum = (int) *sbs_hcal_rowblk;
    //coin_atime_corr = coin_atime - hcal_atcoin_means[hcal_rownum];
    hcal_time = *sbs_hcal_tdctimeblk;
    coin_time = hcal_time - bb_hodotdc_clus_bar_tdc_meantime[0];
    //int hodo_barid = (int) bb_hodotdc_clus_bar_tdc_id[0];
    //coin_time_corr = coin_time - hodo_tcoin_means[hodo_barid]; 
    int nhodo = *Ndata_bb_hodotdc_clus_bar_tdc_id;
    hodo_tmean.clear();
    hodo_tleft.clear();
    hodo_tright.clear();
    hodo_totleft.clear();
    hodo_totright.clear();
    hodo_timehitpos.clear();
    hodo_bar.clear();
    
    for (int i=0; i<nhodo; i++){
      hodo_tmean.push_back(bb_hodotdc_clus_bar_tdc_meantime[i]);
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
    for (int i=0; i<ntdctrig; i++){
      tdctrig_id.push_back(bb_tdctrig_tdcelemID[i]);
      tdctrig.push_back(bb_tdctrig_tdc[i]);
    }
    
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
  
  cout <<"Saving rootfiles" << endl;
  outfile->cd();
  outfile->Write();
  outfile->Close();
  t.Stop();
  cout << "Total Time: ";
  t.Print();
  cout << " s" << endl;
  
}//End of macro
  
