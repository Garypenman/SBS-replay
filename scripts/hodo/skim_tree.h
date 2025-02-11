//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Nov 28 19:15:35 2024 by ROOT version 6.30/04
// from TTree T/T
// found on file: /home/gpenman/vol/data/out_gen2_He3_qe.root
//////////////////////////////////////////////////////////

#ifndef skim_tree_h
#define skim_tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "ROOT/RVec.hxx"
#include "ROOT/RVec.hxx"

class skim_tree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        BPMAx;
   Double_t        BPMAy;
   Double_t        BeamPol;
   Double_t        He3Pol;
   Double_t        He3PolErr;
   Double_t        Q2;
   Double_t        Rasterx;
   Double_t        Rasterx2;
   Double_t        Rastery;
   Double_t        Rastery2;
   Long64_t        UID;
   Double_t        W2;
   Double_t        bb_ntr;
   Double_t        bb_p;
   Double_t        bb_pathl;
   Double_t        bb_ph_fp;
   Double_t        bb_ph_tg;
   Double_t        bb_px;
   Double_t        bb_py;
   Double_t        bb_pz;
   Double_t        bb_rftime;
   Double_t        bb_th_fp;
   Double_t        bb_th_tg;
   Double_t        bb_tr_nhits;
   Double_t        bb_trigtime;
   Double_t        bb_vx;
   Double_t        bb_vy;
   Double_t        bb_vz;
   Double_t        bb_x;
   Double_t        bb_y;
   Double_t        bb_y_tg;
   Double_t        bbcal_xdiff;
   Double_t        bbcal_ydiff;
   Double_t        cal_sum;
   Double_t        coin_atime;
   Double_t        coin_time;
   Double_t        dev;
   Double_t        dx;
   Double_t        dy;
   Double_t        edivp;
   Double_t        eps;
   Int_t           evt_hel;
   UInt_t          evt_len;
   UInt_t          evt_num;
   ULong64_t       evt_time;
   UInt_t          evt_type;
   Double_t        gr_adc;
   Double_t        gr_size;
   Double_t        gr_tmean;
   Double_t        gr_totmean;
   Double_t        gr_trindex;
   Double_t        gr_xmean;
   Double_t        gr_ymean;
   Double_t        hcal_adctime;
   Double_t        hcal_atimeblk;
   ROOT::VecOps::RVec<double> *hcal_blk_atime;
   ROOT::VecOps::RVec<double> *hcal_blk_col;
   ROOT::VecOps::RVec<double> *hcal_blk_e;
   ROOT::VecOps::RVec<double> *hcal_blk_id;
   ROOT::VecOps::RVec<double> *hcal_blk_row;
   ROOT::VecOps::RVec<double> *hcal_blk_tdctime;
   ROOT::VecOps::RVec<double> *hcal_blk_tdctime_c;
   ROOT::VecOps::RVec<double> *hcal_blk_x;
   ROOT::VecOps::RVec<double> *hcal_blk_y;
   Double_t        hcal_col;
   Double_t        hcal_e;
   Double_t        hcal_id;
   ROOT::VecOps::RVec<double> *hcal_ref;
   ROOT::VecOps::RVec<double> *hcal_refid;
   Double_t        hcal_row;
   Double_t        hcal_tdctime;
   Double_t        hcal_tdctime_c;
   Double_t        hcal_tdctimeblk;
   Double_t        hcal_tdctimeblk_c;
   Double_t        hcal_x;
   Double_t        hcal_xdiff;
   Double_t        hcal_y;
   Double_t        hcal_ydiff;
   Double_t        helicity;
   Int_t           helicity_status;
   ROOT::VecOps::RVec<double> *hodo_bar;
   Double_t        hodo_nbar;
   Double_t        hodo_nclus;
   Double_t        hodo_refleft;
   Double_t        hodo_refright;
   ROOT::VecOps::RVec<double> *hodo_tdiff;
   ROOT::VecOps::RVec<double> *hodo_tdiff_c;
   ROOT::VecOps::RVec<double> *hodo_timehitpos;
   ROOT::VecOps::RVec<double> *hodo_tleft;
   ROOT::VecOps::RVec<double> *hodo_tleft_c;
   ROOT::VecOps::RVec<double> *hodo_tmean;
   ROOT::VecOps::RVec<double> *hodo_tmean_c;
   ROOT::VecOps::RVec<double> *hodo_totleft;
   ROOT::VecOps::RVec<double> *hodo_totmean;
   ROOT::VecOps::RVec<double> *hodo_totright;
   ROOT::VecOps::RVec<double> *hodo_tright;
   ROOT::VecOps::RVec<double> *hodo_tright_c;
   Double_t        nu;
   Double_t        p_elas;
   Double_t        pdiff;
   Double_t        pexp_th;
   Double_t        ph_neutron;
   Double_t        ph_proton;
   Double_t        phi_e;
   Double_t        pred_ang_horiz;
   Double_t        pred_ang_vert;
   Double_t        pred_beta;
   Double_t        pred_mom;
   Double_t        pred_x;
   Double_t        pred_y;
   Double_t        proton_deflection;
   Double_t        ps_atime;
   Double_t        ps_e;
   Double_t        ps_ep;
   Double_t        ps_x;
   Double_t        ps_y;
   UInt_t          run_num;
   Double_t        sbs_pathl;
   Double_t        sbs_rftime;
   Double_t        sbs_trigtime;
   Double_t        sh_atime;
   Double_t        sh_e;
   Double_t        sh_ep;
   Double_t        sh_x;
   Double_t        sh_y;
   Double_t        tau;
   ROOT::VecOps::RVec<double> *tdctrig;
   ROOT::VecOps::RVec<int> *tdctrig_id;
   Double_t        th_lab;
   Double_t        th_lab_deg;
   Double_t        th_neutron;
   Double_t        th_pq_neutron;
   Double_t        th_pq_proton;
   Double_t        th_proton;
   Double_t        theta_e;
   Double_t        thetabend;
   UInt_t          trigbits;
   Double_t        trx_sh;
   Double_t        try_sh;
   Double_t        x_bj;

   // List of branches
   TBranch        *b_BPMAx;   //!
   TBranch        *b_BPMAy;   //!
   TBranch        *b_BeamPol;   //!
   TBranch        *b_He3Pol;   //!
   TBranch        *b_He3PolErr;   //!
   TBranch        *b_Q2;   //!
   TBranch        *b_Rasterx;   //!
   TBranch        *b_Rasterx2;   //!
   TBranch        *b_Rastery;   //!
   TBranch        *b_Rastery2;   //!
   TBranch        *b_UID;   //!
   TBranch        *b_W2;   //!
   TBranch        *b_bb_ntr;   //!
   TBranch        *b_bb_p;   //!
   TBranch        *b_bb_pathl;   //!
   TBranch        *b_bb_ph_fp;   //!
   TBranch        *b_bb_ph_tg;   //!
   TBranch        *b_bb_px;   //!
   TBranch        *b_bb_py;   //!
   TBranch        *b_bb_pz;   //!
   TBranch        *b_bb_rftime;   //!
   TBranch        *b_bb_th_fp;   //!
   TBranch        *b_bb_th_tg;   //!
   TBranch        *b_bb_tr_nhits;   //!
   TBranch        *b_bb_trigtime;   //!
   TBranch        *b_bb_vx;   //!
   TBranch        *b_bb_vy;   //!
   TBranch        *b_bb_vz;   //!
   TBranch        *b_bb_x;   //!
   TBranch        *b_bb_y;   //!
   TBranch        *b_bb_y_tg;   //!
   TBranch        *b_bbcal_xdiff;   //!
   TBranch        *b_bbcal_ydiff;   //!
   TBranch        *b_cal_sum;   //!
   TBranch        *b_coin_atime;   //!
   TBranch        *b_coin_time;   //!
   TBranch        *b_dev;   //!
   TBranch        *b_dx;   //!
   TBranch        *b_dy;   //!
   TBranch        *b_edivp;   //!
   TBranch        *b_eps;   //!
   TBranch        *b_evt_hel;   //!
   TBranch        *b_evt_len;   //!
   TBranch        *b_evt_num;   //!
   TBranch        *b_evt_time;   //!
   TBranch        *b_evt_type;   //!
   TBranch        *b_gr_adc;   //!
   TBranch        *b_gr_size;   //!
   TBranch        *b_gr_tmean;   //!
   TBranch        *b_gr_totmean;   //!
   TBranch        *b_gr_trindex;   //!
   TBranch        *b_gr_xmean;   //!
   TBranch        *b_gr_ymean;   //!
   TBranch        *b_hcal_adctime;   //!
   TBranch        *b_hcal_atimeblk;   //!
   TBranch        *b_hcal_blk_atime;   //!
   TBranch        *b_hcal_blk_col;   //!
   TBranch        *b_hcal_blk_e;   //!
   TBranch        *b_hcal_blk_id;   //!
   TBranch        *b_hcal_blk_row;   //!
   TBranch        *b_hcal_blk_tdctime;   //!
   TBranch        *b_hcal_blk_tdctime_c;   //!
   TBranch        *b_hcal_blk_x;   //!
   TBranch        *b_hcal_blk_y;   //!
   TBranch        *b_hcal_col;   //!
   TBranch        *b_hcal_e;   //!
   TBranch        *b_hcal_id;   //!
   TBranch        *b_hcal_ref;   //!
   TBranch        *b_hcal_refid;   //!
   TBranch        *b_hcal_row;   //!
   TBranch        *b_hcal_tdctime;   //!
   TBranch        *b_hcal_tdctime_c;   //!
   TBranch        *b_hcal_tdctimeblk;   //!
   TBranch        *b_hcal_tdctimeblk_c;   //!
   TBranch        *b_hcal_x;   //!
   TBranch        *b_hcal_xdiff;   //!
   TBranch        *b_hcal_y;   //!
   TBranch        *b_hcal_ydiff;   //!
   TBranch        *b_helicity;   //!
   TBranch        *b_helicity_status;   //!
   TBranch        *b_hodo_bar;   //!
   TBranch        *b_hodo_nbar;   //!
   TBranch        *b_hodo_nclus;   //!
   TBranch        *b_hodo_refleft;   //!
   TBranch        *b_hodo_refright;   //!
   TBranch        *b_hodo_tdiff;   //!
   TBranch        *b_hodo_tdiff_c;   //!
   TBranch        *b_hodo_timehitpos;   //!
   TBranch        *b_hodo_tleft;   //!
   TBranch        *b_hodo_tleft_c;   //!
   TBranch        *b_hodo_tmean;   //!
   TBranch        *b_hodo_tmean_c;   //!
   TBranch        *b_hodo_totleft;   //!
   TBranch        *b_hodo_totmean;   //!
   TBranch        *b_hodo_totright;   //!
   TBranch        *b_hodo_tright;   //!
   TBranch        *b_hodo_tright_c;   //!
   TBranch        *b_nu;   //!
   TBranch        *b_p_elas;   //!
   TBranch        *b_pdiff;   //!
   TBranch        *b_pexp_th;   //!
   TBranch        *b_ph_neutron;   //!
   TBranch        *b_ph_proton;   //!
   TBranch        *b_phi_e;   //!
   TBranch        *b_pred_ang_horiz;   //!
   TBranch        *b_pred_ang_vert;   //!
   TBranch        *b_pred_beta;   //!
   TBranch        *b_pred_mom;   //!
   TBranch        *b_pred_x;   //!
   TBranch        *b_pred_y;   //!
   TBranch        *b_proton_deflection;   //!
   TBranch        *b_ps_atime;   //!
   TBranch        *b_ps_e;   //!
   TBranch        *b_ps_ep;   //!
   TBranch        *b_ps_x;   //!
   TBranch        *b_ps_y;   //!
   TBranch        *b_run_num;   //!
   TBranch        *b_sbs_pathl;   //!
   TBranch        *b_sbs_rftime;   //!
   TBranch        *b_sbs_trigtime;   //!
   TBranch        *b_sh_atime;   //!
   TBranch        *b_sh_e;   //!
   TBranch        *b_sh_ep;   //!
   TBranch        *b_sh_x;   //!
   TBranch        *b_sh_y;   //!
   TBranch        *b_tau;   //!
   TBranch        *b_tdctrig;   //!
   TBranch        *b_tdctrig_id;   //!
   TBranch        *b_th_lab;   //!
   TBranch        *b_th_lab_deg;   //!
   TBranch        *b_th_neutron;   //!
   TBranch        *b_th_pq_neutron;   //!
   TBranch        *b_th_pq_proton;   //!
   TBranch        *b_th_proton;   //!
   TBranch        *b_theta_e;   //!
   TBranch        *b_thetabend;   //!
   TBranch        *b_trigbits;   //!
   TBranch        *b_trx_sh;   //!
   TBranch        *b_try_sh;   //!
   TBranch        *b_x_bj;   //!

   skim_tree(TTree *tree=0);
   virtual ~skim_tree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef skim_tree_cxx
skim_tree::skim_tree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/home/gpenman/vol/data/out_gen2_He3_qe.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/home/gpenman/vol/data/out_gen2_He3_qe.root");
      }
      f->GetObject("T",tree);

   }
   Init(tree);
}

skim_tree::~skim_tree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t skim_tree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t skim_tree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void skim_tree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   hcal_blk_atime = 0;
   hcal_blk_col = 0;
   hcal_blk_e = 0;
   hcal_blk_id = 0;
   hcal_blk_row = 0;
   hcal_blk_tdctime = 0;
   hcal_blk_tdctime_c = 0;
   hcal_blk_x = 0;
   hcal_blk_y = 0;
   hcal_ref = 0;
   hcal_refid = 0;
   hodo_bar = 0;
   hodo_tdiff = 0;
   hodo_tdiff_c = 0;
   hodo_timehitpos = 0;
   hodo_tleft = 0;
   hodo_tleft_c = 0;
   hodo_tmean = 0;
   hodo_tmean_c = 0;
   hodo_totleft = 0;
   hodo_totmean = 0;
   hodo_totright = 0;
   hodo_tright = 0;
   hodo_tright_c = 0;
   tdctrig = 0;
   tdctrig_id = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("BPMAx", &BPMAx, &b_BPMAx);
   fChain->SetBranchAddress("BPMAy", &BPMAy, &b_BPMAy);
   fChain->SetBranchAddress("BeamPol", &BeamPol, &b_BeamPol);
   fChain->SetBranchAddress("He3Pol", &He3Pol, &b_He3Pol);
   fChain->SetBranchAddress("He3PolErr", &He3PolErr, &b_He3PolErr);
   fChain->SetBranchAddress("Q2", &Q2, &b_Q2);
   fChain->SetBranchAddress("Rasterx", &Rasterx, &b_Rasterx);
   fChain->SetBranchAddress("Rasterx2", &Rasterx2, &b_Rasterx2);
   fChain->SetBranchAddress("Rastery", &Rastery, &b_Rastery);
   fChain->SetBranchAddress("Rastery2", &Rastery2, &b_Rastery2);
   fChain->SetBranchAddress("UID", &UID, &b_UID);
   fChain->SetBranchAddress("W2", &W2, &b_W2);
   fChain->SetBranchAddress("bb_ntr", &bb_ntr, &b_bb_ntr);
   fChain->SetBranchAddress("bb_p", &bb_p, &b_bb_p);
   fChain->SetBranchAddress("bb_pathl", &bb_pathl, &b_bb_pathl);
   fChain->SetBranchAddress("bb_ph_fp", &bb_ph_fp, &b_bb_ph_fp);
   fChain->SetBranchAddress("bb_ph_tg", &bb_ph_tg, &b_bb_ph_tg);
   fChain->SetBranchAddress("bb_px", &bb_px, &b_bb_px);
   fChain->SetBranchAddress("bb_py", &bb_py, &b_bb_py);
   fChain->SetBranchAddress("bb_pz", &bb_pz, &b_bb_pz);
   fChain->SetBranchAddress("bb_rftime", &bb_rftime, &b_bb_rftime);
   fChain->SetBranchAddress("bb_th_fp", &bb_th_fp, &b_bb_th_fp);
   fChain->SetBranchAddress("bb_th_tg", &bb_th_tg, &b_bb_th_tg);
   fChain->SetBranchAddress("bb_tr_nhits", &bb_tr_nhits, &b_bb_tr_nhits);
   fChain->SetBranchAddress("bb_trigtime", &bb_trigtime, &b_bb_trigtime);
   fChain->SetBranchAddress("bb_vx", &bb_vx, &b_bb_vx);
   fChain->SetBranchAddress("bb_vy", &bb_vy, &b_bb_vy);
   fChain->SetBranchAddress("bb_vz", &bb_vz, &b_bb_vz);
   fChain->SetBranchAddress("bb_x", &bb_x, &b_bb_x);
   fChain->SetBranchAddress("bb_y", &bb_y, &b_bb_y);
   fChain->SetBranchAddress("bb_y_tg", &bb_y_tg, &b_bb_y_tg);
   fChain->SetBranchAddress("bbcal_xdiff", &bbcal_xdiff, &b_bbcal_xdiff);
   fChain->SetBranchAddress("bbcal_ydiff", &bbcal_ydiff, &b_bbcal_ydiff);
   fChain->SetBranchAddress("cal_sum", &cal_sum, &b_cal_sum);
   fChain->SetBranchAddress("coin_atime", &coin_atime, &b_coin_atime);
   fChain->SetBranchAddress("coin_time", &coin_time, &b_coin_time);
   fChain->SetBranchAddress("dev", &dev, &b_dev);
   fChain->SetBranchAddress("dx", &dx, &b_dx);
   fChain->SetBranchAddress("dy", &dy, &b_dy);
   fChain->SetBranchAddress("edivp", &edivp, &b_edivp);
   fChain->SetBranchAddress("eps", &eps, &b_eps);
   fChain->SetBranchAddress("evt_hel", &evt_hel, &b_evt_hel);
   fChain->SetBranchAddress("evt_len", &evt_len, &b_evt_len);
   fChain->SetBranchAddress("evt_num", &evt_num, &b_evt_num);
   fChain->SetBranchAddress("evt_time", &evt_time, &b_evt_time);
   fChain->SetBranchAddress("evt_type", &evt_type, &b_evt_type);
   fChain->SetBranchAddress("gr_adc", &gr_adc, &b_gr_adc);
   fChain->SetBranchAddress("gr_size", &gr_size, &b_gr_size);
   fChain->SetBranchAddress("gr_tmean", &gr_tmean, &b_gr_tmean);
   fChain->SetBranchAddress("gr_totmean", &gr_totmean, &b_gr_totmean);
   fChain->SetBranchAddress("gr_trindex", &gr_trindex, &b_gr_trindex);
   fChain->SetBranchAddress("gr_xmean", &gr_xmean, &b_gr_xmean);
   fChain->SetBranchAddress("gr_ymean", &gr_ymean, &b_gr_ymean);
   fChain->SetBranchAddress("hcal_adctime", &hcal_adctime, &b_hcal_adctime);
   fChain->SetBranchAddress("hcal_atimeblk", &hcal_atimeblk, &b_hcal_atimeblk);
   fChain->SetBranchAddress("hcal_blk_atime", &hcal_blk_atime, &b_hcal_blk_atime);
   fChain->SetBranchAddress("hcal_blk_col", &hcal_blk_col, &b_hcal_blk_col);
   fChain->SetBranchAddress("hcal_blk_e", &hcal_blk_e, &b_hcal_blk_e);
   fChain->SetBranchAddress("hcal_blk_id", &hcal_blk_id, &b_hcal_blk_id);
   fChain->SetBranchAddress("hcal_blk_row", &hcal_blk_row, &b_hcal_blk_row);
   fChain->SetBranchAddress("hcal_blk_tdctime", &hcal_blk_tdctime, &b_hcal_blk_tdctime);
   fChain->SetBranchAddress("hcal_blk_tdctime_c", &hcal_blk_tdctime_c, &b_hcal_blk_tdctime_c);
   fChain->SetBranchAddress("hcal_blk_x", &hcal_blk_x, &b_hcal_blk_x);
   fChain->SetBranchAddress("hcal_blk_y", &hcal_blk_y, &b_hcal_blk_y);
   fChain->SetBranchAddress("hcal_col", &hcal_col, &b_hcal_col);
   fChain->SetBranchAddress("hcal_e", &hcal_e, &b_hcal_e);
   fChain->SetBranchAddress("hcal_id", &hcal_id, &b_hcal_id);
   fChain->SetBranchAddress("hcal_ref", &hcal_ref, &b_hcal_ref);
   fChain->SetBranchAddress("hcal_refid", &hcal_refid, &b_hcal_refid);
   fChain->SetBranchAddress("hcal_row", &hcal_row, &b_hcal_row);
   fChain->SetBranchAddress("hcal_tdctime", &hcal_tdctime, &b_hcal_tdctime);
   fChain->SetBranchAddress("hcal_tdctime_c", &hcal_tdctime_c, &b_hcal_tdctime_c);
   fChain->SetBranchAddress("hcal_tdctimeblk", &hcal_tdctimeblk, &b_hcal_tdctimeblk);
   fChain->SetBranchAddress("hcal_tdctimeblk_c", &hcal_tdctimeblk_c, &b_hcal_tdctimeblk_c);
   fChain->SetBranchAddress("hcal_x", &hcal_x, &b_hcal_x);
   fChain->SetBranchAddress("hcal_xdiff", &hcal_xdiff, &b_hcal_xdiff);
   fChain->SetBranchAddress("hcal_y", &hcal_y, &b_hcal_y);
   fChain->SetBranchAddress("hcal_ydiff", &hcal_ydiff, &b_hcal_ydiff);
   fChain->SetBranchAddress("helicity", &helicity, &b_helicity);
   fChain->SetBranchAddress("helicity_status", &helicity_status, &b_helicity_status);
   fChain->SetBranchAddress("hodo_bar", &hodo_bar, &b_hodo_bar);
   fChain->SetBranchAddress("hodo_nbar", &hodo_nbar, &b_hodo_nbar);
   fChain->SetBranchAddress("hodo_nclus", &hodo_nclus, &b_hodo_nclus);
   fChain->SetBranchAddress("hodo_refleft", &hodo_refleft, &b_hodo_refleft);
   fChain->SetBranchAddress("hodo_refright", &hodo_refright, &b_hodo_refright);
   fChain->SetBranchAddress("hodo_tdiff", &hodo_tdiff, &b_hodo_tdiff);
   fChain->SetBranchAddress("hodo_tdiff_c", &hodo_tdiff_c, &b_hodo_tdiff_c);
   fChain->SetBranchAddress("hodo_timehitpos", &hodo_timehitpos, &b_hodo_timehitpos);
   fChain->SetBranchAddress("hodo_tleft", &hodo_tleft, &b_hodo_tleft);
   fChain->SetBranchAddress("hodo_tleft_c", &hodo_tleft_c, &b_hodo_tleft_c);
   fChain->SetBranchAddress("hodo_tmean", &hodo_tmean, &b_hodo_tmean);
   fChain->SetBranchAddress("hodo_tmean_c", &hodo_tmean_c, &b_hodo_tmean_c);
   fChain->SetBranchAddress("hodo_totleft", &hodo_totleft, &b_hodo_totleft);
   fChain->SetBranchAddress("hodo_totmean", &hodo_totmean, &b_hodo_totmean);
   fChain->SetBranchAddress("hodo_totright", &hodo_totright, &b_hodo_totright);
   fChain->SetBranchAddress("hodo_tright", &hodo_tright, &b_hodo_tright);
   fChain->SetBranchAddress("hodo_tright_c", &hodo_tright_c, &b_hodo_tright_c);
   fChain->SetBranchAddress("nu", &nu, &b_nu);
   fChain->SetBranchAddress("p_elas", &p_elas, &b_p_elas);
   fChain->SetBranchAddress("pdiff", &pdiff, &b_pdiff);
   fChain->SetBranchAddress("pexp_th", &pexp_th, &b_pexp_th);
   fChain->SetBranchAddress("ph_neutron", &ph_neutron, &b_ph_neutron);
   fChain->SetBranchAddress("ph_proton", &ph_proton, &b_ph_proton);
   fChain->SetBranchAddress("phi_e", &phi_e, &b_phi_e);
   fChain->SetBranchAddress("pred_ang_horiz", &pred_ang_horiz, &b_pred_ang_horiz);
   fChain->SetBranchAddress("pred_ang_vert", &pred_ang_vert, &b_pred_ang_vert);
   fChain->SetBranchAddress("pred_beta", &pred_beta, &b_pred_beta);
   fChain->SetBranchAddress("pred_mom", &pred_mom, &b_pred_mom);
   fChain->SetBranchAddress("pred_x", &pred_x, &b_pred_x);
   fChain->SetBranchAddress("pred_y", &pred_y, &b_pred_y);
   fChain->SetBranchAddress("proton_deflection", &proton_deflection, &b_proton_deflection);
   fChain->SetBranchAddress("ps_atime", &ps_atime, &b_ps_atime);
   fChain->SetBranchAddress("ps_e", &ps_e, &b_ps_e);
   fChain->SetBranchAddress("ps_ep", &ps_ep, &b_ps_ep);
   fChain->SetBranchAddress("ps_x", &ps_x, &b_ps_x);
   fChain->SetBranchAddress("ps_y", &ps_y, &b_ps_y);
   fChain->SetBranchAddress("run_num", &run_num, &b_run_num);
   fChain->SetBranchAddress("sbs_pathl", &sbs_pathl, &b_sbs_pathl);
   fChain->SetBranchAddress("sbs_rftime", &sbs_rftime, &b_sbs_rftime);
   fChain->SetBranchAddress("sbs_trigtime", &sbs_trigtime, &b_sbs_trigtime);
   fChain->SetBranchAddress("sh_atime", &sh_atime, &b_sh_atime);
   fChain->SetBranchAddress("sh_e", &sh_e, &b_sh_e);
   fChain->SetBranchAddress("sh_ep", &sh_ep, &b_sh_ep);
   fChain->SetBranchAddress("sh_x", &sh_x, &b_sh_x);
   fChain->SetBranchAddress("sh_y", &sh_y, &b_sh_y);
   fChain->SetBranchAddress("tau", &tau, &b_tau);
   fChain->SetBranchAddress("tdctrig", &tdctrig, &b_tdctrig);
   fChain->SetBranchAddress("tdctrig_id", &tdctrig_id, &b_tdctrig_id);
   fChain->SetBranchAddress("th_lab", &th_lab, &b_th_lab);
   fChain->SetBranchAddress("th_lab_deg", &th_lab_deg, &b_th_lab_deg);
   fChain->SetBranchAddress("th_neutron", &th_neutron, &b_th_neutron);
   fChain->SetBranchAddress("th_pq_neutron", &th_pq_neutron, &b_th_pq_neutron);
   fChain->SetBranchAddress("th_pq_proton", &th_pq_proton, &b_th_pq_proton);
   fChain->SetBranchAddress("th_proton", &th_proton, &b_th_proton);
   fChain->SetBranchAddress("theta_e", &theta_e, &b_theta_e);
   fChain->SetBranchAddress("thetabend", &thetabend, &b_thetabend);
   fChain->SetBranchAddress("trigbits", &trigbits, &b_trigbits);
   fChain->SetBranchAddress("trx_sh", &trx_sh, &b_trx_sh);
   fChain->SetBranchAddress("try_sh", &try_sh, &b_try_sh);
   fChain->SetBranchAddress("x_bj", &x_bj, &b_x_bj);
   Notify();
}

Bool_t skim_tree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void skim_tree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t skim_tree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef skim_tree_cxx
