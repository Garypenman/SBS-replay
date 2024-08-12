//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Aug 12 12:04:14 2024 by ROOT version 6.30/04
// from TTree T/Output of GEN tof Skim
// found on file: /home/gpenman/vol/data/out_tof_gen2_H2.root
//////////////////////////////////////////////////////////

#ifndef skim_tree_h
#define skim_tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class skim_tree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        Q2;
   Double_t        W2;
   Double_t        eps;
   Double_t        nu;
   Double_t        x_bj;
   Double_t        bb_x;
   Double_t        bb_y;
   Double_t        bb_px;
   Double_t        bb_py;
   Double_t        bb_pz;
   Double_t        bb_p;
   Double_t        bb_vx;
   Double_t        bb_vy;
   Double_t        bb_vz;
   Double_t        bb_th_fp;
   Double_t        bb_th_tg;
   Double_t        bb_ph_fp;
   Double_t        bb_ph_tg;
   Double_t        bb_y_tg;
   Double_t        bb_pathl;
   Double_t        bb_tr_nhits;
   Double_t        bb_ntr;
   Double_t        trigtime;
   Double_t        sh_e;
   Double_t        sh_x;
   Double_t        sh_y;
   Double_t        ps_e;
   Double_t        ps_x;
   Double_t        ps_y;
   Double_t        hcal_e;
   Double_t        hcal_eblk;
   Double_t        hcal_x;
   Double_t        hcal_y;
   Double_t        hcal_row;
   Double_t        hcal_col;
   Double_t        hcal_id;
   Double_t        pred_ang_horiz;
   Double_t        pred_ang_vert;
   Double_t        pred_y;
   Double_t        pred_x;
   Double_t        pred_mom;
   Double_t        dx;
   Double_t        dy;
   Double_t        dev;
   Double_t        ps_atime;
   Double_t        sh_atime;
   Double_t        hcal_atime;
   Double_t        coin_atime;
   Double_t        hcal_time;
   Double_t        coin_time;
   vector<double>  *hodo_tmean;
   vector<double>  *hodo_tleft;
   vector<double>  *hodo_tright;
   vector<double>  *hodo_totleft;
   vector<double>  *hodo_totright;
   vector<double>  *hodo_timehitpos;
   vector<double>  *hodo_bar;
   Double_t        gr_adc;
   Double_t        gr_size;
   Double_t        gr_tmean;
   Double_t        gr_totmean;
   Double_t        gr_trindex;
   Double_t        gr_xmean;
   Double_t        gr_ymean;
   Double_t        helicity;
   Int_t           helicity_status;
   Double_t        He3Pol;
   Double_t        BeamPol;
   ULong64_t       evt_time;
   UInt_t          evt_num;
   UInt_t          evt_type;
   UInt_t          evt_len;
   Int_t           evt_hel;
   UInt_t          run_num;
   UInt_t          trigbits;
   vector<int>     *tdctrig_id;
   vector<double>  *tdctrig;
   Long64_t        UID;

   // List of branches
   TBranch        *b_Q2;   //!
   TBranch        *b_W2;   //!
   TBranch        *b_eps;   //!
   TBranch        *b_nu;   //!
   TBranch        *b_x_bj;   //!
   TBranch        *b_bb_x;   //!
   TBranch        *b_bb_y;   //!
   TBranch        *b_bb_px;   //!
   TBranch        *b_bb_py;   //!
   TBranch        *b_bb_pz;   //!
   TBranch        *b_bb_p;   //!
   TBranch        *b_bb_vx;   //!
   TBranch        *b_bb_vy;   //!
   TBranch        *b_bb_vz;   //!
   TBranch        *b_bb_th_fp;   //!
   TBranch        *b_bb_th_tg;   //!
   TBranch        *b_bb_ph_fp;   //!
   TBranch        *b_bb_ph_tg;   //!
   TBranch        *b_bb_y_tg;   //!
   TBranch        *b_bb_pathl;   //!
   TBranch        *b_bb_tr_nhits;   //!
   TBranch        *b_bb_ntr;   //!
   TBranch        *b_trigtime;   //!
   TBranch        *b_sh_e;   //!
   TBranch        *b_sh_x;   //!
   TBranch        *b_sh_y;   //!
   TBranch        *b_ps_e;   //!
   TBranch        *b_ps_x;   //!
   TBranch        *b_ps_y;   //!
   TBranch        *b_hcal_e;   //!
   TBranch        *b_hcal_eblk;   //!
   TBranch        *b_hcal_x;   //!
   TBranch        *b_hcal_y;   //!
   TBranch        *b_hcal_row;   //!
   TBranch        *b_hcal_col;   //!
   TBranch        *b_hcal_id;   //!
   TBranch        *b_pred_ang_horiz;   //!
   TBranch        *b_pred_ang_vert;   //!
   TBranch        *b_pred_y;   //!
   TBranch        *b_pred_x;   //!
   TBranch        *b_pred_mom;   //!
   TBranch        *b_dx;   //!
   TBranch        *b_dy;   //!
   TBranch        *b_dev;   //!
   TBranch        *b_ps_atime;   //!
   TBranch        *b_sh_atime;   //!
   TBranch        *b_hcal_atime;   //!
   TBranch        *b_coin_atime;   //!
   TBranch        *b_hcal_time;   //!
   TBranch        *b_coin_time;   //!
   TBranch        *b_hodo_tmean;   //!
   TBranch        *b_hodo_tleft;   //!
   TBranch        *b_hodo_tright;   //!
   TBranch        *b_hodo_totleft;   //!
   TBranch        *b_hodo_totright;   //!
   TBranch        *b_hodo_timehitpos;   //!
   TBranch        *b_hodo_bar;   //!
   TBranch        *b_gr_adc;   //!
   TBranch        *b_gr_size;   //!
   TBranch        *b_gr_tmean;   //!
   TBranch        *b_gr_totmean;   //!
   TBranch        *b_gr_trindex;   //!
   TBranch        *b_gr_xmean;   //!
   TBranch        *b_gr_ymean;   //!
   TBranch        *b_helicity;   //!
   TBranch        *b_helicity_status;   //!
   TBranch        *b_He3Pol;   //!
   TBranch        *b_BeamPol;   //!
   TBranch        *b_evt_time;   //!
   TBranch        *b_evt_num;   //!
   TBranch        *b_evt_type;   //!
   TBranch        *b_evt_len;   //!
   TBranch        *b_evt_hel;   //!
   TBranch        *b_run_num;   //!
   TBranch        *b_trigbits;   //!
   TBranch        *b_tdctrig_id;   //!
   TBranch        *b_tdctrig;   //!
   TBranch        *b_UID;   //!

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
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/home/gpenman/vol/data/out_tof_gen2_H2.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/home/gpenman/vol/data/out_tof_gen2_H2.root");
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
   hodo_tmean = 0;
   hodo_tleft = 0;
   hodo_tright = 0;
   hodo_totleft = 0;
   hodo_totright = 0;
   hodo_timehitpos = 0;
   hodo_bar = 0;
   tdctrig_id = 0;
   tdctrig = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Q2", &Q2, &b_Q2);
   fChain->SetBranchAddress("W2", &W2, &b_W2);
   fChain->SetBranchAddress("eps", &eps, &b_eps);
   fChain->SetBranchAddress("nu", &nu, &b_nu);
   fChain->SetBranchAddress("x_bj", &x_bj, &b_x_bj);
   fChain->SetBranchAddress("bb_x", &bb_x, &b_bb_x);
   fChain->SetBranchAddress("bb_y", &bb_y, &b_bb_y);
   fChain->SetBranchAddress("bb_px", &bb_px, &b_bb_px);
   fChain->SetBranchAddress("bb_py", &bb_py, &b_bb_py);
   fChain->SetBranchAddress("bb_pz", &bb_pz, &b_bb_pz);
   fChain->SetBranchAddress("bb_p", &bb_p, &b_bb_p);
   fChain->SetBranchAddress("bb_vx", &bb_vx, &b_bb_vx);
   fChain->SetBranchAddress("bb_vy", &bb_vy, &b_bb_vy);
   fChain->SetBranchAddress("bb_vz", &bb_vz, &b_bb_vz);
   fChain->SetBranchAddress("bb_th_fp", &bb_th_fp, &b_bb_th_fp);
   fChain->SetBranchAddress("bb_th_tg", &bb_th_tg, &b_bb_th_tg);
   fChain->SetBranchAddress("bb_ph_fp", &bb_ph_fp, &b_bb_ph_fp);
   fChain->SetBranchAddress("bb_ph_tg", &bb_ph_tg, &b_bb_ph_tg);
   fChain->SetBranchAddress("bb_y_tg", &bb_y_tg, &b_bb_y_tg);
   fChain->SetBranchAddress("bb_pathl", &bb_pathl, &b_bb_pathl);
   fChain->SetBranchAddress("bb_tr_nhits", &bb_tr_nhits, &b_bb_tr_nhits);
   fChain->SetBranchAddress("bb_ntr", &bb_ntr, &b_bb_ntr);
   fChain->SetBranchAddress("trigtime", &trigtime, &b_trigtime);
   fChain->SetBranchAddress("sh_e", &sh_e, &b_sh_e);
   fChain->SetBranchAddress("sh_x", &sh_x, &b_sh_x);
   fChain->SetBranchAddress("sh_y", &sh_y, &b_sh_y);
   fChain->SetBranchAddress("ps_e", &ps_e, &b_ps_e);
   fChain->SetBranchAddress("ps_x", &ps_x, &b_ps_x);
   fChain->SetBranchAddress("ps_y", &ps_y, &b_ps_y);
   fChain->SetBranchAddress("hcal_e", &hcal_e, &b_hcal_e);
   fChain->SetBranchAddress("hcal_eblk", &hcal_eblk, &b_hcal_eblk);
   fChain->SetBranchAddress("hcal_x", &hcal_x, &b_hcal_x);
   fChain->SetBranchAddress("hcal_y", &hcal_y, &b_hcal_y);
   fChain->SetBranchAddress("hcal_row", &hcal_row, &b_hcal_row);
   fChain->SetBranchAddress("hcal_col", &hcal_col, &b_hcal_col);
   fChain->SetBranchAddress("hcal_id", &hcal_id, &b_hcal_id);
   fChain->SetBranchAddress("pred_ang_horiz", &pred_ang_horiz, &b_pred_ang_horiz);
   fChain->SetBranchAddress("pred_ang_vert", &pred_ang_vert, &b_pred_ang_vert);
   fChain->SetBranchAddress("pred_y", &pred_y, &b_pred_y);
   fChain->SetBranchAddress("pred_x", &pred_x, &b_pred_x);
   fChain->SetBranchAddress("pred_mom", &pred_mom, &b_pred_mom);
   fChain->SetBranchAddress("dx", &dx, &b_dx);
   fChain->SetBranchAddress("dy", &dy, &b_dy);
   fChain->SetBranchAddress("dev", &dev, &b_dev);
   fChain->SetBranchAddress("ps_atime", &ps_atime, &b_ps_atime);
   fChain->SetBranchAddress("sh_atime", &sh_atime, &b_sh_atime);
   fChain->SetBranchAddress("hcal_atime", &hcal_atime, &b_hcal_atime);
   fChain->SetBranchAddress("coin_atime", &coin_atime, &b_coin_atime);
   fChain->SetBranchAddress("hcal_time", &hcal_time, &b_hcal_time);
   fChain->SetBranchAddress("coin_time", &coin_time, &b_coin_time);
   fChain->SetBranchAddress("hodo_tmean", &hodo_tmean, &b_hodo_tmean);
   fChain->SetBranchAddress("hodo_tleft", &hodo_tleft, &b_hodo_tleft);
   fChain->SetBranchAddress("hodo_tright", &hodo_tright, &b_hodo_tright);
   fChain->SetBranchAddress("hodo_totleft", &hodo_totleft, &b_hodo_totleft);
   fChain->SetBranchAddress("hodo_totright", &hodo_totright, &b_hodo_totright);
   fChain->SetBranchAddress("hodo_timehitpos", &hodo_timehitpos, &b_hodo_timehitpos);
   fChain->SetBranchAddress("hodo_bar", &hodo_bar, &b_hodo_bar);
   fChain->SetBranchAddress("gr_adc", &gr_adc, &b_gr_adc);
   fChain->SetBranchAddress("gr_size", &gr_size, &b_gr_size);
   fChain->SetBranchAddress("gr_tmean", &gr_tmean, &b_gr_tmean);
   fChain->SetBranchAddress("gr_totmean", &gr_totmean, &b_gr_totmean);
   fChain->SetBranchAddress("gr_trindex", &gr_trindex, &b_gr_trindex);
   fChain->SetBranchAddress("gr_xmean", &gr_xmean, &b_gr_xmean);
   fChain->SetBranchAddress("gr_ymean", &gr_ymean, &b_gr_ymean);
   fChain->SetBranchAddress("helicity", &helicity, &b_helicity);
   fChain->SetBranchAddress("helicity_status", &helicity_status, &b_helicity_status);
   fChain->SetBranchAddress("He3Pol", &He3Pol, &b_He3Pol);
   fChain->SetBranchAddress("BeamPol", &BeamPol, &b_BeamPol);
   fChain->SetBranchAddress("evt_time", &evt_time, &b_evt_time);
   fChain->SetBranchAddress("evt_num", &evt_num, &b_evt_num);
   fChain->SetBranchAddress("evt_type", &evt_type, &b_evt_type);
   fChain->SetBranchAddress("evt_len", &evt_len, &b_evt_len);
   fChain->SetBranchAddress("evt_hel", &evt_hel, &b_evt_hel);
   fChain->SetBranchAddress("run_num", &run_num, &b_run_num);
   fChain->SetBranchAddress("trigbits", &trigbits, &b_trigbits);
   fChain->SetBranchAddress("tdctrig_id", &tdctrig_id, &b_tdctrig_id);
   fChain->SetBranchAddress("tdctrig", &tdctrig, &b_tdctrig);
   fChain->SetBranchAddress("UID", &UID, &b_UID);
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
