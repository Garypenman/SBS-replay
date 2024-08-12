//#include "gmn_tree.C"
#include "skim_tree.C"
#include "TChain.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
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

const double Mp=0.938272;

void TOFcal(const char *inputfilename="testconfig", const char *outputfilename="TOFcal_out.root"){

  TString currentline;

  /*
  ifstream configfile(inputfilename);
  
  if( !configfile ) return;
  
  TChain *C = new TChain("T");
  
  while( currentline.ReadLine(configfile) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      C->Add(currentline.Data());
    }
  }

  TCut globalcut = "";
  
  while( currentline.ReadLine(configfile) && !currentline.BeginsWith("endcut") ){
    if( !currentline.BeginsWith("#") ){
      globalcut += currentline.Data();
    }
  }
  */
  
  TChain *C = new TChain("T");
  C->Add("~/vol/data/out_tof_gen2_H2.root");
  
  //Default reference PMTs for aligning all other PMTs:
  int hodorefID = 43; // reference BAR (not PMT); we'll set the LEFT PMT to the reference by convention:
  int hcalrefID = 199; //row 16 column 6
  
  double zhodo = 1.854454; //meters
  double Lbar_hodo = 0.6; //meters
  double etof0 = (1.96+3.0)/0.299792458; // "central" TOF value:

  double thetaHCAL = 34.7;
  double dHCAL = 17.0;
  double dSBS = 2.25;

  double W2min=0.6, W2max=1.2;
  
  double sbsmaxfield = 1.26;
  double sbsfield = 1;

  const double Dgap = 48.0*2.54/100.0;
  
  double Ebeam = 4.291;
  
  /*
  while( currentline.ReadLine(configfile) && !currentline.BeginsWith("endconfig") ){
    if( !currentline.BeginsWith("#") ){
      TObjArray *currentline_tokens = currentline.Tokenize(" ");
      if( currentline_tokens->GetEntries() >= 2 ){
	TString skey = ( (TObjString*) (*currentline_tokens)[0] )->GetString();

	TString sval = ( (TObjString*) (*currentline_tokens)[1] )->GetString();
	if( skey.BeginsWith("hodorefID") ){
	  hodorefID = sval.Atoi();
	}

	if( skey.BeginsWith("hcalrefID") ){
	  hcalrefID = sval.Atoi();
	}

	if( skey.BeginsWith("zhodo") ){
	  zhodo = sval.Atof();
	}

	if( skey.BeginsWith("Lbar_hodo") ){
	  Lbar_hodo = sval.Atof();
	}

	if( skey.BeginsWith("etof0") ){
	  etof0 = sval.Atof();
	}
      }
    }
 }
  */

  thetaHCAL *= TMath::Pi()/180.0;
  
  TVector3 zaxis_HCAL(-sin(thetaHCAL),0,cos(thetaHCAL));
  TVector3 xaxis_HCAL(0,-1,0);
  TVector3 yaxis_HCAL = (zaxis_HCAL.Cross(xaxis_HCAL)).Unit();
  
  TVector3 HCALorigin = dHCAL*zaxis_HCAL;


  //int nparams_hodo = 3*180; //offset, walk correction, propagation speed for left and right PMTs for 90 bars minus offset fixed at zero for reference PMT
  //actually, vscint should be determined per paddle;
  int nparams_hodo = 180+2*90; //We actually only really want one zero offset (paddle-specific) that applies to the mean time, so we need 90 zero offsets, 90 vscint values, and 180 walk correction slopes.
  int nparams_HCAL = 288*2; //offset and walk correction for 288 modules.

  //read in previously generated config for iterative calibration
  vector<double> hodoparams, hcalparams;
  vector<double> HODOt0, HODOwL, HODOwR, vscint, HCALt0, HCALw;
  TString hodoconfigfile = "hodoparams.txt";
  TString hcalconfigfile = "hcalparams.txt";
  
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
  ifstream hcalconfigin(hodoconfigfile);
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
    //HODOt0R.push_back( hodoparams[i+90] );
    HODOwL.push_back( hodoparams[i+180] );
    HODOwR.push_back( hodoparams[i+270] );
    vscint.push_back( 1.0/hodoparams[i+90] );
  }
  
  for (int i=0; i<288; i++ ){
    HCALt0.push_back( hcalparams[i] );
    HCALw.push_back( hcalparams[i+288] );
  }
  
  //Ordering of params:
  // t0ileft, t0iright, walkileft, walkiright, 1/vileft, 1/viright
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
  TH2D *htmean_hodoID_new = new TH2D("htmean_hodoID_new","NEW ;bar ID; bar t_{mean} [ns]", 90,-0.5,89.5, nbins,-30,30);

  TH2D *htdiff_yhodo_old = new TH2D("htdiff_yhodo_old", "OLD; y_{hodo} (m); bar t_{L}-t_{R} [ns]", nbins,-0.3,0.3,nbins,-15,15);
  TH2D *htdiff_yhodo_new = new TH2D("htdiff_yhodo_new", "NEW; y_{hodo} (m); bar t_{L}-t_{R} [ns]", nbins,-0.3,0.3,nbins,-15,15);

  TH2D *htmean_trigRF_old = new TH2D("htmean_trigRF_old", " OLD; bar ID; bar t_{mean} + bb_{trig} - fmod{rf,2.004) ", 90,-0.5,89.5,nbins,-30,30);
  TH2D *htmean_trigRF_new = new TH2D("htmean_trigRF_new", " NEW; bar ID; bar t_{mean} + bb_{trig} - fmod{rf,2.004) ", 90,-0.5,89.5,nbins,-30,30);
  
  TH2D *htmean_hcalID_old = new TH2D("htmean_hcalID_old","OLD ;hcal ID; clus t_{mean} [ns]", 288,-0.5,287.5,nbins,-30,30);
  TH2D *htmean_hcalID_new = new TH2D("htmean_hcalID_new","NEW ;hcal ID; clus t_{mean} [ns]", 288,-0.5,287.5,nbins,-30,30);
  
  TH2D *htmean_hcale_old = new TH2D("htmean_hcale_old","OLD ;hcal edep [GeV]; clus t_{mean} [ns]", nbins,0,0.5,nbins,-30,30);
  TH2D *htmean_hcale_new = new TH2D("htmean_hcale_new","NEW ;hcal edep [GeV]; clus t_{mean} [ns]", nbins,0,0.5,nbins,-30,30);
  
  //Really want a large quantity of zero-field LH2 data for this purpose:
  TH2D *htdiffHCAL_vs_HCALID_old = new TH2D("htdiffHCAL_vs_HCALID_old", "OLD hodo, OLD HCAL; HCAL ID; t_{HCAL}-t_{HODO,corr} [ns]",224,-0.5,223.5,250,-50,50); 
  TH2D *htdiffHCAL_vs_HCALID_new = new TH2D("htdiffHCAL_vs_HCALID_new", "NEW hodo, NEW HCAL; HCAL ID; t_{HCAL}-t_{HODO,corr} [ns]",224,-0.5,223.5,250,-50,50); 

  TH2D *hdxdy = new TH2D("hdxdy",";#Deltay (m); #Deltax (m)",250,-1.25,1.25,250,-1.25,1.25);

  //event loop stuff
  skim_tree *T = new skim_tree(C);
  int treenum=-1,oldtreenum=-1;
  long nev = C->GetEntries();
  for( long ev=0; ev<nev; ev++ ){
    
    T->GetEntry(ev);
    if( ev % 10000 == 0 ) cout << "nevent = " << ev << " / " << nev << endl;
    
    //define some simple cuts here
    //skim file takes care of exclusive
    //condition of atleast 1 track and 1 hcal cluster
    //if (T->bb_ntr<1) continue;
    //if (T->sbs_hcal_nclus<1) continue;
    if (T->hodo_bar->size()==0) continue;
    
    //3 hits is the minimum anyway.
    //this can be made tighter, 4 or 5 hits, if we want
    if (T->bb_tr_nhits < 3) continue;

    //calorimer energy cuts
    //is hcal energy reconstruction reliable?
    //if (T->hcal_e < 0.05) continue;
    if (T->ps_e<0.2) continue;

    //in principle we can clean up electron selection
    //jack has calibrated grinch timing for all kinematics
    //if (T->gr_trindex ==-1) continue;
    //if (T->gr_size < 3) continue;

    //quasi invariant mass squared cut
    //for GEN2 can leave fairly wide. even 0<W2<2 would probably be fine
    if (T->W2 < 0 || T->W2 > 2.0) continue;
    

    //grab needed track parameters:
    double vz = T->bb_vz;
    double pathl = T->bb_pathl;
    double etof = pathl/0.299792458; //electron TOF from vertex to hodo.
    
    //bbcal trigger time
    double trigtime = T->trigtime; 
    
    //rf time
    double rftime=-999.;
    for(uint i=0; i<T->tdctrig_id->size(); i++){
      if(T->tdctrig_id->at(i)==4) rftime = T->tdctrig->at(i);
    }
    if (rftime==-999) continue;
    
    //relevant hodo cluster parameters
    double tleft = T->hodo_tleft->at(0);
    //tleft += trigtime - fmod(rftime,2.004);
    double tright = T->hodo_tright->at(0);
    //tright += trigtime - fmod(rftime,2.004);
    double totleft = T->hodo_totleft->at(0);
    double totright = T->hodo_totright->at(0);

    int ID = int(T->hodo_bar->at(0));
      
    double yhodo = T->bb_y+zhodo*T->bb_ph_fp;
      
    //relevant hcal parameters
    double xHCAL = T->hcal_x;
    double yHCAL = T->hcal_y;
    int idblkHCAL = (int) T->hcal_id - 1;
    double tHCAL = T->hcal_time;
    double eblkHCAL = T->hcal_e;
      
      
    // how to define our chi2 statistic? We want to correct all hodo PMT times to tvertex = 0:
    // tPMT = etof + t0 - walk * TOT + d/vscint
    // chi2 = sum_{i=1}^N \sum_{j=1}^Nhit (tPMT-etof-t0+walk*TOT - d/vscint)^2
    // dchi2/dt0_k = 2*\sum_i,j (tPMT-etof-t0+walk*TOT-d/vscint)*-1
    // dchi2/dwalk_k = 2*sum_i,j (tPMT-etof-t0+walk*TOT-d/vscint)*TOT
    // dchi2/d(1/v)_k = 2*sum_i,j (tPMT-etof-t0+walk*TOT-d/vscint)*-d

    double dLEFT = std::min(Lbar_hodo,std::max(0.0,Lbar_hodo/2.0 - yhodo));
    double dRIGHT = std::min(Lbar_hodo,std::max(0.0,Lbar_hodo/2.0 + yhodo));
      
    //hcal stuff
    TVector3 ep(T->bb_px, T->bb_py, T->bb_pz);
    TLorentzVector k(0,0,Ebeam,Ebeam);
    TLorentzVector P(0,0,0,Mp);

    TLorentzVector kprime(ep,ep.Mag());
    TLorentzVector q = k-kprime;
    TLorentzVector Pprime = P + q;

    double Q2 = q.M2();
      
    double etheta = ep.Theta();
    double ephi = ep.Phi();
      
    double Eprime_etheta = Ebeam/(1.0+Ebeam/Mp*(1.0-cos(etheta)));
    double Tp_etheta = Ebeam-Eprime_etheta;
    double pp_etheta = sqrt(pow(Tp_etheta,2)+2.0*Mp*Tp_etheta);
    double ptheta_etheta = acos( (Ebeam-Eprime_etheta*cos(etheta))/(pp_etheta));
    double pphi_ephi = ephi + TMath::Pi();

    double ptheta_4vect = q.Vect().Theta();
    double pphi_4vect = q.Vect().Phi();
      
    TVector3 pnhat_expect( sin(ptheta_etheta)*cos(pphi_ephi),
			   sin(ptheta_etheta)*sin(pphi_ephi),
			   cos(ptheta_etheta) );

    TVector3 vertex(0,0,vz);

    double sintersect = (HCALorigin-vertex).Dot(zaxis_HCAL)/(pnhat_expect.Dot(zaxis_HCAL));

    TVector3 HCAL_intersect = vertex + sintersect * pnhat_expect;

    double xHCAL_expect = (HCAL_intersect - HCALorigin).Dot( xaxis_HCAL );
    double yHCAL_expect = (HCAL_intersect - HCALorigin).Dot( yaxis_HCAL );

    double Lpath_HCAL = (HCAL_intersect - vertex).Mag();

    double beta_proton = pp_etheta/sqrt(pow(pp_etheta,2)+pow(Mp,2));
      
    double TOF_HCAL = Lpath_HCAL/(beta_proton*0.299792458);
    double HCALtcent = dHCAL/(beta_proton*0.299792458);
      
    double W2 = T->W2;

    // thetabend = 0.3/p * BdL ;
    double BdL = Dgap * sbsmaxfield * sbsfield;
    double thetabend = 0.3/pp_etheta * BdL;
     
    double protondeflection = tan(thetabend)*(dHCAL-(dSBS+Dgap/2.0)); 
      
    double deltax = xHCAL-(xHCAL_expect-protondeflection);
    double deltay = yHCAL-(yHCAL_expect);
    //zero data in sbs which is wrongly pushed back as a "zero" in the datastream
    if (idblkHCAL == -1) continue;
      
    //if(iter==0){
    // cout << "Best hodo hit (tL,tR, totL, totR, yhodo, dL, dR, pathl, etof)=("
    // 	   << tleft << ", " << tright << ", " << totleft << ", " << totright << ", " << yhodo << ", "
    // 	   << dLEFT << ", " << dRIGHT << ", " << pathl << ", " << etof << ")" << endl;

    // we want to minimize
    // chi2 = sum_{i=1}^Nevent (tleft_corr^2 + tright_corr^2)
    // dchi2/dt0_j = sum_i (2*tleft_corr_i*dtleft_corr_i/dt0_j + 2*tright_corr_i*dtright_corr_i/dt0_j) = 0
    //we COULD just use TMinuit but that's lazy.
    //let parameters 0-89 be the t0L, 90-179 be t0R, 180-269 be wL, 270-359 be wR, and 360-449 be 1/vscint:
    int ipar_t0 = ID;
    int ipar_vinv = ID+90;
    //int ipar_t0R = ID+90;
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
	
    //hcal pmts:
    ipar_t0 = idblkHCAL;
    int ipar_w = idblkHCAL + 288; 
    Mhcal(ipar_t0,ipar_t0) += 1;
    Mhcal(ipar_t0,ipar_w) += -eblkHCAL;
    Mhcal(ipar_w,ipar_t0) +=  -eblkHCAL;
    Mhcal(ipar_w,ipar_w) += pow(eblkHCAL,2);
	
    bhcal(ipar_t0) += tHCAL - (TOF_HCAL - HCALtcent);
    bhcal(ipar_w) += (tHCAL - (TOF_HCAL - HCALtcent))*(-eblkHCAL);
    //if this is the second time running, calculate corrected variables, otherwise
    //these will just come out the same as the "old" variables since all new factors
    //will be zero upon first iteration.
    double tleft_CORR = tleft - (etof-etof0) - HODOt0[ID] + HODOwL[ID]*totleft - dLEFT/vscint[ID];
    double tright_CORR = tright - (etof-etof0) - HODOt0[ID] + HODOwR[ID]*totright - dRIGHT/vscint[ID];

    double tmean_CORR = 0.5 * (tleft_CORR + tright_CORR);
	
    double tmean_old = (tleft + tright)/2;
    double tdiff_old = tleft - tright;
	
    double tHCAL_CORR = tHCAL - HCALt0[idblkHCAL] + HCALw[idblkHCAL]*eblkHCAL;
	
    //calculate a corrected time difference without the propagation correction:
    double tdiff_CORR = (tleft - HODOt0[ID]+ HODOwL[ID]*totleft) - (tright-HODOt0[ID]+HODOwR[ID]*totright);
      
    //Fill histograms
    htmean_hodoID_old->Fill( ID, tmean_old );
    htmean_hodoID_new->Fill( ID, tmean_CORR );
    htdiff_yhodo_old->Fill( yhodo, tdiff_old );
    htdiff_yhodo_new->Fill( yhodo, tdiff_CORR );
    htmean_trigRF_old->Fill( ID, tmean_old + trigtime - fmod(rftime,2.004) );
    htmean_trigRF_new->Fill( ID, tmean_CORR+ trigtime - fmod(rftime,2.004) );
    htmean_hcalID_old->Fill( idblkHCAL, tHCAL );
    htmean_hcalID_new->Fill( idblkHCAL, tHCAL_CORR );
    htmean_hcale_old->Fill( eblkHCAL, tHCAL );
    htmean_hcale_new->Fill( eblkHCAL, tHCAL_CORR );
	
    hdxdy->Fill( deltay,deltax );
      
    //if( W2min < W2 && W2 < W2max && sqrt(pow(deltax,2)+pow(deltay,2))<=0.24 && eblkHCAL>0.02 ){
    if( eblkHCAL>0.02 && W2>W2min && W2<W2max && sqrt(pow(deltax,2)+pow(deltay,2))<=0.24 ){
      // cout << "(Eprime,etheta,ephi)=("
      //      << ep.Mag() << ", " << etheta*57.3 << ", "
      //      << ephi*57.3 << ")" << endl;
      // cout << "(ptheta_etheta,pp_etheta)=("
      //      << ptheta_etheta*57.3 << ", " << pp_etheta << ")" << endl;
	
      // cout << "(ephi,pphi)=(" << ephi*57.3 << ", " << pphi_ephi*57.3 << ")"
      //      << endl;
      
      htdiffHCAL_vs_HCALID_old->Fill( idblkHCAL, tHCAL-tmean_old);
      htdiffHCAL_vs_HCALID_new->Fill( idblkHCAL, tHCAL-tmean_CORR-(TOF_HCAL-HCALtcent));
    }
  }//end of event loop
    
  TDecompSVD Ahodo(Mhodo);
  Ahodo.Solve(bhodo);

  TDecompSVD Ahcal(Mhcal);
  Ahcal.Solve(bhcal);
  
  bhodo.Print();
  bhcal.Print();
	      
  //double HODOt0[90], HODOwL[90], HODOwR[90], vscint[90];
  //double HCALt0[288], HCALw[288];
  
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

