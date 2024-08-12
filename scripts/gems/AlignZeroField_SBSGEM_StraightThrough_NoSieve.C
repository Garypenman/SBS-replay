#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include <iostream>
#include <fstream>
#include "TMinuit.h"
#include "TFitResult.h"
#include "TF1.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TString.h"
#include "TCut.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTreeFormula.h"
#include "TVector3.h"
#include "TRotation.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompSVD.h"
#include <iostream>
#include <fstream>

const double PI = TMath::Pi();

// Align SBS GEM to Hall A coordinate system using single-foil, 
// straight-through data

int NMAX=100000;

int NTRACKS;

vector<double> XTRACK,YTRACK,XPTRACK,YPTRACK,xHCAL,yHCAL;

void CHI2_FCN( int &npar, double *gin, double &f, double *par, int flag ){
  double chi2 = 0.0;
  
  TVector3 GEMPOS(par[0],par[1],par[2]);
  
  double ax = par[3];
  double ay = par[4];
  double az = par[5];

  TVector3 HCALPOS(par[6],par[7],par[8]);
  
  //assuming these angles represent small misalignments, it basically doesn't matter in what order we apply the rotations.
  TRotation Rtot;
  Rtot.RotateZ(az);
  Rtot.RotateX(ax);
  Rtot.RotateY(ay);
  
  TRotation Rinv = Rtot.Inverse();
  
  TVector3 Global_yaxis(0,1,0);
  TVector3 Global_xaxis(1,0,0);
  TVector3 Global_zaxis(0,0,1);
  
  TVector3 GEM_zaxis = Rtot * Global_zaxis;
  //TVector3 GEM_xaxis = (Global_yaxis.Cross(GEM_zaxis)).Unit();
  //TVector3 GEM_yaxis = (GEM_zaxis.Cross(GEM_xaxis)).Unit();

  TVector3 GEM_xaxis = Rtot * Global_xaxis;
  TVector3 GEM_yaxis = Rtot * Global_yaxis;

  //Here our goal is to assume that all tracks come from the origin, 
  // (xtar,ytar,ztar=0)=(0,0,0)

  TVector3 BeamPos(0.0,0,0); //This should not be hard-coded, but for now based on survey results for GEN-RP, we take the beam to be 6 mm below the SBS magnet center and the SBS GEM center.
  
  chi2 = 0.0;
  for( int i=0; i<NTRACKS; i++ ){
    TVector3 TrackDirLocal(XPTRACK[i],YPTRACK[i],1.0);
    TrackDirLocal = TrackDirLocal.Unit();
    
    //    TVector3 TrackDirGlobal = TrackDirLocal.X() * GEM_xaxis +
    //  TrackDirLocal.Y() * GEM_yaxis + 
    //  TrackDirLocal.Z() * GEM_zaxis;

    TVector3 TrackDirGlobal = Rinv * TrackDirLocal;
    //For the direction we need to apply the INVERSE of the total rotation:
    // Why? because it seems to work, but it "feels" wrong
    //                   | RXX RYX RZX |   | X' |   | RXX X' + RYX Y' + RZX Z' |   
    // Rinv * dirlocal = | RXY RYY RZY | * | Y' | = | RXY X' + RYY Y' + RZY Z' | 
    //                   | RXZ RYZ RZZ |   | Z' |   | RXZ X' + RYZ Y' + RZZ Z' |
    //                         | RXX RXY RXZ |   | X |   | RXX X + RXY Y + RXZ Z |   
    // GEMPOS + R * poslocal = | RXY RYY RYZ | * | Y | = | RYX X + RYY Y + RYZ Z | = GEMPOS + X * xaxis + Y * yaxis + Z * zaxis
    //                         | RZX RZY RZZ |   | Z |   | RZX X + RZY Y + RZZ Z |   
    //TVector3 TrackDirGlobal = Rinv * TrackDirLocal; 
    
    //TVector3 TrackDirGlobal = GEM_xaxis * TrackDirLocal.X() + GEM_yaxis * TrackDirLocal.Y() + GEM_zaxis * TrackDirLocal.Z();

    //Let's ASSUME track originates at the origin:
    // (BeamPos + s * nhatg - GEMPOS ) dot GEMzaxis = 0
    double sint = (GEMPOS-BeamPos).Dot( GEM_zaxis ) / TrackDirGlobal.Dot( GEM_zaxis );

    TVector3 TrackIntersect_GEM = BeamPos + sint * TrackDirGlobal;

    double xfp_expect = (TrackIntersect_GEM - GEMPOS).Dot( GEM_xaxis );
    double yfp_expect = (TrackIntersect_GEM - GEMPOS).Dot( GEM_yaxis );
    
    TVector3 TrackPosLocal(XTRACK[i], YTRACK[i], 0.0 );
    TVector3 TrackPosGlobal = GEMPOS + Rtot * TrackPosLocal;
    //  TVector3 TrackPosGlobal = GEMPOS + TrackPosLocal.X() * GEM_xaxis + TrackPosLocal.Y() * GEM_yaxis + TrackPosLocal.Z() * GEM_zaxis;

    double sint_target = -TrackPosGlobal.Z() / TrackDirGlobal.Z();

    TVector3 TrackIntersect_target = TrackPosGlobal + sint_target * TrackDirGlobal;

    double ytarget = TrackIntersect_target.Y();
    double xtarget = TrackIntersect_target.X();
    
    //Here let's not use the measured track direction but define the global track direction according to to the global position,
    //reflecting the assumption that the track starts from the origin!
    //TVector3 TrackDirGlobal_pos = TrackPosGlobal.Unit(); //However, this will FORCE the track to point exactly at the origin. 
    
    //Now calculate track intersection with the xy plane: 

    //double sint_HCAL = (HCALPOS - TrackPosGlobal).Dot( Global_zaxis ) / TrackDirGlobal.Dot( Global_zaxis );

    double sint_HCAL = HCALPOS.Z() / TrackDirGlobal.Z(); 
    
    TVector3 TrackIntersectHCAL = sint_HCAL * TrackDirGlobal;

    double thHCAL = (xHCAL[i]+HCALPOS.X()-BeamPos.X())/HCALPOS.Z();
    double phHCAL = (yHCAL[i]+HCALPOS.Y()-BeamPos.Y())/HCALPOS.Z();
    double thGEM = TrackDirGlobal.X()/TrackDirGlobal.Z();
    double phGEM = TrackDirGlobal.Y()/TrackDirGlobal.Z();

    //double thGEM_pos = TrackDirGlobal_pos.X()/TrackDirGlobal_pos.Z();
    //double phGEM_pos = TrackDirGlobal_pos.Y()/TrackDirGlobal_pos.Z();

    //   chi2 += pow( (TrackPosLocal.X() - xfp_expect)/0.002, 2 ) + pow( (TrackPosLocal.Y() - yfp_expect)/0.002, 2 );
    chi2 += pow( (TrackPosLocal.X() - xfp_expect)/0.002, 2 ) + pow( (TrackPosLocal.Y() - yfp_expect)/0.002, 2 ) +
      pow( (thGEM - thHCAL)/0.006,2 ) + pow( (phGEM - phHCAL)/0.005, 2 ) + pow( (TrackIntersectHCAL.X()-(xHCAL[i]+HCALPOS.X()))/0.055,2) + pow( (TrackIntersectHCAL.Y()-(yHCAL[i]+HCALPOS.Y()))/0.05,2) + pow( (xtarget-BeamPos.X())/.002, 2) + pow( (ytarget-BeamPos.Y())/.002,2 );
    // chi2 += pow( TrackIntersectGlobal.X()/0.002, 2 ) + pow(TrackIntersectGlobal.Y()/0.002, 2 ) + pow( (TrackIntersectHCAL.X()-(xHCAL[i]+HCALPOS.X()))/0.055,2) + pow( (TrackIntersectHCAL.Y()-(yHCAL[i]+HCALPOS.Y()))/0.05,2) +
    //   pow( (thGEM-thHCAL)/0.006, 2 ) + pow( (phGEM-phHCAL)/0.005, 2 ) +
    //   pow( (thGEM-thGEM_pos)/0.002, 2 ) + pow( (phGEM-phGEM_pos)/0.002, 2 );
      
  }
  
  f = chi2;
}

void AlignZeroField_SBSGEM_StraightThrough_NoSieve( const char *configfilename, const char *outfilename = "ZeroFieldAlign_SBS_nosieve_temp.root" ){
  TChain *C = new TChain("T");

  double GEMX0=0.0, GEMY0=0.0, GEMZ0=4.11;
  double GEMtheta = 0.0;
  double GEMphi = 0.0;
  double GEMax=0.0, GEMay=0.0, GEMaz=0.0; //need to work in terms of yaw, pitch, roll angles for sufficient degrees of freedom:

  double HCALX0 = 0.0, HCALY0 = 0.0, HCALZ0 = 9.007;

  TString prefix = "sbs.gem";

  TCut globalcut = "";
  
  ifstream infile(configfilename);
  if( infile ){
    TString currentline;
    
    while( currentline.ReadLine(infile) && !currentline.BeginsWith("endlist") ){
      if( !currentline.BeginsWith("#") ){
	C->Add(currentline.Data());
      }
    }

    while( currentline.ReadLine(infile) && !currentline.BeginsWith("endcut") ){
      if( !currentline.BeginsWith("#") ){
	globalcut += currentline;
      }
    }

    while( currentline.ReadLine(infile) && !currentline.BeginsWith("endconfig") ){
      if( !currentline.BeginsWith("#") ){
	TObjArray *tokens = currentline.Tokenize(" ");
	int ntokens = tokens->GetEntries();
	if( ntokens >= 2 ){
	  TString skey = ( (TObjString*) (*tokens)[0] )->GetString();
	  TString sval = ( (TObjString*) (*tokens)[1] )->GetString();
	  if( skey == "prefix" ){
	    prefix = sval;
	  }

	  if( skey == "GEMX0" ){
	    GEMX0 = sval.Atof();
	  }

	  if( skey == "GEMY0" ){
	    GEMY0 = sval.Atof();
	  }

	  if( skey == "GEMZ0" ){
	    GEMZ0 = sval.Atof();
	  }

	  if( skey == "HCALX0" ){
	    HCALX0 = sval.Atof();
	  }

	  if( skey == "HCALY0" ){
	    HCALY0 = sval.Atof();
	  }

	  if( skey == "HCALZ0" ){
	    HCALZ0 = sval.Atof();
	  }

	  if( skey == "GEMtheta" ){
	    GEMtheta = sval.Atof();
	  }

	  if( skey == "GEMphi" ){
	    GEMphi = sval.Atof();
	  }

	  if( skey == "GEMax" ){
	    GEMax = sval.Atof();
	  }

	  if( skey == "GEMay" ){
	    GEMay = sval.Atof();
	  }

	  if( skey == "GEMaz" ){
	    GEMaz = sval.Atof();
	  }

	  if( skey == "NMAX" ){
	    NMAX = sval.Atoi();
	  }
	}
      }
    }
    
    C->SetBranchStatus("*",0);
    
    C->SetBranchStatus("sbs.x_bcp",1);
    C->SetBranchStatus("sbs.y_bcp",1);
    C->SetBranchStatus("sbs.z_bcp",1);

    C->SetBranchStatus("sbs.x_fcp",1);
    C->SetBranchStatus("sbs.y_fcp",1);
    C->SetBranchStatus("sbs.z_fcp",1);

    C->SetBranchStatus("sbs.hcal.e",1);
    C->SetBranchStatus("sbs.hcal.x",1);
    C->SetBranchStatus("sbs.hcal.y",1);
    C->SetBranchStatus("sbs.hcal.atimeblk",1);

    C->SetBranchStatus("bb.sh.atimeblk",1);
    C->SetBranchStatus("bb.ps.atimeblk",1);
    C->SetBranchStatus("bb.sh.e",1);
    C->SetBranchStatus("bb.ps.e",1);

    C->SetBranchStatus("bb.tr.*",1); //cut on BB track variables if needed

    TString branchname;

    C->SetBranchStatus( branchname.Format( "%s.track.*", prefix.Data() ), 1 );

    C->SetBranchStatus( "sbs.tr.*", 1 );
   
    const int MAXNTRACKS = 100;

    double ntrack;
    vector<double> tracknhits(MAXNTRACKS);
    vector<double> trackChi2NDF(MAXNTRACKS);
    vector<double> trackX(MAXNTRACKS);
    vector<double> trackY(MAXNTRACKS);
    vector<double> trackXp(MAXNTRACKS);
    vector<double> trackYp(MAXNTRACKS);

    C->SetBranchAddress( branchname.Format("%s.track.ntrack",prefix.Data()), &ntrack );
    C->SetBranchAddress( branchname.Format("%s.track.chi2ndf",prefix.Data()), &(trackChi2NDF[0]) );
    C->SetBranchAddress( branchname.Format("%s.track.nhits",prefix.Data()), &(tracknhits[0]) );
    C->SetBranchAddress( branchname.Format("%s.track.x",prefix.Data()), &(trackX[0]) );
    C->SetBranchAddress( branchname.Format("%s.track.y",prefix.Data()), &(trackY[0]) );
    C->SetBranchAddress( branchname.Format("%s.track.xp",prefix.Data()), &(trackXp[0]) );
    C->SetBranchAddress( branchname.Format("%s.track.yp",prefix.Data()), &(trackYp[0]) );

    double xhcal,yhcal;
    C->SetBranchAddress( "sbs.hcal.x",&xhcal);
    C->SetBranchAddress( "sbs.hcal.y",&yhcal);

    TFile *fout = new TFile(outfilename,"RECREATE");

    TString dbfilename = outfilename;
    dbfilename.ReplaceAll(".root",".dat");

    ofstream dbfile(dbfilename);
    
    TH1D *hxtar_old = new TH1D("hxtar_old","Old alignment; Target x (m);",250,-0.2,0.2);
    TH1D *hytar_old = new TH1D("hytar_old","Old alignment; Target y (m);",250,-0.2,0.2);
    TH2D *hxytar_old = new TH2D("hxytar_old", "Old alignment; Target y (m); Target x (m)", 250,-0.2,0.2, 250,-0.2,0.2);

    TH2D *hxtar_x_old = new TH2D("hxtar_x_old","Old alignment; track x (m); target x (m)", 250, -0.8,0.8,250,-0.2,0.2);
    TH2D *hxtar_y_old = new TH2D("hxtar_y_old","Old alignment; track y (m); target x (m)", 250, -0.3,0.3,250,-0.2,0.2);
    TH2D *hxtar_th_old = new TH2D("hxtar_th_old","Old alignment; track dx/dz; target x (m)", 250, -0.5,0.5,250,-0.2,0.2);
    TH2D *hxtar_ph_old = new TH2D("hxtar_ph_old","Old alignment; track dy/dz; target x (m)", 250, -0.125,0.125,250,-0.2,0.2);

    TH2D *hytar_x_old = new TH2D("hytar_x_old","Old alignment; track x (m); target y (m)", 250, -0.8,0.8,250,-0.2,0.2);
    TH2D *hytar_y_old = new TH2D("hytar_y_old","Old alignment; track y (m); target y (m)", 250, -0.3,0.3,250,-0.2,0.2);
    TH2D *hytar_th_old = new TH2D("hytar_th_old","Old alignment; track dx/dz; target y (m)", 250, -0.5,0.5,250,-0.2,0.2);
    TH2D *hytar_ph_old = new TH2D("hytar_ph_old","Old alignment; track dy/dz; target y (m)", 250, -0.125,0.125,250,-0.2,0.2);

    TH1D *hxtar_new = new TH1D("hxtar_new","New alignment; Target x (m);",250,-0.2,0.2);
    TH1D *hytar_new = new TH1D("hytar_new","New alignment; Target y (m);",250,-0.2,0.2);
    TH2D *hxytar_new = new TH2D("hxytar_new", "New alignment; Target y (m); Target x (m)", 250,-0.2,0.2, 250,-0.2,0.2);

    TH2D *hxtar_x_new = new TH2D("hxtar_x_new","New alignment; track x (m); target x (m)", 250, -0.8,0.8,250,-0.2,0.2);
    TH2D *hxtar_y_new = new TH2D("hxtar_y_new","New alignment; track y (m); target x (m)", 250, -0.3,0.3,250,-0.2,0.2);
    TH2D *hxtar_th_new = new TH2D("hxtar_th_new","New alignment; track dx/dz; target x (m)", 250, -0.5,0.5,250,-0.2,0.2);
    TH2D *hxtar_ph_new = new TH2D("hxtar_ph_new","New alignment; track dy/dz; target x (m)", 250, -0.125,0.125,250,-0.2,0.2);

    TH2D *hytar_x_new = new TH2D("hytar_x_new","New alignment; track x (m); target y (m)", 250, -0.8,0.8,250,-0.2,0.2);
    TH2D *hytar_y_new = new TH2D("hytar_y_new","New alignment; track y (m); target y (m)", 250, -0.3,0.3,250,-0.2,0.2);
    TH2D *hytar_th_new = new TH2D("hytar_th_new","New alignment; track dx/dz; target y (m)", 250, -0.5,0.5,250,-0.2,0.2);
    TH2D *hytar_ph_new = new TH2D("hytar_ph_new","New alignment; track dy/dz; target y (m)", 250, -0.125,0.125,250,-0.2,0.2);
    
    TTreeFormula *GlobalCut = new TTreeFormula("GlobalCut",globalcut,C);

    long nevent=0; 
    
    int treenum=-1, currenttreenum=-1;
    
    xHCAL.clear();
    yHCAL.clear();
    XTRACK.clear();
    YTRACK.clear();
    XPTRACK.clear();
    YPTRACK.clear();
    NTRACKS = 0;

    TVector3 GEMPOS(GEMX0,GEMY0,GEMZ0);
    //	TVector3 GEMdir(tan(GEMtheta), tan(GEMphi), 1.0 );
    //GEMdir = GEMdir.Unit();
    
    TVector3 xaxis_g(1,0,0);
    TVector3 yaxis_g(0,1,0);
    TVector3 zaxis_g(0,0,1);
    
    TRotation R;
    R.RotateZ(GEMaz);
    R.RotateX(GEMax);
    R.RotateY(GEMay);
    
    
    TRotation Rinverse = R.Inverse();
    //Rinverse.Invert();
    
    //TVector3 GEMdir = Rinverse * zaxis_g = | Rxx  Ryx  Rzx | * | 0 |   | Rzx |
    //                                       | Rxy  Ryy  Rzy | * | 0 | = | Rzy |, WRONG!
    //                                       | Rxz  Ryz  Rzz | * | 1 |   | Rzz |
    
    
    TVector3 GEMdir = R * zaxis_g;
    //    GEMdir.SetXYZ( R.XZ(), R.YZ(), R.ZZ() );
    
    //   | Rxx Rxy Rxz | * | 0 |   | Rxz | 
    // = | Ryx Ryy Ryz | * | 0 | = | Ryz |, RIGHT! 
    //   | Rzx Rzy Rzz | * | 1 |   | Rzz |
    
    //TVector3 xaxis_gem = yaxis_g.Cross(GEMdir).Unit();
    //TVector3 yaxis_gem = GEMdir.Cross(xaxis_gem).Unit();
    //TVector3 xaxis_gem = R * xaxis_g; // = ( Rxx Ryx Rzx )
    //TVector3 yaxis_gem = R * yaxis_g; // = ( Rxy Ryy Rzy )
    TVector3 xaxis_gem = R * xaxis_g;
    //  xaxis_gem.SetXYZ( R.XX(), R.YX(), R.ZX() );
    TVector3 yaxis_gem = R * yaxis_g;
    //   yaxis_gem.SetXYZ( R.XY(), R.YY(), R.ZY() );
    
    //Let's try a linearized alignment procedure: 
    
    int niter = 10;
    TMatrixD M_LHS(6,6);
    TVectorD b_RHS(6);
    
    //    for( int iter=0; iter<niter; iter++ ){
      
    for( int ipar=0; ipar<6; ipar++ ){
      for( int jpar=0; jpar<6; jpar++ ){
	M_LHS(ipar,jpar) = 0.0;
      }
      b_RHS(ipar) = 0.0;
    }
    
    while( C->GetEntry(nevent++) ){
      if( nevent % 10000 == 0 ) cout << "Event " << nevent << endl;
      
      currenttreenum = C->GetTreeNumber();
      if( currenttreenum != treenum ){
	treenum = currenttreenum;
	GlobalCut->UpdateFormulaLeaves();
      }
      
      bool passed_cut = GlobalCut->EvalInstance(0) != 0;
      
      if( passed_cut ){
	if( NTRACKS < NMAX && int(ntrack)>0 ){
	  XTRACK.push_back( trackX[0] );
	  YTRACK.push_back( trackY[0] );
	  XPTRACK.push_back( trackXp[0] );
	  YPTRACK.push_back( trackYp[0] );
	  xHCAL.push_back( xhcal );
	  yHCAL.push_back( yhcal );
	  NTRACKS++;
	}
	//coefficients of the offset parameters in terms of the xtarget equation:
	//  order is { x0, y0, z0, ax, ay, az }:
	double xcoeff[6] =
	  { -1.0,
	    0.0,
	    trackXp[0],
	    (GEMZ0 * trackYp[0] + trackXp[0]*trackY[0]),
	    -trackXp[0]*(GEMZ0 + trackX[0]),
	    GEMZ0 * trackYp[0] };
	double ycoeff[6] =
	  { 0.0,
	    -1.0,
	    trackYp[0],
	    trackYp[0]*(GEMZ0 + trackY[0]),
	    -(GEMZ0 * trackXp[0] + trackYp[0]*trackX[0]),
	    -GEMZ0 * trackXp[0] };

	double xRHS = trackX[0] - trackXp[0]*GEMZ0 - 0.006;
	double yRHS = trackY[0] - trackYp[0]*GEMZ0;

	for( int ipar=0; ipar<6; ipar++ ){
	  for( int jpar=0; jpar<6; jpar++ ){
	    M_LHS(ipar,jpar) += xcoeff[ipar]*xcoeff[jpar] + ycoeff[ipar]*ycoeff[jpar];
	  }
	  b_RHS(ipar) += xcoeff[ipar]*xRHS + ycoeff[ipar]*yRHS;
	}
	
	TVector3 TrackPosLocal(trackX[0],trackY[0],0.0);
	
	//TVector3 TrackPosGlobal = GEMPOS + xaxis_gem * TrackPosLocal.X() + yaxis_gem * TrackPosLocal.Y() + GEMdir * TrackPosLocal.Z();
	TVector3 TrackPosGlobal = GEMPOS + R * TrackPosLocal;
	
	TVector3 TrackDirLocal(trackXp[0],trackYp[0],1.0);
	TrackDirLocal = TrackDirLocal.Unit();
	  
	TVector3 TrackDirGlobal = Rinverse * TrackDirLocal;
	//TVector3 TrackDirGlobal = TrackDirLocal.X() * xaxis_gem + TrackDirLocal.Y() * yaxis_gem + TrackDirLocal.Z() * GEMdir;
	
	//   | Rxx Ryx Rzx |   | X' |   | Rxx * X' + Ryx * Y' + Rzx * Z' |   | dirlocal dot xaxisglobal |
	// = | Rxy Ryy Rzy | * | Y' | = | Rxy * X' + Ryy * Y' + Rzy * Z' | = | dirlocal dot yaxisglobal | = ? 
	//   | Rxz Ryz RZZ |   | Z' |   | Rxz * X' + Ryz * Y' + Rzz * Z' |   | dirlocal dot zaxisglobal |
	// OTOH;
	//                                                         | Rxx |        | Rxy |        | Rxz | 
	// X' * GEM_xaxis + Y' * GEM_yaxis + Z' * GEM_zaxis = X' * | Ryx | + Y' * | Ryy | + Z' * | Ryz |       
	//                                                         | Rzx |        | Rzy |        | Rzz |
	double sintersect = -TrackPosGlobal.Dot( zaxis_g ) / TrackDirGlobal.Dot( zaxis_g );
	TVector3 TrackIntersectGlobal = TrackPosGlobal + sintersect * TrackDirGlobal;
	
	//if( iter == 0 ){
	hxtar_old->Fill( TrackIntersectGlobal.X() );
	hytar_old->Fill( TrackIntersectGlobal.Y() );
	hxytar_old->Fill( TrackIntersectGlobal.Y(), TrackIntersectGlobal.X() );
	
	double xtar = TrackIntersectGlobal.X();
	double ytar = TrackIntersectGlobal.Y();
	hxtar_x_old->Fill( trackX[0], xtar );
	hxtar_y_old->Fill( trackY[0], xtar );
	hxtar_th_old->Fill( trackXp[0], xtar );
	hxtar_ph_old->Fill( trackYp[0], xtar );
	
	hytar_x_old->Fill( trackX[0], ytar );
	hytar_y_old->Fill( trackY[0], ytar );
	hytar_th_old->Fill( trackXp[0], ytar );
	hytar_ph_old->Fill( trackYp[0], ytar );

	
	
      }
    }
     

    TMinuit *FitFunc = new TMinuit( 9 );
    FitFunc->SetFCN( CHI2_FCN );

    

    int ierflg = 0;
    FitFunc->mnparm( 0, "GEMX0", GEMX0, 0.01,0,0,ierflg );
    FitFunc->mnparm( 1, "GEMY0", GEMY0, 0.02,0,0,ierflg );
    FitFunc->mnparm( 2, "GEMZ0", GEMZ0, 0.05,0,0,ierflg );
    FitFunc->mnparm( 3, "GEMax", GEMax, 0.5*TMath::Pi()/180.0,0,0,ierflg );
    FitFunc->mnparm( 4, "GEMay", GEMay, 0.5*TMath::Pi()/180.0,0,0,ierflg );
    FitFunc->mnparm( 5, "GEMaz", GEMaz, 0.5*TMath::Pi()/180.0,0,0,ierflg );
    FitFunc->mnparm( 6, "HCALX0", HCALX0, 0.001,0,0,ierflg);
    FitFunc->mnparm( 7, "HCALY0", HCALY0, 0.001,0,0,ierflg);
    FitFunc->mnparm( 8, "HCALZ0", HCALZ0, 0.001,0,0,ierflg);

    //Fix HCAL position
    FitFunc->FixParameter(6);
    FitFunc->FixParameter(7);
    FitFunc->FixParameter(8);

    //Fix orientation angles. Comment out to fit rotational offsets
    //FitFunc->FixParameter(3);
    //FitFunc->FixParameter(4);
    //FitFunc->FixParameter(5);
    //FitFunc->FixParameter(8);

    double arglist[10];
    arglist[0] = 1;
    FitFunc->mnexcm( "SET ERR", arglist, 1, ierflg );
    arglist[0] = 5000;
    arglist[1] = 1.;
    
    FitFunc->mnexcm("MIGRAD",arglist,2,ierflg);

    double dGEMax, dGEMay, dGEMaz;
    double dGEMX0,dGEMY0,dGEMZ0,dGEMtheta,dGEMphi;
    double dHCALX0,dHCALY0,dHCALZ0;
    FitFunc->GetParameter( 0, GEMX0, dGEMX0 );
    FitFunc->GetParameter( 1, GEMY0, dGEMY0 );
    FitFunc->GetParameter( 2, GEMZ0, dGEMZ0 );
    FitFunc->GetParameter( 3, GEMax, dGEMax );
    FitFunc->GetParameter( 4, GEMay, dGEMay );
    FitFunc->GetParameter( 5, GEMaz, dGEMaz );
    FitFunc->GetParameter( 6, HCALX0, dHCALX0 );
    FitFunc->GetParameter( 7, HCALY0, dHCALY0 );
    FitFunc->GetParameter( 8, HCALZ0, dHCALZ0 );

    
    std::cout << "After fit: (GEMX0,GEMY0,GEMZ0)=(" << GEMX0 << ", " << GEMY0
	      << ", " << GEMZ0 << ")" << std::endl;
    std::cout << "After fit: (GEMax,GEMay,GEMaz)=(" << GEMax << ", " << GEMay
	      << ", " << GEMaz << ")" << std::endl;

    // TVector3 xaxis_g(1,0,0);
    // TVector3 yaxis_g(0,1,0);
    // TVector3 zaxis_g(0,0,1);
    
    TRotation Rtot; 
    Rtot.RotateZ(GEMaz);
    Rtot.RotateX(GEMax);
    Rtot.RotateY(GEMay);
    

    TRotation Rinv = Rtot.Inverse();
    
    GEMdir = Rtot * zaxis_g;
    //GEMdir.SetXYZ( Rtot.ZX(), Rtot.ZY(), Rtot.ZZ() );
    //	TVector3 xaxis_gem = yaxis_g.Cross(GEMdir).Unit();
    //TVector3 yaxis_gem = GEMdir.Cross(xaxis_gem).Unit();
    
    xaxis_gem = Rtot * xaxis_g;
    yaxis_gem = Rtot * yaxis_g;

    //xaxis_gem.SetXYZ( Rtot.XX(), Rtot.XY(), Rtot.XZ() );
    //yaxis_gem.SetXYZ( Rtot.YX(), Rtot.YY(), Rtot.YZ() );

    std::cout << "Total rotation = " << std::endl;
    std::cout << "| R_xx,  R_xy,  R_xz | = |" << Rtot.XX() << ", " << Rtot.XY() << ", " << Rtot.XZ() << "|" << std::endl
	      << "| R_yx,  R_yy,  R_yz | = |" << Rtot.YX() << ", " << Rtot.YY() << ", " << Rtot.YZ() << "|" << std::endl
	      << "| R_zx,  R_zy,  R_zz | = |" << Rtot.ZX() << ", " << Rtot.ZY() << ", " << Rtot.ZZ() << "|" << std::endl;

    std::cout << "GEM x axis = \n";
    xaxis_gem.Print();
    std::cout << "GEM y axis = \n";
    yaxis_gem.Print(); 
    std::cout << "GEM z axis = \n";
    GEMdir.Print();

    std::cout << "Rtot * xaxis global = \n";
    TVector3 xtest = Rtot * xaxis_g;
    xtest.Print();
    std::cout << "Rtot * yaxis global = \n";
    TVector3 ytest = Rtot * yaxis_g;
    ytest.Print();
    std::cout << "Rtot * zaxis global = \n";
    TVector3 ztest = Rtot * zaxis_g;
    ztest.Print();

    std::cout << "Rinv * xaxis global = \n";
    xtest = Rinv * xaxis_g;
    xtest.Print();
    std::cout << "Rinv * yaxis global = \n";
    ytest = Rinv * yaxis_g;
    ytest.Print();
    std::cout << "Rinv * zaxis global = \n";
    ztest = Rinv * zaxis_g;
    ztest.Print();
    
    std::cout << "Inverse rotation = " << std::endl;
    std::cout << "| Rinv_xx,  Rinv_xy,  Rinv_xz | = |" << Rinv.XX() << ", " << Rinv.XY() << ", " << Rinv.XZ() << "|" << std::endl
	      << "| Rinv_yx,  Rinv_yy,  Rinv_yz | = |" << Rinv.YX() << ", " << Rinv.YY() << ", " << Rinv.YZ() << "|" << std::endl
	      << "| Rinv_zx,  Rinv_zy,  Rinv_zz | = |" << Rinv.ZX() << ", " << Rinv.ZY() << ", " << Rinv.ZZ() << "|" << std::endl;

    TDecompSVD LinearFit( M_LHS );
    LinearFit.Solve( b_RHS );

    std::cout << "Linear fit results: (GEMX0, GEMY0, GEMZ0, ax, ay, az) = (";
    for( int ipar=0; ipar<6; ipar++ ){
      cout << b_RHS(ipar);
      cout << (ipar<5 ? ", " : ")");
    }
    cout << endl;

    //The linear fit results should be interpreted as the offsets and rotations that must be applied to the
    // focal-plane coordinates to recover the target coordinates (and slopes). Therefore: 

    // GEMX0 = b_RHS(0);
    // GEMY0 = b_RHS(1);
    // GEMZ0 += b_RHS(2);
    // GEMax = b_RHS(3);
    // GEMay = b_RHS(4);
    // GEMaz = b_RHS(5);

    TRotation Rlinear;
    Rlinear.RotateZ(GEMaz);
    Rlinear.RotateX(GEMax);
    Rlinear.RotateY(GEMay);

    //Rlinear.Invert();
    
    nevent=0; 

    currenttreenum = -1;
    treenum = -1;

    while( C->GetEntry( nevent++ ) ){
      if( nevent % 10000 == 0 ) cout << "Event " << nevent << endl;

      currenttreenum = C->GetTreeNumber();
      if( currenttreenum != treenum ){
	treenum = currenttreenum;
	GlobalCut->UpdateFormulaLeaves();
      }

      bool passed_cut = GlobalCut->EvalInstance(0) != 0; 
      
      if( passed_cut ){
	TVector3 GEMPOS(GEMX0,GEMY0,GEMZ0);
	//TVector3 GEMdir(tan(GEMtheta), tan(GEMphi), 1.0 );
	//GEMdir = GEMdir.Unit();
	
	TVector3 TrackPosLocal( trackX[0], trackY[0], 0.0 );

	//	TVector3 TrackPosGlobal = GEMPOS + Rtot * TrackPosLocal;
	//TVector3 TrackPosGlobal = GEMPOS + TrackPosLocal.X() * xaxis_gem +
	//  TrackPosLocal.Y() * yaxis_gem + TrackPosLocal.Z() * GEMdir;
	TVector3 TrackPosGlobal = GEMPOS + Rtot * TrackPosLocal;
	//TVector3 TrackPosGlobal = GEMPOS + Rinv * TrackPosLocal;
	TVector3 TrackDirLocal(trackXp[0],trackYp[0],1.0);
	TrackDirLocal = TrackDirLocal.Unit();
	//TVector3 TrackDirGlobal = Rinv * TrackDirLocal;

	//TVector3 TrackDirGlobal = TrackDirLocal.X() * xaxis_gem + TrackDirLocal.Y() * yaxis_gem + TrackDirLocal.Z() * GEMdir;

	TVector3 TrackDirGlobal = Rinv * TrackDirLocal;
	
	// cout << "Local track position:" << endl;
	// TrackPosLocal.Print();
	// cout << "Global track position, standard method:" << endl;
	// TrackPosGlobal.Print();
	// cout << "Global track position, GEMPOS + R * poslocal:" << endl;
	// TVector3 postest1 = GEMPOS + Rtot * TrackPosLocal;
	// postest1.Print();
	// cout << "Global track position, GEMPOS + Rinv * poslocal:" << endl;  
	// TVector3 postest2 = GEMPOS + Rinv * TrackPosLocal;
	// postest2.Print();

	// cout << "Unit vector corresponding to rotated position, standard:" << endl;
	// TrackPosGlobal.Unit().Print();
	// cout << "Unit vector corresponding to rotated position, GEMPOS + R * poslocal:" << endl;
	// postest1.Unit().Print();
	// cout << "Unit vector corresponding to rotated position, GEMPOS + Rinv * poslocal:" << endl;
	// postest2.Unit().Print();

	// cout << "Track direction, local: " << endl;
	// TrackDirLocal.Print();
	// cout << "Track direction, global, Rinv method: " << endl;
	// TrackDirGlobal.Print();
	// TVector3 dirtest1 = GEMdir * TrackDirLocal.Z() + xaxis_gem * TrackDirLocal.X() + yaxis_gem * TrackDirLocal.Y();
	// TVector3 dirtest2 = Rtot * TrackDirLocal;

	// cout << "Track direction, global, standard method: " << endl;
	// dirtest1.Print();
	// cout << "Track direction, global, Rtot method: " << endl;
	// dirtest2.Print();

	// cout << endl << endl;
	
	
	double sintersect = -TrackPosGlobal.Dot( zaxis_g ) / TrackDirGlobal.Dot( zaxis_g );
	TVector3 TrackIntersectGlobal = TrackPosGlobal + sintersect * TrackDirGlobal;

	hxtar_new->Fill( TrackIntersectGlobal.X() );
	hytar_new->Fill( TrackIntersectGlobal.Y() );
	hxytar_new->Fill( TrackIntersectGlobal.Y(), TrackIntersectGlobal.X() );
	
	double xtar = TrackIntersectGlobal.X();
	double ytar = TrackIntersectGlobal.Y();

	hxtar_x_new->Fill( trackX[0], xtar );
	hxtar_y_new->Fill( trackY[0], xtar );
	hxtar_th_new->Fill( trackXp[0], xtar );
	hxtar_ph_new->Fill( trackYp[0], xtar );

	hytar_x_new->Fill( trackX[0], ytar );
	hytar_y_new->Fill( trackY[0], ytar );
	hytar_th_new->Fill( trackXp[0], ytar );
	hytar_ph_new->Fill( trackYp[0], ytar );
	
      }
    }

    fout->Write();

    TString posout;
    posout.Form("sbs.gemorigin_xyz = %12.6g %12.6g %12.6g",GEMX0,GEMY0,GEMZ0);
    dbfile << posout.Data() << endl;
    TString rotout;
    // Output angles in degrees since that seems to be the usual convention
    // for database angles in Podd:
    rotout.Form("sbs.gemangles_xyz = %12.6g %12.6g %12.6g", GEMax*180.0/PI, GEMay*180.0/PI, GEMaz*180.0/PI);
    dbfile << rotout.Data() << endl;
    
  } else {
    cout << "Can't open input file, fix" << endl;
    return;
  }
}
