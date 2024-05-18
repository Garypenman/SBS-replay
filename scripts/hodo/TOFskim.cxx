#include "TH1D.h"
#include "TH2D.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

#include <ROOT/RDataFrame.hxx>

void TOFskim(const char *filebase = "/volatile/halla/sbs/gpenman/GEN_REPLAY/rootfiles/"){

  TChain *C = new TChain("T");
  C->Add(filebase+"e1209016_fullreplay_*.root");

  ROOT::EnableImplicitMT(32);
  ROOT::RDataFrame d0("T",C);
  auto df = d0.Filter("bb.tr.n>0  &&  sbs.hca.e>0.05  &&  bb.ps.e>0.2  &&  bb.grinch_tdc.clus.trackindex==0  &&  bb.grinch_tdc.clus.size>=3  &&  bb.gem.track.nhits[0]>=3  &&  e.kine.W2>0.25  &&  e.kine.W2<1.7");

  ROOT::RDF::RSnapshotOptions opts;
  //opts.fLazy = true;
  opts.fMode = "RECREATE";
  std::vector<std::string> colnames;
  colnames.push_back( "e.kine.W2" );
  colnames.push_back( "bb.tr.p" );
  colnames.push_back( "bb.tr.px" );
  colnames.push_back( "bb.tr.py" );
  colnames.push_back( "bb.tr.pz" );
  colnames.push_back( "bb.tr.vz" );
  colnames.push_back( "bb.tr.th" );
  colnames.push_back( "bb.tr.ph" );
  colnames.push_back( "bb.tr.y" );
  colnames.push_back( "bb.hodotdc.clus.tleft" );
  colnames.push_back( "bb.hodotdc.clus.tright" );
  colnames.push_back( "bb.hodotdc.clus.totleft" );
  colnames.push_back( "bb.hodotdc.clus.totright" );
  colnames.push_back( "bb.hodotdc.clus.id" );
  colnames.push_back( "bb.hodotdc.clus.tmean" );
  colnames.push_back( "bb.hodotdc.clus.tdiff" );
  colnames.push_back( "sbs.hcal.idblk" );
  colnames.push_back( "sbs.hcal.tdctimeblk" );
  colnames.push_back( "sbs.hcal.eblk" );

  auto snap = df.Snapshot("T","TOFskim.root",colnames,opts);
  
}
