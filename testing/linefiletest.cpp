#include "catch.hpp"
#include <fitter_library.h>
#include <TTreeReader.h>
#include <TFile.h>
#include <vector>


int readrun() {
  // input reader
  TFile* ff = new TFile("../testing/idealline1000.root");
  TTreeReader reader("hit_tree", ff);
  // obtain all the required input data from file
  TTreeReaderValue<std::vector<double>> radius(reader, "radius");
  TTreeReaderValue<std::vector<double>> wirex(reader, "wirex");
  TTreeReaderValue<std::vector<double>> wirey(reader, "wirey");
  TTreeReaderValue<std::vector<double>> wirez(reader, "wirez");

  int counter = 0;
  int some = 0;
  GeigerRing ring;
  TrackerHit th;
  std::vector<TrackerHit> rings;
  SNFitter snf;
  int endevent = reader.GetEntries(true);
  for (int i=0; i<endevent; i++) { // event loop
    reader.Next(); // all data available
    for (unsigned int j=0;j<radius->size();j++) {
      ring.rerr   = 0.9;
      ring.zerr = 1.0;
      ring.radius = radius->at(j);
      ring.wirex  = wirex->at(j);
      ring.wirey  = wirey->at(j);
      ring.zcoord = wirez->at(j);
      th.gr = ring;
      rings.push_back(th);
    }
    snf.setData(rings);
    std::vector<LineFit> res = snf.fitline();
    
    for (LineFit entry : res) {
      if (entry.status>0)
	some++;
    }
    if (some==4) {
      counter++; // no valid fit
    }
    some = 0;
    rings.clear();
  }
  return counter;
}


int check_run(){
  return readrun();
}



TEST_CASE( "Line", "[falaise][linefilerun]" ) {
  REQUIRE( check_run() == 29 );
}

