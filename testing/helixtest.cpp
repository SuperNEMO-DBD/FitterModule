#include "catch.hpp"
#include <fitter_library.h>
#include <TTreeReader.h>
#include <TFile.h>
#include <vector>


int readrun() {
  // input reader
  TFile* ff = new TFile("../testing/idealhelix1000.root");
  TTreeReader reader("hit_tree", ff);
  // obtain all the required input data from file
  TTreeReaderValue<std::vector<double>> dirx(reader, "dirx");
  TTreeReaderValue<std::vector<double>> diry(reader, "diry");
  TTreeReaderValue<std::vector<double>> dirz(reader, "dirz");
  TTreeReaderValue<std::vector<double>> pointx(reader, "pointx");
  TTreeReaderValue<std::vector<double>> pointy(reader, "pointy");
  TTreeReaderValue<std::vector<double>> pointz(reader, "pointz");
  TTreeReaderValue<std::vector<double>> radius(reader, "radius");
  TTreeReaderValue<std::vector<double>> wirex(reader, "wirex");
  TTreeReaderValue<std::vector<double>> wirey(reader, "wirey");
  TTreeReaderValue<std::vector<double>> wirez(reader, "wirez");
  TTreeReaderValue<std::vector<int>> grid_id(reader, "grid_id");
  TTreeReaderValue<std::vector<int>> grid_side(reader, "grid_side");
  TTreeReaderValue<std::vector<int>> grid_layer(reader, "grid_layer");
  TTreeReaderValue<std::vector<int>> grid_column(reader, "grid_column");
  TTreeReaderValue<std::vector<int>> calo_id(reader, "calo_id");
  TTreeReaderValue<std::vector<int>> calo_type(reader, "calo_type");
  TTreeReaderValue<std::vector<int>> calo_side(reader, "calo_side");
  TTreeReaderValue<std::vector<int>> calo_wall(reader, "calo_wall");
  TTreeReaderValue<std::vector<int>> calo_column(reader, "calo_column");
  TTreeReaderValue<std::vector<int>> calo_row(reader, "calo_row");

  int counter = 0;
  GeigerRing ring;
  TrackerHit th;
  std::vector<TrackerHit> rings;
  SNFitter snf;
  int endevent = reader.GetEntries(true);
  for (int i=0; i<endevent; i++) { // event loop
    reader.Next(); // all data available
    for (unsigned int j=0;j<radius->size();j++) {
      ring.rerr   = 0.9;
      ring.zerr = 10.0;
      ring.radius = radius->at(j);
      ring.wirex  = wirex->at(j);
      ring.wirey  = wirey->at(j);
      ring.zcoord = wirez->at(j);
      th.gr = ring;
      rings.push_back(th);
    }
    snf.setData(rings);
    if (snf.fithelix().empty())
      counter++;
    rings.clear();
  }    
  return counter;
}


int check_run(){
  return readrun();
}



TEST_CASE( "Helix", "[falaise][helixrun]" ) {
  REQUIRE( check_run() == 2 );
}

