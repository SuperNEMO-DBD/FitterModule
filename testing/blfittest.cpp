#include <TFile.h>
#include <TTreeReader.h>
#include "SNFitter/SNFitter.h"
#include <vector>
#include "catch.hpp"

int readrun() {
  // input reader
  auto* ff = new TFile("../testing/idealline1000.root");
  TTreeReader reader("hit_tree", ff);
  // obtain all the required input data from file
  TTreeReaderValue<std::vector<double>> radius(reader, "radius");
  TTreeReaderValue<std::vector<double>> wirex(reader, "wirex");
  TTreeReaderValue<std::vector<double>> wirey(reader, "wirey");
  TTreeReaderValue<std::vector<double>> wirez(reader, "wirez");
  TTreeReaderValue<std::vector<int>> grid_id(reader, "grid_id");
  TTreeReaderValue<std::vector<int>> grid_side(reader, "grid_side");
  TTreeReaderValue<std::vector<int>> grid_layer(reader, "grid_layer");
  TTreeReaderValue<std::vector<int>> grid_column(reader, "grid_column");

  int counter = 0;
  GeigerRing ring;
  MetaInfo mi;
  TrackerHit th;
  std::vector<TrackerHit> rings;
  SNFitter snf;
  int endevent = reader.GetEntries(true);
  for (int i = 0; i < endevent; i++) {  // event loop
    reader.Next();                      // all data available
    for (unsigned int j = 0; j < radius->size(); j++) {
      ring.rerr = 0.9;
      ring.zerr = 10.0;
      ring.radius = radius->at(j);
      ring.wirex = wirex->at(j);
      ring.wirey = wirey->at(j);
      ring.zcoord = wirez->at(j);
      mi.hitid = grid_id->at(j);
      mi.side = grid_side->at(j);
      mi.row = grid_column->at(j);
      mi.column = grid_layer->at(j);
      th.mi = mi;
      th.gr = ring;
      rings.push_back(th);
    }
    snf.setData(rings);
    if (snf.fitbrokenline().empty()) {
      counter++;
    }
    rings.clear();
  }
  return counter;
}

int check_run() { return readrun(); }

TEST_CASE("BrokenLine", "[falaise][brokenlinerun]") { REQUIRE(check_run() == 2); }
