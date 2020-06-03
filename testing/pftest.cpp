#include <fitter_library.h>
#include <vector>
#include "catch.hpp"

std::vector<PathPointCollection> testPF() {
  // pattern
  // . o o o .
  std::vector<TrackerHit> data;
  GeigerRing gr1;
  MetaInfo mi1;
  TrackerHit TH1;
  gr1.radius = 20.0;  // [mm]
  gr1.wirex = 36.0;
  gr1.wirey = 29.0;
  gr1.zcoord = 0.0;
  gr1.rerr = 0.9;
  gr1.zerr = 9.0;

  mi1.hitid = 0;
  mi1.side = 0;
  mi1.row = 10;
  mi1.column = 1;

  TH1.gr = gr1;
  TH1.mi = mi1;
  data.push_back(TH1);

  GeigerRing gr2;
  MetaInfo mi2;
  TrackerHit TH2;
  gr2.radius = 15.0;  // [mm]
  gr2.wirex = 80.0;
  gr2.wirey = 29.0;
  gr2.zcoord = 0.0;
  gr2.rerr = 0.9;
  gr2.zerr = 9.0;

  mi2.hitid = 1;
  mi2.side = 0;
  mi2.row = 10;
  mi2.column = 2;

  TH2.gr = gr2;
  TH2.mi = mi2;
  data.push_back(TH2);

  GeigerRing gr3;
  MetaInfo mi3;
  TrackerHit TH3;
  gr3.radius = 10.0;  // [mm]
  gr3.wirex = 124.0;
  gr3.wirey = 29.0;
  gr3.zcoord = 0.0;
  gr3.rerr = 0.9;
  gr3.zerr = 9.0;

  mi3.hitid = 2;
  mi3.side = 0;
  mi3.row = 10;
  mi3.column = 3;

  TH3.gr = gr3;
  TH3.mi = mi3;
  data.push_back(TH3);

  PathFinder pf(data);
  pf.create_paths();
  std::vector<PathPointCollection> pathcollection = pf.allpaths();
  return pathcollection;
}

std::vector<PathPointCollection> testPF2() {
  // pattern
  // . o .
  // o o .
  std::vector<TrackerHit> data;
  GeigerRing gr1;
  MetaInfo mi1;
  TrackerHit TH1;
  gr1.radius = 20.0;  // [mm]
  gr1.wirex = 36.0;
  gr1.wirey = 29.0;
  gr1.zcoord = 0.0;
  gr1.rerr = 0.9;
  gr1.zerr = 9.0;

  mi1.hitid = 0;
  mi1.side = 0;
  mi1.row = 10;
  mi1.column = 1;

  TH1.gr = gr1;
  TH1.mi = mi1;
  data.push_back(TH1);

  GeigerRing gr2;
  MetaInfo mi2;
  TrackerHit TH2;
  gr2.radius = 15.0;  // [mm]
  gr2.wirex = 80.0;
  gr2.wirey = 29.0;
  gr2.zcoord = 0.0;
  gr2.rerr = 0.9;
  gr2.zerr = 9.0;

  mi2.hitid = 1;
  mi2.side = 0;
  mi2.row = 10;
  mi2.column = 2;

  TH2.gr = gr2;
  TH2.mi = mi2;
  data.push_back(TH2);

  GeigerRing gr3;
  MetaInfo mi3;
  TrackerHit TH3;
  gr3.radius = 10.0;  // [mm]
  gr3.wirex = 80.0;
  gr3.wirey = 73.0;
  gr3.zcoord = 0.0;
  gr3.rerr = 0.9;
  gr3.zerr = 9.0;

  mi3.hitid = 2;
  mi3.side = 0;
  mi3.row = 11;
  mi3.column = 2;

  TH3.gr = gr3;
  TH3.mi = mi3;
  data.push_back(TH3);

  PathFinder pf(data);
  pf.create_paths();
  std::vector<PathPointCollection> pathcollection = pf.allpaths();
  return pathcollection;
}

std::vector<PathPointCollection> testPFsim() {
  // real line, event 5 in illumination.tsim
  double dataxc[4] = {53.0, 97.0, 141.0, 185.0};
  int col[4] = {1, 2, 3, 4};
  double datayc[4] = {880.0, 880.0, 836.0, 836.0};
  int row[4] = {1, 1, 2, 2};
  double datazc[4] = {365.464, 651.376, 1016.16, 1302.07};
  double rad[4] = {2.76276, 14.4634, 16.2518, 4.55114};

  std::vector<TrackerHit> data;
  GeigerRing gr;
  MetaInfo mi;
  TrackerHit TH;
  for (int i = 0; i < 4; i++) {
    gr.radius = rad[i];
    gr.wirex = dataxc[i];
    gr.wirey = datayc[i];
    gr.zcoord = datazc[i];
    gr.rerr = 0.9;
    gr.zerr = 9.0;

    mi.hitid = i;
    mi.side = 0;
    mi.row = row[i];
    mi.column = col[i];

    TH.gr = gr;
    TH.mi = mi;
    data.push_back(TH);
  }

  PathFinder pf(data);
  pf.create_paths();
  std::vector<PathPointCollection> pathcollection = pf.allpaths();
  return pathcollection;
}

std::vector<PathPointCollection> testPFbl() {
  // real line, event 1 in multiscatter_calo.tsim
  double dataxc[10] = {53.0, 97.0, 141.0, 185.0, 229.0, 273.0, 273.0, 317.0, 361.0, 405.0};
  int col[10] = {0, 1, 2, 3, 4, 5, 5, 6, 7, 8};
  double datayc[10] = {44.0, 44.0, 88.0, 88.0, 132.0, 132.0, 176.0, 176.0, 220.0, 220.0};
  int row[10] = {55, 55, 54, 54, 53, 53, 52, 52, 51, 51};
  double datazc[10] = {5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0};
  double rad[10] = {11.8831, 9.96064, 6.39051, 15.4532, 0.897902,
                    20.9459, 17.2491, 4.5947,  11.7564, 11.2643};

  std::vector<TrackerHit> data;
  GeigerRing gr;
  MetaInfo mi;
  TrackerHit TH;
  for (int i = 0; i < 10; i++) {
    gr.radius = rad[i];
    gr.wirex = dataxc[i];
    gr.wirey = datayc[i];
    gr.zcoord = datazc[i];
    gr.rerr = 0.9;
    gr.zerr = 9.0;

    mi.hitid = i;
    mi.side = 0;
    mi.row = row[i];
    mi.column = col[i];

    TH.gr = gr;
    TH.mi = mi;
    data.push_back(TH);
  }

  PathFinder pf(data);
  pf.create_paths();
  std::vector<PathPointCollection> pathcollection = pf.allpaths();
  return pathcollection;
}

bool check_linepath() {
  bool fails = true;
  std::vector<PathPointCollection> ppc = testPF();
  double l = ppc.at(0).length;
  if (l < 87.44 && l > 87.42) fails = false;
  return fails;
}

unsigned int check_trianglepath() {
  std::vector<PathPointCollection> ppc = testPF2();
  return ppc.size();  // 23
}

unsigned int check_simpath() {
  std::vector<PathPointCollection> ppc = testPFsim();
  return ppc.size() + ppc.at(0).path.size();  // 6
}

bool check_blpath() {
  bool fails = true;
  std::vector<PathPointCollection> ppc = testPFbl();
  double l = ppc.at(1).length;
  if (l < 398.3 && l > 398.2) fails = false;
  return fails;
}

TEST_CASE("PF A", "[falaise][line3p]") { REQUIRE(check_linepath() == false); }

TEST_CASE("PF B", "[falaise][triangle]") { REQUIRE(check_trianglepath() == 23); }

TEST_CASE("PF C", "[falaise][simpath]") { REQUIRE(check_simpath() == 6); }

TEST_CASE("PF D", "[falaise][blpath]") { REQUIRE(check_blpath() == false); }
