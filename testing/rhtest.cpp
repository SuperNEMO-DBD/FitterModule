#include <fitter_library.h>
#include <vector>
#include "catch.hpp"

PathPoint testRH_samey(double r1, double r2, int which) {
  GeigerRing gr1;
  MetaInfo mi1;
  TrackerHit TH1;
  gr1.radius = r1;  // [mm]
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

  GeigerRing gr2;
  MetaInfo mi2;
  TrackerHit TH2;
  gr2.radius = r2;  // [mm]
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

  RelationalHit rh(TH1, TH2);
  std::vector<std::pair<int, int> > iedges = rh.get_edges();
  std::vector<PathPoint> hostpp = rh.nodes_host();
  return hostpp.at(which);
}

PathPoint testRH_samex(double r1, double r2, int which) {
  GeigerRing gr1;
  MetaInfo mi1;
  TrackerHit TH1;
  gr1.radius = r1;  // [mm]
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

  GeigerRing gr2;
  MetaInfo mi2;
  TrackerHit TH2;
  gr2.radius = r2;  // [mm]
  gr2.wirex = 36.0;
  gr2.wirey = 73.0;
  gr2.zcoord = 0.0;
  gr2.rerr = 0.9;
  gr2.zerr = 9.0;

  mi2.hitid = 1;
  mi2.side = 0;
  mi2.row = 11;
  mi2.column = 1;

  TH2.gr = gr2;
  TH2.mi = mi2;

  RelationalHit rh(TH1, TH2);
  std::vector<std::pair<int, int> > iedges = rh.get_edges();
  std::vector<PathPoint> hostpp = rh.nodes_host();
  return hostpp.at(which);
}

PathPoint testRH_sameyminus(double r1, double r2, int which) {
  GeigerRing gr1;
  MetaInfo mi1;
  TrackerHit TH1;
  gr1.radius = r1;  // [mm]
  gr1.wirex = -36.0;
  gr1.wirey = -29.0;
  gr1.zcoord = 0.0;
  gr1.rerr = 0.9;
  gr1.zerr = 9.0;

  mi1.hitid = 0;
  mi1.side = 0;
  mi1.row = 10;
  mi1.column = 1;

  TH1.gr = gr1;
  TH1.mi = mi1;

  GeigerRing gr2;
  MetaInfo mi2;
  TrackerHit TH2;
  gr2.radius = r2;  // [mm]
  gr2.wirex = -80.0;
  gr2.wirey = -29.0;
  gr2.zcoord = 0.0;
  gr2.rerr = 0.9;
  gr2.zerr = 9.0;

  mi2.hitid = 1;
  mi2.side = 0;
  mi2.row = 10;
  mi2.column = 2;

  TH2.gr = gr2;
  TH2.mi = mi2;

  RelationalHit rh(TH1, TH2);
  std::vector<std::pair<int, int> > iedges = rh.get_edges();
  std::vector<PathPoint> hostpp = rh.nodes_host();
  return hostpp.at(which);
}

PathPoint testRH_samexminus(double r1, double r2, int which) {
  GeigerRing gr1;
  MetaInfo mi1;
  TrackerHit TH1;
  gr1.radius = r1;  // [mm]
  gr1.wirex = -36.0;
  gr1.wirey = -29.0;
  gr1.zcoord = 0.0;
  gr1.rerr = 0.9;
  gr1.zerr = 9.0;

  mi1.hitid = 0;
  mi1.side = 0;
  mi1.row = 10;
  mi1.column = 1;

  TH1.gr = gr1;
  TH1.mi = mi1;

  GeigerRing gr2;
  MetaInfo mi2;
  TrackerHit TH2;
  gr2.radius = r2;  // [mm]
  gr2.wirex = -36.0;
  gr2.wirey = -73.0;
  gr2.zcoord = 0.0;
  gr2.rerr = 0.9;
  gr2.zerr = 9.0;

  mi2.hitid = 1;
  mi2.side = 0;
  mi2.row = 11;
  mi2.column = 1;

  TH2.gr = gr2;
  TH2.mi = mi2;

  RelationalHit rh(TH1, TH2);
  std::vector<std::pair<int, int> > iedges = rh.get_edges();
  std::vector<PathPoint> hostpp = rh.nodes_host();
  return hostpp.at(which);
}

PathPoint testRH_diagr(double r1, double r2, int which) {
  GeigerRing gr1;
  MetaInfo mi1;
  TrackerHit TH1;
  gr1.radius = r1;  // [mm]
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

  GeigerRing gr2;
  MetaInfo mi2;
  TrackerHit TH2;
  gr2.radius = r2;  // [mm]
  gr2.wirex = 80.0;
  gr2.wirey = 73.0;
  gr2.zcoord = 0.0;
  gr2.rerr = 0.9;
  gr2.zerr = 9.0;

  mi2.hitid = 1;
  mi2.side = 0;
  mi2.row = 11;
  mi2.column = 2;

  TH2.gr = gr2;
  TH2.mi = mi2;

  RelationalHit rh(TH1, TH2);
  std::vector<std::pair<int, int> > iedges = rh.get_edges();
  std::vector<PathPoint> hostpp = rh.nodes_host();
  return hostpp.at(which);
}

PathPoint testRH_diagl(double r1, double r2, int which) {
  GeigerRing gr1;
  MetaInfo mi1;
  TrackerHit TH1;
  gr1.radius = r1;  // [mm]
  gr1.wirex = 36.0;
  gr1.wirey = 73.0;
  gr1.zcoord = 0.0;
  gr1.rerr = 0.9;
  gr1.zerr = 9.0;

  mi1.hitid = 0;
  mi1.side = 0;
  mi1.row = 11;
  mi1.column = 1;

  TH1.gr = gr1;
  TH1.mi = mi1;

  GeigerRing gr2;
  MetaInfo mi2;
  TrackerHit TH2;
  gr2.radius = r2;  // [mm]
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

  RelationalHit rh(TH1, TH2);
  std::vector<std::pair<int, int> > iedges = rh.get_edges();
  std::vector<PathPoint> hostpp = rh.nodes_host();
  return hostpp.at(which);
}

unsigned int testRH_edges(double r1, double r2) {
  GeigerRing gr1;
  MetaInfo mi1;
  TrackerHit TH1;
  gr1.radius = r1;  // [mm]
  gr1.wirex = 36.0;
  gr1.wirey = 73.0;
  gr1.zcoord = 0.0;
  gr1.rerr = 0.9;
  gr1.zerr = 9.0;

  mi1.hitid = 0;
  mi1.side = 0;
  mi1.row = 11;
  mi1.column = 1;

  TH1.gr = gr1;
  TH1.mi = mi1;

  GeigerRing gr2;
  MetaInfo mi2;
  TrackerHit TH2;
  gr2.radius = r2;  // [mm]
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

  RelationalHit rh(TH1, TH2);
  return rh.get_edges().size();
}

bool check_sx() {
  bool fails = true;
  PathPoint pp = testRH_samex(15.0, 10.0, 0);
  if (pp.xc > 50.9 && pp.xc < 51.0) fails = false;
  return fails;
}

bool check_sxerr() {
  bool fails = true;
  PathPoint pp = testRH_samex(15.0, 10.0, 3);
  if (pp.erry > 1.78 && pp.erry < 1.79) fails = false;
  return fails;
}

bool check_sy() {
  bool fails = true;
  PathPoint pp = testRH_samey(15.0, 10.0, 0);
  if (pp.yc > 43.9 && pp.yc < 44.0) fails = false;
  return fails;
}

bool check_syerr() {
  bool fails = true;
  PathPoint pp = testRH_samey(15.0, 10.0, 3);
  if (pp.erry > 1.58 && pp.erry < 1.59) fails = false;
  return fails;
}

bool check_sxm() {
  bool fails = true;
  PathPoint pp = testRH_samexminus(15.0, 10.0, 1);
  if (pp.xc > -21.1 && pp.xc < -21.0) fails = false;
  return fails;
}

bool check_sxmerr() {
  bool fails = true;
  PathPoint pp = testRH_samexminus(15.0, 10.0, 0);
  if (pp.erry > 1.41 && pp.erry < 1.42) fails = false;
  return fails;
}

bool check_sym() {
  bool fails = true;
  PathPoint pp = testRH_sameyminus(15.0, 10.0, 1);
  if (pp.yc > -14.1 && pp.yc < -14.0) fails = false;
  return fails;
}

bool check_symerr() {
  bool fails = true;
  PathPoint pp = testRH_sameyminus(15.0, 10.0, 0);
  if (pp.errx > 1.21 && pp.errx < 1.22) fails = false;
  return fails;
}

bool check_dr() {
  bool fails = true;
  PathPoint pp = testRH_diagr(15.0, 10.0, 1);
  if (pp.xc > 47.4 && pp.xc < 47.8) fails = false;
  return fails;
}

bool check_drerr() {
  bool fails = true;
  PathPoint pp = testRH_diagr(15.0, 10.0, 2);
  if (pp.erry > 1.17 && pp.erry < 1.18) fails = false;
  return fails;
}

bool check_dl() {
  bool fails = true;
  PathPoint pp = testRH_diagl(15.0, 10.0, 1);
  if (pp.yc > 61.5 && pp.yc < 61.6) fails = false;
  return fails;
}

bool check_dlerr() {
  bool fails = true;
  PathPoint pp = testRH_diagl(15.0, 10.0, 2);
  if (pp.erry > 1.263 && pp.erry < 1.264) fails = false;
  return fails;
}

bool check_drsmall() {
  bool fails = true;
  PathPoint pp = testRH_diagr(20.0, 1.0, 0);
  if (pp.yc > 46.8 && pp.yc < 46.9) fails = false;
  return fails;
}

bool check_drsmallerr() {
  bool fails = true;
  PathPoint pp = testRH_diagr(20.0, 1.0, 2);
  if (pp.errx > 1.09 && pp.errx < 1.1) fails = false;
  return fails;
}

bool check_drsmallrev() {
  bool fails = true;
  PathPoint pp = testRH_diagr(1.0, 20.0, 0);
  if (pp.yc > 28.0 && pp.yc < 30.0) fails = false;
  return fails;
}

bool check_drsmallerrrev() {
  bool fails = true;
  PathPoint pp = testRH_diagr(1.0, 20.0, 2);
  if (pp.errx > 0.99 && pp.errx < 1.01) fails = false;
  return fails;
}

unsigned int check_edges() { return testRH_edges(1.0, 20.0); }

TEST_CASE("RH A", "[falaise][samex]") { REQUIRE(check_sx() == false); }

TEST_CASE("RH B", "[falaise][samexerr]") { REQUIRE(check_sxerr() == false); }

TEST_CASE("RH C", "[falaise][samey]") { REQUIRE(check_sy() == false); }

TEST_CASE("RH D", "[falaise][sameyerr]") { REQUIRE(check_syerr() == false); }

TEST_CASE("RH E", "[falaise][minusx]") { REQUIRE(check_sxm() == false); }

TEST_CASE("RH F", "[falaise][minusxerr]") { REQUIRE(check_sxmerr() == false); }

TEST_CASE("RH G", "[falaise][minusy]") { REQUIRE(check_sym() == false); }

TEST_CASE("RH H", "[falaise][minusyerr]") { REQUIRE(check_symerr() == false); }

TEST_CASE("RH I", "[falaise][diagr]") { REQUIRE(check_dr() == false); }

TEST_CASE("RH J", "[falaise][diagrerr]") { REQUIRE(check_drerr() == false); }

TEST_CASE("RH K", "[falaise][diagl]") { REQUIRE(check_dl() == false); }

TEST_CASE("RH L", "[falaise][diaglerr]") { REQUIRE(check_dlerr() == false); }

TEST_CASE("RH M", "[falaise][diagrsmall]") { REQUIRE(check_drsmall() == false); }

TEST_CASE("RH N", "[falaise][diagrsmallerr]") { REQUIRE(check_drsmallerr() == false); }

TEST_CASE("RH O", "[falaise][smallrev]") { REQUIRE(check_drsmallrev() == false); }

TEST_CASE("RH P", "[falaise][smallreverr]") { REQUIRE(check_drsmallerrrev() == false); }

TEST_CASE("RH Q", "[falaise][nedges]") { REQUIRE(check_edges() == 3); }
