#ifndef SNFITTER_SNFITTER_H
#define SNFITTER_SNFITTER_H

#include <vector>

#include "SNFitter/DataProducts.h"

// General fitter class with line and helix fitting methods
class SNFitter {
 public:
  SNFitter() = default;
  SNFitter(std::vector<TrackerHit> th) : rings(th) {
    for (auto& hit : rings) grings.push_back(hit.gr);
  }  // Constructor

  ~SNFitter() = default;

  void setData(std::vector<TrackerHit> th) {
    rings.clear();
    grings.clear();
    rings = th;
    for (auto& hit : rings) grings.push_back(hit.gr);
  }

  std::vector<LineFit> fitline();  // sets the data and operates
  std::vector<HelixFit> fithelix();
  std::vector<BrokenLineFit> fitbrokenline();

 protected:
  std::vector<double> line_initials(double frad);
  std::vector<double> helix_initials();
  std::vector<double> helixbackup();
  static std::vector<int> kink_finder(std::vector<double> betaangles,
                                      std::vector<double> errangles);
  static LineFit fitline2D(const std::vector<PathPoint>& data);
  static bool peak_alarm(std::vector<double> betaangles);
  static double martingale(double previous, double val, double alpha);

 private:
  std::vector<TrackerHit> rings;
  std::vector<GeigerRing> grings;
};

#endif