#ifndef SNFITTER_DATAPRODUCTS_H
#define SNFITTER_DATAPRODUCTS_H

#include <utility>
#include <vector>

// data storage objects
struct MetaInfo {
  // uniquely tag a geiger hit
  int hitid;
  int side;
  int row;
  int column;
};

// tracker data storage
struct GeigerRing {
  double radius;
  double wirex;
  double wirey;
  double zcoord;
  double rerr;
  double zerr;
};

// path point store
struct PathPoint {
  // just a 3D point with errors for simple fitting
  std::pair<int, size_t> pointid;
  double xc;
  double yc;
  double zc;
  double errx;
  double erry;
  double errz;
};

// complete path storage
struct PathPointCollection {
  double length;  // path length in Euclidean metric
  std::vector<PathPoint> path;
};

struct TrackerHit {
  // full tracker hit info
  // unique to each
  int clid;
  GeigerRing gr;
  MetaInfo mi;
};

// line fit store, 4 parameter with errors
struct LineFit {
  // par
  double ixy;
  double slxy;
  double ixz;
  double slxz;
  double errixy;
  double errslxy;
  double errixz;
  double errslxz;
  // fit diagnostics
  double chi2;
  double prob;
  int status;
  int clid;
};

// helix fit store, 5 parameter with errors
struct HelixFit {
  // par
  double radius;
  double pitch;
  double xc;
  double yc;
  double zc;
  double raderr;
  double errpitch;
  double errxc;
  double erryc;
  double errzc;
  // fit diagnostics
  double chi2;
  double prob;
  int status;
  int clid;
};

// broken line fit store, 4 parameter with errors and break
struct BrokenLineFit {
  LineFit linefit1;              // has its own diagnostics
  LineFit linefit2;              // could be empty
  std::vector<int> breakpoints;  // first and last of interest
  // fit diagnostics for full BL fit
  std::vector<double> angles;
  std::vector<PathPoint> path;
  double length;  // path length in Euclidean metric
  double chi2;
  double prob;
  int status;
  int clid;
};

#endif