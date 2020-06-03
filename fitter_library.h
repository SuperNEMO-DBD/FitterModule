#ifndef YR_SNFitter
#define YR_SNFitter

// std libraries
#include <utility>
#include <vector>

// ROOT includes
#include <Math/Vector2D.h>
#include <Math/Vector3D.h>
#include <TMath.h>
#include <TVector3.h>

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

// function Object to be minimized
class LineDistance2 {
  // data member
 private:
  std::vector<GeigerRing>* frings;

 protected:
  static double linedistance(GeigerRing gr, const double* p);  // calculate distance line-cylinder

 public:
  LineDistance2(std::vector<GeigerRing>* g) : frings(g) {}
  ~LineDistance2() {}

  // implementation of the function to be minimized
  double operator()(const double* par) {
    double sum = 0.0;
    for (GeigerRing entry : *frings) {
      double d = linedistance(entry, par);
      sum += d * d;  // squared weighted distance
    }
    return sum;
  }
};

// function Object to be minimized
class HelixDistance2 {
  // data member
 private:
  std::vector<GeigerRing>* frings;

 protected:
  static double helixdistance(GeigerRing gr, const double* p);  // calculate distance helix-cylinder

 public:
  HelixDistance2(std::vector<GeigerRing>* g) : frings(g) {}
  ~HelixDistance2() {}

  // implementation of the function to be minimized
  double operator()(const double* par) {
    double sum = 0.0;
    for (GeigerRing entry : *frings) {
      double d = helixdistance(entry, par);
      sum += d * d;  // squared weighted distance
    }
    return sum;
  }
};

// *** Broken Line section
// ***********************
// function Object to be minimized
class BLDistance2 {
  // data member
 private:
  std::vector<PathPoint>* fdata;
  std::vector<double> angles;
  std::vector<double> errangles;

  double bldistance(PathPoint pp, const double* p, ROOT::Math::XYZVector& model,
                    ROOT::Math::XYZVector& weight) {
    // distance line ring is D= | (xp-x0) cross  ux | - r
    // where ux is direction of line and x0 is a point on the line (like t = 0)
    // line not parallel to y-z plane, i.e calo plane
    double x = pp.xc;
    double y = pp.yc;
    double z = pp.zc;

    ROOT::Math::XYZVector xp(x, y, z);
    ROOT::Math::XYZVector x0(0., p[0], p[2]);
    ROOT::Math::XYZVector x1(1., p[0] + p[1], p[2] + p[3]);
    ROOT::Math::XYZVector u = (x1 - x0).Unit();
    ROOT::Math::XYZVector dvec = (xp - x0).Cross(u);

    // errors from xp in cross product
    model = xp + dvec;  // to return from the function
    weight.SetXYZ((pp.erry * u.Z() * pp.erry * u.Z()) + (pp.errz * u.Y() * pp.errz * u.Y()),
                  (pp.errx * u.Z() * pp.errx * u.Z()) + (pp.errz * u.X() * pp.errz * u.X()),
                  (pp.errx * u.Y() * pp.errx * u.Y()) +
                      (pp.erry * u.X() * pp.erry * u.X()));  // weight sqr vector

    double d2 = dvec.Mag2();
    double w2 = weight.X() + weight.Y() + weight.Z();  // sum of squares is .Mag2()
    return d2 / w2;
  }

  std::vector<double> calculate_beta(const double* par) {
    // needs TVector3 objects internally for angles
    std::vector<double> beta;
    ROOT::Math::XYZVector um1(0, 0, 0);                             // to be overwritten
    ROOT::Math::XYZVector ui(0, 0, 0);                              // in bldistance
    ROOT::Math::XYZVector up1(0, 0, 0);                             // function
    ROOT::Math::XYZVector dweight(0, 0, 0);                         // dummy, not needed here
    for (unsigned int j = 1; j < fdata->size() - 1; j++) {          // N-2 angles
      double d1 = bldistance(fdata->at(j - 1), par, um1, dweight);  // ordered path points assumed
      double d2 = bldistance(fdata->at(j), par, ui, dweight);
      double d3 = bldistance(fdata->at(j + 1), par, up1, dweight);

      TVector3 v1((ui - um1).X(), (ui - um1).Y(), (ui - um1).Z());  // make a tvector3
      TVector3 v2((up1 - ui).X(), (up1 - ui).Y(), (up1 - ui).Z());  // make a tvector3
      double angle = v2.Angle(v1);
      if (fabs(angle) > TMath::Pi() / 2.0)  // backward to forward hemisphere
        angle -= TMath::Pi();
      beta.push_back(angle);
      // cout << "bestfit angle check: " << endl;
      // cout << "v1 = (" << v1.x() << ", " << v1.y() << ", " << v1.z() << ")" << endl;
      // cout << "v2 = (" << v2.x() << ", " << v2.y() << ", " << v2.z() << ")" << endl;
      // cout << "angle: " << angle << endl;
    }
    return beta;
  }

  std::vector<double> calculate_angle_errors(const double* par) {
    // needs TVector3 objects internally for angles
    std::vector<double> errbeta;
    ROOT::Math::XYZVector um1(0, 0, 0);  // to be overwritten
    ROOT::Math::XYZVector ui(0, 0, 0);   // in bldistance
    ROOT::Math::XYZVector up1(0, 0, 0);  // function
    ROOT::Math::XYZVector w1(0, 0, 0);
    ROOT::Math::XYZVector w2(0, 0, 0);
    ROOT::Math::XYZVector w3(0, 0, 0);
    for (unsigned int j = 1; j < fdata->size() - 1; j++) {  // N-2 errors
      PathPoint pp1 = fdata->at(j - 1);
      PathPoint pp2 = fdata->at(j);
      PathPoint pp3 = fdata->at(j + 1);
      // residual and error square on connection vector
      double d1 = bldistance(pp1, par, um1, w1);
      double d2 = bldistance(pp2, par, ui, w2);
      double d3 = bldistance(pp3, par, up1, w3);
      // error square on connection vectors v1 and v2
      ROOT::Math::XYZVector dv1 = w1 + w2;
      ROOT::Math::XYZVector dv2 = w2 + w3;
      // prepare the error on the single angle
      // from acos(v1.v2/|v1||v2|)
      TVector3 v1((ui - um1).X(), (ui - um1).Y(), (ui - um1).Z());  // make a tvector3
      TVector3 v2((up1 - ui).X(), (up1 - ui).Y(), (up1 - ui).Z());  // make a tvector3
      double angle = v2.Angle(v1);
      if (angle > TMath::Pi() / 2.0)  // backward to forward hemisphere
        angle -= TMath::Pi();
      double constant = fabs(TMath::Sin(angle));
      // denominators
      double denom1 = v1.Mag() * v2.Mag();
      double denom2 = v1.Mag() * v1.Mag() * v1.Mag() * v2.Mag();
      double denom3 = v1.Mag() * v2.Mag() * v2.Mag() * v2.Mag();
      // derivative terms for v1 and v2
      double tv1x = v2.x() / denom1 - v1.x() * (v1.Dot(v2)) / denom2;
      double tv1y = v2.y() / denom1 - v1.y() * (v1.Dot(v2)) / denom2;
      double tv1z = v2.z() / denom1 - v1.z() * (v1.Dot(v2)) / denom2;
      double tv2x = v1.x() / denom1 - v2.x() * (v1.Dot(v2)) / denom3;
      double tv2y = v1.y() / denom1 - v2.y() * (v1.Dot(v2)) / denom3;
      double tv2z = v1.z() / denom1 - v2.z() * (v1.Dot(v2)) / denom3;
      // component error square terms
      double cerrx = constant * (dv1.X() * tv1x * tv1x + dv2.X() * tv2x * tv2x);
      double cerry = constant * (dv1.Y() * tv1y * tv1y + dv2.Y() * tv2y * tv2y);
      double cerrz = constant * (dv1.Z() * tv1z * tv1z + dv2.Z() * tv2z * tv2z);
      double sinerr = TMath::Sqrt(cerrx + cerry + cerrz);
      if (sinerr > 1.0) sinerr = 1.0;
      double angle_error = TMath::ASin(sinerr);
      if (angle_error < 1.0e-5) angle_error = 1.0e-5;
      errbeta.push_back(angle_error);
    }
    return errbeta;
  }

 public:
  // constructor
  BLDistance2(std::vector<PathPoint>* g) : fdata(g) {}

  // implementation of the function to be minimized
  double operator()(const double* par) {
    angles.clear();
    errangles.clear();
    ROOT::Math::XYZVector dummy1(0, 0, 0);  // not needed here
    ROOT::Math::XYZVector dummy2(0, 0, 0);  // but later inside angle functions
    double sum = 0;
    for (PathPoint entry : *fdata) {
      double d = bldistance(entry, par, dummy1, dummy2);
      sum += d;  // squared weighted distance
    }
    // broken line angle contribution
    angles = calculate_beta(par);
    errangles = calculate_angle_errors(par);
    for (unsigned int j = 0; j < angles.size(); j++) {
      if (errangles.at(j) > 0.0) sum += angles.at(j) * angles.at(j) / errangles.at(j);
    }
    return sum;
  }
  std::vector<double> get_angles(const double* par) { return calculate_beta(par); }
  std::vector<double> get_errors(const double* par) { return calculate_angle_errors(par); }
};

// Usage as helper class: Create complex Graph nodes or edges
// in a countable container independent of this graph class. Once finished
// with a static container content of those complex vertex objects,
// all it takes is to use the unique element index in that container
// to initate this graph and use it entirely with integers throughout.
// Weights can be doubles representing space distances between nodes.
// Any result from here can then address complex objects in their
// container by index, again outside this class.
class WeightedGraph {
  // undirected, weighted graph object using integers as Node identifiers
  // for simple usage as a pure helper class.
 private:
  int V;  // number of vertices
  // Our adjacency list.
  std::vector<std::vector<std::pair<int, double> > > adjList;
  std::vector<std::pair<double, int> > dist;  // First is dist, second is the previous node.
  std::vector<int> sPath;                     // result path

 public:
  WeightedGraph(int nv);  // Constructor with number of vertices

  void addEdge(int v, int w, double weight);    // function to add an edge to graph
  bool isReachable(int s, int d);               // returns true if there is a path from s to d
  double dijkstraPaths(int start, int target);  // all paths between s and t
  std::vector<int> path() { return sPath; }
};

class Interval {
  // simple interval mostly for convenience
  // on checking pairs of doubles
  // on rings for overlap
 private:
  double lower;
  double upper;

 protected:
  static int checkforPiHalf(double& angle);  // on ring beyond pi/2 angle, get sign problems

 public:
  Interval();                    // Default Constructor
  Interval(double s, double e);  // Constructor with lower and upper limit

  Interval hull(Interval other) const;  // return interval containing both this and other
  double midinterval() { return 0.5 * (lower + upper); }  // mean interval value
  double angle_midinterval() const;                       // mean angular interval value
  double angle_dphi() const;                              // half width with angles
  double from() { return lower; }                         // boundary return
  double to() { return upper; }                           // boundary return
  bool empty() { return lower == upper; }                 // check for empty interval
  bool overlap(Interval other) const;                     // return true if overlap exists
  bool angle_overlap(Interval other) const;               // for intervals on a circle
};

class RelationalHit {
  // A TrackerHit in relation to another TrackerHit
  // where the relevant relation is expressed in
  // tangent point intervals, i.e. the error interval
  // around each tangent point. Such error intervals
  // can overlap and hence merge.
 private:
  TrackerHit host;
  TrackerHit guest;
  std::vector<std::pair<ROOT::Math::XYVector, ROOT::Math::XYVector> >
      onthis;  // node position and error interval
  std::vector<std::pair<ROOT::Math::XYVector, ROOT::Math::XYVector> > onother;  // on the two ring
  std::vector<std::pair<int, int> > hostnode_map;   // node indices mapping from, to
  std::vector<std::pair<int, int> > guestnode_map;  // node indices mapping from, to
  struct InternalStore {
    double xa;
    double ya;
    double xb;
    double yb;
    double rsmall;
    double ds;
    double rlarge;
    double dl;
    int order;
  } intern;

  std::pair<ROOT::Math::XYVector, ROOT::Math::XYVector> TP_pair;  // point pair
  std::pair<Interval, Interval> TP_error;                         // angle error interval pair
  std::vector<ROOT::Math::XYVector> allhost_points;               // point collection
  std::vector<ROOT::Math::XYVector> allguest_points;              // point collection
  std::vector<Interval> allhost_errors;                           // error in x=radial, y=dphi
  std::vector<Interval> allguest_errors;                          // error in x=radial, y=dphi
  std::vector<PathPoint> hostnodes;   // handing back 3D point with errors
  std::vector<PathPoint> guestnodes;  // handing back 3D point with errors

  // case handling
  void calculate_overlapping(double distance);
  void calculate_four(double distance);
  void calculate_overlapping_errors(double distance);
  void calculate_four_errors(double distance);
  void equal_rings(double distance, bool oflag);
  void intersectioncalc(double distance);
  // initiate case handling
  void size_sorting();  // sets internal store
  // general tangent point calculators
  void pointcalc(double xw, double yw, double rad, double h, double xo, double yo);
  void calc_onvertical(double xw, double yw, double rad, double h, double xo, double yo);
  void calc_onhorizontal(double xw, double yw, double rad, double h, double xo, double yo);
  // helpers
  static int checkforpihalf(double& angle);
  static double checkforpi(double angle);
  static double check_subtraction(double r, double dr);
  static std::vector<Interval> deltaxy_to_deltaphi(
      const std::vector<ROOT::Math::XYVector>& xy, std::vector<ROOT::Math::XYVector> dxdy,
      std::vector<ROOT::Math::XYVector> ctr);  // conversion to error on ring
  static ROOT::Math::XYVector rphitoxy(const ROOT::Math::XYVector& tp,
                                       const ROOT::Math::XYVector& ctr, double r, double dr,
                                       double dphi);

 protected:
  void calculate_tangentpoints();
  void create_maps();

 public:
  RelationalHit() { ; }                             // default constructor
  RelationalHit(TrackerHit hit, TrackerHit other);  // Constructor

  MetaInfo get_hostid() { return host.mi; }    // identify host tracker hit with unique meta info
  MetaInfo get_guestid() { return guest.mi; }  // identify guest tracker hit with unique meta info
  std::vector<PathPoint> nodes_host() { return hostnodes; }
  std::vector<PathPoint> nodes_guest() { return guestnodes; }

  std::vector<std::pair<int, int> > get_edges();  // countable, for the path graph eventually
};

class PathFinder {
  // Input all the TrackerHits in the event
  // creates a RelationalHit collection and sorts through them
  // in order to make the full list of all available edges
  // which then give a graph. Shortest paths result from that
  // graph and that results in data made up of points in 3D
  // with errors that can be fit, e.g. with the broken line fitter.
 private:
  int width;
  int height;
  std::vector<TrackerHit> hits;
  std::vector<PathPointCollection> paths;
  std::vector<std::pair<int, int> > edges;  // node connections with indices
  std::vector<PathPoint> ordered_points;    // node order for reverse translating
  std::vector<std::pair<PathPoint, PathPoint> > related_points;  // edge data points
  std::vector<RelationalHit> allpairs;
  // internal functions
  static bool overlap_x(PathPoint p1, PathPoint p2);
  static bool overlap_y(PathPoint p1, PathPoint p2);
  double edgelength(int e1, int e2);
  static std::vector<int> find_allint(const std::vector<PathPoint>& myvector, int pint);
  static size_t make_ahash(PathPoint pp);
  std::vector<int> column_hits(int col);

 protected:
  void make_edges();
  void find_paths();
  void clean_pointpairs(const std::vector<std::pair<PathPoint, PathPoint> >& im);
  static bool has_overlap(PathPoint p1, PathPoint p2);
  static bool is_neighbour(TrackerHit start, TrackerHit target);

 public:
  PathFinder(std::vector<TrackerHit> th);  // Constructor

  void create_paths();  // initialization for all calculations
  std::vector<PathPointCollection> allpaths() { return paths; }
};

// General fitter class with line and helix fitting methods
class SNFitter {
 private:
  std::vector<TrackerHit> rings;
  std::vector<GeigerRing> grings;

 protected:
  std::vector<double> line_initials(double frad);
  std::vector<double> helix_initials();
  std::vector<double> helixbackup();
  static std::vector<int> kink_finder(std::vector<double> betaangles,
                                      std::vector<double> errangles);
  static LineFit fitline2D(const std::vector<PathPoint>& data);
  static bool peak_alarm(std::vector<double> betaangles);
  static double martingale(double previous, double val, double alpha);

 public:
  SNFitter() {
    rings.clear();
    grings.clear();
  }  // default Constructor
  SNFitter(std::vector<TrackerHit> th) : rings(th) {
    for (auto& hit : rings) grings.push_back(hit.gr);
  }  // Constructor
  ~SNFitter() { rings.clear(); }

  void setData(std::vector<TrackerHit> th) {
    rings.clear();
    grings.clear();
    rings = th;
    for (auto& hit : rings) grings.push_back(hit.gr);
  }
  std::vector<LineFit> fitline();  // sets the data and operates
  std::vector<HelixFit> fithelix();
  std::vector<BrokenLineFit> fitbrokenline();
};

#endif
