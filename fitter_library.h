#ifndef YR_SNFitter
#define YR_SNFitter

// std libraries
#include <vector>


// tracker data storage
struct GeigerRing {
  double radius;
  double wirex;
  double wirey;
  double zcoord;
  double rerr;
  double zerr;
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
};


// function Object to be minimized
class LineDistance2 {
  // data member
 private:
  std::vector<GeigerRing>* frings;
  
 protected: 
  double linedistance(GeigerRing gr, const double *p); // calculate distance line-cylinder
  
 public:
 LineDistance2(std::vector<GeigerRing> * g) : frings(g) { }
  ~Linedistance2() { }
  
  // implementation of the function to be minimized
  double operator() (const double * par) {
    double sum = 0.0;
    for (GeigerRing entry : *frings) {
      double d = linedistance(entry, par);
      sum += d*d; // squared weighted distance
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
  double helixdistance(GeigerRing gr, const double *p); // calculate distance helix-cylinder
  
 public:
 HelixDistance2(std::vector<GeigerRing> * g) : frings(g) {}
  ~Helixdistance2() { }
  
  // implementation of the function to be minimized
  double operator() (const double * par) {
    double sum = 0.0;
    for (GeigerRing entry : *frings) {
      double d = helixdistance(entry, par);
      sum += d*d; // squared weighted distance
    }
    return sum;
  }
  
};



// General fitter class with line and helix fitting methods
class SNFitter {

 private:
  std::vector<GeigerRing> rings;
  
 protected:
  std::vector<double> line_initials(double frad);
  std::vector<double> helix_initials();
  std::vector<double> helixbackup();

 public:
 SNFitter(std::vector<GeigerRing> g) : rings(g) {} // Constructor
  ~SNFitter() {
    rings.clear();
  }
  
  std::vector<LineFit> fitline(); // sets the data and operates
  std::vector<HelixFit> fithelix();
};


#endif
