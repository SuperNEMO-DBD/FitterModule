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

// function Object to be minimized
/* struct LineDistance2 { */
/*   // data member */
/*   SNFitter* snf; */

/*   LineDistance2(SNFitter * snfitter) : snf(snfitter) {} */
  
/*   // implementation of the function to be minimized */
/*   double operator() (const double * par) { */
/*     double sum = 0; */
/*     for (GeigerRing entry : *snf->allRings()) { */
/*       double d = snf->linedistance(entry, par); */
/*       sum += d*d; // squared weighted distance */
/*     } */
/*     return sum; */
/*   } */
  
/* }; */



// General fitter class with line and helix fitting methods
class SNFitter {

 private:
  std::vector<GeigerRing> rings;

  struct LineDistance2 {
    // implementation of the function to be minimized
    double operator() (const double * par) {
      double sum = 0;
      for (GeigerRing entry : rings) {
	double d = linedistance(entry, par);
	sum += d*d; // squared weighted distance
      }
      return sum;
    }
  };
  
  
 protected:
  std::vector<double> line_initials(double frad);
  std::vector<GeigerRing> allRings() {return rings;}
  double linedistance(GeigerRing gr, const double *p); // calculate distance line-cylinder
  //  double helixdistance(GeigerRing gr, const double *p);
  
  
 public:
  SNFitter() {} // Constructor
  ~SNFitter() {
    rings.clear();
  }
  
  std::vector<LineFit> fitline(std::vector<GeigerRing> gr); // sets the data and operates
  //  std::vector<HelixFit> fithelix(std::vector<GeigerRing> rings);
};


#endif
