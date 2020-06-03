// us
#include <fitter_library.h>

// ROOT includes
#include <Fit/Fitter.h>
#include <Math/Functor.h>
#include <TSpectrum.h>
#include "TH1D.h"
#include "TKDTree.h"
#include "TLinearFitter.h"

// std
#include <algorithm>
#include <functional>
#include <iostream>
#include <list>
#include <queue>
#include <set>
#include <valarray>

// *** Line section
// *****************

// calculate distance line-cylinder
double LineDistance2::linedistance(GeigerRing gr, const double* p) {
  // distance line cylinder is D= | (xp-x0) cross  ux | - r
  // where ux is direction of line and x0 is a point on the line (like t = 0)
  // line not parallel to y-z plane, i.e calo plane
  double radius = gr.radius;
  double x = gr.wirex;
  double y = gr.wirey;
  double z = gr.zcoord;
  ROOT::Math::XYZVector weight;

  ROOT::Math::XYZVector xp(x, y, z);
  ROOT::Math::XYZVector x0(0., p[0], p[2]);
  ROOT::Math::XYZVector x1(1., p[0] + p[1], p[2] + p[3]);
  ROOT::Math::XYZVector u = (x1 - x0).Unit();
  double d2 = ((xp - x0).Cross(u)).Mag2();
  weight.SetXYZ((gr.rerr * u.Z() * gr.rerr * u.Z()) + (gr.zerr * u.Y() * gr.zerr * u.Y()),
                (gr.rerr * u.Z() * gr.rerr * u.Z()) + (gr.zerr * u.X() * gr.zerr * u.X()),
                (gr.rerr * u.Y() * gr.rerr * u.Y()) +
                    (gr.rerr * u.X() * gr.rerr * u.X()));  // weight sqr vector
  double w2 = weight.X() + weight.Y() + weight.Z();        // sum of squares is .Mag2()
  return (TMath::Sqrt(d2) - radius) / TMath::Sqrt(w2);     // only this is ring specific
}

// *** Helix section
// *****************

// calculate distance helix-ring
double HelixDistance2::helixdistance(GeigerRing gr, const double* p) {
  double ringr = gr.radius;
  double ringx = gr.wirex;
  double ringy = gr.wirey;
  double ringz = gr.zcoord;

  // weighting squared
  double weight = TMath::Sqrt(gr.rerr * gr.rerr + gr.zerr * gr.zerr);  // weight sqrt

  // p[0]=r, p[1]=h, p[2,3,4]=ref point
  ROOT::Math::XYZVector xp(ringx, ringy, ringz);  // ring centre
  ROOT::Math::XYZVector xref(p[2], p[3], p[4]);   // helix centre point
  ROOT::Math::XYZVector ray = xref - xp;          // connection vector

  double dir = TMath::ATan2(ray.y(), ray.x()) + TMath::Pi();
  ROOT::Math::XYZVector xhel(p[2] + p[0] * TMath::Cos(dir), p[3] + p[0] * TMath::Sin(dir),
                             p[4] + p[1] * dir);  // point on helix
  double d2 = (xhel - xp).Mag2();                 // distance squared

  return (TMath::Sqrt(d2) - ringr) / weight;  // only this is ring specific;
}

// *** SNFitter section
// ********************
// lines in SNFitter
std::vector<double> SNFitter::line_initials(double frad) {
  // provide options for x-y plane initials. X-Z plane not needed
  std::vector<double> results;

  TLinearFitter* lfxy = new TLinearFitter();
  lfxy->SetFormula("pol1");

  unsigned int nd = grings.size();
  double* xx = new double[nd];
  double* yy = new double[nd];
  double* ey = new double[nd];
  for (unsigned int j = 0; j < nd; j++) {
    xx[j] = grings.at(j).wirex;
    yy[j] = grings.at(j).wirey;
    ey[j] = grings.at(j).rerr;
  }
  lfxy->AssignData((int)nd, 1, xx, yy, ey);
  lfxy->Eval();
  // each attempt in dublets for the 2 parameter
  // part 1
  results.push_back(lfxy->GetParameter(0) - frad);  // - shift intercept in y by radius
  double slope = lfxy->GetParameter(1);
  double slopeaddon = 0.01;
  results.push_back(slope + slopeaddon);
  // part 2
  results.push_back(lfxy->GetParameter(0) - frad);  // - shift intercept in y
  results.push_back(slope - slopeaddon);            // minus add-on slope
  // part 3
  results.push_back(lfxy->GetParameter(0) + frad);  // + shift intercept in y
  results.push_back(slope + slopeaddon);            // same slope as before
  // part 4
  results.push_back(lfxy->GetParameter(0) + frad);  // + shift intercept in y
  results.push_back(slope - slopeaddon);

  delete[] xx;
  delete[] yy;
  delete[] ey;
  delete lfxy;

  return results;
}

std::vector<LineFit> SNFitter::fitline() {
  std::vector<LineFit> results;
  LineFit lf;
  ROOT::Fit::Fitter fitter;
  fitter.Config().SetDefaultMinimizer("Minuit2");

  // make the functor object
  LineDistance2 ldist(&grings);
  ROOT::Math::Functor fcn(ldist, 4);

  // first initials
  std::vector<double> dummy;
  std::vector<double>::iterator mit;
  double leftright = grings.at(0).wirex;
  int id = rings.at(0).clid;
  for (GeigerRing gg : grings) dummy.push_back(gg.wirex);
  if (leftright > 0)
    mit = std::min_element(dummy.begin(), dummy.end());
  else
    mit = std::max_element(dummy.begin(), dummy.end());

  int idx = mit - dummy.begin();
  double first_radius = grings.at(idx).radius;

  // initials for x-y plane
  std::vector<double> res = line_initials(first_radius);
  // for x-z plane
  unsigned int nd = grings.size();
  double* xx = new double[nd];
  double* zz = new double[nd];
  double* ez = new double[nd];
  for (unsigned int j = 0; j < nd; j++) {
    xx[j] = grings.at(j).wirex;
    zz[j] = grings.at(j).zcoord;
    ez[j] = grings.at(j).zerr;
  }
  // fit in x-z plane, once only
  TLinearFitter* lfxz = new TLinearFitter();
  lfxz->SetFormula("pol1");
  lfxz->AssignData((int)nd, 1, xx, zz, ez);
  lfxz->Eval();
  //  std::cout << "xz: i " << lfxz->GetParameter(0) << ", sl " << lfxz->GetParameter(1) <<
  //  std::endl;

  // full fit attempts
  double iniline[4];
  for (int i = 0; i < 8; i += 2) {  // 4 candidates for initials in x-y
    //    std::cout << "xy: i " << res[i] << ", sl " << res[i+1] << std::endl;
    iniline[0] = res[i];
    iniline[1] = res[i + 1];
    iniline[2] = lfxz->GetParameter(0);
    iniline[3] = lfxz->GetParameter(1);

    // fit
    fitter.SetFCN(fcn, iniline, (unsigned int)grings.size(), true);

    bool ok = fitter.FitFCN();
    const ROOT::Fit::FitResult& result = fitter.Result();
    lf.status = result.Status();  // fit to ring in 2d status
    lf.chi2 = result.Chi2();
    lf.prob = result.Prob();
    lf.clid = id;
    //    result.Print(std::cout);

    const double* bestfit = result.GetParams();
    const double* errors = result.GetErrors();
    lf.ixy = bestfit[0];
    lf.errixy = errors[0];
    lf.slxy = bestfit[1];
    lf.errslxy = errors[1];
    lf.ixz = bestfit[2];
    lf.errixz = errors[2];
    lf.slxz = bestfit[3];
    lf.errslxz = errors[3];

    results.push_back(lf);
  }
  delete[] ez;
  delete[] zz;
  delete[] xx;
  delete lfxz;
  return results;
}

std::vector<double> SNFitter::helix_initials() {
  std::vector<double> results;

  unsigned int nd = grings.size();
  double* xx = new double[nd];
  double* yy = new double[nd];
  double* zz = new double[nd];
  for (unsigned int j = 0; j < nd; j++) {
    xx[j] = grings.at(j).wirex;
    yy[j] = grings.at(j).wirey;
    zz[j] = grings.at(j).zcoord;
  }

  std::valarray<double> valx(xx, nd);
  std::valarray<double> valy(yy, nd);

  double meanx = valx.sum() / nd;
  double meany = valy.sum() / nd;
  // transform to u,v coordinates
  std::valarray<double> u = valx - meanx;
  std::valarray<double> v = valy - meany;

  // helper sums
  std::valarray<double> u2 = u * u;
  std::valarray<double> v2 = v * v;
  std::valarray<double> uv = u * v;
  double suu = u2.sum();
  double svv = v2.sum();
  double suv = uv.sum();
  std::valarray<double> u3 = u2 * u;
  std::valarray<double> v3 = v2 * v;
  std::valarray<double> uv2 = v2 * u;
  std::valarray<double> u2v = v * u2;
  double suuu = u3.sum();
  double svvv = v3.sum();
  double suvv = uv2.sum();
  double suuv = u2v.sum();

  // solve for circle centre
  double c1 = 0.5 * (suuu + suvv);
  double c2 = 0.5 * (svvv + suuv);
  double denominator = (suv * suv - suu * svv);  // prevent div by zero

  if (fabs(denominator) < 1.0e-16)  // not a  helix but a line in wires
    return results;                 // return empty as signal

  double uc = (c2 * suv - c1 * svv) / denominator;
  double vc = (c1 * suv - c2 * suu) / denominator;
  double rad = uc * uc + vc * vc + (suu + svv) / nd;  // squared radius
  results.push_back(TMath::Sqrt(rad));                // helix parameter 0

  TLinearFitter* lfxz = new TLinearFitter();
  lfxz->SetFormula("pol1");  // straight line x-z plane
  lfxz->AssignData((int)nd, 1, xx, zz);
  lfxz->Eval();
  if (lfxz->GetParameter(1) > 0)  // up or down pitch
    results.push_back(1.0);       // right sign
  else
    results.push_back(-1.0);  // right sign

  results.push_back(uc + meanx);  // centre x
  results.push_back(vc + meany);  // centre y

  std::vector<double> dummy;
  for (GeigerRing gg : grings) dummy.push_back(gg.wirex);
  std::vector<double>::iterator maxit = std::max_element(dummy.begin(), dummy.end());
  std::vector<double>::iterator minit = std::min_element(dummy.begin(), dummy.end());
  int maxidx = maxit - dummy.begin();
  int minidx = minit - dummy.begin();
  results.push_back(0.5 * (grings.at(maxidx).zcoord + grings.at(minidx).zcoord));  // z centre

  delete[] zz;
  delete[] yy;
  delete[] xx;
  delete lfxz;

  return results;
}

std::vector<double> SNFitter::helixbackup() {
  std::vector<double> results;
  std::vector<double> dummy;
  for (GeigerRing gg : grings) dummy.push_back(gg.wirex);
  std::vector<double>::iterator maxit = std::max_element(dummy.begin(), dummy.end());
  std::vector<double>::iterator minit = std::min_element(dummy.begin(), dummy.end());
  int maxidx = maxit - dummy.begin();
  int minidx = minit - dummy.begin();
  double zmean = 0.5 * (grings.at(maxidx).zcoord + grings.at(minidx).zcoord);
  double ymean = 0.5 * (grings.at(maxidx).wirey + grings.at(minidx).wirey);
  // some guess dummy values
  results.push_back(100.0);  // r

  if (zmean >= 0)            // up or down pitch
    results.push_back(1.0);  // h
  else
    results.push_back(-1.0);

  results.push_back(0.0);  // xc
  if (ymean >= 0)
    results.push_back(100.0);  // yc
  else
    results.push_back(-100.0);  // yc
  results.push_back(zmean);     // z centre
  return results;
}

std::vector<HelixFit> SNFitter::fithelix() {
  std::vector<HelixFit> results;
  HelixFit hf;
  ROOT::Fit::Fitter fitter;

  // make the functor object
  HelixDistance2 hdist(&grings);
  ROOT::Math::Functor fcn(hdist, 5);

  int id = rings.at(0).clid;
  std::vector<double> res = helix_initials();
  if (res.size() < 1) {  // not a helix in wires, make some up
    res.clear();
    res = helixbackup();
  }

  double inihelix[5] = {res[0], res[1], res[2], res[3], res[4]};
  //  std::cout << "helix3Dfit initials: " << res[0] << ", " << res[1] << ", " << res[2] << ", " <<
  //  res[3] << ", " << res[4] << std::endl;

  fitter.SetFCN(fcn, inihelix, (unsigned int)grings.size(), true);

  // fit
  bool ok = fitter.FitFCN();
  if (!ok) {
    std::cout << "helix3Dfit: Helix3D Fit failed: 2nd attempt started." << std::endl;
    // 2nd possibility
    res = helixbackup();
    double inihelix2[5] = {res[0], res[1], res[2], res[3], res[4]};

    fitter.SetFCN(fcn, inihelix2, (unsigned int)grings.size(), true);

    // fit
    ok = fitter.FitFCN();

    if (!ok) {
      std::cout << "helix3Dfit: Helix3D Fit failed" << std::endl;
    }
  }
  const ROOT::Fit::FitResult& result = fitter.Result();
  hf.status = result.Status();
  hf.chi2 = result.Chi2();
  hf.prob = result.Prob();
  hf.clid = id;

  const double* bestfit = result.GetParams();
  const double* errors = result.GetErrors();
  hf.radius = bestfit[0];
  hf.raderr = errors[0];
  hf.pitch = bestfit[1];
  hf.errpitch = errors[1];
  hf.xc = bestfit[2];
  hf.errxc = errors[2];
  hf.yc = bestfit[3];
  hf.erryc = errors[3];
  hf.zc = bestfit[4];
  hf.errzc = errors[4];

  results.push_back(hf);
  return results;
}

double SNFitter::martingale(double previous, double val, double alpha) {
  return alpha * previous + (1.0 - alpha) * val;
}

// change detect filter
bool SNFitter::peak_alarm(std::vector<double> betaangles) {
  if (betaangles.size() < 5) return false;  // pointless search
  double alpha = 0.9;
  double threshold = 0.11;
  double prval = betaangles.front();  // first
  unsigned int len = betaangles.size();
  unsigned int count = 0;
  for (double angle : betaangles) {
    double val = martingale(prval, angle, alpha);
    if (fabs(val - prval) > threshold && count < len - 1)
      return true;  // exclude final value kinks - dont' work with TSpectrum
    prval = val;
    count++;
  }
  return false;
}

std::vector<int> SNFitter::kink_finder(std::vector<double> betaangles,
                                       std::vector<double> errangles) {
  // Attempts to locate kinks (usually just one) using the error list
  // of angles along fitted multiple scattering model. Peaks in those
  // errors along fitted model point at kink locations.
  std::vector<int> bins;
  unsigned int n = betaangles.size();
  TH1D* hist = new TH1D("hist", "title", n, 1, n);
  for (unsigned int i = 0; i < n; i++) {
    hist->SetBinContent(i + 1, fabs(betaangles.at(i)));
    hist->SetBinError(i + 1, TMath::Sqrt(errangles.at(i)));
  }
  TSpectrum* sp = new TSpectrum();
  int npeaks = sp->Search(hist, 1, "goff", 0.9);  // 90% threshold for second peak
  if (npeaks == 0)                                // try one more
    npeaks = sp->Search(hist, 2.0, "goff", 0.9);
  if (npeaks > 0) {
    double* peakpos = sp->GetPositionX();
    for (int j = 0; j < npeaks; j++) {
      int bin = hist->FindBin(peakpos[j]);
      bins.push_back(bin);
    }
  }
  delete sp;
  delete hist;
  return bins;
}

LineFit SNFitter::fitline2D(std::vector<PathPoint> data) {
  LineFit lf;
  // separate linear fits in xy and xz
  // good for initial conditions or separate broken line pieces
  TLinearFitter* lfxy = new TLinearFitter();
  TLinearFitter* lfxz = new TLinearFitter();

  int nd = (int)data.size();
  double* xx = new double[nd];
  double* yy = new double[nd];
  double* zz = new double[nd];
  double* ey = new double[nd];
  double* ez = new double[nd];
  int j = 0;
  for (PathPoint pp : data) {
    xx[j] = pp.xc;
    yy[j] = pp.yc;
    zz[j] = pp.zc;
    ey[j] = pp.erry;
    ez[j] = pp.errz;
    j++;
  }

  // fit in x-y plane
  lfxy->SetFormula("pol1");
  lfxy->AssignData(nd, 1, xx, yy, ey);
  lfxy->Eval();
  //  std::cout << "xy: i " << lfxy->GetParameter(0) << ", slxy " << lfxy->GetParameter(1) <<
  //  std::endl;

  // fit in x-z plane
  lfxz->SetFormula("pol1");
  lfxz->AssignData(nd, 1, xx, zz, ez);
  lfxz->Eval();
  //  std::cout << "xz: i " << lfxz->GetParameter(0) << ", slxz " << lfxz->GetParameter(1) <<
  //  std::endl;

  lf.ixy = lfxy->GetParameter(0);
  lf.slxy = lfxy->GetParameter(1);
  lf.ixz = lfxz->GetParameter(0);
  lf.slxz = lfxz->GetParameter(1);

  lf.errixy = lfxy->GetParError(0);  // four best fit errors
  lf.errslxy = lfxy->GetParError(1);
  lf.errixz = lfxz->GetParError(0);
  lf.errslxz = lfxz->GetParError(1);

  lf.status = -1;                  // not available
  lf.chi2 = lfxy->GetChisquare();  // better than nothing
  lf.prob = -1.0;                  // not available

  delete[] xx;
  delete[] yy;
  delete[] zz;
  delete[] ey;
  delete[] ez;
  delete lfxy;
  delete lfxz;

  return lf;
}

std::vector<BrokenLineFit> SNFitter::fitbrokenline() {
  std::vector<BrokenLineFit> results;
  BrokenLineFit blf;
  LineFit lf;
  LineFit backup;
  int id = rings.at(0).clid;

  PathFinder pf(rings);
  pf.create_paths();  // make all the tangent point data paths for fitting
  std::vector<PathPointCollection> pathcollection = pf.allpaths();
  if (pathcollection.empty()) return results;  // return empty container

  ROOT::Fit::Fitter fitter;
  fitter.Config().SetDefaultMinimizer("Minuit2");

  for (auto& ppc : pathcollection) {
    // make the functor object
    BLDistance2 bldist(&ppc.path);
    ROOT::Math::Functor fcn(bldist, 4);

    // store the data
    blf.path = ppc.path;
    blf.length = ppc.length;

    // initials from x-y, x-z plane fits
    backup = fitline2D(ppc.path);

    double iniline[4];
    iniline[0] = backup.ixy;
    iniline[1] = backup.slxy;
    iniline[2] = backup.ixz;
    iniline[3] = backup.slxz;
    //    cout << "initial lines: " << iniline[0] << ", " << iniline[1] << ", " << iniline[2] << ",
    //    " << iniline[3] << endl;

    // do the fit
    fitter.SetFCN(fcn, iniline, ppc.path.size(), true);

    bool ok = fitter.FitFCN();
    const ROOT::Fit::FitResult& result = fitter.Result();

    const double* bestfit = result.GetParams();  // bldist does not store angles or errors
    blf.angles = bldist.get_angles(bestfit);
    std::vector<double> errangles = bldist.get_errors(bestfit);

    //    cout << "n angles " << blf.angles.size() << " and n err angles " << errangles.size() <<
    //    endl;
    // check
    // for (double a : blf.angles)
    //   std::cout << "angle = " << a << std::endl;
    // for (double ea : errangles)
    //   std::cout << "angle error = " << ea << std::endl;

    blf.status = result.Status();  // fit to ring in 2d status
    blf.chi2 = result.Chi2();
    blf.prob = result.Prob();
    blf.clid = id;

    //    result.Print(std::cout);

    // check for break points for at least 5 or more angles
    if (peak_alarm(blf.angles)) {  // otherwise nothing
      blf.breakpoints = kink_finder(blf.angles, errangles);

      std::vector<PathPoint> foilside;
      std::vector<PathPoint> caloside;

      if (blf.breakpoints.empty()) {  //  no kinks, empty
        lf.ixy = bestfit[0];          // four best fit parameter
        lf.slxy = bestfit[1];
        lf.ixz = bestfit[2];
        lf.slxz = bestfit[3];

        const double* err = result.GetErrors();
        lf.errixy = err[0];  // four best fit errors
        lf.errslxy = err[1];
        lf.errixz = err[2];
        lf.errslxz = err[3];

        lf.status = blf.status;  // same as overall
        lf.chi2 = blf.chi2;
        lf.prob = blf.prob;
        blf.linefit1 = lf;  // no breaks, one line fit only
        lf.chi2 = -1.0;
        blf.linefit2 = lf;                      // second line fit not valid chi2 = -1
      } else if (blf.breakpoints.size() > 1) {  // more than one
        std::sort(blf.breakpoints.begin(),
                  blf.breakpoints.end());  // is vector<int>, ascending order
        int bin1 = blf.breakpoints.front();
        if (bin1 < 2)  // have at least 2 points to fit to
          bin1 += 2;
        for (unsigned int i = 0; i < bin1; i++) foilside.push_back(ppc.path.at(i));

        int bin2 = blf.breakpoints.back();
        if (bin2 < 2)  // have at least 2 points to fit to
          bin2 += 2;
        for (unsigned int i = ppc.path.size() - bin2; i < ppc.path.size(); i++)
          caloside.push_back(ppc.path.at(i));

        blf.linefit1 = fitline2D(foilside);
        blf.linefit2 = fitline2D(caloside);
      } else {  // exactly one break point
        int bin1 = blf.breakpoints.front();
        if (bin1 < 2)  // have at least 2 points to fit to
          bin1 += 2;
        //    	cout << "foil side rings: " << endl;
        for (unsigned int i = 0; i < bin1; i++) {
          foilside.push_back(ppc.path.at(i));
        }
        // cout << "calo side rings: " << endl;
        for (unsigned int i = bin1; i < ppc.path.size(); i++) {
          caloside.push_back(ppc.path.at(i));
        }
        blf.linefit1 = fitline2D(foilside);
        blf.linefit2 = fitline2D(caloside);
      }
    } else {
      //      cout << "** No peak alarm. **" << endl;
      LineFit lf;
      lf.ixy = bestfit[0];  // four best fit parameter
      lf.slxy = bestfit[1];
      lf.ixz = bestfit[2];
      lf.slxz = bestfit[3];

      const double* err = result.GetErrors();
      lf.errixy = err[0];  // four best fit errors
      lf.errslxy = err[1];
      lf.errixz = err[2];
      lf.errslxz = err[3];

      lf.status = blf.status;  // same as overall
      lf.chi2 = blf.chi2;
      lf.prob = blf.prob;
      blf.linefit1 = lf;  // no breaks, one line fit only
      lf.chi2 = -1.0;
      blf.linefit2 = lf;  // second line fit not valid chi2 = -1
    }

    // // store results
    results.push_back(blf);
  }
  // sort according to chi2 in ascending order
  std::sort(results.begin(), results.end(),
            [](BrokenLineFit& bl1, BrokenLineFit& bl2) { return (bl1.chi2 < bl2.chi2); });

  return results;
}

// Helper methods

// *****
// ** PathFinder methods
// *****
PathFinder::PathFinder(std::vector<TrackerHit> th) {
  width = 9;     // set tracker columns, not going to change anytime soon
  height = 113;  // tracker rows, needed for start/end point finding
  hits = th;     // copy
  // rest to be filled by methods
  paths.clear();
  allpairs.clear();
  edges.clear();
  related_points.clear();
  ordered_points.clear();
}

void PathFinder::create_paths() {
  make_edges();  // build all the RelationalHits and create edges from collection
  find_paths();
}

bool PathFinder::overlap_x(PathPoint p1, PathPoint p2) {
  // two error bar length overlap, maybe multiples required for safety
  int n = 3;
  if (p1.xc >= p2.xc) {                                // overlap in x check
    if (p1.xc - n * p1.errx <= p2.xc + n * p2.errx) {  // overlap error bars
      return true;
    }
  } else {
    if (p1.xc + n * p1.errx >= p2.xc - n * p2.errx) {  // overlap error bars
      return true;
    }
  }
  return false;
}

bool PathFinder::overlap_y(PathPoint p1, PathPoint p2) {
  int n = 3;
  if (p1.yc >= p2.yc) {                                // overlap in y check
    if (p1.yc - n * p1.erry <= p2.yc + n * p2.erry) {  // overlap error bars
      return true;
    }
  } else {
    if (p1.yc + n * p1.erry >= p2.yc - n * p2.erry) {  // overlap error bars
      return true;
    }
  }
  return false;
}

bool PathFinder::has_overlap(PathPoint p1, PathPoint p2) {
  if (p1.pointid != p2.pointid) {
    //    std::cout << "in overlap: x is " << overlap_x(p1, p2) << " and y is " << overlap_y(p1, p2)
    //    << std::endl;
    if (overlap_x(p1, p2) && overlap_y(p1, p2)) return true;
  }
  return false;
}

bool PathFinder::is_neighbour(TrackerHit start, TrackerHit target) {
  MetaInfo ppstart = start.mi;
  MetaInfo pptarget = target.mi;

  if (ppstart.column == pptarget.column &&
      ppstart.row == pptarget.row)  // exclude itself as neighbour
    return false;

  // check grid x-y
  if (fabs(ppstart.column - pptarget.column) < 2 &&
      fabs(ppstart.row - pptarget.row) < 2)  // one grid place only
    return true;
  return false;
}

size_t PathFinder::make_ahash(PathPoint pp) {
  return (std::hash<double>()(pp.xc) ^ std::hash<double>()(pp.yc));
}

std::vector<int> PathFinder::column_hits(int col) {
  std::vector<int> found;
  std::vector<int> rows;
  for (auto& hit : hits) rows.push_back(hit.mi.row);  // orders row integers, min to max
  std::sort(rows.begin(), rows.end());
  int minrow = rows.front();
  int maxrow = rows.back();
  //  std::cout << "max row = " << maxrow << ", min row = " << minrow << std::endl;
  for (auto& hit : hits) {
    MetaInfo trmi = hit.mi;
    int trhitid = trmi.hitid;
    if (trmi.column == col) {
      // std::cout << "found col = " << col << std::endl;
      // found.push_back(trhitid);
      // if (trmi.row == maxrow || trmi.row == minrow || trmi.row == maxrow-1 || trmi.row ==
      // minrow+1) {
      if (trmi.row == maxrow || trmi.row == minrow) {
        // std::cout << "found row = " << trmi.row << std::endl;
        found.push_back(trhitid);
      }
    }
  }

  std::vector<int> columnpoints;  // store ordered point indices with that hit id
  for (int hitid : found) {       // find all nodes with that hit id
    std::vector<int> temp = find_allint(ordered_points, hitid);
    columnpoints.insert(columnpoints.end(), temp.begin(),
                        temp.end());  // all of temp into columnpoints
    // std::cout << "id=" << hitid << " found times " << temp.size() << std::endl;
  }
  return columnpoints;
}

std::vector<int> PathFinder::find_allint(std::vector<PathPoint> myvector, int pint) {
  std::vector<int> temp;
  for (auto& pp : myvector) temp.push_back(pp.pointid.first);  // collect all int hit ids
  std::vector<int>::iterator it;
  std::vector<int> position;
  int pos;

  // count the values
  int mycount = (int)std::count(temp.begin(), temp.end(), pint);

  it = std::find(temp.begin(), temp.end(), pint);
  for (int i = 0; i < mycount; i++) {
    pos = it - temp.begin();
    position.push_back(pos);
    ++it;
    it = std::find(it, temp.end(), pint);
  }

  return position;
}

double PathFinder::edgelength(int e1, int e2) {
  ROOT::Math::XYZVector x0(ordered_points.at(e1).xc, ordered_points.at(e1).yc,
                           ordered_points.at(e1).zc);
  ROOT::Math::XYZVector x1(ordered_points.at(e2).xc, ordered_points.at(e2).yc,
                           ordered_points.at(e2).zc);
  return TMath::Sqrt((x1 - x0).Mag2());
}

void PathFinder::clean_pointpairs(std::vector<std::pair<PathPoint, PathPoint> > im) {
  // fill store of related_points, check for pathpoint overlaps and clean store
  // check
  // for (auto& rp : im) {
  //   std::cout << "path point from id (" << rp.first.pointid.first << "," <<
  //   rp.first.pointid.second << ") to id (" << rp.second.pointid.first << "," <<
  //   rp.second.pointid.second << ")" << std::endl; std::cout << "path point from tp " <<
  //   rp.first.xc << ", " << rp.first.yc << std::endl; std::cout << "path point to tp " <<
  //   rp.second.xc << ", " << rp.second.yc << std::endl;
  // }
  related_points.clear();
  std::vector<PathPoint> onhost;
  std::vector<PathPoint> onguest;
  for (auto& rp : im) {
    onhost.push_back(rp.first);  // keep order
    onguest.push_back(rp.second);
  }
  std::set<int> temp;
  for (auto& pp : onhost) temp.insert(pp.pointid.first);  // unique hit ids

  std::set<std::pair<int, int> > mergehost;
  std::set<std::pair<int, int> > mergeguest;
  for (int nn : temp) {
    std::vector<int> where = find_allint(onhost, nn);  // ring hit id search
    //    std::cout << " hit id " << nn << " found " << where.size() << " times in onhost." <<
    //    std::endl;
    for (unsigned int which = 0; which < where.size() - 1; which++) {
      for (unsigned int target = which + 1; target < where.size();
           target++) {  // check on all subsequent
        if (has_overlap(onhost.at(where.at(which)),
                        onhost.at(where.at(target)))) {  // test all points
          // std::cout << "merger on host = (" << onhost.at(where.at(which)).pointid.first << "," <<
          // onhost.at(where.at(which)).pointid.second << ") with ("
          //    << onhost.at(where.at(target)).pointid.first << "," <<
          //    onhost.at(where.at(target)).pointid.second << ")" << std::endl;
          mergehost.insert(
              std::make_pair(where.at(which), where.at(target)));  // store, don't modify
        }
      }
    }
  }
  temp.clear();
  for (auto& pp : onguest) temp.insert(pp.pointid.first);  // unique hit ids

  for (int nn : temp) {
    std::vector<int> where = find_allint(onguest, (int)nn);  // ring hit id search
    //    std::cout << " hit id " << nn << " found " << where.size() << " times in onguest." <<
    //    std::endl;
    for (unsigned int which = 0; which < where.size() - 1; which++) {
      for (unsigned int target = which + 1; target < where.size();
           target++) {  // check on all subsequent
        if (has_overlap(onguest.at(where.at(which)),
                        onguest.at(where.at(target)))) {  // test all points
          // std::cout << "merger on guest = (" << onguest.at(where.at(which)).pointid.first << ","
          // << onguest.at(where.at(which)).pointid.second << ") with ("
          //    << onguest.at(where.at(target)).pointid.first << "," <<
          //    onguest.at(where.at(target)).pointid.second << ")" << std::endl;
          mergeguest.insert(
              std::make_pair(where.at(which), where.at(target)));  // store, don't modify
        }
      }
    }
  }
  // now can modify intermediate container content
  for (auto& mpair : mergehost) {
    onhost.at(mpair.second) = onhost.at(mpair.first);  // target = which PathPoint
    //    std::cout << "h: overwrite nr " << mpair.second << " with " << mpair.first << std::endl;
  }
  for (auto& mpair : mergeguest) {
    onguest.at(mpair.second) = onguest.at(mpair.first);  // target = which PathPoint
    //    std::cout << "g: overwrite nr " << mpair.second << " with " << mpair.first << std::endl;
  }

  //  std::cout << "size im " << im.size() << " with on host " << onhost.size() << " and on guest "
  //  << onguest.size() << " sizes " << std::endl;
  std::set<std::pair<int, size_t> > node_cleaner;  // only unique entries
  for (unsigned int nn = 0; nn < im.size(); nn++) {
    related_points.push_back(
        std::make_pair(onhost.at(nn), onguest.at(nn)));  // stored all valid point pairs
    node_cleaner.insert(onhost.at(nn).pointid);          // insert all nodes after overwrite
    node_cleaner.insert(onguest.at(nn).pointid);         // without multiple entries
  }
  std::vector<std::pair<int, size_t> > nodes;  // convenience node container by unique ids
  for (auto& nd : node_cleaner) {
    nodes.push_back(nd);  // make nodes countable with index
    //    std::cout << " Got node : (" << nd.first << ", " << nd.second << ")" << std::endl;
  }

  ordered_points.resize(nodes.size());
  edges.clear();
  std::vector<std::pair<int, size_t> >::iterator itn;
  for (auto& rp : related_points) {  // store edges with indices
    itn = std::find(nodes.begin(), nodes.end(), rp.first.pointid);
    int pos1 = itn - nodes.begin();
    itn = std::find(nodes.begin(), nodes.end(), rp.second.pointid);
    int pos2 = itn - nodes.begin();
    edges.push_back(std::make_pair(pos1, pos2));  // edge integers map to ordered_point indices
    ordered_points.at(pos1) = rp.first;           // insert PathPoint in node order
    ordered_points.at(pos2) = rp.second;  // make sure to capture all incl. overwriting if required

    // std::cout << "edges from id (" << nodes.at(pos1).first << "," << nodes.at(pos1).second << ")
    // to id (" << nodes.at(pos2).first << "," << nodes.at(pos2).second << ")" << std::endl;
    // std::cout << " node index (" << pos1 << "," << pos2 << ")" << std::endl;
    // std::cout << "edges from tp " << rp.first.xc << ", " << rp.first.yc << std::endl;
    // std::cout << "edges to tp " << rp.second.xc << ", " << rp.second.yc << std::endl;
  }
  // check
  // for (auto& op : ordered_points)
  //   std::cout << "ordered points: " << op.pointid.first << ", " << op.xc << ", " << op.yc <<
  //   std::endl;
}

// data is in tracker hits
void PathFinder::make_edges() {
  // fill tree with doubles by construction
  int nentries = (int)hits.size();
  int nneighbours = (nentries >= 8) ? 8 : nentries;
  double* x = new double[nentries];  // kdTree needs arrays as input
  double* y = new double[nentries];  // kdTree needs arrays as input
  double* z = new double[nentries];  // kdTree needs arrays as input

  int counter = 0;
  for (auto& entry : hits) {               // as stored in hits, in order of hits
    MetaInfo pp = entry.mi;                // get metainfo
    x[counter] = pp.column * 44.0 + 22.0;  // turn to [mm] regular centrepoint at 22 mm
    y[counter] = pp.row * 44.0 + 22.0;     // turn to [mm] in units of cell 44 mm
    z[counter] = entry.gr.zcoord;          // is already [mm]
    counter++;
  }

  TKDTreeID* hittree = new TKDTreeID(nentries, 3, 1);
  // fill kdTree here
  hittree->SetData(0, x);
  hittree->SetData(1, y);
  hittree->SetData(2, z);
  hittree->Build();
  // kdTree ready to work, make edges for a graph using this

  // turn the kdTree container content into container of edges
  TrackerHit start;   // to hold the starter hit
  TrackerHit target;  // to hold the target hit

  double point[3];                         // needed for kdTree requests
  double* dist = new double[nneighbours];  // check on nearest 8 neighbours in grid
  int* indx = new int[nneighbours];        // where 8 is the maximum posible in a square grid

  allpairs.clear();
  for (counter = 0; counter < nentries; counter++) {
    start = hits.at(counter);
    point[0] = x[counter];
    point[1] = y[counter];
    point[2] = z[counter];
    hittree->FindNearestNeighbors(point, nneighbours, indx,
                                  dist);  // index needs to correspond to hit index
    for (int j = 0; j < nneighbours; j++) {
      target = hits.at(indx[j]);
      if (is_neighbour(start, target)) {
        RelationalHit rh(start, target);
        //	std::cout << "pair made from " << start.mi.hitid << " and " << target.mi.hitid <<
        //std::endl;
        allpairs.push_back(rh);  // RH has internal edges with own indices
      }
    }
  }

  //  std::cout << "all pairs made, done all relational hits, size = " << allpairs.size() <<
  //  std::endl;

  // turn all RH into all PathPoint pairs as intermediate for global id
  std::vector<std::pair<PathPoint, PathPoint> > intermediate;
  std::vector<PathPoint> hostpp;
  std::vector<PathPoint> guestpp;
  for (auto& rh : allpairs) {
    std::vector<std::pair<int, int> > iedges = rh.get_edges();
    hostpp = rh.nodes_host();  // only defined after get_edges call
    guestpp = rh.nodes_guest();
    for (auto& iedge : iedges) {  // unfold internal edges into global container
      //      std::cout << "edge: " << iedge.first << ", " << iedge.second << std::endl;
      PathPoint pph = hostpp.at(iedge.first);
      pph.pointid.second = make_ahash(pph);  // overwrite local second id for global unique id
      PathPoint ppg = guestpp.at(iedge.second);
      ppg.pointid.second = make_ahash(ppg);  // overwrite local second id for global unique id
      intermediate.push_back(std::make_pair(pph, ppg));
    }
  }
  //  std::cout << "all edges before cleaning, size = " << intermediate.size() << std::endl;
  clean_pointpairs(intermediate);  // set storage of edges of valid PathPoints

  delete[] indx;
  delete[] dist;
  delete[] x;
  delete[] y;
  delete[] z;
  delete hittree;
}

void PathFinder::find_paths() {
  // fill graph with ordered point indices as ints
  WeightedGraph gr((int)ordered_points.size());  // graph object

  for (auto& edge : edges) {  // fill with edge indices
    double weight = edgelength(edge.first, edge.second);
    gr.addEdge(edge.first, edge.second, weight);
  }
  // these int containers hold indices of ordered points container
  std::vector<int> tempstarts = column_hits(0);           // starter Tracker Hit tangent points
  std::vector<int> temptargets = column_hits(width - 1);  // target Tracker Hit tangent points

  // short or curving paths to consider
  if (temptargets.empty()) {
    int whichcol = width - 2;
    std::vector<int> nextcolumn = column_hits(whichcol);  // column width-2 for targets
    while (nextcolumn.empty() && whichcol > 0) {
      whichcol--;
      nextcolumn = column_hits(whichcol);  // column whichcol for targets
    }
    if (!nextcolumn.empty()) temptargets = nextcolumn;
  }
  if (tempstarts.empty()) {
    int whichcol = 1;
    std::vector<int> nextcolumn = column_hits(whichcol);  // column 1 for starts
    while (nextcolumn.empty() && whichcol < width - 1) {
      whichcol++;
      nextcolumn = column_hits(whichcol);  // column whichcol for starts
    }
    if (!nextcolumn.empty()) tempstarts = nextcolumn;
  }

  // permit only unique indices
  std::set<int> pathstarts;
  std::set<int> pathtargets;
  for (int idx : tempstarts) pathstarts.insert(idx);    // get unique starts
  for (int idx : temptargets) pathtargets.insert(idx);  // get unique targets

  // check
  // std::cout << "Starts ids:" << std::endl;
  // for (int idx : pathstarts)
  //   std::cout << idx << std::endl;
  // std::cout << "End ids:" << std::endl;
  // for (int idx : pathtargets)
  //   std::cout << idx << std::endl;
  //  for (int s : pathstarts) {
  //    for (int t : pathtargets) {
  //      if (gr.isReachable(s, t) && s != t)
  //      	std::cout << "Reach " << t << " from " << s << std::endl;
  //    }
  //  }

  PathPointCollection ppc;
  std::vector<PathPoint> path;
  paths.clear();
  for (int s : pathstarts) {
    for (int t : pathtargets) {
      if (gr.isReachable(s, t) && s != t) {
        //        std::cout << "Reach " << t << " from " << s << std::endl;
        ppc.length = gr.dijkstraPaths(s, t);
        // std::cout << "direct length return: " << ppc.length << std::endl;
        for (auto& pid : gr.path()) path.push_back(ordered_points.at(pid));  // fill with PathPoints
        // check
        // for (auto& op : path)
        //    std::cout << "ordered points on path: " << op.xc << ", " << op.yc << std::endl;

        ppc.path = path;
        paths.push_back(ppc);
        path.clear();
      }
      //   std::cout << "Path as ids: " << std::endl;
      //   for (auto& pid : pp)
      //     std::cout << "id=" << pid << endl;
    }
  }
  std::sort(paths.begin(), paths.end(), [](PathPointCollection ppc1, PathPointCollection ppc2) {
    return (ppc1.length < ppc2.length);
  });  // sort according to path length

  // check
  // for (auto& ppc : paths) {
  //   std::cout << "Check length = " << ppc.length << std::endl;
  //   std::cout << "entries in path: " << ppc.path.size() << std::endl;
  // }

  // check pathpointcollection
  //  for (auto& ppc : paths) {
  //    std::cout << "Check path: " << std::endl;
  //    for (auto& pp : ppc.path)
  //      std::cout << "id=" << pp.pointid.first << std::endl;
  //  }
}

// *****
// ** Relational hit methods
// *****
RelationalHit::RelationalHit(TrackerHit hit, TrackerHit other) {
  // steer process to final edges.

  host = hit;     // focus on this geiger hit
  guest = other;  // this geiger hit is the neighbour

  // all other data members empty - to be filled by methods.
  calculate_tangentpoints();  // first process
  create_maps();
}

void RelationalHit::create_maps() {
  onthis.clear();
  onthis.resize(allhost_points.size());  // can slot in tp in place
  onother.clear();
  onother.resize(allguest_points.size());
  std::vector<int> booked;
  std::vector<int> merged;
  std::vector<int>::iterator it;
  ROOT::Math::XYVector hostanode(host.gr.wirex, host.gr.wirey);
  ROOT::Math::XYVector guestanode(guest.gr.wirex, guest.gr.wirey);
  bool block = false;

  // nodes on the focus ring
  for (unsigned int i = 0; i < allhost_errors.size() && !block; i++) {
    it = std::find(merged.begin(), merged.end(), i);
    if (it == merged.end()) {
      Interval ival = allhost_errors.at(i);
      //      std::cout << "create map: host ival[" << ival.from() << ", " << ival.to() << "]" <<
      //      std::endl;
      for (unsigned int j = i + 1; j < allhost_errors.size(); j++) {  //
        Interval trial = allhost_errors.at(j);
        //	std::cout << "create map: host trial[" << trial.from() << ", " << trial.to() << "]"
        //<< std::endl;
        if (ival.angle_overlap(trial)) {
          merged.push_back(j);                           // not using this at later iteration
          booked.push_back(i);                           // both i,j are replaced due to overlap
          hostnode_map.push_back(std::make_pair(i, i));  // i guest tp for both, i,j
          hostnode_map.push_back(std::make_pair(j, i));  // tangent points on host
          // std::cout << "create map: overlap edge (" << i << ", " << i << ") inserted." <<
          // std::endl; std::cout << "create map: and edge (" << j << ", " << i << ") inserted." <<
          // std::endl;
          Interval combined = ival.hull(trial);  // merge tangent points and errors
          //	  std::cout << "Hull interval: [" << combined.from() << ", " << combined.to() << "]"
          //<< std::endl;
          // check for pihalf
          int signx = 1;
          double start = combined.from();
          double end = combined.to();
          int sign1 = checkforpihalf(start);  // changes combined.lower/upper
          int sign2 = checkforpihalf(end);    // if necessary.
          //	  std::cout << "hull bounds after pi half check: [" << start << ", " << end << "]"
          //<< std::endl;
          if (sign1 < 0 && sign2 < 0) {
            combined = Interval(start, end);
            //	    std::cout << "New Hull interval: [" << combined.from() << ", " << combined.to()
            //<< "]" << std::endl;
            signx = -1;
          }
          double combiangle = combined.angle_midinterval();  // into one node
          //	  std::cout << " mid interval angle: " << combiangle << std::endl;
          double rad = host.gr.radius;
          ROOT::Math::XYVector relativetp(signx * rad * TMath::Cos(combiangle),
                                          signx * rad * TMath::Sin(combiangle));  // origin 0.0
          ROOT::Math::XYVector tp = relativetp + hostanode;                       // combined TP
          // now the error
          // capture extremely large arc for tiny rings
          if (combined.to() - combiangle > TMath::Pi() / 4.0) {
            //	    std::cout << " host - extreme arc case " << std::endl;
            ROOT::Math::XYVector dxdy(rad, rad);  // dx, dy as size of ring
            for (unsigned int n = 0; n < onthis.size(); n++)
              onthis.at(n) = std::make_pair(hostanode, dxdy);  // anode position for all nodes
            block = true;
            break;  // done, out
          } else {
            ROOT::Math::XYVector dxdy =
                rphitoxy(tp, ROOT::Math::XYVector(host.gr.wirex, host.gr.wirey), rad, host.gr.rerr,
                         combined.angle_dphi());      // returns dx, dy in vector
            onthis.at(i) = std::make_pair(tp, dxdy);  // same according to map
            onthis.at(j) = std::make_pair(tp, dxdy);
            //	    std::cout << "host - overlap edge (" << j << ", " << i << ") inserted." <<
            //std::endl;
          }
        }
      }
      it = std::find(booked.begin(), booked.end(), i);
      if (it == booked.end()) {  // store only once
        ROOT::Math::XYVector dxdy =
            rphitoxy(allhost_points.at(i), ROOT::Math::XYVector(host.gr.wirex, host.gr.wirey),
                     host.gr.radius, host.gr.rerr, ival.angle_dphi());  // returns dx, dy in vector
        onthis.at(i) = std::make_pair(allhost_points.at(i), dxdy);
        booked.push_back(i);
        hostnode_map.push_back(std::make_pair(i, i));  // host i has partner i as in allhost_pairs
        //	std::cout << "host - create map: no overlap edge (" << i << ", " << i << ")
        //inserted." << std::endl;
      }
    }
  }

  // nodes on the guest ring
  booked.clear();
  merged.clear();  // reset
  block = false;
  for (unsigned int i = 0; i < allguest_errors.size() && !block; i++) {
    it = std::find(merged.begin(), merged.end(), i);
    if (it == merged.end()) {
      Interval ival = allguest_errors.at(i);  // check on only the guest ring interval
      //      std::cout << "create map: guest ival[" << ival.from() << ", " << ival.to() << "]" <<
      //      std::endl;
      for (unsigned int j = i + 1; j < allguest_errors.size(); j++) {  //
        Interval trial = allguest_errors.at(j);
        //	std::cout << "create map: guest trial[" << trial.from() << ", " << trial.to() << "]"
        //<< std::endl;
        if (ival.angle_overlap(trial)) {
          merged.push_back(j);                            // not using this at later iteration
          booked.push_back(i);                            // both i,j are replaced due to overlap
          guestnode_map.push_back(std::make_pair(i, i));  // i host tp for both, i,j
          guestnode_map.push_back(std::make_pair(j, i));  // tangent points on guest
          // std::cout << "create map: overlap edge (" << i << ", " << i << ") inserted." <<
          // std::endl; std::cout << "create map: and edge (" << j << ", " << i << ") inserted." <<
          // std::endl;
          Interval combined = ival.hull(trial);  // merge tangent points and errors
          //	  std::cout << "guest Hull interval: [" << combined.from() << ", " << combined.to()
          //<< "]" << std::endl;
          // check for pihalf
          int signx = 1;
          double start = combined.from();
          double end = combined.to();
          int sign1 = checkforpihalf(start);  // changes start/end values
          int sign2 = checkforpihalf(end);    // if necessary.
          if (sign1 < 0 && sign2 < 0) {
            combined = Interval(start, end);
            //	    std::cout << "guest New Hull interval: [" << combined.from() << ", " <<
            //combined.to() << "]" << std::endl;
            signx = -1;
          }
          double combiangle = combined.angle_midinterval();  // into one node
          //	  std::cout << " mid interval angle: " << combiangle << std::endl;
          double rad = guest.gr.radius;
          ROOT::Math::XYVector relativetp(signx * rad * TMath::Cos(combiangle),
                                          signx * rad * TMath::Sin(combiangle));  // origin 0.0
          ROOT::Math::XYVector tp = relativetp + guestanode;                      // combined TP
          // now the error
          // capture extremely large arc for tiny rings
          if (combined.to() - combiangle > TMath::Pi() / 4.0) {
            ROOT::Math::XYVector dxdy(rad, rad);  // dx, dy as size of ring
            for (unsigned int n = 0; n < onother.size(); n++)
              onother.at(n) = std::make_pair(guestanode, dxdy);  // anode position for all nodes
            block = true;
            break;  // done, out
          } else {
            ROOT::Math::XYVector dxdy =
                rphitoxy(tp, ROOT::Math::XYVector(guest.gr.wirex, guest.gr.wirey), rad,
                         guest.gr.rerr, combined.angle_dphi());  // returns dx, dy in vector
            onother.at(i) = std::make_pair(tp, dxdy);            // same according to map
            onother.at(j) = std::make_pair(tp, dxdy);
            //	    std::cout << "guest - overlap edge (" << j << ", " << i << ") inserted." <<
            //std::endl;
          }
        }
      }
      it = std::find(booked.begin(), booked.end(), i);
      if (it == booked.end()) {  // store only once
        ROOT::Math::XYVector dxdy = rphitoxy(
            allguest_points.at(i), ROOT::Math::XYVector(guest.gr.wirex, guest.gr.wirey),
            guest.gr.radius, guest.gr.rerr, ival.angle_dphi());  // returns dx, dy in vector
        onother.at(i) = std::make_pair(allguest_points.at(i), dxdy);
        booked.push_back(i);
        guestnode_map.push_back(std::make_pair(i, i));  // guest i has partner i as in allhost_pairs
        //	std::cout << "guest - create map: no overlap edge (" << i << ", " << i << ")
        //inserted." << std::endl;
      }
    }
  }
}

void RelationalHit::calculate_tangentpoints() {
  // step 1 in process: all tangent points and errors,
  // safely stored in data members.
  double xwire1 = host.gr.wirex;
  double ywire1 = host.gr.wirey;
  double xwire2 = guest.gr.wirex;
  double ywire2 = guest.gr.wirey;
  double distance =
      TMath::Sqrt((xwire1 - xwire2) * (xwire1 - xwire2) + (ywire1 - ywire2) * (ywire1 - ywire2));
  double ring_overlap = (host.gr.radius + guest.gr.radius) - distance;

  if (ring_overlap >= 0.0) {
    //    std::cout << "tangent points: overlapping case" << std::endl;
    calculate_overlapping(distance);
    calculate_overlapping_errors(distance);
  } else {
    //    std::cout << "tangent points: not overlapping case" << std::endl;
    calculate_four(distance);
    calculate_four_errors(distance);
  }
  return;
}

void RelationalHit::calculate_overlapping(double distance) {
  if (fabs(host.gr.radius - guest.gr.radius) < 0.01) {  // capture equal radius case
    equal_rings(distance, true);
    return;
  }

  size_sorting();  // initiate internal storage, must come before error calculation

  // on smaller ring first
  double h = intern.rsmall * distance / (intern.rlarge - intern.rsmall);
  double xo = (1.0 + h / distance) * intern.xa - h / distance * intern.xb;
  double yo = (1.0 + h / distance) * intern.ya - h / distance * intern.yb;
  pointcalc(intern.xa, intern.ya, intern.rsmall, h, xo, yo);  // pair of points done
  if (intern.order < 1) {                                     // host is host
    allhost_points.push_back(TP_pair.first);
    allhost_points.push_back(TP_pair.second);
  } else {
    allguest_points.push_back(TP_pair.first);
    allguest_points.push_back(TP_pair.second);
  }

  // on larger ring
  h = intern.rlarge * distance / (intern.rlarge - intern.rsmall);
  pointcalc(intern.xb, intern.yb, intern.rlarge, h, xo, yo);  // pair of points done
  if (intern.order < 1) {                                     // host is host
    allguest_points.push_back(TP_pair.first);
    allguest_points.push_back(TP_pair.second);
  } else {
    allhost_points.push_back(TP_pair.first);
    allhost_points.push_back(TP_pair.second);
  }

  // finally the intersection points
  intersectioncalc(distance);
  allhost_points.push_back(TP_pair.first);
  allhost_points.push_back(TP_pair.second);

  allguest_points.push_back(TP_pair.first);
  allguest_points.push_back(TP_pair.second);
}

void RelationalHit::calculate_overlapping_errors(double distance) {
  bool swapped = false;  // origin must be on the correct side.
  std::vector<std::pair<ROOT::Math::XYVector, ROOT::Math::XYVector> > l_pairs;
  std::vector<std::pair<ROOT::Math::XYVector, ROOT::Math::XYVector> > s_pairs;
  double rstemp = intern.rsmall;
  double rltemp = intern.rlarge;
  if ((rstemp + intern.ds) > rltemp - intern.dl) {
    swapped = true;  // needs change
    intern.rsmall = check_subtraction(rltemp, intern.dl);
    intern.rlarge = rstemp + intern.ds;
    double xd = intern.xa;
    double yd = intern.ya;
    intern.xa = intern.xb;  // swap
    intern.xb = xd;
    intern.ya = intern.yb;
    intern.yb = yd;
  } else {
    intern.rsmall = rstemp + intern.ds;
    intern.rlarge = check_subtraction(rltemp, intern.dl);
  }
  // vary the ring radii, get tangent points on those -> errors in x,y
  double h = intern.rsmall * distance / (intern.rlarge - intern.rsmall);
  double xo = (1.0 + h / distance) * intern.xa - h / distance * intern.xb;
  double yo = (1.0 + h / distance) * intern.ya - h / distance * intern.yb;
  pointcalc(intern.xa, intern.ya, intern.rsmall, h, xo, yo);  // pair of points done
  l_pairs.push_back(TP_pair);                                 // l1, l2 in here

  h = intern.rlarge * distance / (intern.rlarge - intern.rsmall);
  pointcalc(intern.xb, intern.yb, intern.rlarge, h, xo, yo);  // pair of points done
  l_pairs.push_back(TP_pair);                                 // l3, l4

  // intersection
  intersectioncalc(distance);
  l_pairs.push_back(TP_pair);  // l5, l6

  // swap centres back to initial
  intern.rsmall = check_subtraction(rstemp, intern.ds);
  intern.rlarge = rltemp + intern.dl;
  if (swapped) {  // swap back to initial
    double xd = intern.xa;
    double yd = intern.ya;
    intern.xa = intern.xb;  // swap
    intern.xb = xd;
    intern.ya = intern.yb;
    intern.yb = yd;
  }
  h = intern.rsmall * distance / (intern.rlarge - intern.rsmall);
  xo = (1.0 + h / distance) * intern.xa - h / distance * intern.xb;
  yo = (1.0 + h / distance) * intern.ya - h / distance * intern.yb;
  pointcalc(intern.xa, intern.ya, intern.rsmall, h, xo, yo);  // pair of points done
  s_pairs.push_back(TP_pair);                                 // s1, s2 in here

  h = intern.rlarge * distance / (intern.rlarge - intern.rsmall);
  pointcalc(intern.xb, intern.yb, intern.rlarge, h, xo, yo);  // pair of points done
  s_pairs.push_back(TP_pair);                                 // s3, s4

  // intersection
  intersectioncalc(distance);
  s_pairs.push_back(TP_pair);  // s5, s6

  // prepare inputs to angle error conversion (x,y err to r, phi err)
  std::vector<ROOT::Math::XYVector> xy;
  std::vector<ROOT::Math::XYVector> dxdy;
  std::vector<ROOT::Math::XYVector> centres;
  if (swapped) {
    xy.push_back(allguest_points.at(0));  // fix 1-6
    xy.push_back(allguest_points.at(1));
    xy.push_back(allhost_points.at(0));
    xy.push_back(allhost_points.at(1));
    xy.push_back(allhost_points.at(3));
    xy.push_back(allhost_points.at(2));
    dxdy.push_back(l_pairs.at(0).first - s_pairs.at(1).second);  // l1-s4 gives dx and dy
    dxdy.push_back(l_pairs.at(0).second -
                   s_pairs.at(1).first);  // l2-s3 since XYVectors can subtract
    dxdy.push_back(l_pairs.at(1).first - s_pairs.at(0).second);  // l3-s2
    dxdy.push_back(l_pairs.at(1).second - s_pairs.at(0).first);  // l4-s1
    dxdy.push_back(l_pairs.at(2).first - s_pairs.at(2).second);  // l5-s6
    dxdy.push_back(l_pairs.at(2).second - s_pairs.at(2).first);  // l6-s5
    centres.push_back(ROOT::Math::XYVector(intern.xb, intern.yb));
    centres.push_back(ROOT::Math::XYVector(intern.xb, intern.yb));
    centres.push_back(ROOT::Math::XYVector(intern.xa, intern.ya));
    centres.push_back(ROOT::Math::XYVector(intern.xa, intern.ya));
    centres.push_back(ROOT::Math::XYVector(intern.xa, intern.ya));
    centres.push_back(ROOT::Math::XYVector(intern.xa, intern.ya));
  } else {
    xy.push_back(allhost_points.at(0));  // fix 1-6
    xy.push_back(allhost_points.at(1));
    xy.push_back(allguest_points.at(0));
    xy.push_back(allguest_points.at(1));
    xy.push_back(allhost_points.at(2));
    xy.push_back(allhost_points.at(3));
    dxdy.push_back(l_pairs.at(0).first - s_pairs.at(0).first);  // l1-s1 gives dx and dy
    dxdy.push_back(l_pairs.at(0).second -
                   s_pairs.at(0).second);  // l2-s2 since XYVectors can subtract
    dxdy.push_back(l_pairs.at(1).first - s_pairs.at(1).first);    // l3-s3
    dxdy.push_back(l_pairs.at(1).second - s_pairs.at(1).second);  // l4-s4
    dxdy.push_back(l_pairs.at(2).first - s_pairs.at(2).first);    // l5-s5
    dxdy.push_back(l_pairs.at(2).second - s_pairs.at(2).second);  // l6-s6
    centres.push_back(ROOT::Math::XYVector(intern.xa, intern.ya));
    centres.push_back(ROOT::Math::XYVector(intern.xa, intern.ya));
    centres.push_back(ROOT::Math::XYVector(intern.xb, intern.yb));
    centres.push_back(ROOT::Math::XYVector(intern.xb, intern.yb));
    centres.push_back(ROOT::Math::XYVector(intern.xa, intern.ya));
    centres.push_back(ROOT::Math::XYVector(intern.xa, intern.ya));
  }
  std::vector<Interval> dphi = deltaxy_to_deltaphi(xy, dxdy, centres);
  if (intern.order < 1) {
    TP_error = std::make_pair(dphi.at(0), dphi.at(1));  // result set, ready to store
    allhost_errors.push_back(TP_error.first);
    allhost_errors.push_back(TP_error.second);
    TP_error = std::make_pair(dphi.at(4), dphi.at(5));  // result set, ready to store
    allhost_errors.push_back(TP_error.first);
    allhost_errors.push_back(TP_error.second);
    TP_error = std::make_pair(dphi.at(2), dphi.at(3));  // result set, ready to store
    allguest_errors.push_back(TP_error.first);
    allguest_errors.push_back(TP_error.second);
    TP_error = std::make_pair(dphi.at(4), dphi.at(5));  // result set, ready to store
    allguest_errors.push_back(TP_error.first);
    allguest_errors.push_back(TP_error.second);
  } else {
    TP_error = std::make_pair(dphi.at(0), dphi.at(1));  // result set, ready to store
    allguest_errors.push_back(TP_error.first);
    allguest_errors.push_back(TP_error.second);
    TP_error = std::make_pair(dphi.at(4), dphi.at(5));  // result set, ready to store
    allguest_errors.push_back(TP_error.first);
    allguest_errors.push_back(TP_error.second);
    TP_error = std::make_pair(dphi.at(2), dphi.at(3));  // result set, ready to store
    allhost_errors.push_back(TP_error.first);
    allhost_errors.push_back(TP_error.second);
    TP_error = std::make_pair(dphi.at(4), dphi.at(5));  // result set, ready to store
    allhost_errors.push_back(TP_error.first);
    allhost_errors.push_back(TP_error.second);
  }
}

void RelationalHit::calculate_four(double distance) {
  //  std::cout << "tangent points: calc four case" << std::endl;
  if (fabs(host.gr.radius - guest.gr.radius) < 0.01) {  // capture equal radius case
    //    std::cout << "tangent points: equal rings case" << std::endl;
    equal_rings(distance, false);
    return;
  }

  size_sorting();  // initiate internal storage, must come before error calculation

  // on smaller ring first
  double h = intern.rsmall * distance / (intern.rlarge - intern.rsmall);
  double xo = (1.0 + h / distance) * intern.xa - h / distance * intern.xb;
  double yo = (1.0 + h / distance) * intern.ya - h / distance * intern.yb;
  pointcalc(intern.xa, intern.ya, intern.rsmall, h, xo, yo);  // pair of points done
  // std::cout << "tangent points: TP pair smaller ring" << std::endl;
  // std::cout << "x1=" << TP_pair.first.x() << " y1=" << TP_pair.first.y() << " x2=" <<
  // TP_pair.second.x() << " y2=" << TP_pair.second.y() << std::endl;
  if (intern.order < 1) {  // host is host
    allhost_points.push_back(TP_pair.first);
    allhost_points.push_back(TP_pair.second);
  } else {
    allguest_points.push_back(TP_pair.first);
    allguest_points.push_back(TP_pair.second);
  }

  // on larger ring
  h = intern.rlarge * distance / (intern.rlarge - intern.rsmall);
  pointcalc(intern.xb, intern.yb, intern.rlarge, h, xo, yo);  // pair of points done
  // std::cout << "tangent points: TP pair larger ring" << std::endl;
  // std::cout << "x1=" << TP_pair.first.x() << " y1=" << TP_pair.first.y() << " x2=" <<
  // TP_pair.second.x() << " y2=" << TP_pair.second.y() << std::endl;
  if (intern.order < 1) {  // host is host
    allguest_points.push_back(TP_pair.first);
    allguest_points.push_back(TP_pair.second);
  } else {
    allhost_points.push_back(TP_pair.first);
    allhost_points.push_back(TP_pair.second);
  }

  // origin between rings
  double ha = intern.rsmall * distance / (intern.rlarge + intern.rsmall);
  xo = (1.0 - ha / distance) * intern.xa + ha / distance * intern.xb;
  yo = (1.0 - ha / distance) * intern.ya + ha / distance * intern.yb;
  pointcalc(intern.xa, intern.ya, intern.rsmall, ha, xo, yo);  // pair of points done
  // std::cout << "tangent points: TP pair between rings" << std::endl;
  // std::cout << "x1=" << TP_pair.first.x() << " y1=" << TP_pair.first.y() << " x2=" <<
  // TP_pair.second.x() << " y2=" << TP_pair.second.y() << std::endl;
  if (intern.order < 1) {  // host is host
    allhost_points.push_back(TP_pair.first);
    allhost_points.push_back(TP_pair.second);
  } else {
    allguest_points.push_back(TP_pair.first);
    allguest_points.push_back(TP_pair.second);
  }

  // same origin between rings for larger ring
  double hb = intern.rlarge * distance / (intern.rlarge + intern.rsmall);
  pointcalc(intern.xb, intern.yb, intern.rlarge, hb, xo, yo);  // pair of points done
  // std::cout << "tangent points: TP pair between rings, larger ring" << std::endl;
  // std::cout << "x1=" << TP_pair.first.x() << " y1=" << TP_pair.first.y() << " x2=" <<
  // TP_pair.second.x() << " y2=" << TP_pair.second.y() << std::endl;
  if (intern.order < 1) {  // host is host
    allguest_points.push_back(TP_pair.first);
    allguest_points.push_back(TP_pair.second);
  } else {
    allhost_points.push_back(TP_pair.first);
    allhost_points.push_back(TP_pair.second);
  }
}

void RelationalHit::calculate_four_errors(double distance) {
  //  std::cout << "tangent points: calc four errors" << std::endl;
  bool swapped = false;  // origin must be on the correct side.
  std::vector<std::pair<ROOT::Math::XYVector, ROOT::Math::XYVector> > l_pairs;
  std::vector<std::pair<ROOT::Math::XYVector, ROOT::Math::XYVector> > s_pairs;
  double rstemp = intern.rsmall;
  double rltemp = intern.rlarge;
  if ((rstemp + intern.ds) > rltemp - intern.dl) {
    swapped = true;  // needs change
    intern.rsmall = check_subtraction(rltemp, intern.dl);
    intern.rlarge = rstemp + intern.ds;
    double xd = intern.xa;
    double yd = intern.ya;
    intern.xa = intern.xb;  // swap
    intern.xb = xd;
    intern.ya = intern.yb;
    intern.yb = yd;
  } else if ((rstemp + intern.ds) == rltemp - intern.dl) {  // div by zero risk
    intern.rsmall = rstemp + intern.ds - 0.1;               // allow fixed small separation
    intern.rlarge = check_subtraction(rltemp, intern.dl - 0.1);
  } else {
    intern.rsmall = rstemp + intern.ds;
    intern.rlarge = check_subtraction(rltemp, intern.dl);
  }
  // std::cout << "four errors: swap=" << swapped << std::endl;
  // std::cout << "four errors: rsmall=" << intern.rsmall << " rlarge=" << intern.rlarge <<
  // std::endl;

  // vary the ring radii, get tangent points on those -> errors in x,y
  double h = intern.rsmall * distance / (intern.rlarge - intern.rsmall);  // potential div by zero
  double xo = (1.0 + h / distance) * intern.xa - h / distance * intern.xb;
  double yo = (1.0 + h / distance) * intern.ya - h / distance * intern.yb;
  pointcalc(intern.xa, intern.ya, intern.rsmall, h, xo, yo);  // pair of points done
  l_pairs.push_back(TP_pair);                                 // l1, l2 in here

  h = intern.rlarge * distance / (intern.rlarge - intern.rsmall);
  pointcalc(intern.xb, intern.yb, intern.rlarge, h, xo, yo);  // pair of points done
  l_pairs.push_back(TP_pair);                                 // l3, l4

  // swap centres back to initial
  intern.rsmall = check_subtraction(rstemp, intern.ds);
  intern.rlarge = rltemp + intern.dl;
  if (swapped) {  // swap back to initial
    double xd = intern.xa;
    double yd = intern.ya;
    intern.xa = intern.xb;  // swap
    intern.xb = xd;
    intern.ya = intern.yb;
    intern.yb = yd;
  }
  h = intern.rsmall * distance / (intern.rlarge - intern.rsmall);
  xo = (1.0 + h / distance) * intern.xa - h / distance * intern.xb;
  yo = (1.0 + h / distance) * intern.ya - h / distance * intern.yb;
  pointcalc(intern.xa, intern.ya, intern.rsmall, h, xo, yo);  // pair of points done
  s_pairs.push_back(TP_pair);                                 // s1, s2 in here

  h = intern.rlarge * distance / (intern.rlarge - intern.rsmall);
  pointcalc(intern.xb, intern.yb, intern.rlarge, h, xo, yo);  // pair of points done
  s_pairs.push_back(TP_pair);                                 // s3, s4

  // prepare inputs to angle error conversion (x,y err to r, phi err)
  std::vector<ROOT::Math::XYVector> xy;
  std::vector<ROOT::Math::XYVector> dxdy;
  std::vector<ROOT::Math::XYVector> centres;
  if (intern.order) {
    xy.push_back(allguest_points.at(0));  // fix 1-4
    xy.push_back(allguest_points.at(1));
    xy.push_back(allhost_points.at(0));
    xy.push_back(allhost_points.at(1));
    if (swapped) {
      dxdy.push_back(l_pairs.at(0).first - s_pairs.at(1).second);  // l1-s4 gives dx and dy
      dxdy.push_back(l_pairs.at(0).second -
                     s_pairs.at(1).first);  // l2-s3 since XYVectors can subtract
      dxdy.push_back(l_pairs.at(1).first - s_pairs.at(0).second);  // l3-s2 gives dx and dy
      dxdy.push_back(l_pairs.at(1).second -
                     s_pairs.at(0).first);  // l4-s1 since XYVectors can subtract
    } else {
      dxdy.push_back(l_pairs.at(0).first - s_pairs.at(0).first);  // l1-s1 gives dx and dy
      dxdy.push_back(l_pairs.at(0).second -
                     s_pairs.at(0).second);  // l2-s2 since XYVectors can subtract
      dxdy.push_back(l_pairs.at(1).first - s_pairs.at(1).first);  // l3-s3 gives dx and dy
      dxdy.push_back(l_pairs.at(1).second -
                     s_pairs.at(1).second);  // l4-s4 since XYVectors can subtract
    }
  } else {
    xy.push_back(allhost_points.at(0));  // order reversed
    xy.push_back(allhost_points.at(1));
    xy.push_back(allguest_points.at(0));
    xy.push_back(allguest_points.at(1));
    if (swapped) {
      dxdy.push_back(l_pairs.at(0).first - s_pairs.at(1).second);  // l1-s4 gives dx and dy
      dxdy.push_back(l_pairs.at(0).second -
                     s_pairs.at(1).first);  // l2-s3 since XYVectors can subtract
      dxdy.push_back(l_pairs.at(1).first - s_pairs.at(0).second);  // l3-s2 gives dx and dy
      dxdy.push_back(l_pairs.at(1).second -
                     s_pairs.at(0).first);  // l4-s1 since XYVectors can subtract
    } else {
      dxdy.push_back(l_pairs.at(0).first - s_pairs.at(0).first);  // l1-s1 gives dx and dy
      dxdy.push_back(l_pairs.at(0).second -
                     s_pairs.at(0).second);  // l2-s2 since XYVectors can subtract
      dxdy.push_back(l_pairs.at(1).first - s_pairs.at(1).first);  // l3-s3 gives dx and dy
      dxdy.push_back(l_pairs.at(1).second -
                     s_pairs.at(1).second);  // l4-s4 since XYVectors can subtract
    }
  }
  centres.push_back(ROOT::Math::XYVector(intern.xa, intern.ya));
  centres.push_back(ROOT::Math::XYVector(intern.xa, intern.ya));
  centres.push_back(ROOT::Math::XYVector(intern.xb, intern.yb));
  centres.push_back(ROOT::Math::XYVector(intern.xb, intern.yb));

  std::vector<Interval> dphi = deltaxy_to_deltaphi(xy, dxdy, centres);
  //  std::cout << "four errors: done dphi, store TP_errors" << std::endl;
  if (intern.order < 1) {
    TP_error = std::make_pair(dphi.at(0), dphi.at(1));  // result set, ready to store
    allhost_errors.push_back(TP_error.first);
    allhost_errors.push_back(TP_error.second);
    TP_error = std::make_pair(dphi.at(2), dphi.at(3));  // result set, ready to store
    allguest_errors.push_back(TP_error.first);
    allguest_errors.push_back(TP_error.second);
  } else {
    TP_error = std::make_pair(dphi.at(0), dphi.at(1));  // result set, ready to store
    allguest_errors.push_back(TP_error.first);
    allguest_errors.push_back(TP_error.second);
    TP_error = std::make_pair(dphi.at(2), dphi.at(3));  // result set, ready to store
    allhost_errors.push_back(TP_error.first);
    allhost_errors.push_back(TP_error.second);
  }

  // between rings
  xy.clear();
  dxdy.clear();  // centres unchanged
  l_pairs.clear();
  s_pairs.clear();
  // origin between rings - here vary radii in the same
  // direction to get the error variation on the rings
  // change of symmetry compared to outer lines.
  // increase radii by errors
  intern.rsmall = rstemp + intern.ds;
  intern.rlarge = rltemp + intern.dl;
  double ring_overlap = (intern.rlarge + intern.rsmall) - distance;
  if (ring_overlap >= 0.0) {  // no origin between rings if they overlap
    intersectioncalc(distance);
    l_pairs.push_back(TP_pair);  // l1, l2 in here
    l_pairs.push_back(
        std::make_pair(TP_pair.second, TP_pair.first));  // swapped order of tp, l3, l4
  } else {
    // not overlapping on larger errors
    double ha = intern.rsmall * distance / (intern.rlarge + intern.rsmall);
    xo = (1.0 - ha / distance) * intern.xa + ha / distance * intern.xb;
    yo = (1.0 - ha / distance) * intern.ya + ha / distance * intern.yb;
    pointcalc(intern.xa, intern.ya, intern.rsmall, ha, xo, yo);  // pair of points done
    l_pairs.push_back(TP_pair);                                  // l1, l2 in here

    // same origin between rings for larger ring
    double hb = intern.rlarge * distance / (intern.rlarge + intern.rsmall);
    pointcalc(intern.xb, intern.yb, intern.rlarge, hb, xo, yo);  // pair of points done
    l_pairs.push_back(TP_pair);                                  // l3, l4 in here
  }

  // decrease radii by errors
  // no need here to check for overlaps
  intern.rsmall = check_subtraction(rstemp, intern.ds);
  intern.rlarge = check_subtraction(rltemp, intern.dl);
  double ha = intern.rsmall * distance / (intern.rlarge + intern.rsmall);
  xo = (1.0 - ha / distance) * intern.xa + ha / distance * intern.xb;
  yo = (1.0 - ha / distance) * intern.ya + ha / distance * intern.yb;
  pointcalc(intern.xa, intern.ya, intern.rsmall, ha, xo, yo);  // pair of points done
  s_pairs.push_back(TP_pair);                                  // s1, s2
  double hb = intern.rlarge * distance / (intern.rlarge + intern.rsmall);
  pointcalc(intern.xb, intern.yb, intern.rlarge, hb, xo, yo);  // pair of points done
  s_pairs.push_back(TP_pair);                                  // s3, s4

  if (intern.order) {
    xy.push_back(allguest_points.at(2));  // fix 1-4
    xy.push_back(allguest_points.at(3));
    xy.push_back(allhost_points.at(2));
    xy.push_back(allhost_points.at(3));
  } else {
    xy.push_back(allhost_points.at(2));
    xy.push_back(allhost_points.at(3));
    xy.push_back(allguest_points.at(2));  // fix 1-4
    xy.push_back(allguest_points.at(3));
  }
  dxdy.push_back(l_pairs.at(0).first - s_pairs.at(0).first);  // l1-s1 gives dx and dy
  dxdy.push_back(l_pairs.at(0).second -
                 s_pairs.at(0).second);                       // l2-s2 since XYVectors can subtract
  dxdy.push_back(l_pairs.at(1).first - s_pairs.at(1).first);  // l3-s3 gives dx and dy
  dxdy.push_back(l_pairs.at(1).second -
                 s_pairs.at(1).second);  // l4-s4 since XYVectors can subtract

  dphi = deltaxy_to_deltaphi(xy, dxdy, centres);
  //  std::cout << "four errors: done dphi, between rings" << std::endl;
  if (intern.order < 1) {
    TP_error = std::make_pair(dphi.at(0), dphi.at(1));  // result set, ready to store
    allhost_errors.push_back(TP_error.first);
    allhost_errors.push_back(TP_error.second);
    TP_error = std::make_pair(dphi.at(2), dphi.at(3));  // result set, ready to store
    allguest_errors.push_back(TP_error.first);
    allguest_errors.push_back(TP_error.second);
  } else {
    TP_error = std::make_pair(dphi.at(0), dphi.at(1));  // result set, ready to store
    allguest_errors.push_back(TP_error.first);
    allguest_errors.push_back(TP_error.second);
    TP_error = std::make_pair(dphi.at(2), dphi.at(3));  // result set, ready to store
    allhost_errors.push_back(TP_error.first);
    allhost_errors.push_back(TP_error.second);
  }
}

void RelationalHit::equal_rings(double distance, bool oflag) {
  size_sorting();  // initiate internal storage, must come before error calculation

  double alpha = TMath::ATan2(intern.yb - intern.ya, intern.xb - intern.xa);
  ROOT::Math::XYVector tp1(intern.xa + intern.rsmall * TMath::Sin(alpha),
                           intern.ya - intern.rsmall * TMath::Cos(alpha));
  ROOT::Math::XYVector tp2(intern.xa - intern.rsmall * TMath::Sin(alpha),
                           intern.ya + intern.rsmall * TMath::Cos(alpha));
  TP_pair = std::make_pair(tp1, tp2);  // result set, ready to store
  if (intern.order < 1) {              // host is host
    allhost_points.push_back(TP_pair.first);
    allhost_points.push_back(TP_pair.second);
  } else {
    allguest_points.push_back(TP_pair.first);
    allguest_points.push_back(TP_pair.second);
  }

  ROOT::Math::XYVector tp3(intern.xb + intern.rlarge * TMath::Sin(alpha),
                           intern.yb - intern.rlarge * TMath::Cos(alpha));
  ROOT::Math::XYVector tp4(intern.xb - intern.rlarge * TMath::Sin(alpha),
                           intern.yb + intern.rlarge * TMath::Cos(alpha));
  TP_pair = std::make_pair(tp3, tp4);  // result set, ready to store
  if (intern.order < 1) {              // host is host
    allguest_points.push_back(TP_pair.first);
    allguest_points.push_back(TP_pair.second);
  } else {
    allhost_points.push_back(TP_pair.first);
    allhost_points.push_back(TP_pair.second);
  }
  if (oflag) {
    intersectioncalc(distance);
    allhost_points.push_back(TP_pair.first);
    allhost_points.push_back(TP_pair.second);
    allguest_points.push_back(TP_pair.first);
    allguest_points.push_back(TP_pair.second);
  } else {
    // origin between rings
    double ha = intern.rsmall * distance / (intern.rlarge + intern.rsmall);
    double xo = (1.0 - ha / distance) * intern.xa + ha / distance * intern.xb;
    double yo = (1.0 - ha / distance) * intern.ya + ha / distance * intern.yb;
    pointcalc(intern.xa, intern.ya, intern.rsmall, ha, xo, yo);  // pair of points done
    if (intern.order < 1) {                                      // host is host
      allhost_points.push_back(TP_pair.first);
      allhost_points.push_back(TP_pair.second);
    } else {
      allguest_points.push_back(TP_pair.first);
      allguest_points.push_back(TP_pair.second);
    }
    double hb = intern.rlarge * distance / (intern.rlarge + intern.rsmall);
    pointcalc(intern.xb, intern.yb, intern.rlarge, hb, xo, yo);  // pair of points done
    if (intern.order < 1) {                                      // host is host
      allguest_points.push_back(TP_pair.first);
      allguest_points.push_back(TP_pair.second);
    } else {
      allhost_points.push_back(TP_pair.first);
      allhost_points.push_back(TP_pair.second);
    }
  }
}

void RelationalHit::intersectioncalc(double distance) {
  // Formulae to calculate the intersection points of two rings
  // in case they overlap.
  double dummy = (intern.rlarge * intern.rlarge - intern.rsmall * intern.rsmall);
  double d = TMath::Sqrt(fabs(distance * distance *
                                  (2.0 * intern.rlarge * intern.rlarge +
                                   2.0 * intern.rsmall * intern.rsmall - distance * distance) -
                              dummy * dummy));
  double p1 = (intern.yb - intern.ya) * d / (2 * distance * distance);
  double p2 = (intern.xb - intern.xa) * d / (2 * distance * distance);
  double x_1 = 0.5 * (intern.xa + intern.xb) -
               ((intern.xb - intern.xa) * dummy) / (2.0 * distance * distance) + p1;
  double x_2 = x_1 - 2.0 * p1;
  double y_1 = 0.5 * (intern.ya + intern.yb) -
               ((intern.yb - intern.ya) * dummy) / (2.0 * distance * distance) - p2;
  double y_2 = y_1 + 2.0 * p2;

  ROOT::Math::XYVector tp1(x_1, y_1);
  ROOT::Math::XYVector tp2(x_2, y_2);
  TP_pair = std::make_pair(tp1, tp2);  // result set, ready to store
}

double RelationalHit::checkforpi(double angle) {
  double pi = TMath::Pi();
  if (angle < -pi)
    angle += 2.0 * pi;
  else if (angle > pi)
    angle -= 2.0 * pi;
  return angle;
}

int RelationalHit::checkforpihalf(double& angle) {
  int sign = 1;
  if (angle < -TMath::Pi() / 2.0) {
    angle += TMath::Pi();
    sign = -1;
  } else if (angle > TMath::Pi() / 2.0) {
    angle -= TMath::Pi();
    sign = -1;
  }
  return sign;
}

double RelationalHit::check_subtraction(double r, double dr) {
  if ((r - dr) <= 0.0)
    return r;
  else
    return r - dr;
}

void RelationalHit::size_sorting() {
  // initiate internal store
  if (host.gr.radius <= guest.gr.radius) {
    intern.rsmall = host.gr.radius;
    intern.ds = host.gr.rerr;
    intern.rlarge = guest.gr.radius;
    intern.dl = guest.gr.rerr;
    intern.xa = host.gr.wirex;
    intern.ya = host.gr.wirey;
    intern.xb = guest.gr.wirex;
    intern.yb = guest.gr.wirey;
    intern.order = 0;
  } else {
    intern.rsmall = guest.gr.radius;
    intern.ds = guest.gr.rerr;
    intern.rlarge = host.gr.radius;
    intern.dl = host.gr.rerr;
    intern.xa = guest.gr.wirex;
    intern.ya = guest.gr.wirey;
    intern.xb = host.gr.wirex;
    intern.yb = host.gr.wirey;
    intern.order = 1;
  }
}

void RelationalHit::calc_onvertical(double xw, double yw, double rad, double h, double xo,
                                    double yo) {
  double y = yw - rad * rad / (h * h) * (yw - yo);
  double p2 = (yw - yo) * rad / (h * h) * TMath::Sqrt(h * h - rad * rad);
  double x_1 = xw - rad * rad / (h * h) * (xw - xo) - p2;
  double x_2 = xw - rad * rad / (h * h) * (xw - xo) + p2;

  ROOT::Math::XYVector tp1(x_1, y);
  ROOT::Math::XYVector tp2(x_2, y);
  TP_pair = std::make_pair(tp1, tp2);  // result set, ready to store
}

void RelationalHit::calc_onhorizontal(double xw, double yw, double rad, double h, double xo,
                                      double yo) {
  double x = xw - rad * rad / (h * h) * (xw - xo);
  double p2 = (xw - xo) * rad / (h * h) * TMath::Sqrt(h * h - rad * rad);
  double y_1 = yw - rad * rad / (h * h) * (yw - yo) - p2;
  double y_2 = yw - rad * rad / (h * h) * (yw - yo) + p2;

  ROOT::Math::XYVector tp1(x, y_1);
  ROOT::Math::XYVector tp2(x, y_2);
  TP_pair = std::make_pair(tp1, tp2);  // result set, ready to store
}

void RelationalHit::pointcalc(double xw, double yw, double rad, double h, double xo, double yo) {
  // Formulae to calculate the tangent points, all four from
  // this little function. Just needs different input:
  // which wire and ring radius is targeted, the
  // distance: centre to 'line connection crossing the
  // centre line between wire centres' = h, and the
  // origin or point of the tangent lines crossing: xo, yo.
  // This routine should cope with singularities like
  // vertical and horizontal lines.
  double EPS = 1.0e-3;  // define as small
  // special cases first
  if (fabs(yw - yo) <= EPS) {
    calc_onhorizontal(xw, yw, rad, h, xo, yo);
    return;
  } else if (fabs(xw - xo) <= EPS) {
    calc_onvertical(xw, yw, rad, h, xo, yo);
    return;
  } else {  // normal calculation
    double p2 = (yw - yo) * rad / (h * h) * TMath::Sqrt(h * h - rad * rad);
    double x_1 = xw - rad * rad / (h * h) * (xw - xo) + p2;
    double x_2 = xw - rad * rad / (h * h) * (xw - xo) - p2;

    double denom1 = rad * rad + (x_1 - xw) * (xw - xo);
    double denom2 = rad * rad + (x_2 - xw) * (xw - xo);
    double y_1 = (xw * yw * xo - x_1 * (xo * yw + xw * (yw - 2.0 * yo)) + x_1 * x_1 * (yw - yo) +
                  rad * rad * yo - xw * xw * yo) /
                 denom1;
    double y_2 = (xw * yw * xo - x_2 * (xo * yw + xw * (yw - 2.0 * yo)) + x_2 * x_2 * (yw - yo) +
                  rad * rad * yo - xw * xw * yo) /
                 denom2;

    ROOT::Math::XYVector tp1(x_1, y_1);
    ROOT::Math::XYVector tp2(x_2, y_2);
    TP_pair = std::make_pair(tp1, tp2);  // result set, ready to store
  }
}

ROOT::Math::XYVector RelationalHit::rphitoxy(ROOT::Math::XYVector tp, ROOT::Math::XYVector ctr,
                                             double r, double dr, double dphi) {
  double pih = TMath::Pi() / 2.0;
  double dx = 0.0;
  double dy = 0.0;
  double rel_x = tp.x() - ctr.x();
  double rel_y = tp.y() - ctr.y();
  double theta = TMath::ATan2(fabs(rel_y), fabs(rel_x));  // all prepared for top right quadrant
  dy = fabs((r + dr) * TMath::Sin(theta + dphi) - r * TMath::Sin(theta));
  dx = fabs((r + dr) * TMath::Cos(theta + dphi) - r * TMath::Cos(theta));

  if (dx < dr) dx = dr;  // set absolute minimal error in [mm]
  if (dy < dr) dy = dr;  // set absolute minimal error in [mm]

  // std::cout << " dx = " << dx << std::endl;
  // std::cout << " dy = " << dy << std::endl;
  ROOT::Math::XYVector err(dx, dy);
  return err;
}

std::vector<Interval> RelationalHit::deltaxy_to_deltaphi(
    std::vector<ROOT::Math::XYVector> xy, std::vector<ROOT::Math::XYVector> dxdy,
    std::vector<ROOT::Math::XYVector> ctr) {  // conversion to error on ring
  std::vector<Interval> dphi;
  int counter = 0;
  for (auto& dublet : xy) {
    double p = TMath::ATan2(dublet.y() - ctr.at(counter).y(), dublet.x() - ctr.at(counter).x());
    // std::cout << "dxy to dphi: point and centre " << dublet.x() << " " << dublet.y() << "; " <<
    // ctr.at(counter).x() << " " << ctr.at(counter).y() << std::endl; std::cout << "dxy to dphi:
    // angle to point " << p << std::endl;
    ROOT::Math::XYVector up = dublet + dxdy.at(counter);
    double pplus1 = TMath::ATan2(up.y() - ctr.at(counter).y(), up.x() - ctr.at(counter).x());
    //    std::cout << "dxy to dphi: angle to point + err " << pplus1 << std::endl;
    ROOT::Math::XYVector down = dublet - dxdy.at(counter);
    double pplus2 = TMath::ATan2(down.y() - ctr.at(counter).y(), down.x() - ctr.at(counter).x());
    //    std::cout << "dxy to dphi: angle to point - err " << pplus2 << std::endl;
    double p1 = checkforpi(pplus1);  // from
    double p2 = checkforpi(pplus2);  // to
    //    std::cout << "dxy to dphi: after pi check " << p1 << " " << p2 << std::endl;
    Interval ival(p1, p2);
    dphi.push_back(ival);
    //    std::cout << "dxy to dphi: interval[" << ival.from() << ", " << ival.to() << "]" <<
    //    std::endl;
    counter++;
  }

  return dphi;
}

// public method
std::vector<std::pair<int, int> > RelationalHit::get_edges() {
  std::vector<int> firsts;
  std::vector<int>::iterator it;
  for (auto& dublet : guestnode_map) firsts.push_back(dublet.first);

  int from;
  int to;
  int trial, idx;
  std::set<std::pair<int, int> > connections;
  for (unsigned int i = 0; i < hostnode_map.size(); i++) {  // same size as guest map
    from = hostnode_map.at(i).second;                       // starter host node
    trial = hostnode_map.at(i).first;                       // intermediate
    it = std::find(firsts.begin(), firsts.end(), trial);
    if (it != firsts.end())       // if found
      idx = it - firsts.begin();  // get the index
    else
      continue;
    to = guestnode_map.at(idx).second;
    connections.insert(std::make_pair(from, to));  // store edge as unique index pair
  }

  firsts.clear();
  for (auto& dublet : hostnode_map) firsts.push_back(dublet.first);
  for (unsigned int i = 0; i < guestnode_map.size(); i++) {  // same size as guest map
    from = guestnode_map.at(i).second;                       // starter host node
    trial = guestnode_map.at(i).first;                       // intermediate
    it = std::find(firsts.begin(), firsts.end(), trial);
    if (it != firsts.end())       // if found
      idx = it - firsts.begin();  // get the index
    else
      continue;
    to = hostnode_map.at(idx).second;
    connections.insert(std::make_pair(from, to));  // store edge as unique index pair
  }

  std::vector<std::pair<int, int> > edges;  // countable container
  for (auto& cn : connections) edges.push_back(cn);

  // and transform nodes into node container, same order vector -> vector
  PathPoint pp;
  size_t counter = 0;
  for (auto& nn : onthis) {
    pp.pointid = std::make_pair(host.mi.hitid, counter);  // unique id
    pp.xc = nn.first.x();
    pp.yc = nn.first.y();
    pp.errx = nn.second.x();
    pp.erry = nn.second.y();
    pp.zc = host.gr.zcoord;
    pp.errz = host.gr.zerr;
    hostnodes.push_back(pp);
    counter++;
  }
  counter = 0;
  for (auto& nn : onother) {
    pp.pointid = std::make_pair(guest.mi.hitid, counter);  // unique id
    pp.xc = nn.first.x();
    pp.yc = nn.first.y();
    pp.errx = nn.second.x();
    pp.erry = nn.second.y();
    pp.zc = guest.gr.zcoord;
    pp.errz = guest.gr.zerr;
    guestnodes.push_back(pp);
    counter++;
  }

  return edges;
}

// *****
// ** Utility Interval methods
// *****
Interval::Interval() {
  lower = 0.0;
  upper = 0.0;
}

Interval::Interval(double s, double e) {
  if (s <= e) {  // start should be smaller than end
    lower = s;
    upper = e;
  } else {  // swap order
    lower = e;
    upper = s;
  }
}

Interval Interval::hull(Interval other) {
  double left;
  double right;
  if (lower < other.from())
    left = lower;
  else
    left = other.from();
  if (upper > other.to())
    right = upper;
  else
    right = other.to();
  return Interval(left, right);
}

bool Interval::overlap(Interval other) {
  if (lower < other.from())
    return upper > other.from();
  else if (upper > other.to())
    return lower < other.to();
  else
    return true;  // this a subset
}

int Interval::checkforPiHalf(double& angle) {
  int sign = 1;
  if (angle < -TMath::Pi() / 2.0) {
    angle += TMath::Pi();
    sign = -1;
  } else if (angle > TMath::Pi() / 2.0) {
    angle -= TMath::Pi();
    sign = -1;
  }
  return sign;
}

double Interval::angle_midinterval() {
  // taking care of border at pi
  Interval temp(0, 0);
  double left = lower;  // subject to change
  double right = upper;
  int signleft = checkforPiHalf(left);
  int signright = checkforPiHalf(right);
  if (signleft < 0 && signright < 0) {
    temp = Interval(left, right);
    //    std::cout << "midint - after check for pi half: [" << left << ", " << right<< "]" <<
    //    std::endl;
    return temp.midinterval();
  }
  return 0.5 * (lower + upper);
}

double Interval::angle_dphi() {
  // taking care of border at pi
  double midpoint = angle_midinterval();
  bool swapped = false;
  double left = lower;                    // subject to change
  double right = upper;                   // one side is enough for dphi
  int signright = checkforPiHalf(right);  // check upper
  int signleft = checkforPiHalf(left);    // check upper
  if (signright < 0 && signleft < 0) {
    //    std::cout << "dphi - after check for pi half: " << right << " - " << midpoint <<
    //    std::endl;
    return fabs(right - midpoint);
  }
  return fabs(upper - midpoint);
}

bool Interval::angle_overlap(Interval other) {
  // std::cout << "Input interval: [" << lower << ", " << upper << "]" << std::endl;
  // std::cout << "compared to: [" << other.from() << ", " << other.to() << "]" << std::endl;
  Interval temp(0, 0);
  Interval tempo(0, 0);
  // bring to right hemisphere
  bool swapped1 = false;
  bool swapped2 = false;
  double left = lower;  // subject to change
  double right = upper;
  int signleft = checkforPiHalf(left);
  int signright = checkforPiHalf(right);
  if (signleft < 0 && signright < 0) {
    temp = Interval(left, right);
    swapped1 = true;
  }
  double otherleft = other.from();  // subject to change
  double otherright = other.to();
  int signoleft = checkforPiHalf(otherleft);
  int signoright = checkforPiHalf(otherright);
  if (signoleft < 0 && signoright < 0) {
    tempo = Interval(otherleft, otherright);
    swapped2 = true;
  }
  // if (swapped1)
  //   std::cout << "after check 1 for pi half: [" << left << ", " << right<< "]" << std::endl;
  // if (swapped2)
  //   std::cout << "after check 2 for pi half: [" << otherleft << ", " << otherright << "]" <<
  //   std::endl;
  if (swapped1 && swapped2)
    return temp.overlap(tempo);
  else if (swapped1 && !swapped2)
    return false;
  else if (!swapped1 && swapped2)
    return false;

  return overlap(other);
}

// *****
// ** Utility WeightedGraph methods
// *****
WeightedGraph::WeightedGraph(int nv) {
  V = nv;  // v vertices
  dist.clear();
  sPath.clear();
  adjList.clear();
  for (int i = 0; i < nv; i++) {
    // Create a vector to represent a row, and add it to the adjList.
    std::vector<std::pair<int, double> > row;
    adjList.push_back(row);
  }
}

void WeightedGraph::addEdge(int v, int w, double weight) {
  // for an undirect graph, vertices of an edge are each accessible from another
  adjList[v].push_back(std::make_pair(w, weight));  // Add w to v’s vector.
  adjList[w].push_back(std::make_pair(v, weight));  // and vice versa
}

// A BFS based function to check whether d is reachable from s.
bool WeightedGraph::isReachable(int s, int d) {
  // Base case
  if (s == d) return true;

  // Mark all the vertices as not visited
  bool* visited = new bool[V];

  for (int i = 0; i < V; i++) visited[i] = false;

  // Create a queue for BFS
  std::list<int> queue;

  // Mark the current node as visited and enqueue it
  visited[s] = true;

  queue.push_back(s);

  while (!queue.empty()) {
    // Dequeue a vertex from queue and print it
    s = queue.front();

    queue.pop_front();

    // Get all adjacent vertices of the dequeued vertex s

    // If a adjacent has not been visited, then mark it visited

    // and enqueue it
    for (auto& entry : adjList[s])  // pair
    {
      // If this adjacent node is the destination node, then return true
      if (entry.first == d) {
        delete[] visited;
        return true;
      }

      // Else, continue to do BFS
      if (!visited[entry.first]) {
        visited[entry.first] = true;

        queue.push_back(entry.first);
      }
    }
  }
  delete[] visited;
  return false;
}

// Given an Adjacency List, find all shortest paths from "start" to target vertex.
double WeightedGraph::dijkstraPaths(int start, int target) {
  // Initialize all source->vertex as infinite.
  dist.clear();
  sPath.clear();
  int n = adjList.size();
  for (int i = 0; i < n; i++)
    dist.push_back(
        std::make_pair(1000000007, i));  // Define "infinity" as necessary by constraints.

  // Create a PQ.
  std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double> >,
                      std::greater<std::pair<int, double> > >
      pq;

  // Add source to pq, where distance is 0.
  pq.push(std::make_pair(start, 0.0));
  dist[start] = std::make_pair(0.0, start);
  ;

  // While pq isn't empty...
  while (pq.empty() == false) {
    // Get min distance vertex from pq. (Call it u.)
    int u = pq.top().first;
    pq.pop();

    // Visit all of u's friends. For each one (called v)....
    for (int i = 0; i < adjList[u].size(); i++) {
      int v = adjList[u][i].first;
      double weight = adjList[u][i].second;

      // If the distance to v is shorter by going through u...
      if (dist[v].first > dist[u].first + weight) {
        // Update the distance of v.
        dist[v].first = dist[u].first + weight;
        // Update the previous node of v.
        dist[v].second = u;
        // Insert v into the pq.
        pq.push(std::make_pair(v, dist[v].first));
      }
    }
  }
  int currenttarget = target;
  while (currenttarget != start) {
    sPath.push_back(currenttarget);
    currenttarget = dist[currenttarget].second;
  }
  sPath.push_back(start);  // finish with the start index
  return dist[target].first;
}
