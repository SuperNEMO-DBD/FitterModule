#include "SNFitter.h"

#include "SNFitterImpl.h"

#include <algorithm>
#include <valarray>

// ROOT
#include "Fit/Fitter.h"
#include "Math/Functor.h"
#include "TLinearFitter.h"
#include "TSpectrum.h"
#include "TH1D.h"

// *** SNFitter section
// ********************
// lines in SNFitter
std::vector<double> SNFitter::line_initials(double frad) {
  // provide options for x-y plane initials. X-Z plane not needed
  std::vector<double> results;

  auto* lfxy = new TLinearFitter();
  lfxy->SetFormula("pol1");

  unsigned int nd = grings.size();
  auto* xx = new double[nd];
  auto* yy = new double[nd];
  auto* ey = new double[nd];
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
  for (GeigerRing gg : grings) {
    dummy.push_back(gg.wirex);
  }
  if (leftright > 0) {
    mit = std::min_element(dummy.begin(), dummy.end());
  } else {
    mit = std::max_element(dummy.begin(), dummy.end());
  }

  int idx = mit - dummy.begin();
  double first_radius = grings.at(idx).radius;

  // initials for x-y plane
  std::vector<double> res = line_initials(first_radius);
  // for x-z plane
  unsigned int nd = grings.size();
  auto* xx = new double[nd];
  auto* zz = new double[nd];
  auto* ez = new double[nd];
  for (unsigned int j = 0; j < nd; j++) {
    xx[j] = grings.at(j).wirex;
    zz[j] = grings.at(j).zcoord;
    ez[j] = grings.at(j).zerr;
  }
  // fit in x-z plane, once only
  auto* lfxz = new TLinearFitter();
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
  auto* xx = new double[nd];
  auto* yy = new double[nd];
  auto* zz = new double[nd];
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

  if (fabs(denominator) < 1.0e-16) {  // not a  helix but a line in wires
    return results;                   // return empty as signal
  }

  double uc = (c2 * suv - c1 * svv) / denominator;
  double vc = (c1 * suv - c2 * suu) / denominator;
  double rad = uc * uc + vc * vc + (suu + svv) / nd;  // squared radius
  results.push_back(TMath::Sqrt(rad));                // helix parameter 0

  auto* lfxz = new TLinearFitter();
  lfxz->SetFormula("pol1");  // straight line x-z plane
  lfxz->AssignData((int)nd, 1, xx, zz);
  lfxz->Eval();
  if (lfxz->GetParameter(1) > 0) {  // up or down pitch
    results.push_back(1.0);         // right sign
  } else {
    results.push_back(-1.0);  // right sign
  }

  results.push_back(uc + meanx);  // centre x
  results.push_back(vc + meany);  // centre y

  std::vector<double> dummy;
  for (GeigerRing gg : grings) {
    dummy.push_back(gg.wirex);
  }
  auto maxit = std::max_element(dummy.begin(), dummy.end());
  auto minit = std::min_element(dummy.begin(), dummy.end());
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
  for (GeigerRing gg : grings) {
    dummy.push_back(gg.wirex);
  }
  auto maxit = std::max_element(dummy.begin(), dummy.end());
  auto minit = std::min_element(dummy.begin(), dummy.end());
  int maxidx = maxit - dummy.begin();
  int minidx = minit - dummy.begin();
  double zmean = 0.5 * (grings.at(maxidx).zcoord + grings.at(minidx).zcoord);
  double ymean = 0.5 * (grings.at(maxidx).wirey + grings.at(minidx).wirey);
  // some guess dummy values
  results.push_back(100.0);  // r

  if (zmean >= 0) {          // up or down pitch
    results.push_back(1.0);  // h
  } else {
    results.push_back(-1.0);
  }

  results.push_back(0.0);  // xc
  if (ymean >= 0) {
    results.push_back(100.0);  // yc
  } else {
    results.push_back(-100.0);  // yc
  }
  results.push_back(zmean);  // z centre
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
  if (res.empty()) {  // not a helix in wires, make some up
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
  if (betaangles.size() < 5) {
    return false;  // pointless search
  }
  double alpha = 0.9;
  double threshold = 0.11;
  double prval = betaangles.front();  // first
  unsigned int len = betaangles.size();
  unsigned int count = 0;
  for (double angle : betaangles) {
    double val = martingale(prval, angle, alpha);
    if (fabs(val - prval) > threshold && count < len - 1) {
      return true;  // exclude final value kinks - dont' work with TSpectrum
    }
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
  auto* sp = new TSpectrum();
  int npeaks = sp->Search(hist, 1, "goff", 0.9);  // 90% threshold for second peak
  if (npeaks == 0) {                              // try one more
    npeaks = sp->Search(hist, 2.0, "goff", 0.9);
  }
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

LineFit SNFitter::fitline2D(const std::vector<PathPoint>& data) {
  LineFit lf;
  // separate linear fits in xy and xz
  // good for initial conditions or separate broken line pieces
  auto* lfxy = new TLinearFitter();
  auto* lfxz = new TLinearFitter();

  int nd = (int)data.size();
  auto* xx = new double[nd];
  auto* yy = new double[nd];
  auto* zz = new double[nd];
  auto* ey = new double[nd];
  auto* ez = new double[nd];
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
  if (pathcollection.empty()) {
    return results;  // return empty container
  }

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
        if (bin1 < 2) {  // have at least 2 points to fit to
          bin1 += 2;
        }
        for (unsigned int i = 0; i < bin1; i++) {
          foilside.push_back(ppc.path.at(i));
        }

        int bin2 = blf.breakpoints.back();
        if (bin2 < 2) {  // have at least 2 points to fit to
          bin2 += 2;
        }
        for (unsigned int i = ppc.path.size() - bin2; i < ppc.path.size(); i++) {
          caloside.push_back(ppc.path.at(i));
        }

        blf.linefit1 = fitline2D(foilside);
        blf.linefit2 = fitline2D(caloside);
      } else {  // exactly one break point
        int bin1 = blf.breakpoints.front();
        if (bin1 < 2) {  // have at least 2 points to fit to
          bin1 += 2;
        }
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