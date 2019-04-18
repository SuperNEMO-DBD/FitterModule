// us
#include <fitter_library.h>

// ROOT includes
#include <Math/Vector3D.h>
#include <Math/Functor.h>
#include <Fit/Fitter.h>

// std
#include <valarray>
#include <algorithm>
#include <iostream>



// *** Line section
// *****************

// calculate distance line-cylinder
double LineDistance2::ldistance(GeigerRing gr, const double *p) {
  // distance line cylinder is D= | (xp-x0) cross  ux | - r
  // where ux is direction of line and x0 is a point on the line (like t = 0)
  // line not parallel to y-z plane, i.e calo plane
  double radius = gr.radius;
  double x = gr.wirex;
  double y = gr.wirey;
  double z = gr.zcoord;
  ROOT::Math::XYZVector weight;

  ROOT::Math::XYZVector xp(x,y,z);
  ROOT::Math::XYZVector x0(0., p[0], p[2]);
  ROOT::Math::XYZVector x1(1., p[0] + p[1], p[2] + p[3]);
  ROOT::Math::XYZVector u = (x1-x0).Unit();
  double d2 = ((xp-x0).Cross(u)).Mag2();
  weight.SetXYZ((gr.rerr*u.Z()*gr.rerr*u.Z())+(gr.zerr*u.Y()*gr.zerr*u.Y()),
		(gr.rerr*u.Z()*gr.rerr*u.Z())+(gr.zerr*u.X()*gr.zerr*u.X()),
		(gr.rerr*u.Y()*gr.rerr*u.Y())+(gr.rerr*u.X()*gr.rerr*u.X())); // weight sqr vector
  double w2 = weight.X() + weight.Y() + weight.Z(); // sum of squares is .Mag2()
  return (TMath::Sqrt(d2) - radius) / TMath::Sqrt(w2); // only this is ring specific
}


// lines in SNFitter
std::vector<double> SNFitter::line_initials(double frad) {
  // provide options for x-y plane initials. X-Z plane not needed
  std::vector<double> results;

  TLinearFitter *lfxy = new TLinearFitter();
  lfxy->SetFormula("pol1");

  unsigned int nd = rings.size();
  double* xx = new double [nd];
  double* yy = new double [nd];
  double* ey = new double [nd];
  for (unsigned int j=0;j<nd;j++) {
    xx[j] = rings.at(j).wirex;
    yy[j] = rings.at(j).wirey;
    ey[j] = rings.at(j).rerr;
  }
  lfxy->AssignData((int)nd,1,xx,yy,ey);
  lfxy->Eval();
  // each attempt in dublets for the 2 parameter
  // part 1
  results.push_back(lfxy->GetParameter(0)-frad); // - shift intercept in y by radius
  double slope = lfxy->GetParameter(1); 
  double slopeaddon = 0.01;
  results.push_back(slope+slopeaddon);
  // part 2
  results.push_back(lfxy->GetParameter(0)-frad); // - shift intercept in y
  results.push_back(slope-slopeaddon); // minus add-on slope
  // part 3
  results.push_back(lfxy->GetParameter(0)+frad); // + shift intercept in y
  results.push_back(slope+slopeaddon); // same slope as before
  // part 4
  results.push_back(lfxy->GetParameter(0)+frad); // + shift intercept in y
  results.push_back(slope-slopeaddon);

  delete [] xx;
  delete [] yy;
  delete [] ey;
  delete lfxy;

  return results;
}

std::vector<LineFit> SNFitter::fitline() {
  std::vector<LineFit> results;
  LineFit lf;
  ROOT::Fit::Fitter fitter;
  fitter.Config().SetDefaultMinimizer("Minuit2");
  
  // make the functor object
  LineDistance2 ldist(&rings);
  ROOT::Math::Functor fcn(ldist,4);

  // first initials
  std::vector<double> dummy;
  std::vector<double>::iterator mit;
  double leftright = rings[0].wirex;
  for (GeigerRing gg  :rings)
    dummy.push_back(gg.wirex);
  if (leftright>0)
    mit = std::min_element(dummy.begin(), dummy.end());
  else
    mit = std::max_element(dummy.begin(), dummy.end());

  int idx = mit - dummy.begin();
  double first_radius = rings.at(idx).radius;

  // initials for x-y plane
  std::vector<double> res = line_initials(first_radius);
  // for x-z plane
  unsigned int nd = rings.size();
  double* xx = new double [nd];
  double* zz = new double [nd];
  double* ez = new double [nd];
  for (unsigned int j=0;j<nd;j++) {
    xx[j] = rings.at(j).wirex;
    zz[j] = rings.at(j).zcoord;
    ez[j] = rings.at(j).zerr;
  }
  // fit in x-z plane, once only
  TLinearFitter *lfxz = new TLinearFitter();
  lfxz->SetFormula("pol1");
  lfxz->AssignData((int)nd,1,xx,zz,ez);
  lfxz->Eval();
  //  cout << "xz: i " << lfxz->GetParameter(0) << ", sl " << lfxz->GetParameter(1) << endl;

  // full fit attempts
  double iniline[4];
  for (int i=0;i<8;i+=2) { // 4 candidates for initials in x-y
    //    cout << "xy: i " << res[i] << ", sl " << res[i+1] << endl;
    iniline[0] = res[i];
    iniline[1] = res[i+1];
    iniline[2] = lfxz->GetParameter(0);
    iniline[3] = lfxz->GetParameter(1);

    // fit
    fitter.SetFCN(fcn, iniline, (unsigned int)rings.size(), true);
    
    bool ok = fitter.FitFCN();
    const ROOT::Fit::FitResult & result = fitter.Result();
    lf.status = result.Status(); // fit to ring in 2d status
    lf.chi2 = result.Chi2();
    lf.prob = result.Prob();
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
  delete [] ez;
  delete [] zz;
  delete [] xx;
  delete lfxz;
  return results;
}


// *** Helix section
// *****************

// calculate distance helix-ring
double HelixDistance2::helixdistance(GeigerRing gr, const double *p) {
  double ringr = gr.radius;
  double ringx = gr.wirex;
  double ringy = gr.wirey;
  double ringz = gr.zcoord;
  
  // weighting squared
  double weight = TMath::Sqrt(gr.rerr * gr.rerr + gr.zerr * gr.zerr); // weight sqrt
  
  // p[0]=r, p[1]=h, p[2,3,4]=ref point
  ROOT::Math::XYZVector xp(ringx, ringy, ringz); // ring centre
  ROOT::Math::XYZVector xref(p[2], p[3], p[4]); // helix centre point
  ROOT::Math::XYZVector ray = xref - xp; // connection vector
  
  double dir = TMath::ATan2(ray.y(), ray.x())+TMath::Pi();
  ROOT::Math::XYZVector xhel(p[2] + p[0]*TMath::Cos(dir), p[3] + p[0]*TMath::Sin(dir), p[4] + p[1]*dir); // point on helix
  double d2 = (xhel - xp).Mag2(); // distance squared

  return (TMath::Sqrt(d2) - ringr) / weight; // only this is ring specific;
}


std::vector<double> SNFitter::helix_initials() {
  std::vector<double> results;

  unsigned int nd = rings.size();
  double* xx = new double [nd];
  double* yy = new double [nd];
  double* zz = new double [nd];
  for (unsigned int j=0;j<nd;j++) {
    xx[j] = rings.at(j).wirex;
    yy[j] = rings.at(j).wirey;
    zz[j] = rings.at(j).zcoord;
  }

  std::valarray<double> valx (xx, nd);
  std::valarray<double> valy (yy, nd);
  
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
  double denominator = (suv*suv - suu*svv); // prevent div by zero

  if (fabs(denominator) < 1.0e-16) // not a  helix but a line in wires
    return results; // return empty as signal

  double uc = (c2 * suv - c1 * svv) / denominator;
  double vc = (c1 * suv - c2 * suu) / denominator;
  double rad = uc*uc + vc*vc + (suu + svv)/ nd; // squared radius
  results.push_back(TMath::Sqrt(rad)); // helix parameter 0

  TLinearFitter *lfxz = new TLinearFitter();
  lfxz->SetFormula("pol1"); // straight line x-z plane
  lfxz->AssignData((int)nd,1,xx,zz);
  lfxz->Eval();
  if (lfxz->GetParameter(1)>0) // up or down pitch
    results.push_back(1.0); // right sign 
  else
    results.push_back(-1.0); // right sign 

  results.push_back(uc + meanx); // centre x
  results.push_back(vc + meany); // centre y

  std::vector<double> dummy;
  for (GeigerRing gg  :rings)
    dummy.push_back(gg.wirex);
  std::vector<double>::iterator maxit = std::max_element(dummy.begin(), dummy.end());
  std::vector<double>::iterator minit = std::min_element(dummy.begin(), dummy.end());
  int maxidx = maxit - dummy.begin();
  int minidx = minit - dummy.begin();
  results.push_back(0.5*(rings.at(maxidx).zcoord+rings.at(minidx).zcoord)); // z centre

  delete [] zz;
  delete [] yy;
  delete [] xx;
  delete lfxz;
  
  return results;
}


std::vector<double> SNFitter::helixbackup() {
  std::vector<double> results;
  std::vector<double> dummy;
  for (GeigerRing gg  : rings)
    dummy.push_back(gg.wirex);
  std::vector<double>::iterator maxit = std::max_element(dummy.begin(), dummy.end());
  std::vector<double>::iterator minit = std::min_element(dummy.begin(), dummy.end());
  int maxidx = maxit - dummy.begin();
  int minidx = minit - dummy.begin();
  double zmean = 0.5*(rings.at(maxidx).zcoord+rings.at(minidx).zcoord);
  double ymean = 0.5*(rings.at(maxidx).wirey+rings.at(minidx).wirey);
  // some guess dummy values
  results.push_back(100.0); // r

  if (zmean>=0) // up or down pitch
    results.push_back(1.0); // h
  else
    results.push_back(-1.0);

  results.push_back(0.0); // xc
  if (ymean>=0)
    results.push_back(100.0); // yc
  else
    results.push_back(-100.0); // yc
  results.push_back(zmean); // z centre
  return results;
}


std::vector<HelixFit> SNFitter::fithelix() {
  std::vector<HelixFit> results;
  HelixFit hf;
  bool store = true;
  ROOT::Fit::Fitter fitter;
  
  // make the functor object
  HelixDistance2 hdist(&rings);
  ROOT::Math::Functor fcn(hdist,5);

  std::vector<double> res = helix_initials();
  if (res.size()<1) {// not a helix in wires, make some up
    res.clear();
    res = helixbackup();
  }

  double inihelix[5] = {res[0], res[1], res[2], res[3], res[4]};

  fitter.SetFCN(fcn, inihelix, (unsigned int)rings.size(), true);

  // fit
  bool ok = fitter.FitFCN();
  if (!ok) {
    std::cout << "helix3Dfit: Helix3D Fit failed: 2nd attempt started." << std::endl;
    // 2nd possibility
    res = helixbackup();
    double inihelix2[5] = {res[0], res[1], res[2], res[3], res[4]};

    fitter.SetFCN(fcn, inihelix2, (unsigned int)rings.size(), true);

    // fit
    ok = fitter.FitFCN();

    if (!ok) {
      std::cout << "helix3Dfit: Helix3D Fit failed" << std::endl;
      std::cout << "chi square " << result.Chi2() << " is valid: " << result.Status() << " cov matrix: " << result.CovMatrixStatus() << std::endl;
      store = false;
    }
  }
  if (store) {
    const ROOT::Fit::FitResult & result = fitter.Result();
    hf.status = result.Status(); // fit to ring in 2d status
    hf.chi2 = result.Chi2();
    hf.prob = result.Prob();
    //    result.Print(std::cout);
    
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
  }
  return results; // could be empty - for no fit
}

