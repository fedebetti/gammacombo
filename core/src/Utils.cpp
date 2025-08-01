#include <boost/algorithm/string.hpp>

#include <algorithm>
#include <format>
#include <memory>
#include <ranges>
#include <string>
#include <vector>

#include <TColor.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphSmooth.h>
#include <TMatrixDSymEigen.h>
#include <TPaveText.h>
#include <TVectorD.h>

#include <RooMinimizer.h>
#include <RooProdPdf.h>
#include <RooRandom.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>

#include <Utils.h>

using namespace std;
using namespace RooFit;

int Utils::countFitBringBackAngle;     ///< counts how many times an angle needed to be brought back
int Utils::countAllFitBringBackAngle;  ///< counts how many times fitBringBackAngle() was called

///
/// Fit PDF to minimum.
/// \param pdf The PDF.
/// \param thorough Activate Hesse and Minos
/// \param printLevel -1 = no output, 1 verbose output
///
std::unique_ptr<RooFitResult> Utils::fitToMin(RooAbsPdf* pdf, bool thorough, int printLevel) {
  RooMsgService::instance().setGlobalKillBelow(ERROR);

  // pdf->Print("v");
  // const RooProdPdf *prod = (RooProdPdf*)pdf;
  // const RooArgList& pdflist = prod->pdfList();
  // double chi2=0;
  // for (int i=0; i<pdflist.getSize(); i++){
  // const RooAbsPdf *ipdf = (RooAbsPdf*)pdflist.at(i);
  // chi2 += -2*TMath::Log(ipdf->getVal());
  // cout << ipdf->GetName() << " " << -2*TMath::Log(ipdf->getVal()) << " " << chi2 << endl;
  //////ipdf->Print("v");
  //}
  // cout << chi2 << endl;

  RooFormulaVar ll("ll", "ll", "-2*log(@0)", RooArgSet(*pdf));
  // cout << ll.getVal() << endl;
  // string a;
  // cin >> a;
  bool quiet = printLevel < 0;
  RooMinimizer m(ll);
  if (quiet) {
    m.setPrintLevel(-2);
  } else
    m.setPrintLevel(1);
  // if (quiet) m.setLogFile();
  m.setErrorLevel(1.0);
  m.setStrategy(2);
  m.setProfile(0);  // 1 enables migrad timer
  unsigned long long start = rdtsc();
  int status = m.migrad();
  // m.simplex();
  // m.migrad();
  // m.simplex();
  if (thorough) {
    m.hesse();
    // MINOS seems to fail in more complicated scenarios
    // Can't just use m.minos() because there's a cout that cannot be turned off (root 5-34-03).
    // It's not there when we run minos only on selected parameters- so select them all!
    // //m.minos();
    // RooArgSet floatingPars;
    // for ( const auto p : *pdf->getVariables()) {
    //   if ( !p->isConstant() ) floatingPars.add(*p);
    // }
    // m.minos(floatingPars);
    // IMPROVE doesn't really improve much
    // //m.improve();
  }
  unsigned long long stop = rdtsc();
  if (!quiet) std::printf("Fit took %llu clock cycles.\n", stop - start);
  std::unique_ptr<RooFitResult> r(m.save());
  // if (!quiet) r->Print("v");
  RooMsgService::instance().setGlobalKillBelow(INFO);
  return r;
}

///
/// Return an equivalent angle between 0 and 2pi.
/// \param angle Angle that is possibly smaller than 0 or larger than 2pi.
///
double Utils::bringBackAngle(double angle) {
  double val = fmod(angle, 2. * TMath::Pi());
  if (val < 0.0) val = val + 2. * TMath::Pi();
  return val;
}

///
/// Compute the difference between 2 angles. This will never be
/// larger than pi because they wrap around!
///
/// \param angle1 - first angle
/// \param angle2 - second angle
/// \return difference
///
double Utils::angularDifference(double angle1, double angle2) {
  double angleSmaller = max(bringBackAngle(angle1), bringBackAngle(angle2));
  double angleLarger = min(bringBackAngle(angle1), bringBackAngle(angle2));
  double diff1 = angleLarger - angleSmaller;
  double diff2 = (2. * TMath::Pi() - angleLarger) + angleSmaller;
  return min(diff1, diff2);
}

///
/// Fit a pdf to the minimum, but keep angular parameters in a range of
/// [0,2pi]. If after an initial fit, a parameter has walked outside this
/// interval, add multiples of 2pi to bring it back. Then, refit.
/// All variables that have unit 'rad' are taken to be angles.
///
std::unique_ptr<RooFitResult> Utils::fitToMinBringBackAngles(RooAbsPdf* pdf, bool thorough, int printLevel) {
  countAllFitBringBackAngle++;
  auto r = fitToMin(pdf, thorough, printLevel);
  bool refit = false;
  for (const auto pAbs : r->floatParsFinal()) {
    const auto p = static_cast<RooRealVar*>(pAbs);
    if (!isAngle(p)) continue;
    if (p->getVal() < 0.0 || p->getVal() > 2. * TMath::Pi()) {
      std::unique_ptr<RooArgSet> pdfPars(pdf->getParameters(RooArgSet()));
      auto pdfPar = static_cast<RooRealVar*>(pdfPars->find(p->GetName()));
      pdfPar->setVal(bringBackAngle(p->getVal()));
      refit = true;
    }
  }
  if (refit) {
    countFitBringBackAngle++;
    r = fitToMin(pdf, thorough, printLevel);
  }
  return r;
}

///
/// Find the global minimum in a more thorough way.
/// First fit with external start parameters, then
/// for each parameter that starts with "d" or "r" (typically angles and ratios):
///   - at upper scan range, rest at start parameters
///   - at lower scan range, rest at start parameters
/// This amounts to a maximum of 1+2^n fits, where n is the number
/// of parameters to be varied.
///
/// \param w Workspace holding the pdf.
/// \param name Name of the pdf without leading "pdf_".
/// \param forceVariables Apply the force method for these variables only. Format
/// "var1,var2,var3," (list must end with comma). Default is to apply for all angles,
/// all ratios except rD_k3pi and rD_kpi, and the k3pi coherence factor.
///
std::unique_ptr<RooFitResult> Utils::fitToMinForce(RooWorkspace* w, TString name, TString forceVariables, bool debug) {
  TString parsName = "par_" + name;
  TString obsName = "obs_" + name;
  TString pdfName = "pdf_" + name;
  int printlevel = -1;
  RooMsgService::instance().setGlobalKillBelow(ERROR);

  // save start parameters
  if (!w->set(parsName)) {
    cout << "MethodProbScan::scan2d() : ERROR : parsName not found: " << parsName << endl;
    exit(1);
  }
  RooDataSet startPars("startParsForce", "startParsForce", *w->set(parsName));
  startPars.add(*w->set(parsName));

  // set up parameters and ranges
  auto varyPars = std::make_unique<RooArgList>();
  for (const auto p : *w->set(parsName)) {
    if (p->isConstant()) continue;
    if (forceVariables == "" &&
        (false ||
         TString(p->GetName()).BeginsWith("d")  ///< use these variables
         // || TString(p->GetName()).BeginsWith("r")
         || TString(p->GetName()).BeginsWith("k") || TString(p->GetName()) == "g") &&
        !(TString(p->GetName()) == "rD_k3pi"  ///< don't use these
          || TString(p->GetName()) == "rD_kpi"
          // || TString(p->GetName()) == "dD_kpi"
          || TString(p->GetName()) == "d_dk" || TString(p->GetName()) == "d_dsk")) {
      varyPars->add(*p);
    } else if (forceVariables.Contains(TString(p->GetName()) + ",")) {
      varyPars->add(*p);
    }
  }
  int nPars = varyPars->getSize();
  if (debug) cout << "Utils::fitToMinForce() : nPars = " << nPars << " => " << pow(2., nPars) << " fits" << endl;
  if (debug) cout << "Utils::fitToMinForce() : varying ";
  if (debug) varyPars->Print();

  //////////

  auto r = fitToMinBringBackAngles(w->pdf(pdfName), false, printlevel);

  //////////

  int nErrors = 0;

  // We define a binary mask where each bit corresponds
  // to parameter at max or at min.
  for (int i = 0; i < pow(2., nPars); i++) {
    if (debug) cout << "Utils::fitToMinForce() : fit " << i << "        \r" << flush;
    setParameters(w, parsName, startPars.get(0));

    for (int ip = 0; ip < nPars; ip++) {
      auto p = (RooRealVar*)varyPars->at(ip);
      double oldMin = p->getMin();
      double oldMax = p->getMax();
      setLimit(w, p->GetName(), "force");
      if (i / (int)pow(2., ip) % 2 == 0) { p->setVal(p->getMin()); }
      if (i / (int)pow(2., ip) % 2 == 1) { p->setVal(p->getMax()); }
      p->setRange(oldMin, oldMax);
    }

    // check if start parameters are sensible, skip if they're not
    double startParChi2 = getChi2(w->pdf(pdfName));
    if (startParChi2 > 2000) {
      nErrors += 1;
      continue;
    }

    // refit
    auto r2 = fitToMinBringBackAngles(w->pdf(pdfName), false, printlevel);

    // In case the initial fit failed, accept the second one.
    // If both failed, still select the second one and hope the
    // next fit succeeds.
    if (!(r->edm() < 1 && r->covQual() == 3)) {
      r = std::move(r2);
    } else if (r2->edm() < 1 && r2->covQual() == 3 && r2->minNll() < r->minNll()) {
      // better minimum found!
      r = std::move(r2);
    }
  }

  if (debug) cout << endl;
  if (debug) cout << "Utils::fitToMinForce() : nErrors = " << nErrors << endl;

  RooMsgService::instance().setGlobalKillBelow(INFO);

  // (re)set to best parameters
  setParameters(w, parsName, r.get());

  return r;
}

///
/// Fit to the minimum using an improve method.
/// The idea is to kill known minima. For this, we add a multivariate
/// Gaussian to the chi2 at the position of the minimum. Ideally it should
/// be a parabola shaped function, but the Gaussian is readily available.
/// A true parabola like in the original IMPROVE algorithm doesn't work
/// because apparently it affects the other regions of
/// the chi2 too much. There are two parameters to tune: the width of the
/// Gaussian, which is controlled by the error definition of the first fit,
/// and the height of the Gaussian, which is set to 16 (=4 sigma).
/// Three fits are performed: 1) initial fit, 2) fit with improved FCN,
/// 3) fit of the non-improved FCN using the result from step 2 as the start
/// parameters.
/// So far it is only available for the Prob method, via the probimprove
/// command line flag.
///
std::unique_ptr<RooFitResult> Utils::fitToMinImprove(RooWorkspace* w, TString name) {
  TString parsName = "par_" + name;
  TString obsName = "obs_" + name;
  TString pdfName = "pdf_" + name;
  int printlevel = -1;
  RooMsgService::instance().setGlobalKillBelow(ERROR);

  // step 1: find a minimum to start with
  std::unique_ptr<RooFitResult> r1 = nullptr;
  {
    RooFormulaVar ll("ll", "ll", "-2*log(@0)", RooArgSet(*w->pdf(pdfName)));
    // auto r1 = fitToMin(&ll, printlevel);
    RooMinimizer m(ll);
    m.setPrintLevel(-2);
    m.setErrorLevel(4.0);  ///< define 2 sigma errors. This will make the hesse PDF 2 sigma wide!
    int status = m.migrad();
    r1 = std::unique_ptr<RooFitResult>(m.save());
    // if ( 102<RadToDeg(w->var("g")->getVal())&&RadToDeg(w->var("g")->getVal())<103 )
    // {
    //   cout << "step 1" << endl;
    //   r1->Print("v");
    //   gStyle->SetPalette(1);
    //   double xmin = 0.;
    //   double xmax = 3.14;
    //   double ymin = 0.;
    //   double ymax = 0.2;
    //   auto histo = new TH2F("histo", "histo", 100, xmin, xmax, 100, ymin, ymax);
    //   for ( int ix=0; ix<100; ix++ )
    //   for ( int iy=0; iy<100; iy++ )
    //   {
    //     double x = xmin + (xmax-xmin)*(double)ix/(double)100;
    //     double y = ymin + (ymax-ymin)*(double)iy/(double)100;
    //     w->var("d_dk")->setVal(x);
    //     w->var("r_dk")->setVal(y);
    //     histo->SetBinContent(ix+1,iy+1,ll.getVal());
    //   }
    //   newNoWarnTCanvas("c2");
    //   histo->GetZaxis()->SetRangeUser(0,20);
    //   histo->Draw("colz");
    //   setParameters(w, parsName, r1);
    // }
  }

  // step 2: build and fit the improved fcn
  std::unique_ptr<RooFitResult> r2 = nullptr;
  {
    // create a Hesse PDF, import both PDFs into a new workspace,
    // so that their parameters are linked
    RooAbsPdf* hessePdf = r1->createHessePdf(*w->set(parsName));
    if (!hessePdf) return r1;
    auto wImprove = std::make_unique<RooWorkspace>();
    wImprove->import(*w->pdf(pdfName));
    wImprove->import(*hessePdf);
    hessePdf = wImprove->pdf(hessePdf->GetName());
    RooAbsPdf* fullPdf = wImprove->pdf(pdfName);

    RooFormulaVar ll("ll", "ll", "-2*log(@0) +16*@1", RooArgSet(*fullPdf, *hessePdf));
    // auto r2 = fitToMin(&ll, printlevel);
    RooMinimizer m(ll);
    m.setPrintLevel(-2);
    m.setErrorLevel(1.0);
    int status = m.migrad();
    r2 = std::unique_ptr<RooFitResult>(m.save());

    // if ( 102<RadToDeg(w->var("g")->getVal())&&RadToDeg(w->var("g")->getVal())<103 )
    // {
    //   cout << "step 3" << endl;
    //   r2->Print("v");
    //
    //   gStyle->SetPalette(1);
    //   double xmin = 0.;
    //   double xmax = 3.14;
    //   double ymin = 0.;
    //   double ymax = 0.2;
    //   auto histo = new TH2F("histo", "histo", 100, xmin, xmax, 100, ymin, ymax);
    //   for ( int ix=0; ix<100; ix++ )
    //   for ( int iy=0; iy<100; iy++ )
    //   {
    //     double x = xmin + (xmax-xmin)*(double)ix/(double)100;
    //     double y = ymin + (ymax-ymin)*(double)iy/(double)100;
    //     wImprove->var("d_dk")->setVal(x);
    //     wImprove->var("r_dk")->setVal(y);
    //     histo->SetBinContent(ix+1,iy+1,ll.getVal());
    //   }
    //   newNoWarnTCanvas("c7");
    //   histo->GetZaxis()->SetRangeUser(0,20);
    //   histo->Draw("colz");
    //   // setParameters(wImprove, parsName, r1);
    // }
  }

  // step 3: use the fit result of the improved fit as
  // start parameters for the nominal fcn
  std::unique_ptr<RooFitResult> r3 = nullptr;
  {
    setParameters(w, parsName, r2.get());
    RooFormulaVar ll("ll", "ll", "-2*log(@0)", RooArgSet(*w->pdf(pdfName)));
    RooMinimizer m(ll);
    m.setPrintLevel(-2);
    m.setErrorLevel(1.0);
    int status = m.migrad();
    r3 = std::unique_ptr<RooFitResult>(m.save());
    // if ( 102<RadToDeg(w->var("g")->getVal())&&RadToDeg(w->var("g")->getVal())<103 )
    // {
    //   cout << "step 3" << endl;
    //   r3->Print("v");
    // }
  }

  // step 5: chose better minimum
  // cout << r1->minNll() << " " << r3->minNll() << endl;
  std::unique_ptr<RooFitResult> r;
  if (r1->minNll() < r3->minNll()) {
    r = std::move(r1);
  } else {
    r = std::move(r3);
    // cout << "Utils::fitToMinImprove() : improved fit is better!" << endl;
  }

  RooMsgService::instance().setGlobalKillBelow(INFO);

  // set to best parameters
  setParameters(w, parsName, r.get());
  return r;
}

double Utils::getChi2(RooAbsPdf* pdf) {
  RooFormulaVar ll("ll", "ll", "-2*log(@0)", RooArgSet(*pdf));
  return ll.getVal();
}

//
// Randomize all parameters of a set defined in a given
// workspace.
//
void Utils::randomizeParameters(RooWorkspace* w, TString setname) {
  for (const auto pAbs : *w->set(setname)) {
    const auto p = static_cast<RooRealVar*>(pAbs);
    if (p->isConstant()) continue;
    // sample from uniform distribution
    p->randomize();
  }
}

//
// Randomize all parameters according to a Gaussian centered on their
// best fit value with a width of their uncertainty
//
void Utils::randomizeParametersGaussian(RooWorkspace* w, TString setname, RooSlimFitResult* r) {
  RooArgList list = r->floatParsFinal();
  for (const auto pAbs : list) {
    const auto p = static_cast<RooRealVar*>(pAbs);
    auto var = (RooRealVar*)w->var(p->GetName());
    if (w->set(setname)) {
      if (!w->set(setname)->find(p->GetName())) continue;
    }
    if (p->isConstant()) {
      var->setVal(p->getVal());
    } else {
      double randnumb = RooRandom::randomGenerator()->Gaus(p->getVal(), p->getError());
      // make sure it's in range (if not then throw again)
      while (randnumb < var->getMin() || randnumb > var->getMax()) {
        randnumb = RooRandom::randomGenerator()->Gaus(p->getVal(), p->getError());
      }
      var->setVal(randnumb);
    }
  }
}

//
// Randomize all parameters according to a flat distribution within
// some number of sigma of the best fit value
//
void Utils::randomizeParametersUniform(RooWorkspace* w, TString setname, RooSlimFitResult* r, double sigmaRange) {
  RooArgList list = r->floatParsFinal();
  for (const auto pAbs : list) {
    const auto p = static_cast<RooRealVar*>(pAbs);
    auto var = (RooRealVar*)w->var(p->GetName());
    if (p->isConstant()) {
      var->setVal(p->getVal());
    } else {
      double randnumb = RooRandom::randomGenerator()->Uniform(p->getVal() - sigmaRange * p->getError(),
                                                              p->getVal() + sigmaRange * p->getError());
      // make sure it's in range (if not then throw again)
      while (randnumb < var->getMin() || randnumb > var->getMax()) {
        randnumb = RooRandom::randomGenerator()->Gaus(p->getVal(), p->getError());
      }
      var->setVal(randnumb);
    }
  }
}

/// Replace all occurrences of a substring in a string.
std::string Utils::replaceAll(const std::string& input, const std::string& toReplace, const std::string& replaceWith) {
  std::string output = input;
  if (toReplace == replaceWith) return input;
  size_t pos = 0;
  size_t start_pos = 0;
  do {
    pos = output.find(toReplace, start_pos);
    if (pos != std::string::npos) {
      output.replace(pos, toReplace.length(), replaceWith);
      start_pos = pos + toReplace.length();
    }
  } while (pos != std::string::npos);
  return output;
}

///
/// Set each parameter in workspace to the values found
/// in the fit result
///
void Utils::setParameters(RooWorkspace* w, const RooFitResult* values) {
  RooArgList list = values->floatParsFinal();
  list.add(values->constPars());
  for (const auto pAbs : list) {
    const auto p = static_cast<RooRealVar*>(pAbs);
    auto var = dynamic_cast<RooRealVar*>(w->allVars().find(p->GetName()));
    if (!(var)) {
      std::cout << "WARNING in Utils::setParameters(RooWorkspace,RooFitResult) -- no Var found with name "
                << p->GetName() << " in Workspace!" << endl;
    } else {
      var->setVal(p->getVal());
    }
  }
  return;
};

///
/// Set each parameter in setMe to the value found in values.
/// Do nothing if parameter is not found in values.
///
void Utils::setParameters(const RooAbsCollection* setMe, const RooAbsCollection* values) {
  for (const auto pAbs : *setMe) {
    const auto p = static_cast<RooRealVar*>(pAbs);
    auto var = (RooRealVar*)values->find(p->GetName());
    if (var) p->setVal(var->getVal());
  }
}

///
/// Set each floating parameter in setMe to the value found in values.
/// Do nothing if parameter is not found in values.
///
void Utils::setParametersFloating(const RooAbsCollection* setMe, const RooAbsCollection* values) {
  for (const auto pAbs : *setMe) {
    const auto p = static_cast<RooRealVar*>(pAbs);
    if (p->isConstant()) continue;
    auto var = (RooRealVar*)values->find(p->GetName());
    if (var) p->setVal(var->getVal());
  }
}

///
/// Set each parameter in the named set parname inside workspace w
/// to the value found in set. Do nothing if a parameter is present
/// in the parname set, but not found in set.
/// \param w workspace containing a parameter set of name parname
/// \param parname Name of the parameter set containing the "destination" parameters.
/// \param set parameter set holding the "from" parameters.
///
void Utils::setParameters(RooWorkspace* w, TString parname, const RooAbsCollection* set) {
  if (!w->set(parname)) {
    cout << "Utils::setParameters() : ERROR : set not found in workspace: " << parname << endl;
    assert(0);
  }
  setParameters(w->set(parname), set);
}

///
/// Set each floating parameter in the named set parname inside workspace w
/// to the value found in set. Do nothing if a parameter is present
/// in the parname set, but not found in set.
///
/// \param w workspace containing a parameter set of name parname
/// \param parname Name of the parameter set containing the "destination" parameters.
/// \param set parameter set holding the "from" parameters.
///
void Utils::setParametersFloating(RooWorkspace* w, TString parname, const RooAbsCollection* set) {
  setParametersFloating(w->set(parname), set);
}

///
/// Set each parameter in the named set parname inside workspace w
/// to the value found in the final set of floating (or floating and constant)
/// fit parameters in r.
///
/// \param w - workspace containing a parameter set of name parname
/// \param parname - Name of the parameter set containing the "destination" parameters.
/// \param r - a fit result holding the parameter values to be set
/// \param constAndFloat - If set to true, parameter values will be copied from both
///                        constant and floating fit parameters in the RooFitResult.
///                        Default is false.
///
void Utils::setParameters(RooWorkspace* w, TString parname, const RooFitResult* r, bool constAndFloat) {
  if (constAndFloat) {
    RooArgList list = r->floatParsFinal();
    list.add(r->constPars());
    setParameters(w, parname, &list);
    return;
  }
  setParameters(w, parname, &(r->floatParsFinal()));
}

void Utils::setParameters(RooWorkspace* w, TString parname, const RooSlimFitResult* r, bool constAndFloat) {
  // avoid calls to floatParsFinal on a RooSlimFitResult - errgh!
  const auto names = r->getParsNames();
  for (int i = 0; i < names.size(); i++) {
    auto var = (RooRealVar*)w->var(names[i].c_str());
    if (var) var->setVal(r->getParVal(i));
  }
}

///
/// Set each floating parameter in the named set parname inside workspace w
/// to the value found in the final set of floating fit parameters
/// in r.
///
void Utils::setParametersFloating(RooWorkspace* w, TString parname, const RooFitResult* r) {
  setParametersFloating(w, parname, &(r->floatParsFinal()));
}

///
/// Set each parameter in the named set parname inside workspace w
/// to the value found in the first row of the provided dataset.
///
void Utils::setParameters(RooWorkspace* w, TString parname, const RooDataSet* d) {
  setParameters(w, parname, d->get(0));
}

///
/// Set each floating parameter in the named set parname inside workspace w
/// to the value found in the first row of the provided dataset.
///
void Utils::setParametersFloating(RooWorkspace* w, TString parname, const RooDataSet* d) {
  setParametersFloating(w, parname, d->get(0));
}

///
/// Fix each parameter in the named set "parname" inside workspace "w".
///
void Utils::fixParameters(RooWorkspace* w, TString parname) {
  if (w->set(parname)) fixParameters(w->set(parname));
}

void Utils::fixParameters(const RooAbsCollection* set) {
  for (const auto pAbs : *set) {
    const auto p = static_cast<RooRealVar*>(pAbs);
    p->setConstant(true);
  }
}

///
/// Float each parameter in the named set "parname" inside workspace "w".
///
void Utils::floatParameters(RooWorkspace* w, TString parname) {
  for (const auto pAbs : *w->set(parname)) {
    const auto p = static_cast<RooRealVar*>(pAbs);
    p->setConstant(false);
  }
}

void Utils::floatParameters(const RooAbsCollection* set) {
  for (const auto pAbs : *set) {
    const auto p = static_cast<RooRealVar*>(pAbs);
    p->setConstant(false);
  }
}

///
/// Load a named parameter range for a certain parameter.
///
/// \param v - The parameter which will get the limit set.
/// \param limitname - Name of the limit to set.
///
void Utils::setLimit(RooRealVar* v, TString limitname) {
  RooMsgService::instance().setGlobalKillBelow(ERROR);
  v->setRange(v->getMin(limitname), v->getMax(limitname));
  RooMsgService::instance().setGlobalKillBelow(INFO);
}

///
/// Load a named parameter range for a certain parameter,
/// which is found inside a workspace.
///
/// \param w - The workspace holding the parameter.
/// \param parname - The name of the parameter.
/// \param limitname - Name of the limit to set.
///
void Utils::setLimit(RooWorkspace* w, TString parname, TString limitname) {
  RooMsgService::instance().setGlobalKillBelow(ERROR);
  w->var(parname)->setRange(w->var(parname)->getMin(limitname), w->var(parname)->getMax(limitname));
  RooMsgService::instance().setGlobalKillBelow(INFO);
}

///
/// Load a named parameter range for a list of parameters.
///
/// \param set - The list holding the parameters.
/// \param limitname - Name of the limit to set.
///
void Utils::setLimit(const RooAbsCollection* set, TString limitname) {
  RooMsgService::instance().setGlobalKillBelow(ERROR);
  for (const auto pAbs : *set) {
    const auto p = static_cast<RooRealVar*>(pAbs);
    p->setRange(p->getMin(limitname), p->getMax(limitname));
  }
  RooMsgService::instance().setGlobalKillBelow(INFO);
}

///
/// Build a full correlation matrix.
///
/// First, check that the input matrix has the right format.
/// Allowed formats are:
///
///   1)  1.  x   y
///       x   1.  z
///       y   z   1.
///
///   2)  1.  x   y
///       0.  1.  z
///       0.  0.  1.
///
///   3)  1.  0.  0.
///       x   1.  0.
///       y   z   1.
///
///   4)  0.  x   y
///       0.  0.  z
///       0.  0.  0.
///
///   5)  0.  0.  0.
///       x   0.  0.
///       y   z   0.
///
/// If none of the formats are matched, an error is raised and the program
/// execution is stopped.
///
/// In cases 2 to 5, the off-diagonal terms in the lower part of the matrix
/// are set to the value of their symmetric element.
/// In cases 4 and 5, the diagonal terms are set to 1.
///
/// @return `true` if the input matrix was well defined, otherwise `false`.
///
bool Utils::buildCorMatrix(TMatrixDSym& cor) {
  const int n = cor.GetNcols();
  const double tol = 1e-6;  // tolerance for double precision
  auto ill_formed = [cor](const std::string reason) {
    std::cerr << "FAILURE: Input correlation matrix is ill-formed:\n"
                 "         "
              << reason << std::endl;
    cor.Print();
  };

  // --- Check the diagonal
  const bool first_el_one = std::abs(cor[0][0] - 1.) < tol;
  const bool first_el_zero = std::abs(cor[0][0]) < tol;
  if (!(first_el_one || first_el_zero)) {
    ill_formed("first element is not 0. nor 1.");
    return false;
  }
  for (int i = 1; i < n; ++i) {
    if (first_el_one && std::abs(cor[i][i] - 1.) > tol) {
      ill_formed("First element is 1, but diagonal has elements which differ from 1 (index " + std::to_string(i) + ")");
      return false;
    }
    if (first_el_zero && std::abs(cor[i][i]) > tol) {
      ill_formed("First element is 0, but diagonal has elements which differ from 0 (index " + std::to_string(i) + ")");
      return false;
    }
  }
  // Fill the diagonal
  for (int i = 0; i < n; ++i) cor[i][i] = 1.;

  // --- Check the off-diagonal elements
  // Check if at least one of the two off-diagonal triangles is zero
  bool upper_triangle_zero = true;
  bool lower_triangle_zero = true;
  for (int i = 0; i < n - 1; ++i) {
    for (int j = i + 1; j < n; ++j) {
      if (std::abs(cor[i][j]) > tol) upper_triangle_zero = false;
      if (std::abs(cor[j][i]) > tol) lower_triangle_zero = false;
    }
  }
  // If the diagonal is zero, at least one of the two triangles should be zero
  if (first_el_zero && !(upper_triangle_zero || lower_triangle_zero)) {
    ill_formed("Diagonal elements are 0, but the off-diagonal elements in the lower (or upper) part of the matrix are "
               "not all 0");
    return false;
  }
  // If both triangles are nonzero, the matrix should be symmetric
  if (!upper_triangle_zero && !lower_triangle_zero) {
    for (int i = 0; i < n - 1; ++i) {
      for (int j = i + 1; j < n; ++j) {
        if (std::abs(cor[i][j] - cor[j][i]) > tol) {
          ill_formed("The matrix is not symmetric (element " + std::to_string(i) + "," + std::to_string(j) + ")");
          return false;
        }
      }
    }
  }
  // Fill the zero triangle
  for (int i = 0; i < n - 1; ++i) {
    for (int j = i + 1; j < n; ++j) {
      if (std::abs(cor[i][j]) < tol)
        cor[i][j] = cor[j][i];
      else
        cor[j][i] = cor[i][j];
    }
  }
  return true;
}

///
/// Build a covariance matrix
/// from a correlation matrix and error vectors.
///
TMatrixDSym* Utils::buildCovMatrix(TMatrixDSym& cor, double* err) {
  int n = cor.GetNcols();
  TMatrixDSym cov(n);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) { cov[i][j] = err[i] * cor[i][j] * err[j]; }
  return new TMatrixDSym(cov);
}

///
/// Build a covariance matrix
/// from a correlation matrix and error vectors.
///
TMatrixDSym* Utils::buildCovMatrix(TMatrixDSym& cor, vector<double>& err) {
  int n = cor.GetNcols();
  TMatrixDSym cov(n);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) { cov[i][j] = err[i] * cor[i][j] * err[j]; }
  return new TMatrixDSym(cov);
}

///
/// Make a theory var in which the parmaeter list gets slimmed down to only
/// contain the relevant dependents.
/// This is a workaround for changes in RooFormulaVar that break things when
/// making subset pdfs
RooFormulaVar* Utils::makeTheoryVar(TString name, TString title, TString formula, RooArgList* pars) {
  auto explicitDependents = new RooArgList();
  for (int i = 0; i < pars->getSize(); i++) {
    if (formula.Contains(TString(pars->at(i)->GetName()))) { explicitDependents->add(*(pars->at(i))); }
  }
  return new RooFormulaVar(name, title, formula, *explicitDependents);
}

void Utils::addSetNamesToList(vector<string>& list, RooWorkspace* w, TString setName) {
  for (const auto p : *w->set(setName)) list.push_back(p->GetName());
}

///
/// Make a named set from a list of strings of object names
/// Duplicates will only be contained once
void Utils::makeNamedSet(RooWorkspace* w, TString mergedSet, vector<string>& names) {

  // 1. remove duplicates
  std::ranges::sort(names);
  vector<string> vars;
  vars.push_back(names[0]);
  string previous = names[0];
  for (int i = 1; i < names.size(); i++) {
    if (previous == names[i]) continue;
    vars.push_back(names[i]);
    previous = names[i];
  }

  // 2. make new, combined set on the workspace
  TString varsCommaList = "";
  for (int i = 0; i < vars.size(); i++) {
    varsCommaList.Append(vars[i]);
    if (i < vars.size() - 1) varsCommaList.Append(",");
  }
  w->defineSet(mergedSet, varsCommaList);
}

///
/// Merge two named sets of variables inside a RooWorkspace.
/// Duplicate variables will only be contained once.
///
void Utils::mergeNamedSets(RooWorkspace* w, TString mergedSet, TString set1, TString set2) {
  // 1. fill all variables into a vector
  vector<string> varsAll;
  for (const auto p : *w->set(set1)) varsAll.push_back(p->GetName());
  for (const auto p : *w->set(set2)) varsAll.push_back(p->GetName());

  // 2. remove duplicates
  std::ranges::sort(varsAll);
  vector<string> vars;
  vars.push_back(varsAll[0]);
  string previous = varsAll[0];
  for (int i = 1; i < varsAll.size(); i++) {
    if (previous == varsAll[i]) continue;
    vars.push_back(varsAll[i]);
    previous = varsAll[i];
  }

  // 3. make new, combined set on the workspace
  TString varsCommaList = "";
  for (int i = 0; i < vars.size(); i++) {
    varsCommaList.Append(vars[i]);
    if (i < vars.size() - 1) varsCommaList.Append(",");
  }
  w->defineSet(mergedSet, varsCommaList);
}

/*
 * doesn't work with 3GB files.
 *
 */
bool Utils::FileExists(TString strFilename) {
  struct stat stFileInfo;
  bool blnReturn;
  int intStat;
  // Attempt to get the file attributes
  intStat = stat(strFilename, &stFileInfo);
  if (intStat == 0) {
    // We were able to get the file attributes
    // so the file obviously exists.
    blnReturn = true;
  } else {
    // We were not able to get the file attributes.
    // This may mean that we don't have permission to
    // access the folder which contains this file. If you
    // need to do that level of checking, lookup the
    // return values of stat which will give you
    // more details on why stat failed.
    blnReturn = false;
  }
  return (blnReturn);
}

void Utils::savePlot(const TCanvas* c1, const TString name, const std::vector<std::string> extensions) {
  cout << "saving plot (pdf and other formats) to: plots/pdf/" + name + ".pdf" << endl;
  gErrorIgnoreLevel = kWarning;
  for (const auto ext : extensions) c1->Print("plots/" + ext + "/" + name + "." + ext);
  gErrorIgnoreLevel = kInfo;
}

///
/// Round a number to a certain number of
/// decimal points.
///
double Utils::Round(double value, int digits) { return TString(Form("%.*f", digits, value)).Atof(); }

///
/// Compute number of digits needed behind the decimal
/// point to achieve a certain number of significant digits.
/// The result can be used in the printf() function.
/// Examples:
///   23.243: returns 0 with sigdigits=2
///   2.2634: returns 1 with sigdigits=2
///
int Utils::calcNsubdigits(double value, int sigdigits) {
  if (value == 0) return 1;
  double myvalue = value;
  int count = 0;
  for (; fabs(value) < pow(10., sigdigits - 1); value *= 10) {
    if (count == 10) break;
    count++;
  }

  // do it again to catch a rounding issue: 9.97 at 2 digit
  // precision gives count=1 above, but we want count=0 to
  // get "10" instead of "10.0"
  value = TString(Form("%.*f", count, myvalue)).Atof();
  if (value == 0) return 1;
  count = 0;
  for (; fabs(value) < pow(10., sigdigits - 1); value *= 10) {
    if (count == 10) break;
    count++;
  }
  return count;
}

///
/// Converts a RooDataSet to a TTree which then can be
/// browsed.
///
TTree* Utils::convertRooDatasetToTTree(RooDataSet* d) {
  // set up the TTree based on the content of the first
  // row of the dataset
  map<string, double> variables;  ///< the proxy variables
  auto t = new TTree("tree", "tree");
  for (const auto pAbs : *d->get(0)) {
    const auto p = static_cast<RooRealVar*>(pAbs);
    variables.insert(pair<string, double>(p->GetName(), p->getVal()));
    t->Branch(p->GetName(), &variables[p->GetName()], TString(p->GetName()) + "/F");
  }

  // loop over the dataset, filling the tree
  int nEntries = d->sumEntries();
  for (int i = 0; i < nEntries; i++) {
    for (const auto pAbs : *d->get(i)) {
      const auto p = static_cast<RooRealVar*>(pAbs);
      variables[p->GetName()] = p->getVal();
    }
    t->Fill();
  }

  return t;
}

/// Convert a TH1 to a TGraph.
std::unique_ptr<TGraph> Utils::convertTH1ToTGraph(TH1* h, bool withErrors) {
  std::unique_ptr<TGraph> g = nullptr;
  if (withErrors)
    g = std::make_unique<TGraphErrors>(h->GetNbinsX());
  else
    g = std::make_unique<TGraph>(h->GetNbinsX());
  g->SetName(getUniqueRootName());
  for (int i = 0; i < h->GetNbinsX(); i++) {
    g->SetPoint(i, h->GetBinCenter(i + 1), h->GetBinContent(i + 1));
    if (withErrors) (static_cast<TGraphErrors*>(g.get()))->SetPointError(i, 0., h->GetBinError(i + 1));
  }
  return g;
}

/// Smooths a graph
TGraph* Utils::smoothGraph(TGraph* g, int option) {
  auto smoother = TGraphSmooth();
  TGraph* gr;
  if (option == 0)
    gr = (TGraph*)smoother.SmoothSuper(g)->Clone(Form("sm%s", g->GetName()));
  else if (option == 1)
    gr = (TGraph*)smoother.Approx(g)->Clone(Form("sm%s", g->GetName()));
  else {
    cout << "Utils::smoothGraph() : ERROR - no such option " << option << endl;
    exit(1);
  }
  return gr;
}

// Smooths a histogram
TGraph* Utils::smoothHist(TH1* h, int option) {
  auto g = convertTH1ToTGraph(h);
  return smoothGraph(g.get());
}

/**
 * Adds a point to a TGraph at the first position where the x value is larger than the new x value.
 *
 * If no such position is found, the point is added at the end.
 */
std::unique_ptr<TGraph> Utils::addPointToGraphAtFirstMatchingX(const TGraph* g, const double xNew, const double yNew) {
  // Get x and y coordinates as vectors- the TGraph interface is just not suited to what we want to do
  std::vector<double> xVec;
  std::vector<double> yVec;
  for (int i = 0; i < g->GetN(); ++i) {
    double xOld, yOld;
    g->GetPoint(i, xOld, yOld);
    xVec.push_back(xOld);
    yVec.push_back(yOld);
  }

  const auto it = std::ranges::find_if(xVec, [xNew](auto x) { return x >= xNew; });
  const int iPos = it - xVec.begin();
  xVec.insert(xVec.begin() + iPos, xNew);
  yVec.insert(yVec.begin() + iPos, yNew);

  // create a new graph of the right kind
  const bool isTGraphErrors = TString(g->ClassName()).EqualTo("TGraphErrors");
  std::unique_ptr<TGraph> gNew;
  if (isTGraphErrors)
    gNew = std::make_unique<TGraphErrors>(g->GetN() + 1);
  else
    gNew = std::make_unique<TGraph>(g->GetN() + 1);
  gNew->SetName(getUniqueRootName());
  for (int i = 0; i < xVec.size(); ++i) { gNew->SetPoint(i, xVec[i], yVec[i]); }
  if (isTGraphErrors) {
    for (int i = 0; i < xVec.size(); ++i) {
      double yErr = (i == iPos) ? 0. : g->GetErrorY(i - (i < iPos ? 0 : 1));
      static_cast<TGraphErrors*>(gNew.get())->SetPointError(i, 0., yErr);
    }
  }
  return gNew;
}
///
/// Creates a fresh, independent copy of the input histogram.
/// We cannot use Root's Clone() or the like, because that
/// crap always copies Draw options and whatnot.
///
/// \param h - the input histogram
/// \param copyContent - true: also copy content. false: initialize with zeroes
/// \param uniqueName - true: append a unique string to the histogram name
/// \return a new histogram. Caller assumes ownership.
///
std::unique_ptr<TH1> Utils::histHardCopy(const TH1* h, bool copyContent, bool uniqueName, TString specName) {
  TString name = h->GetTitle();
  if (specName != "") name = specName;
  if (uniqueName) name += getUniqueRootName();
  auto hNew =
      std::make_unique<TH1F>(name, h->GetTitle(), h->GetNbinsX(), h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
  hNew->SetDirectory(0);
  for (int l = 1; l <= h->GetNbinsX(); l++) {
    if (copyContent)
      hNew->SetBinContent(l, h->GetBinContent(l));
    else
      hNew->SetBinContent(l, 0);
  }
  return hNew;
}

///
/// Creates a fresh, independent copy of the input histogram.
/// 2d version of TH1F* Utils::histHardCopy().
///
std::unique_ptr<TH2> Utils::histHardCopy(const TH2* h, bool copyContent, bool uniqueName, TString specName) {
  TString name = h->GetTitle();
  if (specName != "") name = specName;
  if (uniqueName) name += getUniqueRootName();
  auto hNew =
      std::make_unique<TH2F>(name, h->GetTitle(), h->GetNbinsX(), h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax(),
                             h->GetNbinsY(), h->GetYaxis()->GetXmin(), h->GetYaxis()->GetXmax());
  hNew->SetDirectory(0);
  for (int k = 1; k <= h->GetNbinsX(); k++)
    for (int l = 1; l <= h->GetNbinsY(); l++) {
      if (copyContent)
        hNew->SetBinContent(k, l, h->GetBinContent(k, l));
      else
        hNew->SetBinContent(k, l, 0);
    }
  return hNew;
}

///
/// Check if a matrix is positive definite.
///
bool Utils::isPosDef(TMatrixDSym* c) {
  TMatrixDSymEigen eigen(*c);
  TVectorD eigenvalues = eigen.GetEigenValues();
  Double_t minEigenVal = eigenvalues[0];
  for (int k = 0; k < c->GetNcols(); k++) minEigenVal = TMath::Min(minEigenVal, eigenvalues[k]);
  if (minEigenVal < 0) {
    cout << "isPosDef() : ERROR : Matrix not pos. def." << endl;
    return false;
  }
  return true;
}

///
/// Check if a RooRealVar is an angle.
///
bool Utils::isAngle(RooRealVar* v) { return v->getUnit() == TString("Rad") || v->getUnit() == TString("rad"); }

int Utils::makeNewColor(string hex) {
  int ci = TColor::GetFreeColorIndex();
  int ri, gi, bi;
  sscanf(hex.c_str(), "#%02x%02x%02x", &ri, &gi, &bi);
  double r = double(ri) / 255.;
  double g = double(gi) / 255.;
  double b = double(bi) / 255.;
  TColor col(ci, r, g, b);
  cout << ci << " " << hex << " " << r << " " << g << " " << b << endl;
  return col.GetNumber();
}

///
/// function filling a RooArgList with parameters within a Workspace
/// parameter names given by an vector with TStrings
///
void Utils::fillArgList(RooArgList* list, RooWorkspace* w, std::vector<TString> names) {
  for (auto name : names) {
    if (!list->add(*w->var(name), kTRUE)) {  //> add silent
      std::cout << "WARNING: Utils::fillArgList - Var either already in ArgList: " << list->GetName()
                << " or the List does own its vars" << endl;
    };
  }
}

///
/// Fills vector with floating pars names
///
void Utils::getParameters(const RooFitResult& result, std::vector<TString>& names) {
  RooArgList pars = result.floatParsFinal();
  for (const auto p : pars) { names.push_back(TString(p->GetName())); }
};

///
/// finds a variable within a RooArgSet with a specific name sub-string
/// return all variables containing the sub-string in a std::vector
///
std::vector<TString> Utils::getParsWithName(const TString& subString, const RooArgSet& set) {
  std::string vars = set.contentsString();
  std::vector<std::string> _names;
  std::vector<TString> _results;
  boost::split(_names, vars, boost::is_any_of(","));
  for (std::string str : _names) {
    std::size_t found = str.find(subString.Data());
    if (found != std::string::npos) { _results.push_back(TString(str)); }
  }
  return _results;
}
///
/// searches for unphysical values in a fit result, fills problematic var names in vectors
/// returns true if everything is within boundaries
///
bool Utils::checkBoundary(const RooSlimFitResult& r, std::vector<TString> lowProb, std::vector<TString> highProb) {
  // make sure vectors are empty
  lowProb.clear();
  highProb.clear();
  RooArgList floats = r.floatParsFinal();
  for (const auto varAbs : floats) {
    const auto var = static_cast<RooRealVar*>(varAbs);
    // check lower bound
    if (var->getVal() < var->getMin("phys")) { lowProb.push_back(var->GetName()); }
    // check higher bound
    if (var->getVal() > var->getMax("phys")) { highProb.push_back(var->GetName()); }
  }
  return (lowProb.size() != 0 || highProb.size() != 0) ? false : true;
}

void Utils::setParsConstToBound(RooWorkspace* w, std::vector<TString> namesLow, std::vector<TString> namesHigh) {
  setParsConstToBound(w, namesLow, true);
  setParsConstToBound(w, namesHigh, false);
};

void Utils::setParsConstToBound(RooWorkspace* w, std::vector<TString> names, bool low) {
  if (names.size() == 0) return;
  for (auto& n : names) {
    double_t valueToSet = (low) ? w->var(n)->getMin("phys") : w->var(n)->getMax("phys");
    w->var(n)->setVal(valueToSet);
    w->var(n)->setConstant(kTRUE);
  }
}

void Utils::setParametersFloating(RooWorkspace* w, const std::vector<TString>& names,
                                  const std::vector<TString>& names2) {
  setParametersFloating(w, names);
  setParametersFloating(w, names2);
};

void Utils::setParametersFloating(RooWorkspace* w, const std::vector<TString>& names) {
  if (names.size() == 0) return;
  for (auto& n : names) { w->var(n)->setConstant(kFALSE); }
};

///
/// Debug tools: print the content of a vector to stdout.
///
void Utils::dump_vector(const std::vector<int>& l) {
  for (auto const& el : l) { cout << el << endl; }
}
void Utils::dump_vector(const std::vector<double>& l) {
  for (auto const& el : l) { cout << el << endl; }
}
void Utils::dump_matrix(const std::vector<std::vector<int>>& l) {
  for (int ix = 0; ix < l.size(); ix++) {
    for (int iy = 0; iy < l[0].size(); iy++) { cout << printf("%5i", l[ix][iy]) << " "; }
    cout << endl;
  }
}

///
/// Debug tools: print the content of a 2d map to stdout.
///
void Utils::dump_map(const std::map<int, std::vector<int>>& map) {
  for (auto const& [key, val] : map) {
    cout << "Key: " << key << endl;
    cout << "Values" << endl;
    dump_vector(val);
  }
}

std::vector<std::vector<int>> Utils::transpose(std::vector<std::vector<int>>& v) {
  std::vector<std::vector<int>> newVector;
  int oldNx = v.size();
  if (oldNx == 0) {
    cout << "Utils::transpose() : ERROR : x dimension is 0" << endl;
    return newVector;
  }
  int oldNy = v[0].size();
  if (oldNy == 0) {
    cout << "Utils::transpose() : ERROR : y dimension is 0" << endl;
    return newVector;
  }
  // check if rectangular
  for (int j = 1; j < oldNx; j++) {
    if (v[j].size() != oldNy) {
      cout << "Utils::transpose() : ERROR : vector not rectangular" << endl;
      return newVector;
    }
  }
  // transpose
  for (int iy = 0; iy < oldNy; iy++) {
    std::vector<int> tmp;
    for (int ix = 0; ix < oldNx; ix++) { tmp.push_back(v[ix][iy]); }
    newVector.push_back(tmp);
  }
  return newVector;
}

std::unique_ptr<TCanvas> Utils::newNoWarnTCanvas(TString name, TString title, int width, int height) {
  gErrorIgnoreLevel = kError;
  auto c = std::make_unique<TCanvas>(name, title, width, height);
  gErrorIgnoreLevel = kInfo;
  return c;
}

std::unique_ptr<TCanvas> Utils::newNoWarnTCanvas(TString name, TString title, int x, int y, int width, int height) {
  gErrorIgnoreLevel = kError;
  auto c = std::make_unique<TCanvas>(name, title, x, y, width, height);
  gErrorIgnoreLevel = kInfo;
  return c;
}

void Utils::HFAGLabel(const TString& label, Double_t xpos, Double_t ypos, Double_t scale) {
  TVirtualPad* thePad;

  if (!(thePad = TVirtualPad::Pad())) return;

  UInt_t pad_width(thePad->XtoPixel(thePad->GetX2()));
  UInt_t pad_height(thePad->YtoPixel(thePad->GetY1()));

  Double_t ysiz_pixel(25);
  Double_t ysiz(Double_t(ysiz_pixel) / Double_t(pad_height));
  Double_t xsiz(4.8 * ysiz * Double_t(pad_height) / Double_t(pad_width));

  Double_t x1, x2, y1, y2;
  xsiz = scale * xsiz;
  ysiz = scale * ysiz;

  if (xpos >= 0) {
    x1 = xpos;
    x2 = xpos + xsiz;
  } else {
    x1 = 1 + xpos - xsiz;
    x2 = 1 + xpos;
  }

  if (ypos >= 0) {
    y1 = ypos + 0.9 * ysiz;
    y2 = ypos + 0.9 * ysiz + ysiz;
  } else {
    y1 = 1 + ypos - ysiz;
    y2 = 1 + ypos;
  }

  auto tbox1 = new TPaveText(x1, y1, x2, y2, "BRNDC");
  // tbox1->SetLineColor(1);
  // tbox1->SetLineStyle(1);
  // tbox1->SetLineWidth(2);
  tbox1->SetFillColor(kBlack);
  tbox1->SetFillStyle(1001);
  // tbox1->SetBorderSize(1);
  tbox1->SetShadowColor(kWhite);
  tbox1->SetTextColor(kWhite);
  tbox1->SetTextFont(76);
  tbox1->SetTextSize(24 * scale);
  tbox1->SetTextAlign(22);  // center-adjusted and vertically centered
  tbox1->AddText(TString("HFLAV"));
  tbox1->Draw();
  //
  auto tbox2 = new TPaveText(x1, y1 - 0.9 * ysiz, x2, y2 - ysiz, "BRNDC");
  // tbox2->SetLineColor(1);
  // tbox2->SetLineStyle(1);
  // tbox2->SetLineWidth(2);
  tbox2->SetFillColor(kWhite);
  tbox2->SetFillStyle(1001);
  // tbox2->SetBorderSize(1);
  tbox2->SetShadowColor(kWhite);
  tbox2->SetTextColor(kBlack);
  tbox2->SetTextFont(76);
  tbox2->SetTextSize(18 * scale);
  tbox2->SetTextAlign(22);  // center-adjusted and vertically centered
  tbox2->AddText(label);
  tbox2->Draw();
  return;
}

void Utils::assertFileExists(TString strFilename) {
  if (!FileExists(strFilename)) {
    cout << "ERROR : File not found: " + strFilename << endl;
    exit(EXIT_FAILURE);
  }
}

std::vector<double> Utils::computeNormalQuantiles(std::vector<double>& values, int nsigma) {

  // std::sort( values.begin(), values.end() );
  std::vector<double> probs;  // = { TMath::Prob(4,1), TMath::Prob(1,1), 0.5, 1.-TMath::Prob(1,1), 1.-TMath::Prob(4,1)
                              // };
  for (int i = nsigma; i > 0; i--) probs.push_back(TMath::Prob(sq(i), 1));
  probs.push_back(0.5);
  for (int i = 0; i < nsigma; i++) probs.push_back(1. - TMath::Prob(sq(i + 1), 1));

  vector<double> quants(nsigma * 2 + 1);

  TMath::Quantiles(values.size(), probs.size(), &values[0], &quants[0], &probs[0], false);

  std::vector<double> quantiles;
  for (auto const& quant : quants) { quantiles.push_back(quant); }
  return quantiles;
}

double Utils::getCorrelationFactor(const vector<double>& a, const vector<double>& b) {

  assert(a.size() == b.size());

  double mean_a = mean(a);
  double mean_b = mean(b);

  double s = 0.;
  for (int i = 0; i < a.size(); i++) { s += (a[i] - mean_a) * (b[i] - mean_b); }

  return s / (a.size() * stddev(a) * stddev(b));
}
