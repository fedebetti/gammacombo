#include <CLIntervalMaker.h>
#include <CLIntervalPrinter.h>
#include <FileNameBuilder.h>
#include <FitResultCache.h>
#include <MethodAbsScan.h>
#include <OneMinusClPlotAbs.h>
#include <PullPlotter.h>
#include <Utils.h>

#include <algorithm>
#include <array>
#include <cstdlib>
#include <format>
#include <limits>
#include <optional>
#include <vector>

#include <RooRealVar.h>

#include <TDatime.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TSpline.h>
#include <TStyle.h>

using namespace std;
using namespace RooFit;
using namespace Utils;

namespace {
  auto msgBase = [](const std::string& prefix, const std::string& msg, std::ostream& stream = std::cout) {
    auto msgOut = Utils::replaceAll(msg, "\n", "\n" + std::string(prefix.size(), ' '));
    stream << prefix << msgOut << std::endl;
  };

  auto errBase = [](const std::string& prefix, const std::string& msg) {
    msgBase(prefix, msg + ". Exit...", std::cerr);
    exit(1);
  };
}  // namespace

MethodAbsScan::MethodAbsScan(Combiner* c) : MethodAbsScan(c->getArg()) {
  combiner = c;
  w = c->getWorkspace();
  name = c->getName();
  title = c->getTitle();
  pdfName = "pdf_" + combiner->getPdfName();
  obsName = "obs_" + combiner->getPdfName();
  parsName = "par_" + combiner->getPdfName();
  thName = "th_" + combiner->getPdfName();

  // check workspace content
  if (!w->pdf(pdfName)) {
    cout << "MethodAbsScan::MethodAbsScan() : ERROR : not found in workspace : " << pdfName << endl;
    exit(1);
  }
  if (!w->set(obsName)) {
    cout << "MethodAbsScan::MethodAbsScan() : ERROR : not found in workspace : " << obsName << endl;
    exit(1);
  }
  if (!w->set(parsName)) {
    cout << "MethodAbsScan::MethodAbsScan() : ERROR : not found in workspace : " << parsName << endl;
    exit(1);
  }
  if (!w->set(thName)) {
    cout << "MethodAbsScan::MethodAbsScan() : ERROR : not found in workspace : " << thName << endl;
    exit(1);
  }
}

/// Constructor without combiner. This is still needed for the datasets stuff
MethodAbsScan::MethodAbsScan(const OptParser* opt)
    : arg(opt), scanVar1(opt->var[0]), verbose(opt->verbose), nPoints1d(opt->npoints1d), nPoints2dx(opt->npoints2dx),
      nPoints2dy(opt->npoints2dy) {
  if (opt->var.size() > 1) scanVar2 = opt->var[1];
  if (opt->CL.size() > 0) {
    for (auto level : opt->CL) { ConfidenceLevels.push_back(level / 100.); }
  } else {
    ConfidenceLevels.push_back(0.6827);  // 1sigma
    ConfidenceLevels.push_back(0.9545);  // 2sigma
    ConfidenceLevels.push_back(0.9973);  // 3sigma
  }
}

///
/// Try to find global mininum of the PDF.
/// Despite its name this often finds a local minimum. It's merely
/// used as a starting point. When the scans stumbles upon a better
/// minimum, we'll keep that one.
///
/// Resets parameters to the values they had at function call.
/// Save the RooFitResult of the global minimum (or whatever minimum it found...)
/// into the globalMin member.
///
/// \param force If set to true it fits again, even it the fit was already run before.
///
void MethodAbsScan::doInitialFit(bool force) {
  RooMsgService::instance().setGlobalKillBelow(ERROR);
  if (arg->debug) {
    cout << "\n============================================================" << endl;
    cout << "MethodAbsScan::doInitialFit() : MAKE FIRST FIT ..." << endl;
    cout << "MethodAbsScan::doInitialFit() : PDF " << pdfName << endl << endl;
  }
  if (!force && chi2minGlobalFound) {
    if (arg->debug) {
      cout << "MethodAbsScan::doInitialFit() : Already found previously: chi2minGlobal = " << chi2minGlobal << endl;
      cout << "\n============================================================" << endl;
    }
    return;
  }

  // Save parameter values that were active at function call.
  startPars = std::make_unique<RooDataSet>("startPars", "startPars", *w->set(parsName));
  startPars->add(*w->set(parsName));

  // load parameter range
  combiner->loadParameterLimits();

  fixParameters(w, obsName);     ///< fix observables
  floatParameters(w, parsName);  ///< physics parameters need to be floating to find global minimum

  // load again the parameter values that were specified on the command line -
  // loading a set of start parameters might have changed them
  combiner->setParametersConstant();

  // fix parameters we decided to keep constant (to keep a parameter constant
  // add them manually to the workspace set 'const')
  fixParameters(w, "const");

  // check choice of start parameters
  double nsigma = 10.;
  PullPlotter p(this);
  if (p.hasPullsAboveNsigma(nsigma)) {
    cout << "MethodAbsScan::doInitialFit() : WARNING : Chosen start parameter values result in pulls larger\n"
            "                                WARNING : than "
         << nsigma
         << " sigma. Check the values in your\n"
            "                                WARNING : ParametersAbs class!\n"
            "Offending pulls:"
         << endl;
    p.printPulls(nsigma);
    cout << endl;
  }

  // print init parameters
  if (arg->debug) {
    cout << "MethodAbsScan::doInitialFit() : init parameters:" << endl;
    w->set(parsName)->Print("v");
    cout << "MethodAbsScan::doInitialFit() : init pulls:" << endl;
    p.printPulls(0.);
    cout << "MethodAbsScan::doInitialFit() : PDF evaluated at init parameters: ";
    cout << w->pdf(pdfName)->getVal() << endl;
    RooFormulaVar ll("ll", "ll", "-2*log(@0)", RooArgSet(*w->pdf(pdfName)));
    cout << "MethodAbsScan::doInitialFit() : Chi2 at init parameters: ";
    cout << ll.getVal() << endl;
  }

  int quiet = arg->debug ? 1 : -1;
  auto r = fitToMinBringBackAngles(w->pdf(pdfName), true, quiet);
  if (arg->debug) r->Print("v");
  globalMin = std::move(r);
  chi2minGlobal = globalMin->minNll();
  chi2minGlobalFound = true;

  // reset parameters to their values at function call
  setParameters(w, parsName, startPars->get(0));

  if (arg->debug) cout << "============================================================\n" << endl;
  RooMsgService::instance().setGlobalKillBelow(INFO);
}

///
/// Set the global minimum manually.
///
void MethodAbsScan::setChi2minGlobal(double x) {
  chi2minGlobalFound = true;
  chi2minGlobal = x;
}

void MethodAbsScan::initScan() {
  if (arg->debug) cout << "MethodAbsScan::initScan() : initializing ..." << endl;
  if (m_initialized) {
    cout << "MethodAbsScan::initScan() : already initialized." << endl;
    exit(1);
  }

  // Init the 1-CL histograms. Range is taken from the scan range defined in
  // the ParameterAbs class (and derived ones), unless the --scanrange command
  // line argument is set.
  auto par1 = w->var(scanVar1);
  if (!par1) {
    if (arg->debug) cout << "MethodAbsScan::initScan() : ";
    cout << "ERROR : No such scan parameter: " << scanVar1 << endl;
    cout << "        Choose an existing one using: --var par" << endl << endl;
    cout << "  Available parameters:" << endl;
    cout << "  ---------------------" << endl << endl;
    for (const auto name : combiner->getParameterNames()) { cout << "    " << name << endl; }
    cout << endl;
    exit(1);
  }
  if (!m_xrangeset && arg->scanrangeMin != arg->scanrangeMax) { setXscanRange(arg->scanrangeMin, arg->scanrangeMax); }
  setLimit(w, scanVar1, "scan");
  double min1 = par1->getMin();
  double max1 = par1->getMax();
  hCL = std::make_unique<TH1F>("hCL" + getUniqueRootName(), "hCL" + pdfName, nPoints1d, min1, max1);
  hChi2min = std::make_unique<TH1F>("hChi2min" + getUniqueRootName(), "hChi2min" + pdfName, nPoints1d, min1, max1);
  hCLs = std::make_unique<TH1F>("hCLs" + getUniqueRootName(), "hCLs" + pdfName, nPoints1d, min1, max1);
  hCLsFreq = std::make_unique<TH1F>("hCLsFreq" + getUniqueRootName(), "hCLsFreq" + pdfName, nPoints1d, min1, max1);
  hCLsExp = std::make_unique<TH1F>("hCLsExp" + getUniqueRootName(), "hCLsExp" + pdfName, nPoints1d, min1, max1);
  hCLsErr1Up =
      std::make_unique<TH1F>("hCLsErr1Up" + getUniqueRootName(), "hCLsErr1Up" + pdfName, nPoints1d, min1, max1);
  hCLsErr1Dn =
      std::make_unique<TH1F>("hCLsErr1Dn" + getUniqueRootName(), "hCLsErr1Dn" + pdfName, nPoints1d, min1, max1);
  hCLsErr2Up =
      std::make_unique<TH1F>("hCLsErr2Up" + getUniqueRootName(), "hCLsErr2Up" + pdfName, nPoints1d, min1, max1);
  hCLsErr2Dn =
      std::make_unique<TH1F>("hCLsErr2Dn" + getUniqueRootName(), "hCLsErr2Dn" + pdfName, nPoints1d, min1, max1);

  // fill the chi2 histogram with very unlikely values such
  // that inside scan1d() the if clauses work correctly
  for (int i = 1; i <= nPoints1d; i++) hChi2min->SetBinContent(i, 1e6);

  if (scanVar2 != "") {
    RooRealVar* par2 = w->var(scanVar2);
    if (!par2) {
      if (arg->debug) cout << "MethodAbsScan::initScan() : ";
      cout << "ERROR : No such scan parameter: " << scanVar2 << endl;
      cout << "        Choose an existing one using: --var par" << endl << endl;
      cout << "  Available parameters:" << endl;
      cout << "  ---------------------" << endl << endl;
      for (const auto name : combiner->getParameterNames()) { cout << "    " << name << endl; }
      cout << endl;
      exit(1);
    }
    if (!m_yrangeset && arg->scanrangeyMin != arg->scanrangeyMax) {
      setYscanRange(arg->scanrangeyMin, arg->scanrangeyMax);
    }
    setLimit(w, scanVar2, "scan");
    double min2 = par2->getMin();
    double max2 = par2->getMax();
    hCL2d = std::make_unique<TH2F>("hCL2d" + getUniqueRootName(), "hCL2d" + pdfName, nPoints2dx, min1, max1, nPoints2dy,
                                   min2, max2);
    hChi2min2d = std::make_unique<TH2F>("hChi2min2d" + getUniqueRootName(), "hChi2min", nPoints2dx, min1, max1,
                                        nPoints2dy, min2, max2);
    for (int i = 1; i <= nPoints2dx; i++)
      for (int j = 1; j <= nPoints2dy; j++) hChi2min2d->SetBinContent(i, j, 1e6);
  }

  // Set up storage for the fit results.
  // Clear before so we can call initScan() multiple times.
  // Note that allResults still needs to hold all results, so don't delete the RooFitResults.

  // 1d:
  curveResults.clear();
  curveResults.resize(nPoints1d, nullptr);

  // 2d:
  curveResults2d.clear();
  curveResults2d.resize(nPoints2dx, vector<RooSlimFitResult*>(nPoints2dy, nullptr));

  // global minimum
  doInitialFit();

  // turn off some messages
  RooMsgService::instance().setStreamStatus(0, kFALSE);
  RooMsgService::instance().setStreamStatus(1, kFALSE);
  m_initialized = true;
}

///
/// Save this scanner to a root file placed into plots/scanner.
/// It contains the 1-CL histograms and the solutions.
///
void MethodAbsScan::saveScanner(TString fName) const {
  auto error = [](const std::string& msg) { errBase("MethodAbsScan::saveScanner : ERROR : ", msg); };

  if (fName == "") {
    FileNameBuilder fb(arg);
    fName = fb.getFileNameScanner(this);
  }
  if (arg->debug) cout << "MethodAbsScan::saveScanner() : saving scanner to file " << fName << endl;
  TFile f(fName, "recreate");
  if (f.IsZombie()) error(std::format("Could not open file {:s}", std::string(fName)));
  // save 1-CL histograms
  if (scanVar2 != "") {
    hCL2d->Write("hCL");
    if (hCLs2d) hCLs2d->Write("hCLs");
  } else {
    hCL->Write("hCL");
    if (hCLs) hCLs->Write("hCLs");
    if (hCLsFreq) hCLsFreq->Write("hCLsFreq");
    if (hCLsExp) hCLsExp->Write("hCLsExp");
    if (hCLsErr1Up) hCLsErr1Up->Write("hCLsErr1Up");
    if (hCLsErr1Dn) hCLsErr1Dn->Write("hCLsErr1Dn");
    if (hCLsErr2Up) hCLsErr2Up->Write("hCLsErr2Up");
    if (hCLsErr2Dn) hCLsErr2Dn->Write("hCLsErr2Dn");
  }
  // save chi2 histograms
  if (scanVar2 != "")
    hChi2min2d->Write("hChi2min");
  else
    hChi2min->Write("hChi2min");
  // save solutions
  for (int i = 0; i < solutions.size(); i++) { f.WriteObject(solutions[i].get(), Form("sol%i", i)); }
  f.Close();
}

///
/// Save a scanner from plots/scanner.
/// It contains the 1-CL histograms and the solutions.
///
bool MethodAbsScan::loadScanner(TString fName) {
  auto error = [](const std::string& msg) { errBase("MethodAbsScan::loadScanner : ERROR : ", msg); };

  if (fName == "") {
    FileNameBuilder fb(arg);
    fName = fb.getFileNameScanner(this);
  }
  if (arg->debug) cout << "MethodAbsScan::loadScanner() : ";
  cout << "loading scanner: " << fName << endl;
  if (!FileExists(fName))
    error(std::format("File not found: {:s}\n"
                      "Run first without the '-a plot' option to produce the missing file.",
                      std::string(fName)));
  TFile f(fName, "read");
  if (f.IsZombie()) error(std::format("File {:s} is corrupted", std::string(fName)));
  // load 1-CL histograms
  TObject* obj = f.Get("hCL");
  if (!obj) {
    cout << "MethodAbsScan::loadScanner() : ERROR : 'hCL' not found in root file " << fName << endl;
    exit(1);
  }
  if (scanVar2 != "") {
    hCL2d = Utils::clone<TH2>(obj, "hCL2d" + getUniqueRootName());
  } else {
    hCL = Utils::clone<TH1>(obj, "hCL" + getUniqueRootName());
  }
  // load chi2 histograms
  obj = f.Get("hChi2min");
  if (!obj) {
    cout << "MethodAbsScan::loadScanner() : ERROR : 'hChi2min' not found in root file " << fName << endl;
    // exit(1);
    // return false;
  }
  if (scanVar2 != "") {
    hChi2min2d = Utils::clone<TH2>(obj, "hChi2min2d" + getUniqueRootName());
  } else {
    hChi2min = Utils::clone<TH1>(obj, "hChi2min" + getUniqueRootName());
  }
  // load CLs histograms
  if (std::ranges::find(arg->cls, 1) != arg->cls.end()) {
    obj = f.Get("hCLs");
    if (!obj) {
      cout << "MethodAbsScan::loadScanner() : WARNING : 'hCLs' not found in root file - you can ignore this if you're "
              "not running in dataset mode "
           << fName << endl;
    }
    if (scanVar2 != "") {
      hCLs2d = Utils::clone<TH2>(obj, "hCLs2d" + getUniqueRootName());
    } else {
      hCLs = Utils::clone<TH1>(obj, "hCLs" + getUniqueRootName());
    }
  }
  // load CLs histograms
  bool lookForMixedCLs = std::ranges::find(arg->cls, 2) != arg->cls.end() && !methodName.Contains("Prob");
  if (lookForMixedCLs) {
    obj = f.Get("hCLsFreq");
    if (!obj) {
      cout << "MethodAbsScan::loadScanner() : WARNING : 'hCLsFreq' not found in root file - you can ignore this if "
              "you're not running in dataset mode "
           << fName << endl;
    } else if (scanVar2 == "") {
      hCLsFreq = Utils::clone<TH1>(obj, "hCLsFreq" + getUniqueRootName());
    }
    obj = f.Get("hCLsExp");
    if (!obj) {
      cout << "MethodAbsScan::loadScanner() : WARNING : 'hCLsExp' not found in root file - you can ignore this if "
              "you're not running in dataset mode "
           << fName << endl;
    } else if (scanVar2 == "") {
      hCLsExp = Utils::clone<TH1>(obj, "hCLsExp" + getUniqueRootName());
    }
    obj = f.Get("hCLsErr1Up");
    if (!obj) {
      cout << "MethodAbsScan::loadScanner() : WARNING : 'hCLsErr1Up' not found in root file - you can ignore this if "
              "you're not running in dataset mode "
           << fName << endl;
    } else if (scanVar2 == "") {
      hCLsErr1Up = Utils::clone<TH1>(obj, "hCLsErr1Up" + getUniqueRootName());
    }
    obj = f.Get("hCLsErr1Dn");
    if (!obj) {
      cout << "MethodAbsScan::loadScanner() : WARNING : 'hCLsErr1Dn' not found in root file - you can ignore this if "
              "you're not running in dataset mode "
           << fName << endl;
    } else if (scanVar2 == "") {
      hCLsErr1Dn = Utils::clone<TH1>(obj, "hCLsErr1Dn" + getUniqueRootName());
    }
    obj = f.Get("hCLsErr2Up");
    if (!obj) {
      cout << "MethodAbsScan::loadScanner() : WARNING : 'hCLsErr2Up' not found in root file - you can ignore this if "
              "you're not running in dataset mode "
           << fName << endl;
    } else if (scanVar2 == "") {
      hCLsErr2Up = Utils::clone<TH1>(obj, "hCLsErr2Up" + getUniqueRootName());
    }
    obj = f.Get("hCLsErr2Dn");
    if (!obj) {
      cout << "MethodAbsScan::loadScanner() : WARNING : 'hCLsErr2Dn' not found in root file - you can ignore this if "
              "you're not running in dataset mode "
           << fName << endl;
    } else if (scanVar2 == "") {
      hCLsErr2Dn = Utils::clone<TH1>(obj, "hCLsErr2Dn" + getUniqueRootName());
    }
  }

  // load solutions
  solutions.clear();
  int nSol = 100;
  for (int i = 0; i < nSol; i++) {
    if (auto sfr = dynamic_cast<RooSlimFitResult*>(f.Get(Form("sol%i", i)))) solutions.push_back(sfr->Clone());
  }
  if (f.Get(Form("sol%i", nSol))) {
    cout << "MethodAbsScan::loadScanner() : WARNING : Only the first 100 solutions read from: " << fName << endl;
  }
  f.Close();

  return true;
}

int MethodAbsScan::scan1d() const {
  cout << "MethodAbsScan::scan1d() : not implemented." << endl;
  return 0;
}

int MethodAbsScan::scan2d() const {
  cout << "MethodAbsScan::scan2d() : not implemented." << endl;
  return 0;
}

/**
 * Find the value of x for which a function takes on the value y, given a histogram that approximates the function.
 *
 * The interpolation of x that satisfies h(x) = y is done by means of a linear function between the y values of two
 * neighbouring bins i and i+1, such that y takes on a value in between h[i] and h[i+1].
 *
 * @param h The histogram to be interpolated.
 * @param i Interpolate around this bin (i and i+1 must be are above and below the interpolated value).
 * @param y The y position we want to find the interpolated x for.
 *
 * @return Interpolated x position, or std::nullopt if this value could not be calculated.
 */
std::optional<double> MethodAbsScan::interpolateLinear(const TH1* h, const int i, const double y) const {
  if (!(1 <= i && i <= h->GetNbinsX() - 1)) return {};
  double p1x = h->GetBinCenter(i);
  double p1y = h->GetBinContent(i);
  double p2x = h->GetBinCenter(i + 1);
  double p2y = h->GetBinContent(i + 1);
  if (!((p1y < y && y < p2y) || (p2y < y && y < p1y))) [[unlikely]] {
    std::cerr << "MethodAbsScan::interpolateLinear : ERROR : There is a problem in GammaCombo core" << std::endl;
    exit(1);
  }
  return p2x + (y - p2y) / (p1y - p2y) * (p1x - p2x);
}

///
/// Solve a quadratic equation by means of a modified pq formula:
/// @f[x^2 + \frac{p_1}{p_2} x + \frac{p_0-y}{p2} = 0@f]
///
double MethodAbsScan::pq(double p0, double p1, double p2, double y, int whichSol) const {
  if (whichSol == 0)
    return -p1 / 2. / p2 + sqrt(sq(p1 / 2. / p2) - (p0 - y) / p2);
  else
    return -p1 / 2. / p2 - sqrt(sq(p1 / 2. / p2) - (p0 - y) / p2);
}

/**
 * Find an interpolated x value near a certain bin position of a histogram that is the best estimate for h(x)=y.
 *
 * Interpolates by means of fitting a second grade polynomial to up to five adjacent points.
 * Because that's giving two solutions, we use the central value and knowledge about if it is supposed to be an upper
 * or lower boundary to pick one.
 *
 * @param h       The histogram to be interpolated.
 * @param i       Interpolate around this bin (must be such that i and i+1 are above and below the interpolated value).
 * @param y       The y position we want to find the interpolated x for.
 * @param central Central value of the solution we're trying to get the CL interval for.
 * @param upper   Set to true if we're computing an upper interval boundary.
 *
 * @return        std::nullopt if the interpolation failed, otherwise a pair {val, err} where val is the interpolated
 *                x position, and err the estimated interpolation error.
 */
std::optional<std::pair<double, double>> MethodAbsScan::interpolate(TH1* h, const int i, const double y,
                                                                    const double central, const bool upper) const {
  // cout << "MethodAbsScan::interpolate(): i=" << i << " y=" << y << " central=" << central << endl;
  if (i > h->GetNbinsX() - 2 || i < 3) return {};

  // if method Prob, don't interpolate (no proper error estimate)
  if (methodName.Contains("Prob")) {
    for (int k = 0; k < h->GetNbinsX(); k++) { h->SetBinError(k + 1, 0.); }
  }

  // compute pol2 fit interpolation
  auto g = std::make_unique<TGraphErrors>(3);
  g->SetPoint(0, h->GetBinCenter(i - 1), h->GetBinContent(i - 1));
  g->SetPointError(0, h->GetBinWidth(i - 1) / 2., h->GetBinError(i - 1));
  g->SetPoint(1, h->GetBinCenter(i), h->GetBinContent(i));
  g->SetPointError(1, h->GetBinWidth(i) / 2., h->GetBinError(i));
  g->SetPoint(2, h->GetBinCenter(i + 1), h->GetBinContent(i + 1));
  g->SetPointError(2, h->GetBinWidth(i + 1) / 2., h->GetBinError(i + 1));

  // see if we can add a 4th and 5th point
  if ((h->GetBinContent(i - 2) - h->GetBinError(i - 2) < h->GetBinContent(i - 1) + h->GetBinError(i - 1) &&
       h->GetBinContent(i - 1) < h->GetBinContent(i)) ||
      (h->GetBinContent(i - 2) + h->GetBinError(i - 2) > h->GetBinContent(i - 1) - h->GetBinError(i - 1) &&
       h->GetBinContent(i - 1) > h->GetBinContent(i))) {
    if ((upper && h->FindBin(central) < i - 2) ||
        !upper)  // don't use for upper limit calculation if point is equal or below central value
    {
      // add to the beginning
      auto gNew = std::make_unique<TGraphErrors>(g->GetN() + 1);
      gNew->SetPoint(0, h->GetBinCenter(i - 2), h->GetBinContent(i - 2));
      gNew->SetPointError(0, h->GetBinWidth(i - 2) / 2., h->GetBinError(i - 2));
      Double_t pointx, pointy;
      Double_t pointxerr, pointyerr;
      for (int i = 0; i < g->GetN(); i++) {
        g->GetPoint(i, pointx, pointy);
        pointxerr = g->GetErrorX(i);
        pointyerr = g->GetErrorY(i);
        gNew->SetPoint(i + 1, pointx, pointy);
        gNew->SetPointError(i + 1, pointxerr, pointyerr);
      }
      g = std::move(gNew);
    }
  }

  if ((h->GetBinContent(i + 2) - h->GetBinError(i + 2) < h->GetBinContent(i + 1) + h->GetBinError(i + 1) &&
       h->GetBinContent(i + 1) < h->GetBinContent(i)) ||
      (h->GetBinContent(i + 2) + h->GetBinError(i + 2) > h->GetBinContent(i + 1) - h->GetBinError(i + 1) &&
       h->GetBinContent(i + 1) > h->GetBinContent(i))) {
    if ((!upper && h->FindBin(central) > i + 2) ||
        upper)  // don't use for lower limit calculation if point is equal or above central value
    {
      // add to the end
      g->Set(g->GetN() + 1);
      g->SetPoint(g->GetN() - 1, h->GetBinCenter(i + 2), h->GetBinContent(i + 2));
      g->SetPointError(g->GetN() - 1, h->GetBinWidth(i + 2) / 2., h->GetBinError(i + 2));
    }
  }

  auto f1 = std::make_unique<TF1>("f1", "[0]+[1]*(x-[2])", h->GetBinCenter(i - 2), h->GetBinCenter(i + 2));
  auto f2 =
      std::make_unique<TF1>("f2", "[0]+[1]*(x-[3])+[2]*(x-[3])**2", h->GetBinCenter(i - 2), h->GetBinCenter(i + 2));
  f1->FixParameter(2, h->GetBinCenter(i));
  f2->FixParameter(3, h->GetBinCenter(i));

  f1->SetParameter(1, (h->GetBinContent(i + 1) - h->GetBinContent(i)) / h->GetBinWidth(i));
  g->Fit("f1", "q");  // fit linear to get decent start parameters
  f2->SetParameter(0, f1->GetParameter(0));
  f2->SetParameter(1, f1->GetParameter(1));
  g->Fit("f2", "qf+");  // refit with minuit to get more correct errors (TGraph fit errors bug)
  array<double, 3> p;
  // double e[3];
  // for ( int ii=0; ii<3; ii++ )
  // {
  //  p[ii] = f2->GetParameter(ii);
  //  e[ii] = f2->GetParError(ii);
  // }
  p[0] = f2->GetParameter(2) * (f2->GetParameter(3) * f2->GetParameter(3)) - f2->GetParameter(1) * f2->GetParameter(3) +
         f2->GetParameter(0);
  p[1] = f2->GetParameter(1) - 2 * f2->GetParameter(2) * f2->GetParameter(3);
  p[2] = f2->GetParameter(2);

  double sol0 = pq(p[0], p[1], p[2], y, 0);
  double sol1 = pq(p[0], p[1], p[2], y, 1);
  // cout << upper << " ";
  // printf("%f %f %f\n", central, sol0, sol1);

  // std::cout << central << "\t" << sol0 << "\t" <<sol1 << std::endl;

  // debug: show fitted 1-CL histogram
  if (arg->controlplot) {
    TString debugTitle = methodName + Form(" y=%.2f ", y);
    debugTitle += upper ? Form("%f upper", central) : Form("%f lower", central);
    auto c = newNoWarnTCanvas(getUniqueRootName(), debugTitle);
    g->SetMarkerStyle(3);
    auto th1f = dynamic_cast<TH1F*>(h);
    if (!th1f) {
      std::cerr << "Cannot cast to TH1F!" << std::endl;
      exit(1);
    }
    g->SetHistogram(th1f);
    h->Draw();
    g->Draw("p");
    f2->Draw("SAME");
    savePlot(c.get(), TString(name + "_" + scanVar1 + "_boundary_interpolation_" + methodName + "_" +
                              TString(h->GetName()) + "_" + std::to_string(y)));
  }

  if ((h->GetBinCenter(i - 2) > sol0 || sol0 > h->GetBinCenter(i + 2)) &&
      (h->GetBinCenter(i - 2) > sol1 || sol1 > h->GetBinCenter(i + 2))) {
    if (arg->verbose || arg->debug) {
      cout << "MethodAbsScan::interpolate(): Quadratic interpolation out of bounds [" << h->GetBinCenter(i - 2) << ", "
           << h->GetBinCenter(i + 2) << "]:" << std::endl;
      std::cout << "Solutions are " << central << "(free fit result)\t" << sol0 << "(bound solution 0) \t" << sol1
                << "(bound solution 1)." << std::endl;
    }
    return {};
  } else if (sol0 != sol0 || sol1 != sol1) {
    if (arg->verbose || arg->debug) {
      cout << "MethodAbsScan::interpolate(): Quadratic interpolation leads to NaN:" << std::endl;
      std::cout << "Solutions are " << central << "(free fit result)\t" << sol0 << "(bound solution 0) \t" << sol1
                << "(bound solution 1)." << std::endl;
    }
    return {};
  }

  int useSol = 0;
  if ((sol0 < central && sol1 > central) || (sol1 < central && sol0 > central)) {
    if (upper) {
      useSol = (sol0 < sol1) ? 1 : 0;
    } else {
      useSol = (sol0 < sol1) ? 0 : 1;
    }
  } else {
    useSol = (std::abs(h->GetBinCenter(i) - sol0) < std::abs(h->GetBinCenter(i) - sol1)) ? 0 : 1;
  }

  double val = (useSol == 0) ? sol0 : sol1;

  // try error propagation: sth is wrong in the formulae
  /*
  double err0 = TMath::Max(sq(val-pq(p[0]+e[0], p[1], p[2], y, useSol)), sq(val-pq(p[0]-e[0], p[1], p[2], y, useSol)));
  double err1 = TMath::Max(sq(val-pq(p[0], p[1]+e[1], p[2], y, useSol)), sq(val-pq(p[0], p[1]-e[1], p[2], y, useSol)));
  double err2 = TMath::Max(sq(val-pq(p[0], p[1], p[2]+e[2], y, useSol)), sq(val-pq(p[0], p[1], p[2]-e[2], y, useSol)));
  err = sqrt(err0+err1+err2);
  printf("%f %f %f\n", val, pq(p[0]+e[0], p[1], p[2], y, useSol), pq(p[0]-e[0], p[1], p[2], y, useSol));
  printf("%f %f %f\n", val, pq(p[0], p[1]+e[1], p[2], y, useSol), pq(p[0], p[1]-e[1], p[2], y, useSol));
  printf("%f %f %f\n", val, pq(p[0], p[1], p[2]+e[2], y, useSol), pq(p[0], p[1], p[2]-e[2], y, useSol));
  */
  auto err = std::numeric_limits<double>::quiet_NaN();

  return std::make_pair(val, err);
}

/**
 * Calculate the CL intervals from the CL curve.
 *
 * Start from known local minima and scan upwards and downwards to find the interval boundaries. Then scan again from
 * the boundaries of the scan range to cover the case where an CL interval is not closed yet at the boundary.
 * Use a fit-based interpolation (@see interpolate) if we have more than 25 bins, else revert to a straight line
 * interpolation (@see interpolateLinear).
 */
void MethodAbsScan::calcCLintervals(const int CLsType, const bool calc_expected, const bool quiet) {

  // Messaging
  auto debug = [](const std::string& msg) { msgBase("MethodAbsScan::calcCLintervals() : DEBUG : ", msg); };
  auto info = [](const std::string& msg) { msgBase("MethodAbsScan::calcCLintervals() : ", msg); };
  auto warning = [](const std::string& msg) { msgBase("MethodAbsScan::calcCLintervals() : WARNING : ", msg); };
  auto error = [](const std::string& msg) { errBase("MethodAbsScan::calcCLintervals() : ERROR : ", msg); };
  if (arg->debug) debug(std::format("Calling arguments: {:d}, {:s}, {:s}", CLsType, calc_expected, quiet));

  // TODO
  auto histogramCL = this->getHCL();
  // calc CL intervals with CLs method
  if (CLsType == 1 && this->getHCLs()) {
    histogramCL = this->getHCLs();
  } else if (CLsType == 2 && this->getHCLsFreq()) {
    histogramCL = this->getHCLsFreq();
  }
  if (CLsType == 2 && calc_expected && hCLsExp) {
    histogramCL = this->getHCLsExp();
    std::cout << "Determine expected upper limit:" << std::endl;
  }
  if (!histogramCL) {
    std::cerr << "ERROR : Could not retrieve the histogram. Will not calculate the CLs" << endl;
    return;
  }

  if (CLsType != 0) { std::cout << std::format("Confidence Intervals for CLs method {:d}:", CLsType) << std::endl; }
  if (arg->isQuickhack(8)) {
    // TODO Switch to the new CLIntervalMaker mechanism. It can be activated already using --qh 8,
    //      but it really is in beta stage still
    // TODO Add user specific CL interval
    cout << "\n";
    info(std::format("USING NEW CLIntervalMaker for {:s}\n", std::string(name)));
    CLIntervalMaker clm(arg, histogramCL);
    clm.findMaxima(0.04);  // ignore maxima under pvalue=0.04
    for (int iSol = 0; iSol < solutions.size(); iSol++) {
      auto sol = getScanVar1Solution(iSol);
      clm.provideMorePreciseMaximum(sol, "max PLH");
    }
    clm.calcCLintervals();
    // print
    TString unit = w->var(scanVar1)->getUnit();
    CLIntervalPrinter clp(arg, name, scanVar1, unit,
                          std::string(methodName) + (calc_expected ? "_expected_standardCLs" : ""));
    clp.setDegrees(isAngle(w->var(scanVar1)));
    clp.addIntervals(clm.getClintervals1sigma());
    clp.addIntervals(clm.getClintervals2sigma());
    clp.print();
  }

  // TODO
  if (solutions.empty()) {
    info("Solutions vector empty. Using simple method with linear splines");
    this->calcCLintervalsSimple(CLsType, calc_expected);
    return;
  } else if ((CLsType == 1 || CLsType == 2) && !this->getHCLs()) {
    info("Using simple method with linear splines");
    this->calcCLintervalsSimple(CLsType, calc_expected);
  }

  if (!quiet) cout << std::format("\nCONFIDENCE INTERVALS for combination `{:s}`\n", std::string(name)) << endl;
  clintervals.clear();
  clintervals.resize(ConfidenceLevels.size());

  // Vector containing pairs of starting points and relative bins
  std::vector<std::pair<double, int>> starts;
  const auto n = histogramCL->GetNbinsX();
  for (int i = 0; i < solutions.size(); i++) {
    const auto sol = getScanVar1Solution(i);
    int bin = histogramCL->FindBin(sol);
    if (histogramCL->IsBinOverflow(bin) || histogramCL->IsBinUnderflow(bin)) {
      warning(std::format(
          "Solution {:d} is outside of the scanrange, I will not try to find the relative confidence interval", i));
      starts.emplace_back(std::numeric_limits<double>::quiet_NaN(), -1);
    } else {
      if (bin == 1 || bin == n)
        warning(std::format("Solution {:d} lies at the border of the scanrange, you should increase it", i));
      starts.emplace_back(sol, bin);
    }
  }
  starts.emplace_back(histogramCL->GetXaxis()->GetXmin(), 1);
  starts.emplace_back(histogramCL->GetXaxis()->GetXmax(), n);

  const int minBinsForInterpolation = 25;
  if (n <= minBinsForInterpolation) info("Low number of scan points. Will use linear interpolation");

  // Find and save a confidence interval for each solution within the scan range.
  for (const auto [start, sBin] : starts) {
    if (arg->debug) {
      if (sBin == 1)
        debug("Start scan of low boundary");
      else if (sBin == n)
        debug("Start scan of up boundary");
    }

    for (int c = 0; c < ConfidenceLevels.size(); c++) {
      const double y = 1. - ConfidenceLevels[c];
      auto CLmin = std::numeric_limits<double>::quiet_NaN();
      auto CLmax = std::numeric_limits<double>::quiet_NaN();
      auto CLminErr = std::numeric_limits<double>::quiet_NaN();
      auto CLmaxErr = std::numeric_limits<double>::quiet_NaN();
      bool CLminClosed = false;
      bool CLmaxClosed = false;

      // Case that the solution is outside of scan region, or starting point does not fall into the CL region
      // (e.g. because we are doing a border scan or because the solution is a shallow local minimum)
      if (std::isnan(start) || histogramCL->GetBinContent(sBin) < y) {
        clintervals[c].push_back(nullptr);
        continue;
      }

      // Find lower interval bound
      for (int i = sBin; i > 0; i--) {
        if (i == 1) {
          CLmin = histogramCL->GetXaxis()->GetXmin();
          if (sBin != 1)
            warning(std::format(
                "I am using the lowest bin of histogramCL to calculate the lower end of the CL interval #{:d}.\n"
                "This will lead to wrong results - you need to decrease the minimum of the scan range",
                c));
        }
        if (histogramCL->GetBinContent(i) < y) {
          bool linearInterpolation = true;
          if (n > minBinsForInterpolation) {
            if (auto pair = interpolate(histogramCL, i, y, start, false)) {
              CLmin = pair->first;
              CLminErr = pair->second;
              CLminClosed = true;
              linearInterpolation = false;
            }
          }
          if (linearInterpolation) {
            if (n > minBinsForInterpolation && (arg->verbose || arg->debug)) info("Reverting to linear interpolation.");
            if (auto val = interpolateLinear(histogramCL, i, y)) {
              CLmin = *val;
              CLminClosed = true;
            } else
              warning(std::format("Could not interpolate for bin {:d}", i));
          }
          break;
        }
      }

      // Find upper interval bound
      for (int i = sBin; i <= n; i++) {
        if (i == n) {
          CLmax = histogramCL->GetXaxis()->GetXmax();
          if (sBin != n)
            warning(std::format(
                "I am using the highest bin of histogramCL to calculate the upper end of the CL interval #{:d}.\n"
                "This will lead to wrong results - you need to increase the maximum of the scan range",
                c));
        }
        if (histogramCL->GetBinContent(i) < y) {
          bool linearInterpolation = true;
          if (n > minBinsForInterpolation) {
            if (auto pair = interpolate(histogramCL, i - 1, y, start, true)) {
              CLmax = pair->first;
              CLmaxErr = pair->second;
              CLmaxClosed = true;
              linearInterpolation = false;
            }
          }
          if (linearInterpolation) {
            if (n > minBinsForInterpolation && (arg->verbose || arg->debug)) info("Reverting to linear interpolation.");
            if (auto val = interpolateLinear(histogramCL, i - 1, y)) {
              CLmax = *val;
              CLmaxClosed = true;
            } else
              warning(std::format("Could not interpolate for bin {:d}", i - 1));
          }
          break;
        }
      }

      // Save the interval
      auto cli = std::make_unique<CLInterval>();
      cli->pvalue = y;
      cli->central = (sBin != 1 && sBin != n) ? start : std::numeric_limits<double>::quiet_NaN();
      cli->min = CLmin;
      cli->max = CLmax;
      cli->minclosed = CLminClosed;
      cli->maxclosed = CLmaxClosed;
      if (arg->debug) cli->print();

      clintervals[c].push_back(std::move(cli));
    }
  }

  // Compute the cover of all 1sigma intervals
  // TODO this does not make much sense (especially for angles)
  if (arg->largest) {
    for (int k = 0; k < clintervals[0].size(); k++) {
      auto i = std::make_unique<CLInterval>();
      i->central = clintervals[0][k]->central;
      i->pvalue = clintervals[0][k]->pvalue;
      i->minmethod = "largest";
      i->maxmethod = "largest";
      i->min = clintervals[0][0]->min;
      i->max = clintervals[0][0]->max;
      i->minclosed = clintervals[0][0]->minclosed;
      i->maxclosed = clintervals[0][0]->maxclosed;
      for (const auto& cli : clintervals[0]) {
        if (!cli) continue;
        if (cli->min < i->min) {
          i->min = cli->min;
          i->minclosed = cli->minclosed;
        }
        if (cli->max > i->max) {
          i->max = cli->max;
          i->maxclosed = cli->maxclosed;
        }
      }
      clintervals[0].push_back(std::move(i));
    }
  }
  if (!quiet) printCLintervals(CLsType, calc_expected);

  // Print fit chi2 etc. (not done for datasets running)
  if (!combiner || !combiner->isCombined()) return;
  const auto chi2 = this->getSolution(0)->minNll();
  const auto nObs = combiner->getObservables()->getSize();
  const auto nPar = combiner->getParameters()->getSize();
  if (nObs == nPar) {
    if (std::abs(chi2) > 1e-3)
      cerr << "ERROR : Chi2 is not zero ({:.4f}), even if the number of degrees of freedom is zero\n" << endl;
    else
      cout << std::format(
                  "Fit quality is meaningless for zero degrees of freedom (chi2 = 0): (nObs, nPar) = ({:d}, {:d})\n",
                  nObs, nPar)
           << endl;
  } else {
    const auto prob = TMath::Prob(chi2, nObs - nPar);
    cout << std::format("Fit quality: chi2/(nObs-nPar) = {:.2f}/({:d}-{:d}), P = {:4.1f}%\n", chi2, nObs, nPar,
                        prob * 100.)
         << endl;
  }
}

/**
 * Print the CL intervals.
 */
void MethodAbsScan::printCLintervals(const int CLsType, const bool calc_expected) {
  const auto unit = w->var(scanVar1)->getUnit();
  CLIntervalPrinter clp(arg, name, scanVar1, unit, methodName, CLsType);
  if (calc_expected) {
    clp = CLIntervalPrinter(arg, name, scanVar1, unit, methodName + TString("_expected_standardCLs"));
  }
  clp.setDegrees(isAngle(w->var(scanVar1)));
  clp.addIntervals(clintervals);
  clp.print();
  clp.savePython();
  cout << endl;

  // Print solutions not contained in any of the intervals
  for (int i = 0; i < solutions.size(); i++) {
    double sol = getScanVar1Solution(i);
    bool cont = false;
    for (const auto& clis : clintervals) {
      for (const auto& cli : clis) {
        if (!cli) continue;
        if (cli->min < sol && sol < cli->max) cont = true;
      }
    }
    if (cont) continue;
    if (w->var(scanVar1)->getUnit() == TString("Rad")) sol = RadToDeg(sol);
    int d = arg->digits;
    if (d <= 0) d = 3;
    printf("%s = %7.*f", w->var(scanVar1)->GetName(), d, sol);
    if (unit != "") cout << " [" << unit << "]";
    cout << endl;
  }
}

///
/// Get the CL interval that includes the best-fit value.
/// \param sigma 1,2
///
const CLInterval* MethodAbsScan::getCLintervalCentral(int sigma, bool quiet) { return getCLinterval(0, sigma, quiet); }

/**
 * Get the CL interval that includes a given solution.
 *
 * @param iSol  Index of the interval (by default, solutions 0, 1, 2... followed by border scans).
 * @param index Index identifying the CL (by default, i = 0,1,2 correspond to 1,2,3 sigma).
 * @bool quiet  Sets whether the calculation of CL intervals (in case they have not been calculated yet) should be
 *              verbose or not.
 * @return      A copy of the desired CLInterval.
 */
const CLInterval* MethodAbsScan::getCLinterval(const int iSol, const int index, const bool quiet) {
  auto error = [](const std::string& msg) { errBase("MethodAbsScan::getCLinterval : ERROR : ", msg); };

  if (clintervals.empty()) calcCLintervals(0, false, quiet);
  if (clintervals.empty()) error("This should never happen");

  if (index < 0 || index >= clintervals.size())
    error(std::format("There are no CL intervals with CL identified by index {:d}", index));
  if (iSol < 0 || iSol >= clintervals[index].size())
    error(
        std::format("There are no CL intervals with CL identified by index {:d} and solution identified by index {:d}",
                    index, iSol));

  // compute largest interval
  // TODO this does not make sense at all, in case there are multiple solutions
  // if (arg->largest) {
  //   CLInterval i;
  //   i.pvalue = intervals[iSol].pvalue;
  //   i.min = intervals[iSol].min;
  //   for (const auto interval : intervals) i.min = TMath::Min(i.min, interval.min);
  //   i.max = intervals[iSol].max;
  //   for (const auto interval : intervals) i.max = TMath::Max(i.max, interval.max);
  //   return i;
  // }

  return clintervals[index][iSol].get();
}

double MethodAbsScan::getCL(double val) const { return 1. - hCL->Interpolate(val); }

/**
 * @param[in] \CLsType: 0 (off), 1 (naive CLs t_s+b - t_b), 2 (freq CLs)
 */
void MethodAbsScan::plotOn(OneMinusClPlotAbs* plot, int CLsType) { plot->addScanner(this, CLsType); }

RooRealVar* MethodAbsScan::getScanVar1() { return w->var(scanVar1); }

RooRealVar* MethodAbsScan::getScanVar2() { return w->var(scanVar2); }

void MethodAbsScan::print() const {
  cout << "MethodAbsScan::print() : Method: " << methodName;
  cout << ", Scanner: " << name << endl;
  w->set(parsName)->Print("v");
}

///
/// Make a 1d plot of the NLL in var
///
void MethodAbsScan::plot1d(TString var) {
  cout << "MethodAbsScan::plot1d() : Method: " << methodName;
  cout << ", Scanner: " << name << endl;

  //   RooRealVar* vx = w->var(var);
  //   assert(vx);
  // setLimit(w, var, "plot");
  //
  //   // cout << "MethodAbsScan::plot1d() : loading global minimum ..." << endl;
  //   // if ( !globalMin ){ cout << "MethodAbsScan::plot1d() : no global minimum. Call doInitialFit() first!" << endl;
  //   exit(1); }
  //   // setParameters(w, parsName, globalMinP);
  //   // print();
  //
  //   RooNLLVar nll("nll", "nll", *(w->pdf(pdfName)), *(w->data(dataName))) ;
  //
  //   TString plotName = "plot1d_"+name+"_"+var;
  //   auto c1 = newNoWarnTCanvas();
  //   RooPlot *frame = vx->frame();
  //   // w->pdf(pdfName)->plotOn(frame);
  //   nll.plotOn(frame);
  //   frame->Draw();
  //
  //   savePlot(c1.get(), plotName);
}

///
/// Make a 2d plot of the PDF in varx and vary.
///
void MethodAbsScan::plot2d(TString varx, TString vary) {
  cout << "MethodAbsScan::plot2d() : Method: " << methodName;
  cout << ", scanner: " << name << endl;

  RooRealVar* vx = w->var(varx);
  RooRealVar* vy = w->var(vary);
  assert(vx);
  assert(vy);
  setLimit(w, varx, "plot");
  setLimit(w, vary, "plot");

  cout << "MethodAbsScan::plot2d() : loading global minimum ..." << endl;
  if (!globalMin) {
    cout << "MethodAbsScan::plot2d() : no global minimum. Call doInitialFit() first!" << endl;
    exit(1);
  }

  setParameters(w, parsName, globalMin.get());
  print();

  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPalette(1);

  TString plotName = "plot2d_" + name + "_" + varx + "_" + vary;
  auto c1 = newNoWarnTCanvas(plotName, plotName);
  auto h = w->pdf(pdfName)->createHistogram(plotName, *vx, YVar(*vy));
  h->Draw("colz");

  savePlot(c1.get(), plotName + arg->plotext);
}

/**
 * Load the values at a specific minimum into the workspace.
 *
 * This way we can use it for goodness of fit, start points, etc.
 *
 * @param i Index of the solution, i=0 corresponds to the best one.
 */
bool MethodAbsScan::loadSolution(const int i) {
  if (arg->debug) cout << "MethodAbsScan::loadSolution() : loading solution " << i << endl;
  if (i < 0 || i >= solutions.size()) {
    cout << "MethodAbsScan::loadSolution() : ERROR : solution ID out of range." << endl;
    return false;
  }
  auto tmp = std::make_unique<RooArgSet>();
  tmp->add(solutions[i]->floatParsFinal());
  tmp->add(solutions[i]->constPars());
  setParameters(w, parsName, tmp.get());
  return true;
}

///
/// Load the values given by an (external) fit result.
///
void MethodAbsScan::loadParameters(const RooSlimFitResult* r) {
  if (arg->debug) cout << "MethodAbsScan::loadParameters() : loading a RooSlimFitResult " << endl;
  auto tmp = std::make_unique<RooArgSet>();
  tmp->add(r->floatParsFinal());
  tmp->add(r->constPars());
  setParameters(w, parsName, tmp.get());
}

///
/// Print local minima solutions.
///
void MethodAbsScan::printLocalMinima() const {
  TDatime date;  // lets also print the current date
  if (arg->debug) {
    cout << "MethodAbsScan::printLocalMinima() : LOCAL MINIMA for " << title << endl;
    cout << endl;
  }
  for (int i = 0; i < solutions.size(); i++) {
    cout << "SOLUTION " << i << ":\n" << endl;
    cout << "  combination: " << name << endl;
    cout << "  title:       " << title << endl;
    cout << "  date:        " << date.AsString() << endl;
    solutions[i]->Print(arg->verbose, arg->printcor);
  }
}

///
/// Save local minima solutions.
///
void MethodAbsScan::saveLocalMinima(TString fName) const {
  TDatime date;  // lets also print the current date
  if (arg->debug) {
    cout << "MethodAbsScan::saveLocalMinima() : LOCAL MINIMA for " << title << endl;
    cout << endl;
  }
  ofstream outfile;
  outfile.open(fName.Data());

  for (int i = 0; i < solutions.size(); i++) {
    outfile << "\%SOLUTION " << i << ":\n" << endl;
    outfile << "\%  combination: " << name << endl;
    outfile << "\%  title:       " << title << endl;
    outfile << "\%  date:        " << date.AsString() << endl;
    solutions[i]->SaveLatex(outfile, arg->verbose, arg->printcor);
  }
  outfile.close();
}

/**
 * Get value of scan parameter at a certain solution.
 *
 * @param iVar Index of scan variable, 1 or 2.
 * @param iSol Index of solution. 0 corresponds to the best one, indices increase in order of chi2.
 *
 * @return     Central value of the solution
 *             -999 no solutions available
 *             -99 solution not found
 *             -9999 no such variable
 */
double MethodAbsScan::getScanVarSolution(const int iVar, const int iSol) {
  auto error = [](const std::string& msg, const double err) {
    std::cerr << "MethodAbsScan::getScanVarSolution() : ERROR : " << msg << endl;
    return err;
  };
  if (solutions.empty()) error("The vector of solutions is empty", -999.);
  if (iSol >= solutions.size()) error(std::format("No solution with id {:d}", iSol), -99.);
  auto r = getSolution(iSol);
  assert(r);
  TString varName;
  if (iVar == 1)
    varName = getScanVar1Name();
  else if (iVar == 2)
    varName = getScanVar2Name();
  else
    error(std::format("No such variable {:d}", iVar), -9999.);
  if (r->isConfirmed()) {
    return r->getFloatParFinalVal(varName);
  } else {
    if (nWarnings == 0) cout << "MethodAbsScan::getScanVarSolution() : WARNING : Using unconfirmed solution." << endl;
    nWarnings += 1;
    return r->getConstParVal(varName);
  }
}

/**
 * Get value of scan parameter 1 a certain solution.
 *
 * @param iSol Index of solution (0 corresponds to the best one, indices increase in order of chi2).
 */
double MethodAbsScan::getScanVar1Solution(int iSol) { return getScanVarSolution(1, iSol); }

/**
 * Get value of scan parameter 2 a certain solution.
 *
 * @param iSol Index of solution (0 corresponds to the best one, indices increase in order of chi2).
 */
double MethodAbsScan::getScanVar2Solution(int iSol) { return getScanVarSolution(2, iSol); }

///
/// Sort solutions in order of increasing chi2.
///
void MethodAbsScan::sortSolutions() {
  if (arg->debug) cout << "MethodAbsScan::sortSolutions() : sorting solutions ..." << endl;
  std::ranges::sort(solutions, [](const std::unique_ptr<RooSlimFitResult>& a,
                                  const std::unique_ptr<RooSlimFitResult>& b) { return a->minNll() < b->minNll(); });
}

/**
 * Refit all possible solutions with the scan parameter left free to confirm the solutions.
 * We will reject solutions as fake if the free fit using them as the starting point will move too far away,
 * or if their Delta chi2 value is above 25.
 */
void MethodAbsScan::confirmSolutions() {
  if (arg->debug) cout << "MethodAbsScan::confirmSolutions() : Confirming solutions ..." << endl;
  FitResultCache frCache(arg);
  frCache.storeParsAtFunctionCall(w->set(parsName));

  vector<std::unique_ptr<RooSlimFitResult>> confirmedSolutions;
  auto par1 = w->var(scanVar1);
  auto par2 = w->var(scanVar2);
  if (par1) par1->setConstant(false);
  if (par2) par2->setConstant(false);
  for (int i = 0; i < solutions.size(); i++) {
    if (!loadSolution(i)) {
      cout << std::format(
                  "MethodAbsScan::confirmSolutions() : WARNING Could not load solution {:d}, so I will not confirm", i)
           << endl;
      continue;
    }
    if (arg->debug) {
      cout << "MethodAbsScan::confirmSolutions() : solution " << i;
      cout << " " << par1->GetName() << "=" << par1->getVal();
      if (par2) cout << " " << par2->GetName() << "=" << par2->getVal();
      cout << endl;
    }

    // Refit the solution. `true` uses thorough fit with HESSE, -1 silences output
    auto r = fitToMinBringBackAngles(w->pdf(pdfName), true, -1);

    // Check scan parameter shift.
    // We'll allow for a shift equivalent to 3 step sizes.
    // Express the scan step size in terms of sigmas of the fitted parameters.
    double allowedSigma;
    if (arg->var.size() == 1) {
      // 1d scan
      double par1stepsize = (par1->getMax("scan") - par1->getMin("scan")) / arg->npoints1d;
      auto par1New = (RooRealVar*)r->floatParsFinal().find(par1->GetName());
      double par1stepsizeInSigma = par1New->getError() > 0 ? par1stepsize / par1New->getError() : 0.2;
      allowedSigma = 3. * par1stepsizeInSigma;
    } else if (arg->var.size() == 2) {
      // 2d scan
      double par1stepsize = (par1->getMax("scan") - par1->getMin("scan")) / arg->npoints2dx;
      double par2stepsize = (par2->getMax("scan") - par2->getMin("scan")) / arg->npoints2dy;
      auto par1New = static_cast<RooRealVar*>(r->floatParsFinal().find(par1->GetName()));
      auto par2New = static_cast<RooRealVar*>(r->floatParsFinal().find(par2->GetName()));
      double par1stepsizeInSigma = par1New->getError() > 0 ? par1stepsize / par1New->getError() : 1.;
      double par2stepsizeInSigma = par2New->getError() > 0 ? par2stepsize / par2New->getError() : 1.;
      allowedSigma = std::max(3. * par1stepsizeInSigma, 3. * par2stepsizeInSigma);
    }

    // Warn if a parameter is close to its limit
    for (const auto pAbs : r->floatParsFinal()) {
      auto p = static_cast<RooRealVar*>(pAbs);
      if (p->getMax() - p->getVal() < p->getError() || p->getVal() - p->getMin() < p->getError()) {
        cout << "\nMethodAbsScan::confirmSolutions() : WARNING : " << p->GetName() << " is close to its limit!" << endl;
        cout << "                                  : ";
        p->Print();
        cout << endl;
      }
    }

    // check migration of the parameters
    RooArgList listOld = solutions[i]->floatParsFinal();
    listOld.add(solutions[i]->constPars());
    RooArgList listNew = r->floatParsFinal();
    listNew.add(r->constPars());
    bool isConfirmed = true;
    TString rejectReason = "";
    for (const auto p : *w->set(parsName)) {
      auto pOld = (RooRealVar*)listOld.find(p->GetName());
      auto pNew = (RooRealVar*)listNew.find(p->GetName());
      if (!pOld && !pNew) {
        cout << "MethodAbsScan::confirmSolutions() : ERROR : parameter not found: " << p->GetName() << endl;
        continue;
      }
      if (pNew->getError() > 0) {
        double shift = fabs(pOld->getVal() - pNew->getVal());
        if (isAngle(pOld)) shift = angularDifference(pOld->getVal(), pNew->getVal());
        if (shift / pNew->getError() > allowedSigma) {
          if (arg->debug) {
            cout << "MethodAbsScan::confirmSolutions() : solution " << i << ", too large parameter shift:" << endl;
            pOld->Print();
            pNew->Print();
          }
          isConfirmed = false;
          rejectReason = TString("too large shift in ") + pNew->GetName();
        }
      }
    }
    if (r->minNll() - chi2minGlobal > 25.) {
      cout << "MethodAbsScan::confirmSolutions() : WARNING : local minimum has DeltaChi2>25." << endl;
      isConfirmed = false;
      rejectReason =
          Form("too large chi2: DeltaChi2>25 - chi2minGlobal: %e and confirmed NLL: %e", chi2minGlobal, r->minNll());
    }
    if (isConfirmed) {
      if (arg->debug) cout << "MethodAbsScan::confirmSolutions() : solution " << i << " accepted." << endl;
      auto sr = make_unique<RooSlimFitResult>(r.get(), true);  // true saves correlation matrix
      sr->setConfirmed(true);
      confirmedSolutions.push_back(std::move(sr));
    } else {
      cout << std::format("MethodAbsScan::confirmSolutions() : WARNING : solution {:d} rejected ({:s})", i,
                          std::string(rejectReason))
           << endl;
    }
  }
  // do NOT delete the old solutions! They are still in allResults and curveResults.
  solutions = std::move(confirmedSolutions);
  sortSolutions();
  if (arg->debug) printLocalMinima();
  removeDuplicateSolutions();
  // reset parameters
  setParameters(w, parsName, frCache.getParsAtFunctionCall());
}

/**
 * Remove duplicate solutions from the common solutions storage ('solutions' vector).
 *
 * Duplicate solutions can occur when two unconfirmed solutions converge to the same true local minimum when refitted
 * by confirmSolutions(). No solutions will be removed if --qh 9 is given.
 */
void MethodAbsScan::removeDuplicateSolutions() {
  // TODO upgrade the quickhack to a proper option
  if (arg->isQuickhack(9)) return;
  vector<std::unique_ptr<RooSlimFitResult>> solutionsNoDup;
  for (int i = 0; i < solutions.size(); i++) {
    bool found = false;
    for (int j = i + 1; j < solutions.size(); j++) {
      if (compareSolutions(solutions[i].get(), solutions[j].get())) found = true;
      if (found == true) continue;
    }
    if (!found)
      solutionsNoDup.push_back(solutions[i]->Clone());
    else {
      if (arg->debug) cout << "MethodAbsScan::removeDuplicateSolutions() : removing duplicate solution " << i << endl;
    }
  }
  if (solutions.size() != solutionsNoDup.size()) {
    cout << endl;
    if (arg->debug) cout << "MethodAbsScan::removeDuplicateSolutions() : ";
    cout << "INFO : some equivalent solutions were removed. In case of 2D scans" << endl;
    cout << "       many equivalent solutions may lay on a contour of constant chi2, in" << endl;
    cout << "       that case removing them is perhaps not desired. You can keep all solutions" << endl;
    cout << "       using --qh 9\n" << endl;
  }
  solutions = std::move(solutionsNoDup);
}

///
/// Compare two solutions.
/// \param r1 First solution
/// \param r2 Second solution
/// \return true, if both are equal inside a certain margin
///
bool MethodAbsScan::compareSolutions(const RooSlimFitResult* r1, const RooSlimFitResult* r2) const {
  // compare chi2
  if (fabs(r1->minNll() - r2->minNll()) > 0.05) return false;
  // construct parameter lists
  RooArgList list1 = r1->floatParsFinal();
  list1.add(r1->constPars());
  RooArgList list2 = r2->floatParsFinal();
  list2.add(r2->constPars());
  // compare each parameter
  for (const auto p : *w->set(parsName)) {
    auto p1 = (RooRealVar*)list1.find(p->GetName());
    auto p2 = (RooRealVar*)list2.find(p->GetName());
    if (!p1 && !p2) {
      cout << "MethodAbsScan::compareSolutions() : ERROR : parameter not found: " << p->GetName() << endl;
      continue;
    }
    // We accept two parameters to be equal if they agree within 0.1 sigma.
    double sigma1 = p1->getError() > 0 ? p1->getError() : p1->getVal() / 10.;
    double sigma2 = p2->getError() > 0 ? p2->getError() : p2->getVal() / 10.;
    if (fabs(p1->getVal() - p2->getVal()) / (sqrt(sq(sigma1) + sq(sigma2))) > 0.1) return false;
  }
  return true;
}

///
/// Return a solution corresponding to a minimum of the profile
/// likelihoood.
/// \param i Index of the solution, they are orderd after increasing chi2,
///         i=0 is that with the smallest chi2.
///
RooSlimFitResult* MethodAbsScan::getSolution(const int i) {
  if (i < 0 || i >= solutions.size() || !solutions[i]) {
    cout << Form("MethodAbsScan::getSolution() : ERROR : No solution with id %i.", i) << endl;
    return nullptr;
  }
  return solutions[i].get();
}

///
/// Helper function to copy over solutions from another
/// scanner. Clears the solutions vector and sets the one
/// given.
///
void MethodAbsScan::setSolutions(const vector<std::unique_ptr<RooSlimFitResult>>& sols) {
  solutions.clear();
  for (auto&& sol : sols) solutions.emplace_back(sol->Clone());
}

///
/// Make a pull plot of observables corresponding
/// to the given solution.
///
void MethodAbsScan::plotPulls(int nSolution) {
  PullPlotter p(this);
  p.loadParsFromSolution(nSolution);
  p.savePulls();
  p.plotPulls();
}

void MethodAbsScan::setXscanRange(double min, double max) {
  if (min == max) return;
  RooRealVar* par1 = w->var(scanVar1);
  assert(par1);
  RooMsgService::instance().setGlobalKillBelow(ERROR);
  par1->setRange("scan", min, max);
  RooMsgService::instance().setGlobalKillBelow(INFO);
  if (arg->debug)
    std::cout << "DEBUG in MethodAbsScan::setXscanRange(): setting range for " << scanVar1 << ": " << min << ": " << max
              << std::endl;
  m_xrangeset = true;
}

void MethodAbsScan::setYscanRange(double min, double max) {
  if (min == max) return;
  RooRealVar* par2 = w->var(scanVar2);
  assert(par2);
  RooMsgService::instance().setGlobalKillBelow(ERROR);
  par2->setRange("scan", min, max);
  RooMsgService::instance().setGlobalKillBelow(INFO);
  m_yrangeset = true;
}

void MethodAbsScan::calcCLintervalsSimple(int CLsType, bool calc_expected) {
  clintervals.clear();
  const int nc = 3;
  clintervals.clear();
  clintervals.resize(nc);

  auto histogramCL = this->getHCL();
  if (this->hCLs && CLsType == 1) {
    histogramCL = this->getHCLs();
  } else if (this->hCLsFreq && CLsType == 2) {
    histogramCL = this->getHCLsFreq();
  }
  if (CLsType == 2 && calc_expected && hCLsExp) { histogramCL = this->getHCLsExp(); }
  if (CLsType == 0 || (this->hCLs && CLsType == 1) || (this->hCLsFreq && CLsType == 2)) {
    for (int c = 0; c < nc; c++) {
      const std::pair<double, double> borders = getBorders(TGraph(histogramCL), ConfidenceLevels[c]);
      auto cli = std::make_unique<CLInterval>();
      cli->pvalue = 1. - ConfidenceLevels[c];
      cli->min = borders.first;
      cli->max = borders.second;
      cli->central = -1;
      clintervals[c].push_back(std::move(cli));
      if (CLsType == 1) std::cout << "Simplified CL_s ";
      if (CLsType == 2) std::cout << "Standard CL_s";
      std::cout << "borders at " << ConfidenceLevels[c] << "    [ " << borders.first << " : " << borders.second << "]";
      cout << ", " << methodName << " (simple boundary scan)" << endl;
    }
  }
  ////////////////////////////////////////////////////////////////////////////////////////////
  //// Add a hacky calculation of the CL_s intervals
  //// TODO: Do it properly from the very start by introducing a bkg model and propagate it to the entire framework.
  /// TODO: I think this can now disappear

  if ((!this->hCLs && CLsType == 1) || (!this->hCLsFreq && CLsType == 2)) {
    std::cout << std::string(120, '*') << "\n"
              << "WARNING: hCLs is empty! Will calculate CLs intervals by normalising the p values to the p value of "
                 "the first bin\n"
                 "WARNING: This is only an approximate solution and MIGHT EVEN BE WRONG, if the first bin does not "
                 "represent the background expectation!\n"
              << std::string(120, '*') << std::endl;
    const int nc = 3;
    clintervals.clear();
    clintervals.resize(nc);

    for (int c = 0; c < nc; c++) {
      const std::pair<double, double> borders_CLs = getBorders_CLs(TGraph(histogramCL), ConfidenceLevels[c]);
      auto cli = std::make_unique<CLInterval>();
      cli->pvalue = 1. - ConfidenceLevels[c];
      cli->min = borders_CLs.first;
      cli->max = borders_CLs.second;
      cli->central = -1;
      clintervals[c].push_back(std::move(cli));
      std::cout << "CL_s borders at " << ConfidenceLevels[c] << "  [ " << borders_CLs.first << " : "
                << borders_CLs.second << "]";
      cout << ", " << methodName << " (simple boundary scan)" << endl;
    }
  }
}

/*!
\brief determines the borders of the confidence interval by linear or qubic interpolation.
\param graph The graph holding the p-values.
\param confidence_level The confidence level at which the interval is to be determined.
\param qubic Optional parameter. False by default. If true, qubic interpolation is used.
*/
const std::pair<double, double> MethodAbsScan::getBorders(const TGraph& graph, const double confidence_level,
                                                          bool qubic) const {

  const double p_val = 1 - confidence_level;
  TSpline* splines = nullptr;
  if (qubic) splines = new TSpline3();

  double min_edge = graph.GetX()[0];
  // will never return smaller edge than min_edge
  double max_edge = graph.GetX()[graph.GetN() - 1];
  // will never return higher edge than max_edge
  int scan_steps = 1000;
  double lower_edge = min_edge;
  double upper_edge = max_edge;

  for (double point = min_edge; point < max_edge; point += (max_edge - min_edge) / scan_steps) {

    if (graph.Eval(point, splines) > p_val) {
      lower_edge = point;
      break;
    }
  }
  for (double point = max_edge; point > min_edge; point -= (max_edge - min_edge) / scan_steps) {
    if (graph.Eval(point, splines) > p_val) {
      upper_edge = point;
      break;
    }
  }
  return std::pair<double, double>(lower_edge, upper_edge);
}

////////////////////////////////////////////////////////////////////////////////////////////
//// Do a hacky calculation of the CL_s intervals, where essentially the pValue is normalized to the pValue with
/// n_sig=0. / Let's first assume that the parameter of interest is ALWAYS a parameter correlated with n_sig, so that
/// parameter=0 means n_sig=0. / Therefore the pValue(CL_s) is given by the ratio of the pValue at scanpointand the
/// pValue of the lowest bin. / \todo Do it properly from the very start by introducing a bkg model and propagate it to
/// the entire framework.

const std::pair<double, double> MethodAbsScan::getBorders_CLs(const TGraph& graph, const double confidence_level,
                                                              bool qubic) const {

  const double p_val = 1 - confidence_level;
  TSpline* splines = nullptr;
  if (qubic) splines = new TSpline3();

  double min_edge = graph.GetX()[0];
  // will never return smaller edge than min_edge
  double max_edge = graph.GetX()[graph.GetN() - 1];
  // will never return higher edge than max_edge
  int scan_steps = 1000;
  double lower_edge = min_edge;
  double upper_edge = max_edge;

  for (double point = min_edge; point < max_edge; point += (max_edge - min_edge) / scan_steps) {

    // for CL_s normalize pVal to the pVal at 0 (which has to be the background model)
    if (graph.Eval(point, splines) / graph.Eval(min_edge, splines) > p_val) {
      lower_edge = point;
      break;
    }
  }
  for (double point = max_edge; point > min_edge; point -= (max_edge - min_edge) / scan_steps) {

    // for CL_s normalize pVal to the pVal at 0 (which has to be the background model)
    if (graph.Eval(point, splines) / graph.Eval(min_edge, splines) > p_val) {
      upper_edge = point;
      break;
    }
  }
  return std::pair<double, double>(lower_edge, upper_edge);
}

bool MethodAbsScan::checkCLs() const {
  if (!hCLsExp || !hCLsErr1Up || !hCLsErr1Dn || !hCLsErr2Up || !hCLsErr2Dn) {
    std::cout << "ERROR: ***************************************************" << std::endl;
    std::cout << "ERROR: MethodAbsScan::checkCLs() : No CLs plot available!!" << std::endl;
    std::cout << "ERROR: ***************************************************" << std::endl;
    return false;
  }
  assert(hCLsExp->GetNbinsX() == hCLsErr1Up->GetNbinsX());
  assert(hCLsExp->GetNbinsX() == hCLsErr2Up->GetNbinsX());
  assert(hCLsExp->GetNbinsX() == hCLsErr1Dn->GetNbinsX());
  assert(hCLsExp->GetNbinsX() == hCLsErr2Dn->GetNbinsX());

  // correct for low stats in the lower error
  for (int i = 1; i <= hCLsExp->GetNbinsX(); i++) {
    if (hCLsErr1Dn->GetBinContent(i) >= hCLsExp->GetBinContent(i)) {
      hCLsErr1Dn->SetBinContent(i, hCLsExp->GetBinContent(i) -
                                       (hCLsErr1Up->GetBinContent(i) - hCLsErr1Dn->GetBinContent(i)) / 2.);
    }
    if (hCLsErr1Dn->GetBinContent(i) >= hCLsExp->GetBinContent(i)) {
      hCLsErr1Dn->SetBinContent(i,
                                hCLsExp->GetBinContent(i) - (hCLsErr1Up->GetBinContent(i) - hCLsExp->GetBinContent(i)));
    }
    if (((hCLsExp->GetBinContent(i) - hCLsErr1Dn->GetBinContent(i)) / hCLsExp->GetBinContent(i)) < 0.05) {
      hCLsErr1Dn->SetBinContent(i, hCLsExp->GetBinContent(i) -
                                       (hCLsErr1Up->GetBinContent(i) - hCLsErr1Dn->GetBinContent(i)) / 2.);
    }
    if (hCLsErr2Dn->GetBinContent(i) >= hCLsExp->GetBinContent(i)) {
      hCLsErr2Dn->SetBinContent(i, hCLsExp->GetBinContent(i) -
                                       (hCLsErr2Up->GetBinContent(i) - hCLsErr2Dn->GetBinContent(i)) / 2.);
    }
    if (hCLsErr2Dn->GetBinContent(i) >= hCLsErr1Dn->GetBinContent(i)) {
      hCLsErr2Dn->SetBinContent(i, hCLsExp->GetBinContent(i) -
                                       (hCLsExp->GetBinContent(i) - hCLsErr1Dn->GetBinContent(i)) * 2.);
    }
  }
  return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//// end of CL_s part
////////////////////////////////////////////////////////////////////////////////////////////////////
