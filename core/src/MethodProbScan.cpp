#include <MethodProbScan.h>
#include <RooSlimFitResult.h>
#include <Utils.h>

#include <RooAddition.h>
#include <RooArgSet.h>
#include <RooDataHist.h>
#include <RooDataSet.h>
#include <RooGaussian.h>
#include <RooGlobalFunc.h>
#include <RooMultiVarGaussian.h>
#include <RooPlot.h>
#include <RooPoisson.h>
#include <RooProdPdf.h>
#include <RooRandom.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>

#include <TMarker.h>
#include <TMath.h>
#include <TStopwatch.h>
#include <TStyle.h>
#include <TSystem.h>

#include <format>
#include <iostream>
#include <limits>

using namespace std;
using namespace Utils;

MethodProbScan::MethodProbScan(Combiner* comb) : MethodAbsScan(comb) { methodName = "Prob"; }

MethodProbScan::MethodProbScan(OptParser* opt) : MethodAbsScan(opt) { methodName = "Prob"; }

///
/// Default constructor
///
MethodProbScan::MethodProbScan() { methodName = "Prob"; }

/**
 * Perform a 1d Prob scan.
 *
 * - Scan range defined through the range "scan".
 * - Fill the hCL histogram with the 1-CL curve.
 * - Start at a scan value that is in the middle of the allowed range (preferably a solution), and scan up and down
 *   from there.
 * - use the "probforce" command line flag to enable force minimum finding
 *
 * \param fast    This will scan each scanpoint only once.
 * \param reverse This will scan in reverse direction. When using the drag mode, this can sometimes make a difference.
 * \return The scan status (2 = new global minimum found, 1 = error, 0 = success)
 */
int MethodProbScan::scan1d(bool fast, const bool reverse, const bool quiet) {
  auto debug = [](const std::string& msg) { msgBase("MethodProbScan::scan1d() : DEBUG : ", msg); };
  auto warning = [](const std::string& msg) { msgBase("MethodProbScan::scan1d() : WARNING : ", msg); };
  auto error = [](const std::string& msg) {
    msgBase("MethodProbScan::scan1d() : ERROR : ", msg);
    std::exit(1);
  };

  if (arg->debug) debug(std::format("Call with arguments ({:s}, {:s}, {:s})", fast, reverse, quiet));

  nScansDone++;

  // The "improve" method doesn't need multiple scans.
  if (arg->probforce || arg->probimprove) fast = true;
  if (arg->probforce) scanDisableDragMode = true;

  // Save parameter values that were active at function call.
  startPars = std::make_unique<RooDataSet>("startPars", "startPars", *w->set(parsName));
  startPars->add(*w->set(parsName));

  // load scan parameter and scan range
  setLimit(w, scanVar1, "scan");
  auto par = w->var(scanVar1);
  if (!par) {
    std::cerr << std::format("Could not find the variable {:s} in the RooWorkspace. Exit", std::string(scanVar1))
              << std::endl;
    std::exit(1);
  }
  const auto min = hCL->GetXaxis()->GetXmin();
  const auto max = hCL->GetXaxis()->GetXmax();
  if (std::abs(par->getMin() - min) > 1e-6 || std::abs(par->getMax() - max) > 1e-6)
    warning("Scan range was changed after initScan() was called so the old range will be used.");
  if (arg->verbose) {
    std::cout << "\nProb configuration:" << std::endl;
    std::cout << "  combination : " << title << std::endl;
    std::cout << "  scan variable : " << scanVar1 << std::endl;
    std::cout << "  scan range : " << min << " ... " << max << std::endl;
    std::cout << "  scan steps : " << nPoints1d << std::endl;
    std::cout << "  fast mode : " << fast << std::endl;
    std::cout << std::endl;
  }

  // Set limit to all parameters.
  combiner->loadParameterLimits();

  // fix scan parameter
  par->setConstant(true);

  // Check if there is at least one parameter to be fit
  auto likelihood = w->pdf(pdfName);
  RooFormulaVar nll("nll", "nll", "-2*log(@0)", RooArgSet(*likelihood));
  auto floatPars = std::unique_ptr<RooAbsCollection>(likelihood->getVariables()->selectByAttrib("Constant", false));
  const bool hasFreePars = floatPars->getSize() > 0;
  if (!hasFreePars)
    std::cout << "MethodProbScan::scan1d() : INFO : There are no free parameters. I will scan without fitting"
              << std::endl;

  // Report on the smallest new minimum we come across while scanning.
  // Sometimes the scan doesn't find the minimum that was found before. Warn if this happens.
  const auto bestMinOld = chi2minGlobal;
  auto bestMinFoundInScan = std::numeric_limits<double>::max();
  curveResults.clear();
  curveResults.resize(nPoints1d, nullptr);

  // Scan in the following directions j = 0 : start value -> upper limit
  //                                      1 : upper limit -> start value
  //                                      2 : lower limit -> start value
  //                                      3 : start value -> lower limit
  for (int j = 0; j < 4; j++) {
    if (fast && (j == 1 || j == 2)) continue;
    if (reverse) j = 3 - j;

    double startValue = par->getVal();
    bool scanUp = (j == 0 || j == 2);

    // Improve method doesn't work with drag mode as parameter run at their limits
    if (j == 0 || j == 3 || scanDisableDragMode) setParameters(w, parsName, startPars->get(0));
    double scanStart, scanStop;
    switch (j) {
    case 0:
      scanStart = startValue;
      scanStop = par->getMax();
      break;
    case 1:
      scanStart = par->getMax();
      scanStop = startValue;
      break;
    case 2:
      scanStart = par->getMin();
      scanStop = startValue;
      break;
    case 3:
      scanStart = startValue;
      scanStop = par->getMin();
      break;
    }

    // for the status bar
    auto nTotalSteps = nPoints1d;
    nTotalSteps *= fast ? 1 : 2;
    const double printFreq = nTotalSteps > 15 ? 10. : static_cast<double>(nTotalSteps);

    for (int i = 0; i < nPoints1d; ++i) {
      double scanvalue;
      if (scanUp) {
        scanvalue = min + (max - min) * i / nPoints1d + hCL->GetBinWidth(1) / 2.;
        if (scanvalue < scanStart) continue;
        if (scanvalue > scanStop) break;
      } else {
        scanvalue = max - (max - min) * i / nPoints1d - hCL->GetBinWidth(1) / 2.;
        if (scanvalue > scanStart) continue;
        if (scanvalue < scanStop) break;
      }

      // don't scan in unphysical region
      if (scanvalue < par->getMin() || scanvalue > par->getMax()) continue;

      // set the parameter of interest to the scan point
      par->setVal(scanvalue);

      // status bar
      if (!quiet && (i % static_cast<int>(nTotalSteps / printFreq)) == 0)
        std::cout << "MethodProbScan::scan1d() : scanning " << 100. * i / nTotalSteps << "%   \r" << flush;

      RooSlimFitResult* sfr = nullptr;
      auto chi2minScan = std::numeric_limits<double>::max();
      auto performFit = hasFreePars;
      if (!hasFreePars) {
        // There are no parameters to fit, just calculate NLL
        chi2minScan = nll.getVal();
        // If we found indications of a new global minimum, perform the fit to find it!
        if (chi2minScan < chi2minGlobal) performFit = true;
      }
      if (performFit) {
        if (!hasFreePars) par->setConstant(false);
        std::unique_ptr<RooFitResult> fr = nullptr;
        if (arg->probforce)
          fr = fitToMinForce(w, combiner->getPdfName());
        else if (arg->probimprove)
          fr = fitToMinImprove(w, combiner->getPdfName());
        else
          fr = fitToMinBringBackAngles(likelihood, false, -1);
        chi2minScan = fr->minNll();
        if (std::isinf(chi2minScan))
          chi2minScan = std::numeric_limits<double>::max();  // else the toys in PDF_testConstraint don't work
        allResults.push_back(std::make_unique<RooSlimFitResult>(fr.get()));  // save memory by using the slim fit result
        sfr = allResults.back().get();
        if (!hasFreePars) par->setConstant(true);
      }
      if (chi2minScan < 0) {
        error("I have found a new minimum with negative chi2. This should never happen");
        // TODO are there reasons why we shouldn't exit in the previous command?
        // double newChi2minScan = chi2minGlobal + 25.;  // 5sigma more than best point
        // std::cout << "MethodProbScan::scan1d() : WARNING : " << title
        //      << std::format(" chi2 negative for scan point {:d}: {:.2f} setting to: {:.2f}", i, chi2minScan,
        //                     newChi2minScan)
        //      << std::endl;
        // chi2minScan = newChi2minScan;
      }
      bestMinFoundInScan = std::min(chi2minScan, bestMinFoundInScan);

      // If we find a minimum smaller than the old "global" minimum, this means that all previous 1-CL values are too
      // high.
      if (chi2minScan < chi2minGlobal) {
        if (arg->verbose)
          warning(std::format("'{:s}' new global minimum found! chi2minScan={:.3f}", std::string(title), chi2minScan));
        chi2minGlobal = chi2minScan;
        // recompute previous 1-CL values
        for (int k = 1; k <= hCL->GetNbinsX(); k++) {
          hCL->SetBinContent(k, TMath::Prob(hChi2min->GetBinContent(k) - chi2minGlobal, 1));
        }
      }

      const auto deltaChi2 = chi2minScan - chi2minGlobal;
      const auto oneMinusCL = TMath::Prob(deltaChi2, 1);
      auto deltaChi2Bkg = std::max(chi2minScan - hChi2min->GetBinContent(1), 0.);
      if (i == 0) deltaChi2Bkg = 0.;
      const auto oneMinusCLBkg = TMath::Prob(deltaChi2Bkg, 1);
      hCLs->SetBinContent(hCLs->FindBin(scanvalue), oneMinusCLBkg);

      if (i == 0) chi2minBkg = chi2minScan;

      // Save the 1-CL value and the corresponding fit result. But only if better than before!
      if (const auto iBin = hCL->FindBin(scanvalue); hCL->GetBinContent(iBin) <= oneMinusCL) {
        if (arg->debug)
          debug(std::format("Setting bin {:3d} to {:.4f} (chi2minScan: {:.4f})", iBin, oneMinusCL, chi2minScan));
        hCL->SetBinContent(iBin, oneMinusCL);
        hChi2min->SetBinContent(iBin, chi2minScan);
        if (sfr) { curveResults[iBin - 1] = sfr; }
      }
    }
  }
  std::cout << "MethodProbScan::scan1d() : scan done.           " << std::endl;

  if (bestMinFoundInScan - bestMinOld > 1e-2) {
    warning(std::format("Scan didn't find similar minimum to what was found before!\n"
                        "Too strict parameter limits? Too coarse scan steps? Didn't load global minimum?\n"
                        "chi2 bestMinFoundInScan={:.3f}, bestMinOld={:.3f}",
                        bestMinFoundInScan, bestMinOld));
  }

  // attempt to correct for undercoverage TODO
  if (pvalueCorrector) {
    for (int k = 1; k <= hCL->GetNbinsX(); k++) {
      double pvalueProb = hCL->GetBinContent(k);
      pvalueProb = pvalueCorrector->transform(pvalueProb);
      hCL->SetBinContent(k, pvalueProb);
    }
  }

  setParameters(w, parsName, startPars->get(0));

  saveSolutions();

  if (arg->confirmsols) confirmSolutions();

  if ((bestMinFoundInScan - bestMinOld) / bestMinOld > 1e-2) return 1;
  return 0;
}

///
/// Modify the CL histograms according to the defined test statistics
/// \return status 0 ->potentially can encode debug information here
///
int MethodProbScan::computeCLvalues() {
  auto info = [](const std::string& msg) { msgBase("MethodProbScan::computeCLvalues() : INFO : ", msg); };
  info("Computing CL values based on test statistic decision\n" +
       std::format("using {:d}-sided test statistic", arg->teststatistic));

  if (!globalMin && !this->getSolution()) {
    std::cout << "Could not find a solution so can't redefine the test statistics appropriately" << std::endl;
    return 1;
  }

  double bestfitpoint;
  double bestfitpointerr;

  if (globalMin) {
    bestfitpoint = ((RooRealVar*)globalMin->floatParsFinal().find(scanVar1))->getVal();
    bestfitpointerr = ((RooRealVar*)globalMin->floatParsFinal().find(scanVar1))->getError();
  }

  if (this->getSolution()) {
    if (this->getSolution()->floatParsFinal().find(scanVar1)) {
      bestfitpoint = ((RooRealVar*)this->getSolution()->floatParsFinal().find(scanVar1))->getVal();
      bestfitpointerr = ((RooRealVar*)this->getSolution()->floatParsFinal().find(scanVar1))->getError();
    }
  }

  for (int k = 1; k <= hCL->GetNbinsX(); k++) {
    double scanvalue = hChi2min->GetBinCenter(k);
    double teststat_measured = hChi2min->GetBinContent(k) - chi2minGlobal;
    double CLb = 1. - (normal_cdf(TMath::Sqrt(teststat_measured) + ((scanvalue - 0.) / bestfitpointerr)) +
                       normal_cdf(TMath::Sqrt(teststat_measured) - ((scanvalue - 0.) / bestfitpointerr)) - 1.);
    if (arg->teststatistic == 1) {                                             // use one-sided test statistic
      teststat_measured = bestfitpoint <= scanvalue ? teststat_measured : 0.;  // if mu < muhat then q_mu = 0
      hCL->SetBinContent(k, 1. - normal_cdf(TMath::Sqrt(teststat_measured)));
      // if (scanvalue < bestfitpoint) hCL->SetBinContent(k,1.0);    // should not be here, but looks ugly if solution
      // is drawn as p=1
      CLb = 1. - normal_cdf(TMath::Sqrt(teststat_measured) - ((scanvalue - 0.) / bestfitpointerr));
    }
    if (!arg->cls.empty()) hCLs->SetBinContent(k, min(1., hCL->GetBinContent(k) / CLb));
    // std::cout << "MethodProbScan::" << k << "\t" << hCL->GetBinContent(k) << "\t" << CLb << "\t" <<
    // hCLs->GetBinContent(k) <<std::endl;
  }

  return 0;
}

/**
 * Delete a pointer if it is not included in the curveResults2d std::vector.
 * Also removes it from the allResults std::vector by setting the entry to 0.
 *
 * :return: true if r was deleted, or if it is 0
 */
bool MethodProbScan::deleteIfNotInCurveResults2d(const RooSlimFitResult* r) {
  if (!r) return true;
  bool del = true;
  for (int j = 0; j < hCL2d->GetNbinsX(); j++)
    for (int k = 0; k < hCL2d->GetNbinsY(); k++) {
      if (r == curveResults2d[j][k]) {
        del = false;
        break;
      }
    }
  if (del) {
    for (auto&& result : allResults) {
      if (r == result.get()) result = nullptr;
    }
  }
  return del;
}

///
/// Compute coordinates that can be used as coordinates
/// of that fit result, that we want to use as start parameters
/// in the drag mode. The i,j coordinates spiral out. This function
/// computes coordinates from inner turns of the spiral.
/// If the inner turn doesn't exist, the start coordinates are set
/// as the result, and false is returned.
/// \param iStart center of the spiral
/// \param jStart center of the spiral
/// \param i current coordinates on the outmost turn of the spiral
/// \param j current coordinates on the outmost turn of the spiral
/// \param iResult resulting coordinates of an inner turn
/// \param jResult resulting coordinates of an inner turn
/// \param nTurn number of inner turn to jump to, 1=next-to-outer turn, 2=second-next, etc.
///
bool MethodProbScan::computeInnerTurnCoords(const int iStart, const int jStart, const int i, const int j, int& iResult,
                                            int& jResult, int nTurn) const {
  // compute bin coordinates of start parameters: connect center of
  // the spiral to the scan point with a straight line, go back by sqrt(2)
  // units, take bin this ends us in
  iResult = iStart;
  jResult = jStart;
  if (sq(i - iStart) + sq(j - jStart) > 0) {
    iResult = round((double)i - double(nTurn) * 1.41 * double(i - iStart) / sqrt(sq(i - iStart) + sq(j - jStart)));
    jResult = round((double)j - double(nTurn) * 1.41 * double(j - jStart) / sqrt(sq(i - iStart) + sq(j - jStart)));
  }
  if (iResult - 1 >= curveResults2d.size()) iResult = iStart;
  if (jResult - 1 >= curveResults2d[0].size()) jResult = jStart;
  // check result
  if (iResult - 1 >= curveResults2d.size() || jResult - 1 >= curveResults2d[0].size() || iResult - 1 < 0 ||
      jResult - 1 < 0) {
    std::cout << "MethodProbScan::computeInnerTurnCoords() : ERROR : resulting coordinates out of range! "
              << iResult - 1 << " " << jResult - 1 << std::endl;
  }
  if (iResult == iStart && jResult == jStart) return false;
  return true;
}

void MethodProbScan::sanityChecks() const {
  if (!w->set(parsName)) {
    std::cout << "MethodProbScan::sanityChecks() : ERROR : parsName not found: " << parsName << std::endl;
    std::exit(1);
  }
  if (!w->var(scanVar1)) {
    std::cout << "MethodProbScan::sanityChecks() : ERROR : scanVar1 not found: " << scanVar1 << std::endl;
    std::exit(1);
  }
  if (!w->var(scanVar2)) {
    std::cout << "MethodProbScan::sanityChecks() : ERROR : scanVar2 not found: " << scanVar2 << std::endl;
    std::exit(1);
  }
}

/**
 * Perform a 2d Prob scan.
 *
 * - Scan range is defined by the range "scan".
 * - Fill the hCL2d histogram with the 1-CL curve.
 * - Save all encountered fit results to allResults.
 * - Save the fit results that make it into the 1-CL curve into curveResults2d.
 * - Scan strategy: Spiral out!
 *
 * \return The scan status (2 = new global minimum found, 1 = error, 0 = success)
 */
int MethodProbScan::scan2d() {
  if (arg->debug) std::cout << "MethodProbScan::scan2d() : starting ..." << std::endl;
  nScansDone++;
  sanityChecks();

  // Define whether the 2d contours in hCL are "1D sigma" (ndof=1) or "2D sigma" (ndof=2).
  // Leave this at 1 for now, as the "2D sigma" contours are computed from hChi2min2d, not hCL.
  const int ndof = 1;

  // Set up storage for fit results of this particular scan. This is used for the drag start parameters.
  // We cannot use the curveResults2d member because that only holds better results.
  std::vector<vector<RooSlimFitResult*>> mycurveResults2d(nPoints2dx,
                                                          std::vector<RooSlimFitResult*>(nPoints2dy, nullptr));

  // store start parameters so we can reset them later
  startPars = std::make_unique<RooDataSet>("startPars", "startPars", *w->set(parsName));
  startPars->add(*w->set(parsName));

  auto par1 = w->var(scanVar1);
  auto par2 = w->var(scanVar2);
  if (!par1 || !par2) {
    std::cerr << std::format("Could not find the variables ({:s}, {:s}) in the RooWorkspace. Exit",
                             std::string(scanVar1), std::string(scanVar2))
              << std::endl;
    std::exit(1);
  }

  // Set limit to all parameters.
  combiner->loadParameterLimits();

  // fix scan parameters
  par1->setConstant(true);
  par2->setConstant(true);

  // Check if there is at least one parameter to be fit
  auto likelihood = w->pdf(pdfName);
  RooFormulaVar nll("nll", "nll", "-2*log(@0)", RooArgSet(*likelihood));
  auto floatPars = std::unique_ptr<RooAbsCollection>(likelihood->getVariables()->selectByAttrib("Constant", false));
  const bool hasFreePars = floatPars->getSize() > 0;
  if (!hasFreePars)
    std::cout << "MethodProbScan::scan1d() : INFO : There are not free parameters. I will scan without fitting"
              << std::endl;

  // initialize some control plots
  gStyle->SetOptTitle(1);
  auto cDbg = newNoWarnTCanvas(getUniqueRootName(), Form("DeltaChi2 for 2D scan %i", nScansDone));
  cDbg->SetMargin(0.1, 0.15, 0.1, 0.1);
  double hChi2min2dMin = hChi2min2d->GetMinimum();
  bool firstScanDone = hChi2min2dMin < std::numeric_limits<double>::max();
  auto hDbgChi2min2d =
      histHardCopy(hChi2min2d.get(), firstScanDone, true, TString(hChi2min2d->GetName()) + TString("_Dbg"));
  hDbgChi2min2d->SetTitle(Form("#Delta#chi^{2} for scan %i, %s", nScansDone, title.Data()));
  if (firstScanDone) hDbgChi2min2d->GetZaxis()->SetRangeUser(hChi2min2dMin, hChi2min2dMin + 81);
  hDbgChi2min2d->GetXaxis()->SetTitle(par1->GetTitle());
  hDbgChi2min2d->GetYaxis()->SetTitle(par2->GetTitle());
  hDbgChi2min2d->GetZaxis()->SetTitle("#Delta#chi^{2}");
  auto hDbgStart = histHardCopy(hChi2min2d.get(), false, true, TString(hChi2min2d->GetName()) + TString("_DbgSt"));

  // Report on the smallest new minimum we come across while scanning.
  // Sometimes the scan doesn't find the minimum that was found before. Warn if this happens.
  const auto bestMinOld = chi2minGlobal;
  auto bestMinFoundInScan = std::numeric_limits<double>::max();

  // for the status bar
  int nSteps = 0;
  double nTotalSteps = nPoints2dx * nPoints2dy;
  const double printFreq = nTotalSteps > 100 && !arg->probforce ? 100 : nTotalSteps;

  // Start coordinates. Don't allow the under/overflow bins
  int iStart = min(hCL2d->GetXaxis()->FindBin(par1->getVal()), hCL2d->GetNbinsX());
  int jStart = min(hCL2d->GetYaxis()->FindBin(par2->getVal()), hCL2d->GetNbinsY());
  iStart = max(iStart, 1);
  jStart = max(jStart, 1);
  hDbgStart->SetBinContent(iStart, jStart, 500.);
  TMarker startpointmark(par1->getVal(), par2->getVal(), 3);

  // timer
  TStopwatch tFit;
  TStopwatch tSlimResult;
  TStopwatch tScan;
  TStopwatch tMemory;

  // set up the scan spiral
  int X = 2 * nPoints2dx;
  int Y = 2 * nPoints2dy;
  int x = 0;
  int y = 0;
  int dx = 0;
  int dy = -1;
  auto t = std::max(X, Y);
  auto maxI = t * t;
  for (int spiralstep = 0; spiralstep < maxI; spiralstep++) {
    if ((-X / 2 <= x) && (x <= X / 2) && (-Y / 2 <= y) && (y <= Y / 2)) {
      int i = x + iStart;
      int j = y + jStart;
      if (i > 0 && i <= nPoints2dx && j > 0 && j <= nPoints2dy) {
        tScan.Start(false);

        if (nSteps % static_cast<int>(nTotalSteps / printFreq) == 0)
          std::cout << std::format("MethodProbScan::scan2d() : scanning {:3.0f}%\r", 100. * nSteps / nTotalSteps)
                    << flush;
        nSteps++;

        // status histogram
        if (spiralstep > 0) hDbgStart->SetBinContent(i, j, 500. /*firstScan ? 1. : hChi2min2dMin+36*/);

        // set start parameters from inner turn of the spiral
        int xStartPars, yStartPars;
        computeInnerTurnCoords(iStart, jStart, i, j, xStartPars, yStartPars, 1);
        auto rStartPars = mycurveResults2d[xStartPars - 1][yStartPars - 1];
        if (rStartPars) setParameters(w, parsName, rStartPars);

        // memory management:
        tMemory.Start(false);
        // delete old, inner fit results, that we don't need for start parameters anymore
        // for this we take the second-inner-most turn.
        int iOld, jOld;
        bool innerTurnExists = computeInnerTurnCoords(iStart, jStart, i, j, iOld, jOld, 2);
        if (innerTurnExists) {
          deleteIfNotInCurveResults2d(mycurveResults2d[iOld - 1][jOld - 1]);
          mycurveResults2d[iOld - 1][jOld - 1] = nullptr;
        }
        tMemory.Stop();

        // alternative choice for start parameters: always from what we found at function call
        // setParameters(w, parsName, startPars->get(0));

        // set scan point
        auto scanvalue1 = hCL2d->GetXaxis()->GetBinCenter(i);
        auto scanvalue2 = hCL2d->GetYaxis()->GetBinCenter(j);
        par1->setVal(scanvalue1);
        par2->setVal(scanvalue2);

        RooSlimFitResult* sfr = nullptr;
        auto chi2minScan = std::numeric_limits<double>::max();
        auto performFit = hasFreePars;
        if (!hasFreePars) {
          // There are no parameters to fit, just calculate NLL
          chi2minScan = nll.getVal();
          // If we found another minimum, perform the fit!
          if (chi2minScan < bestMinFoundInScan) performFit = true;
        }
        if (performFit) {
          if (!hasFreePars) {
            par1->setConstant(false);
            par2->setConstant(false);
          }
          tFit.Start(false);
          std::unique_ptr<RooFitResult> fr;
          if (arg->probforce)
            fr = fitToMinForce(w, combiner->getPdfName());
          else
            fr = fitToMinBringBackAngles(likelihood, false, -1);
          chi2minScan = fr->minNll();
          if (std::isinf(chi2minScan))
            chi2minScan = std::numeric_limits<double>::max();  // else the toys in PDF_testConstraint don't work
          tFit.Stop();
          tSlimResult.Start(false);
          allResults.push_back(std::make_unique<RooSlimFitResult>(fr.get()));  // save memory by using the slim fit
                                                                               // result
          sfr = allResults.back().get();
          mycurveResults2d[i - 1][j - 1] = sfr;
          tSlimResult.Stop();
          if (!hasFreePars) {
            par1->setConstant(true);
            par2->setConstant(true);
          }
        }
        bestMinFoundInScan = std::min(chi2minScan, bestMinFoundInScan);

        // If we find a new global minumum, this means that all previous 1-CL values are too high.
        // We'll save the new possible solution, adjust the global minimum, return a status code, and stop. TODO
        if (chi2minScan > -500 && chi2minScan < chi2minGlobal) {
          if (arg->debug || chi2minScan < chi2minGlobal - 1e-2) {
            if (arg->verbose)
              std::cout << "MethodProbScan::scan2d() : WARNING : '" << title
                        << "' new global minimum found! chi2minGlobal=" << chi2minGlobal
                        << " chi2minScan=" << chi2minScan << std::endl;
          }
          chi2minGlobal = chi2minScan;
          // recompute previous 1-CL values
          for (int k = 1; k <= hCL2d->GetNbinsX(); k++)
            for (int l = 1; l <= hCL2d->GetNbinsY(); l++) {
              hCL2d->SetBinContent(k, l, TMath::Prob(hChi2min2d->GetBinContent(k, l) - chi2minGlobal, ndof));
            }
        }

        const double deltaChi2 = chi2minScan - chi2minGlobal;
        const double oneMinusCL = TMath::Prob(deltaChi2, ndof);

        // Save the 1-CL value. But only if better than before!
        if (hCL2d->GetBinContent(i, j) < oneMinusCL) {
          hCL2d->SetBinContent(i, j, oneMinusCL);
          hChi2min2d->SetBinContent(i, j, chi2minScan);
          hDbgChi2min2d->SetBinContent(i, j, chi2minScan);
          curveResults2d[i - 1][j - 1] = sfr;
        }

        // Draw/update histograms. Do that only every nth update depending on value of updateFreq.
        // This saves a lot of time for small combinations
        if ((arg->interactive && ((int)nSteps % arg->updateFreq == 0)) || nSteps == nTotalSteps) {
          hDbgChi2min2d->Draw("colz");
          hDbgStart->Draw("boxsame");
          startpointmark.Draw();
          cDbg->Update();
          cDbg->Modified();
          gSystem->ProcessEvents();
        }
        tScan.Stop();
      }
    }
    // spiral stuff:
    if ((x == y) || ((x < 0) && (x == -y)) || ((x > 0) && (x == 1 - y))) {
      t = dx;
      dx = -dy;
      dy = t;
    }
    x += dx;
    y += dy;
  }
  std::cout << "MethodProbScan::scan2d() : scan done.            " << std::endl;
  if (arg->debug) {
    std::cout << "MethodProbScan::scan2d() : full scan time:             ";
    tScan.Print();
    std::cout << "MethodProbScan::scan2d() : - fitting:                  ";
    tFit.Print();
    std::cout << "MethodProbScan::scan2d() : - create RooSlimFitResults: ";
    tSlimResult.Print();
    std::cout << "MethodProbScan::scan2d() : - memory management:        ";
    tMemory.Print();
  }
  setParameters(w, parsName, startPars->get(0));
  saveSolutions2d();
  if (arg->debug) printLocalMinima();
  if (arg->confirmsols) confirmSolutions();

  // clean all fit results that didn't make it into the final result
  for (auto&& result : allResults) { deleteIfNotInCurveResults2d(result.get()); }

  if (bestMinFoundInScan - bestMinOld > 0.1) {
    std::cout << "MethodProbScan::scan2d() : WARNING: Scan didn't find minimum that was found before!" << std::endl;
    std::cout << "MethodProbScan::scan2d() :          Are you using too strict parameter limits?" << std::endl;
    std::cout << "MethodProbScan::scan2d() :          min chi2 found in scan: " << bestMinFoundInScan
              << ", old min chi2: " << bestMinOld << std::endl;
    return 1;
  }
  return 0;
}

/**
 * Find the RooSlimFitResults corresponding to all local minima from the curveResults std::vector and save them into the
 * solutions std::vector. It will be sorted for minChi2. Index 0 will correspond to the least chi2.
 */
void MethodProbScan::saveSolutions() {
  auto warning = [](const std::string& msg) { msgBase("MethodProbScan::saveSolutions() : WARNING : ", msg); };

  if (arg->debug) std::cout << "MethodProbScan::saveSolutions() : searching for minima in hChi2min ..." << std::endl;

  if (std::ranges::any_of(curveResults, [](const RooSlimFitResult* sfr) { return sfr != nullptr; })) {
    solutions.clear();

    // loop over chi2 histogram to locate local maxima
    for (int i = 2; i < hChi2min->GetNbinsX() - 1; i++) {
      bool oneBinMax = hChi2min->GetBinContent(i - 1) > hChi2min->GetBinContent(i) &&
                       hChi2min->GetBinContent(i + 1) > hChi2min->GetBinContent(i);
      bool twoBinMax = hChi2min->GetBinContent(i - 1) > hChi2min->GetBinContent(i) &&
                       hChi2min->GetBinContent(i) == hChi2min->GetBinContent(i + 1) &&
                       hChi2min->GetBinContent(i + 2) > hChi2min->GetBinContent(i + 1);
      if (!(oneBinMax || twoBinMax)) continue;

      // loop over fit results to find those that produced it
      for (int j = 0; j < curveResults.size(); j++) {
        if (!curveResults[j]) {
          if (arg->debug) warning(std::format("Empty solution at index {:d}", j));
          continue;
        }

        if (hChi2min->FindBin(curveResults[j]->getConstParVal(scanVar1)) == i) {
          if (arg->debug) {
            std::cout << "MethodProbScan::saveSolutions() : saving solution " << j << ":" << std::endl;
            curveResults[j]->Print();
          }
          solutions.push_back(curveResults[j]->Clone());
        }
      }
    }
  } else {
    // This could be due to the fact that there was only one parameter.
    warning("No solutions found during the scan\n"
            "This is normal if there are only 1(2) free parameters in the baseline fit and the scan is 1d (2d)");
  }
  sortSolutions();

  if (solutions.empty()) {
    if (!globalMin) {
      std::cerr << "MethodProbScan::saveSolutions() : ERROR : No solution was found at all. Exit" << std::endl;
      std::exit(1);
    }
    solutions.emplace_back(std::make_unique<RooSlimFitResult>(globalMin.get()));
  }
}

///
/// Find the RooFitResults corresponding to all local
/// minima from the curveResults2d std::vector
/// and save them into the solutions std::vector.
/// It will be sorted for minChi2. Index 0
/// will correspond to the least chi2.
///
/// We use a brute force minimum finding over the histogram,
/// as it is not large: at each bin, inpsect the 8 surrounding
/// ones.
///
void MethodProbScan::saveSolutions2d() {
  if (arg->debug)
    std::cout << "MethodProbScan::saveSolutions2d() : searching for minima in hChi2min2d ..." << std::endl;

  if (std::ranges::any_of(curveResults2d, [](const std::vector<RooSlimFitResult*> v) {
        return std::ranges::any_of(v, [](const RooSlimFitResult* sfr) { return sfr != nullptr; });
      })) {
    solutions.clear();

    // loop over chi2 histogram to locate local minima
    for (int i = 2; i < hChi2min2d->GetNbinsX(); i++) {
      for (int j = 2; j < hChi2min2d->GetNbinsY(); j++) {
        if (!(hChi2min2d->GetBinContent(i - 1, j) > hChi2min2d->GetBinContent(i, j) &&
              hChi2min2d->GetBinContent(i + 1, j) > hChi2min2d->GetBinContent(i, j) &&
              hChi2min2d->GetBinContent(i, j - 1) > hChi2min2d->GetBinContent(i, j) &&
              hChi2min2d->GetBinContent(i, j + 1) > hChi2min2d->GetBinContent(i, j) &&
              hChi2min2d->GetBinContent(i - 1, j - 1) > hChi2min2d->GetBinContent(i, j) &&
              hChi2min2d->GetBinContent(i + 1, j + 1) > hChi2min2d->GetBinContent(i, j) &&
              hChi2min2d->GetBinContent(i + 1, j - 1) > hChi2min2d->GetBinContent(i, j) &&
              hChi2min2d->GetBinContent(i - 1, j + 1) > hChi2min2d->GetBinContent(i, j)))
          continue;

        const auto r = curveResults2d[i - 1][j - 1];  // -1 because it starts counting at 0, but histograms at 1
        if (!r) {
          std::cout << std::format("MethodProbScan::saveSolutions2d() : ERROR : No corresponding RooFitResult found! "
                                   "Skipping (i,j)=({:d},{:d})",
                                   i, j)
                    << std::endl;
          continue;
        }
        if (arg->debug)
          std::cout << std::format("MethodProbScan::saveSolutions2d() : saving solution of bin ({:d},{:d})...", i, j)
                    << std::endl;
        solutions.push_back(curveResults2d[i - 1][j - 1]->Clone());
      }
    }
  } else {
    std::cout << "MethodProbScan::saveSolutions2d() : WARNING : No solutions found in 2D scan!\n"
                 "  This can happen when a solution is too close to the plot boundary.\n"
                 "  In this case, either change the scan range using --scanrange and --scanrangey,\n"
                 "  or increase the number of scan points, using --npoints or --npoints2dx, --npoints2dy"
              << std::endl;
  }

  if (solutions.empty()) {
    if (!globalMin) {
      std::cerr << "MethodProbScan::saveSolutions2d() : ERROR : No solution was found at all. Exit" << std::endl;
      std::exit(1);
    }
    solutions.emplace_back(std::make_unique<RooSlimFitResult>(globalMin.get()));
  } else {
    sortSolutions();
  }
}

///
/// Get the chi2 value of the profile likelihood at a given scan point.
///
/// Requires that the scan was performed before by scan1d().
///
double MethodProbScan::getChi2min(double scanpoint) const {
  assert(hChi2min);
  int iBin = hChi2min->FindBin(scanpoint);
  return hChi2min->GetBinContent(iBin);
}
