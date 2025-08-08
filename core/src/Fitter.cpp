#include <Fitter.h>

#include <OptParser.h>
#include <Utils.h>

#include <RooFitResult.h>
#include <RooWorkspace.h>

#include <TString.h>

#include <cassert>
#include <iostream>

Fitter::Fitter(OptParser* arg, RooWorkspace* w, TString name) {
  this->w = w;
  this->name = name;
  this->arg = arg;

  pdfName = "pdf_" + name;
  obsName = "obs_" + name;
  parsName = "par_" + name;
}

///
/// Perform two fits, each time using different start parameters,
/// retain the smallest chi2. Note: To debug the start paramter
/// issues, remember that fitToMinBringBackAngles() may perform a
/// second fit, which then uses the start parameters of the first.
/// This will show up in the RooFitResult.
///
void Fitter::fitTwice() {
  // first fit
  Utils::setParametersFloating(w, parsName, startparsFirstFit);
  RooFitResult* r1 = Utils::fitToMinBringBackAngles(w->pdf(pdfName), false, -1);
  bool f1failed = !(r1->edm() < 1 && r1->covQual() == 3);

  // second fit
  Utils::setParametersFloating(w, parsName, startparsSecondFit);
  RooFitResult* r2 = Utils::fitToMinBringBackAngles(w->pdf(pdfName), false, -1);
  bool f2failed = !(r2->edm() < 1 && r2->covQual() == 3);

  if (f1failed && f2failed) {
    theResult = r1;
    delete r2;
  } else if (f1failed) {
    nFit2Best++;
    theResult = r2;
    delete r1;
  } else if (f2failed) {
    nFit1Best++;
    theResult = r1;
    delete r2;
  } else if (r1->minNll() < r2->minNll()) {
    nFit1Best++;
    theResult = r1;
    delete r2;
  } else {
    nFit2Best++;
    theResult = r2;
    delete r1;
  }

  Utils::setParametersFloating(w, parsName, theResult);
}

///
/// Force minimum finding. Will use the start parameters set by
/// setStartparsFirstFit().
///
void Fitter::fitForce() {
  Utils::setParametersFloating(w, parsName, startparsFirstFit);
  theResult = Utils::fitToMinForce(w, name);
  Utils::setParametersFloating(w, parsName, theResult);
}

///
/// Returns minimum chi2 value obtained by fit().
/// \return min chi2; 1e6 if fit wasn't performed yet
///
float Fitter::getChi2() const {
  if (!theResult) assert(0);
  if (theResult->minNll() < -10) return -10;  ///< else we have many entries at -1e27 in the ToyTree
  return theResult->minNll();
}

///
/// Check the fit result. We require a good covariance
/// matrix and a reasonable EDM for a status "ok".
/// \return Status code: 0=ok, 1=error, -1=fit() didn't run
///
int Fitter::getStatus() const {
  if (!theResult) return -1;
  if (theResult->floatParsFinal().getSize() == 0) return 0;
  if (theResult->edm() < 1 && theResult->status() == 0 && theResult->covQual() == 3) return 0;
  // theResult->Print("v");
  return 1;
}

///
/// Run the fit. Select a fit method here.
///
void Fitter::fit() {
  if (theResult) delete theResult;
  if (arg->scanforce)
    fitForce();
  else
    fitTwice();
}

void Fitter::print() const {
  std::cout << "Fitter: nFit1Best=" << nFit1Best << " nFit2Best=" << nFit2Best << std::endl;
}
