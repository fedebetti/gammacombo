/**
 * Gamma Combination
 * Author: Till Moritz Karbach, moritz.karbach@cern.ch
 * Date: August 2012
 *
 **/

#ifndef Fitter_h
#define Fitter_h

#include <RooArgSet.h>
#include <RooFitResult.h>
#include <RooWorkspace.h>

#include <TString.h>

class OptParser;

class Fitter {
 public:
  Fitter(OptParser* arg, RooWorkspace* w, TString name);

  void fit();
  void fitForce();
  void fitImprove();
  void fitOnce();
  void fitTwice();
  float getChi2() const;
  int getStatus() const;
  void print() const;
  inline void setStartpars(const RooArgSet* pars) { setStartparsFirstFit(pars); };
  inline void setStartparsFirstFit(const RooArgSet* pars) { startparsFirstFit = pars; };
  inline void setStartparsSecondFit(const RooArgSet* pars) { startparsSecondFit = pars; };

  OptParser* arg = nullptr;   ///< command line arguments
  RooWorkspace* w = nullptr;  ///< holds all input pdfs, parameters, and observables, as well as the combination
  TString name;               ///< Name of the pdf. Call combine() first.

  /// Start parameters to be used by all fit routines that run one fit, and by the first fit of fitTwice()
  RooArgSet const* startparsFirstFit = nullptr;
  /// Start parameters to be used by the second fit of fitTwice()
  RooArgSet const* startparsSecondFit = nullptr;

  int nFit1Best = 0;                  ///< counter, how many times did fit 1 of fitTwice() give smaller chi2
  int nFit2Best = 0;                  ///< counter, how many times did fit 2 of fitTwice() give smaller chi2
  TString pdfName;                    ///< PDF name in workspace, derived from name
  TString obsName;                    ///< dataset name of observables
  TString parsName;                   ///< set name of physics parameters
  RooFitResult* theResult = nullptr;  ///< the final result
};

#endif
