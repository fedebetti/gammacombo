/*
 * Gamma Combination
 * Author: Till Moritz Karbach, moritz.karbach@cern.ch
 * Date: August 2012
 *
 */

#ifndef MethodPluginScan_h
#define MethodPluginScan_h

#include <iostream>

#include <RooDataSet.h>

#include <TH1F.h>
#include <TString.h>

#include "Fitter.h"
#include "MethodAbsScan.h"
#include "MethodProbScan.h"
#include "PDF_Datasets.h"
#include "ProgressBar.h"
#include "ToyTree.h"

class MethodPluginScan : public MethodAbsScan {
 public:
  MethodPluginScan(MethodProbScan* s);
  MethodPluginScan(MethodProbScan* s, PDF_Datasets* pdf, const OptParser* opt);
  MethodPluginScan(Combiner* comb);

  inline void setNtoysPerPoint(int n) { nToys = n; };
  void setParevolPLH(MethodProbScan* s);
  virtual int scan1d(int nRun = 1);
  virtual void scan2d(int nRun = 1);
  virtual void readScan1dTrees(int runMin = 1, int runMax = 1, TString fName = "default");
  void readScan2dTrees(int runMin = 1, int runMax = 1);
  int getNtoys() const { return nToys; };
  double getPvalue1d(RooSlimFitResult* plhScan, double chi2minGlobal, ToyTree* t = 0, int id = 0, bool quiet = false);
  void makeControlPlotsCLs(std::map<int, std::vector<double>> bVals, std::map<int, std::vector<double>> sbVals) const;

 protected:
  TH1F* analyseToys(ToyTree* t, int id = -1, bool quiet = false);
  void computePvalue1d(RooSlimFitResult* plhScan, double chi2minGlobal, ToyTree* t, int id, Fitter* f, ProgressBar* pb);
  RooDataSet* generateToys(int nToys);
  double importance(double pvalue) const;
  RooSlimFitResult* getParevolPoint(float scanpoint);

  int nToys;                   ///< number of toys to be generated at each scan point
  MethodProbScan* profileLH;   ///< external scanner holding the profile likelihood: DeltaChi2 of the scan PDF on data
  MethodProbScan* parevolPLH;  ///< external scanner defining the parameter evolution: set to profileLH unless for the
                               ///< Hybrid Plugin
  RooDataSet* BkgToys = nullptr;
  std::vector<double> chi2minBkgBkgToysvector;     ///< saving the fits of the bkg-only pdf to the bkg-only toy
  std::vector<double> chi2minGlobalBkgToysvector;  ///< saving the fits of the global pdf to the bkg-only toy
};

#endif
