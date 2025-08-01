/*
 * Gamma Combination
 * Author: Matthew Kenzie, matthew.kenzie@cern.ch
 * Date: January 2016
 *
 */

#ifndef MethodCoverageScan_h
#define MethodCoverageScan_h

#include <TH1F.h>
#include <TString.h>
#include <TTree.h>

#include "MethodAbsScan.h"
#include "ParameterCache.h"

#include <memory>
#include <utility>
#include <vector>

class MethodCoverageScan : public MethodAbsScan {
 public:
  MethodCoverageScan(Combiner* comb);
  MethodCoverageScan() = delete;

  void setParameterCache(std::unique_ptr<ParameterCache> _pCache) { pCache = std::move(_pCache); }
  virtual int scan1d(int nRun = 1);
  virtual void readScan1dTrees(int runMin, int runMax);
  virtual void plot() const;
  int getNtoys() const { return nToys; };
  void saveScanner(TString fName = "");
  bool loadScanner(TString fName = "");

 protected:
  std::unique_ptr<ParameterCache> pCache;
  int nToys;  ///< number of toys to be generated at each scan point

  // functions
  std::vector<double> fitHist(TH1* h, TString fitfunc = "p1+exp", bool draw = true) const;
  double transform(std::vector<double> fitParams, TString transFunc, double x) const;
  void printLatexLine(double eta, double finProb, double finProbErr, double finPlug, double finPlugErr) const;

  // result histograms
  TH1F* h_sol = nullptr;
  TH1F* h_pvalue_plugin = nullptr;
  TH1F* h_pvalue_prob = nullptr;
  TH1F* h_pvalue_plugin_notransf = nullptr;
  TH1F* h_pvalue_prob_notransf = nullptr;
  TTree* t_res = nullptr;

  // result values
  Long64_t nentries;
  Long64_t nfailed;
  double n68plugin;
  double n95plugin;
  double n99plugin;
  double n68prob;
  double n95prob;
  double n99prob;
};

#endif
