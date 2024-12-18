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

class MethodCoverageScan : public MethodAbsScan {
 public:
  MethodCoverageScan(Combiner* comb);
  MethodCoverageScan() = delete;

  void setParameterCache(ParameterCache* _pCache) { pCache = _pCache; }
  virtual int scan1d(int nRun = 1);
  virtual void readScan1dTrees(int runMin, int runMax);
  virtual void plot() const;
  int getNtoys() const { return nToys; };
  void saveScanner(TString fName = "");
  bool loadScanner(TString fName = "");

 protected:
  ParameterCache* pCache = nullptr;
  int nToys;  ///< number of toys to be generated at each scan point

  // functions
  std::vector<double> fitHist(TH1* h, TString fitfunc = "p1+exp", bool draw = true) const;
  double transform(std::vector<double> fitParams, TString transFunc, double x) const;
  void printLatexLine(float eta, float finProb, float finProbErr, float finPlug, float finPlugErr) const;

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
  float n68plugin;
  float n95plugin;
  float n99plugin;
  float n68prob;
  float n95prob;
  float n99prob;
};

#endif
