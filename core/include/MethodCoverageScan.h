/*
 * Gamma Combination
 * Author: Matthew Kenzie, matthew.kenzie@cern.ch
 * Date: January 2016
 *
 */

#ifndef MethodCoverageScan_h
#define MethodCoverageScan_h

#include "MethodAbsScan.h"

#include <TString.h>

#include <vector>

class Combiner;
class ParameterCache;

class TH1;
class TH1F;
class TTree;

class MethodCoverageScan : public MethodAbsScan {
 public:
  MethodCoverageScan(Combiner* comb);
  MethodCoverageScan() = delete;

  void setParameterCache(ParameterCache* _pCache) { pCache = _pCache; }
  virtual int scan1d(int nRun = 1);
  virtual void readScan1dTrees(int runMin, int runMax);
  virtual void plot();
  int getNtoys() const { return nToys; };
  void saveScanner(TString fName = "");
  bool loadScanner(TString fName = "");

 protected:
  ParameterCache* pCache = nullptr;
  int nToys = -1;  ///< number of toys to be generated at each scan point

  // functions
  std::vector<double> fitHist(TH1* h, TString fitfunc = "p1+exp", bool draw = true);
  double transform(std::vector<double> fitParams, TString transFunc, double x);
  void printLatexLine(float eta, float finProb, float finProbErr, float finPlug, float finPlugErr);

  // result histograms
  TH1F* h_sol = nullptr;
  TH1F* h_pvalue_plugin = nullptr;
  TH1F* h_pvalue_prob = nullptr;
  TH1F* h_pvalue_plugin_notransf = nullptr;
  TH1F* h_pvalue_prob_notransf = nullptr;
  TTree* t_res = nullptr;

  // result values
  Long64_t nentries = 0ULL;
  Long64_t nfailed = 0ULL;
  double n68plugin = 0.;
  double n95plugin = 0.;
  double n99plugin = 0.;
  double n68prob = 0.;
  double n95prob = 0.;
  double n99prob = 0.;
};

#endif
