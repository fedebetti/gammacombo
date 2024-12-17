/*
 * Gamma Combination
 * Author: Maximilian Schlupp, maxschlupp@gmail.com
 * Author: Konstantin Schubert, schubert.konstantin@gmail.com
 * Date: October 2016
 */

#ifndef MethodDatasetsProbScan_h
#define MethodDatasetsProbScan_h

#include "MethodProbScan.h"
#include "PDF_Datasets.h"
#include "RooSlimFitResult.h"
#include "ToyTree.h"

class MethodDatasetsProbScan : public MethodProbScan {
 public:
  MethodDatasetsProbScan(PDF_Datasets* PDF, OptParser* opt);

  virtual void initScan();
  void loadScanFromFile(TString fileNameBaseIn = "default");
  void loadFitResults(TString file);
  void loadParameterLimits();
  virtual void print();
  virtual int scan1d(bool fast = false, bool reverse = false, bool quiet = false);
  virtual int scan2d();
  virtual bool loadScanner(TString fName);
  inline void setInputFile(TString name) {
    inputFiles.push_back(name);
    explicitInputFile = true;
  };
  inline void addFile(TString name) { inputFiles.push_back(name); };
  void plotFitRes(TString fName);
  int computeCLvalues() const;

  PDF_Datasets* pdf = nullptr;
  TH1F* probPValues = nullptr;
  bool drawPlots = false;
  bool explicitInputFile = false;
  std::vector<TString> inputFiles;
  std::vector<double> bootstrapPVals;
  TChain* chain = nullptr;
  RooFitResult* bkgOnlyFitResult = nullptr;
  ToyTree* probScanTree = nullptr;

 protected:
 private:
  TChain* readFiles(TString fileNameBaseIn = "default");
  void readScan1dTrees(TString fileNameBaseIn = "default");
  RooFitResult* loadAndFit(PDF_Datasets* pdf);
  double getPValueTTestStatistic(double test_statistic_value, bool isCLs = false) const;
  void sanityChecks() const;
  void setAndPrintFitStatusConstrainedToys(const ToyTree& toyTree);
  void setAndPrintFitStatusFreeToys(const ToyTree& toyTree);
  void sethCLFromProbScanTree();
};

#endif
