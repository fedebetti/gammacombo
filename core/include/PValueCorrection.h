#ifndef PValueCorrection_h
#define PValueCorrection_h

#include <memory>
#include <vector>

#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TString.h>

class PValueCorrection {

 public:
  PValueCorrection(TString _transFunc = "none", bool _verbose = false);
  PValueCorrection(int id, bool _verbose = false);

  void setTransFunc(TString tf) { transFunc = tf; }
  void setFitParams(std::vector<double> fP) { fitParams = fP; }
  void setFitParam(int i, double val);
  void readFiles(TString name, int id = 0, bool isPlugin = false);
  void fitHist(TH1* h);
  double transform(double x);

  void checkValid();
  void checkParams();
  void printCoverage(double, double, double, double, TString name = "");

  void write(TString fname);
  void write(TFile* const f);

 private:
  TString transFunc;
  bool verbose;
  TString fitString;
  TF1 fitFunc;
  std::vector<double> fitParams;
  std::vector<TString> allowedFuncs;
  std::unique_ptr<TH1> h_pvalue_before = nullptr;
  std::unique_ptr<TH1> h_pvalue_after = nullptr;
};

#endif
