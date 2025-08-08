#ifndef PValueCorrection_h
#define PValueCorrection_h

#include <TF1.h>
#include <TString.h>

#include <vector>

class TFile;
class TH1;
class TH1F;

class PValueCorrection {

 public:
  PValueCorrection(TString _transFunc = "none", bool _verbose = false);
  PValueCorrection(int id, bool _verbose = false);
  ~PValueCorrection();

  PValueCorrection(PValueCorrection&) = delete;
  PValueCorrection& operator=(const PValueCorrection&) = delete;

  void setTransFunc(TString tf) { transFunc = tf; }
  void setFitParams(std::vector<double> fP) { fitParams = fP; }
  void setFitParam(int i, double val);
  void readFiles(TString name, int id = 0, bool isPlugin = false);
  void fitHist(TH1* h);
  double transform(double x);

  void checkValid() const;
  void checkParams() const;
  void printCoverage(float, float, float, float, TString name = "") const;

  void write(TString fname);
  void write(TFile* f);

 private:
  TString transFunc;
  bool verbose = false;
  TString fitString;
  TF1 fitFunc;
  std::vector<double> fitParams;
  std::vector<TString> allowedFuncs;
  TH1F* h_pvalue_before = nullptr;
  TH1F* h_pvalue_after = nullptr;
};

#endif
