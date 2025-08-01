#ifndef MethodProb_h
#define MethodProb_h

#include <TH1.h>

#include "MethodAbsScan.h"

class MethodProbScan : public MethodAbsScan {
 public:
  MethodProbScan(Combiner* comb);
  MethodProbScan(OptParser* opt);
  MethodProbScan();

  virtual int computeCLvalues();  // compute CL histograms depending on desired test statistic
  double getChi2min(double scanpoint) const;
  inline TH1* getHChi2min() const { return hChi2min.get(); };
  void saveSolutions();
  void saveSolutions2d();
  int scan1d() { return scan1d(false, false, false); }
  virtual int scan1d(bool fast, bool const reverse, const bool quiet);
  virtual int scan2d();
  inline void setScanDisableDragMode(bool f = true) { scanDisableDragMode = f; };

 protected:
  bool computeInnerTurnCoords(const int iStart, const int jStart, const int i, const int j, int& iResult, int& jResult,
                              int nTurn) const;
  bool deleteIfNotInCurveResults2d(const RooSlimFitResult* r);
  void sanityChecks() const;

  int nScansDone = 0;  // count the number of times a scan was done

 private:
  bool scanDisableDragMode = false;
};

#endif
