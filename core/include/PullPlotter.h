/**
 * Gamma Combination
 * Author: Till Moritz Karbach, moritz.karbach@cern.ch
 * Date: August 2013
 *
 * Class to plot pulls.
 *
 **/

#ifndef PullPlotter_h
#define PullPlotter_h

#include <TString.h>

#include <vector>

class MethodAbsScan;
class OptParser;

class PullPlotter {
 public:
  PullPlotter(MethodAbsScan* cmb);

  bool hasPullsAboveNsigma(float nsigma);
  void loadParsFromSolution(int n);
  void savePulls();
  void plotPulls();
  void printPulls(float aboveNsigma = -1.);

 private:
  void defineOrder();
  void plotPullsCanvas(std::vector<TString>& observables, int currentid, int maxid, int nObs);

  MethodAbsScan* cmb = nullptr;   // the scanner to plot pulls for
  OptParser* arg = nullptr;       // command line arguments
  std::vector<TString> obsOrder;  // contains observable names in the desired plot order
  int nSolution = 0;              // index of the solution wrt which the pulls are computed
};

#endif
