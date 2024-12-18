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

#include "MethodAbsScan.h"
#include "OptParser.h"

class MethodAbsScan;

class PullPlotter {
 public:
  PullPlotter(MethodAbsScan* cmb);

  bool hasPullsAboveNsigma(float nsigma) const;
  void loadParsFromSolution(int n);
  void savePulls() const;
  void plotPulls();
  void printPulls(float aboveNsigma = -1.) const;

 private:
  void defineOrder();
  void plotPullsCanvas(std::vector<TString>& observables, int currentid, int maxid, int nObs) const;

  MethodAbsScan* cmb;             // the scanner to plot pulls for
  const OptParser* arg;           // command line arguments
  std::vector<TString> obsOrder;  // contains observable names in the desired plot order
  int nSolution = 0;              // index of the solution wrt which the pulls are computed
};

#endif
