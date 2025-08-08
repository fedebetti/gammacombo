/**
 * Gamma Combination
 * Author: Till Moritz Karbach, moritz.karbach@cern.ch
 * Date: August 2014
 *
 **/

#ifndef CLIntervalMaker_h
#define CLIntervalMaker_h

#include "CLInterval.h"

#include <vector>

class OptParser;

class TH1F;

///
/// Class that makes CL intervals from 1-CL histograms.
///
class CLIntervalMaker {
 public:
  CLIntervalMaker(OptParser* arg, const TH1F& pvalues);
  void calcCLintervals();
  void findMaxima(double pValueThreshold);
  inline std::vector<CLInterval>& getClintervals1sigma() { return _clintervals1sigma; };
  inline std::vector<CLInterval>& getClintervals2sigma() { return _clintervals2sigma; };
  void print();
  void provideMorePreciseMaximum(double value, TString method);

 private:
  int checkNeighboringBins(int i, double y) const;
  bool binsOnSameSide(int i, double y) const;
  double binToValue(int bin) const;
  void findRawIntervals(double pvalue, std::vector<CLInterval>& clis);
  void findRawIntervalsForCentralValues(double pvalue, std::vector<CLInterval>& clis);
  bool interpolateLine(const TH1F* h, int i, double y, double& val) const;
  bool interpolatePol2fit(const TH1F* h, int i, double y, double central, bool upper, double& val, double& err) const;
  bool isInInterval(int binid, double pvalue) const;
  void improveIntervalsLine(std::vector<CLInterval>& clis) const;
  void improveIntervalsPol2fit(std::vector<CLInterval>& clis) const;
  double pq(double p0, double p1, double p2, double y, int whichSol) const;
  void removeBadIntervals();
  bool similarMaximumExists(double value) const;
  void storeRawInterval(int binidLo, int binidHi, double pvalue, std::vector<CLInterval>& clis);
  int valueToBin(double val) const;

  OptParser* _arg = nullptr;                   ///< command line arguments
  const TH1F& _pvalues;                        ///< the pvalue histogram
  std::vector<CLInterval> _clintervals1sigma;  ///< 1 sigma intervals
  std::vector<CLInterval> _clintervals2sigma;  ///< 2 sigma intervals
};

#endif
