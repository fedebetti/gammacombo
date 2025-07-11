/**
 * Gamma Combination
 * Author: Till Moritz Karbach, moritz.karbach@cern.ch
 * Date: August 2014
 *
 **/

#ifndef CLInterval_h
#define CLInterval_h

#include <TString.h>

///
/// Class that represents a confidence interval
///
class CLInterval {
 public:
  void print() const;

  double pvalue = -1.;            // pvalue corresponding to this interval
  double pvalueAtCentral = -1.;   // pvalue at the central value
  double min = -1.;               // lower interval border
  double max = -1.;               // upper interval border
  double central = -1.;           // central value
  bool minclosed = false;         // true if the interval was not closed by limited scan range
  bool maxclosed = false;         // true if the interval was not closed by limited scan range
  TString minmethod = "n/a";      // details on the algorithm that found this interval
  TString maxmethod = "n/a";      // details on the algorithm that found this interval
  TString centralmethod = "n/a";  // details on the algorithm that found the central value
};

#endif
