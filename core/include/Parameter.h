/**
 * Author: Till Moritz Karbach, moritz.karbach@cern.ch
 * Date: April 2013
 *
 * Class representing a (nuisance) parameter.
 *
 **/

#ifndef Parameter_h
#define Parameter_h

#include <TString.h>

class Parameter {
 public:
  inline void setVal(double v) { startvalue = v; };
  void Print() const;

  struct Range {
    float min;
    float max;
  };

 private:
  const static Range default_range;

 public:
  TString name = "not initialized";
  TString title = "not initialized";
  TString unit = "not initialized";
  double startvalue = 0.;
  Range phys = default_range;
  Range scan = default_range;
  Range force = default_range;
  Range bboos = default_range;
  [[deprecated]] Range free = default_range;
};

#endif
