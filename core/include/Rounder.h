/**
 * Gamma Combination
 * Author: Till Moritz Karbach, moritz.karbach@cern.ch
 * Date: April 2013
 *
 * Class to round confidence intervals and central values.
 *
 **/

#ifndef Rounder_h
#define Rounder_h

#include <limits>

class OptParser;

class Rounder {
 public:
  Rounder(OptParser* arg, float cllo, float clhi, float central);

  int getNsubdigits() const;
  float CLlo();
  float CLhi();
  float central();
  float errNeg();
  float errPos();

 private:
  OptParser* arg = nullptr;
  float m_cllo = std::numeric_limits<float>::quiet_NaN();
  float m_clhi = std::numeric_limits<float>::quiet_NaN();
  float m_central = std::numeric_limits<float>::quiet_NaN();
};

#endif
