#ifndef Rounder_h
#define Rounder_h

#include "OptParser.h"

#include <limits>

/**
 * Class to round confidence intervals and their central values.
 */
class Rounder {
 public:
  Rounder(const OptParser* arg, const double clMin, const double clMax, const double central);

  int getNsubdigits() const { return m_nSubDigits; }
  double CLlo() const;
  double CLhi() const;
  double central() const;
  double errNeg() const;
  double errPos() const;

 private:
  const OptParser* arg;  ///< command line arguments
  double m_clMin = std::numeric_limits<double>::quiet_NaN();
  double m_clMax = std::numeric_limits<double>::quiet_NaN();
  double m_central = std::numeric_limits<double>::quiet_NaN();
  int m_nSubDigits = -1;
};

#endif
