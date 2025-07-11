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

#include "OptParser.h"

class Rounder {
 public:
  Rounder(const OptParser* arg, double cllo, double clhi, double central);

  int getNsubdigits() const;
  double CLlo() const;
  double CLhi() const;
  double central() const;
  double errNeg() const;
  double errPos() const;

 private:
  const OptParser* arg;  ///< command line arguments
  double m_cllo;
  double m_clhi;
  double m_central;
};

#endif
