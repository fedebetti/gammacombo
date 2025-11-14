/**
 * Author: Till Moritz Karbach, moritz.karbach@cern.ch
 * Date: April 2013
 *
 * Abstract class to define the (nuisance) parameters.
 * Project: Lb2pktaul
 *
 **/

#ifndef ParametersLb2pktaul_h
#define ParametersLb2pktaul_h

#include <ParametersAbs.h>

class ParametersLb2pktaul : public ParametersAbs {
 public:
  ParametersLb2pktaul();

 protected:
  void defineParameters();
};

#endif
