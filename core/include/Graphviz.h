/**
 * Gamma Combination
 * Author: Till Moritz Karbach, moritz.karbach@cern.ch
 * Date: March 2014
 *
 * Class to print structures into *.dot files that can be
 * proceesed by graphviz. For example the Measurement/parameter
 * structure of the combiners.
 *
 **/

#ifndef Graphviz_h
#define Graphviz_h

#include <TString.h>

class Combiner;
class OptParser;

class Graphviz {
 public:
  Graphviz(OptParser* arg);

  void printCombiner(Combiner* cmb);
  void printCombinerLayer(Combiner* cmb);

 private:
  OptParser* arg = nullptr;  ///< command line arguments
};

#endif
