#ifndef Graphviz_h
#define Graphviz_h

#include "Combiner.h"
#include "OptParser.h"

/**
 * Class to print structures into *.dot files that can be proceesed by graphviz.
 *
 * For example the measurement/parameter structure of the combiners.
 */
class Graphviz {
 public:
  Graphviz(const OptParser* arg);

  void printCombiner(Combiner* cmb) const;
  void printCombinerLayer(Combiner* cmb) const;

 private:
  const OptParser* arg;  ///< command line arguments
};

#endif
