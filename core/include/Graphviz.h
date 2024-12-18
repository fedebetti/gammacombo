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

#include "Combiner.h"
#include "OptParser.h"

#include <fstream>

class Graphviz {
 public:
  Graphviz(const OptParser* arg);
  ~Graphviz();

  void printCombiner(Combiner* cmb) const;
  void printCombinerLayer(Combiner* cmb) const;

 private:
  TString graphvizString(TString s) const;
  // TString   graphvizString(string s);
  bool isDmixingParameter(TString s) const;
  std::ofstream& openFile(TString name) const;
  const OptParser* arg;  ///< command line arguments
};

#endif
