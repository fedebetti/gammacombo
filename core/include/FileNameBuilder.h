/**
 * Gamma Combination
 * Author: Till Moritz Karbach, moritz.karbach@cern.ch
 * Date: August 2014
 *
 **/

#ifndef FileNameBuilder_h
#define FileNameBuilder_h

#include <vector>

#include <TString.h>

#include "OptParser.h"

// forward declarations - include headers only in
// cpp file to deal with circular dependencies
class MethodAbsScan;
class GammaComboEngine;
class Combiner;

///
/// Class that defines the file names of various files being
/// written across the project.
///
class FileNameBuilder {
 public:
  FileNameBuilder(const OptParser* arg, TString name = "gammacombo");
  TString getBaseName();
  TString getFileBaseName(const Combiner* c);
  TString getFileBaseName(const MethodAbsScan* s);
  TString getFileNameScanner(const MethodAbsScan* s);
  TString getFileNameSolution(const MethodAbsScan* s);
  TString getFileNamePar(const Combiner* c);
  TString getFileNamePar(const MethodAbsScan* s);
  TString getFileNamePlot(const std::vector<Combiner*>& cmb);
  TString getFileNamePlotSingle(const std::vector<Combiner*>& cmb, int cId);
  TString getFileNameStartPar(const Combiner* c);
  TString getFileNameStartPar(const MethodAbsScan* s);
  TString getFileNameAsimovPar(const Combiner* c);
  TString getFileNameAsimovPar(const MethodAbsScan* s);
  TString getAsimovCombinerNameAddition(int id);
  TString getPluginNameAddition();
  TString getPluginOnlyNameAddition();
  TString getPreliminaryNameAddition();
  TString getCLsNameAddition();

 private:
  const OptParser* m_arg;       ///< command line arguments
  TString m_basename;           ///< the base name, e.g. "gammacombo"
  TString m_asimov = "Asimov";  ///< literal naming Asimov combiners
};

#endif
