/**
 * Gamma Combination
 * Author: Matthew Kenzie matthew.kenzie@cern.ch
 * Date: Apr 2015
 *
 **/

#ifndef BatchScriptWriter_h
#define BatchScriptWriter_h

#include "Combiner.h"
#include "OptParser.h"
#include "PDF_Abs.h"

class BatchScriptWriter {
 public:
  BatchScriptWriter(int argc, char* argv[]);

  void writeScripts(const OptParser* arg, std::vector<Combiner*>* cmb);
  void writeScripts_datasets(const OptParser* arg, PDF_Abs* pdf);
  void writeScript(TString fname, TString outfloc, int jobn, const OptParser* arg);
  void writeCondorScript(TString fname, const OptParser* arg);
  std::string exec;
  std::string subpkg;
};

#endif
