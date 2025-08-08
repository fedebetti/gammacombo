/**
 * Gamma Combination
 * Author: Matthew Kenzie matthew.kenzie@cern.ch
 * Date: Apr 2015
 *
 **/

#ifndef BatchScriptWriter_h
#define BatchScriptWriter_h

#include <TString.h>

#include <string>
#include <vector>

class Combiner;
class OptParser;
class PDF_Abs;

class BatchScriptWriter {
 public:
  BatchScriptWriter(int argc, char* argv[]);

  void writeScripts(OptParser* arg, std::vector<Combiner*>* cmb);
  void writeScripts_datasets(OptParser* arg, PDF_Abs* pdf);
  void writeScript(TString fname, TString outfloc, int jobn, OptParser* arg);
  void writeCondorScript(TString fname, OptParser* arg);
  std::string exec;
  std::string subpkg;
};

#endif
