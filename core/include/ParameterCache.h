#ifndef ParameterCache_h
#define ParameterCache_h

#include <fstream>
#include <vector>

#include "Combiner.h"
#include "MethodAbsScan.h"
#include "RooSlimFitResult.h"
#include "Utils.h"

#include <TString.h>

class ParameterCache {

 public:
  ParameterCache(OptParser* arg);
  ~ParameterCache();

  void cacheParameters(MethodAbsScan* scanner, TString fileName);
  bool loadPoints(TString fileName);
  void printFitResultToOutStream(std::ofstream& out, RooSlimFitResult* slimFitResult) const;
  void printPoint() const;
  int getNPoints() const;
  void setPoint(Combiner* cmb, int i);
  void setPoint(MethodAbsScan* scanner, int i);
  std::vector<TString> getFixedNames(std::vector<Utils::FixPar> fixPar) const;
  std::vector<std::map<TString, double>> startingValues;

 private:
  bool m_parametersLoaded = false;
  OptParser* m_arg;
};

#endif
