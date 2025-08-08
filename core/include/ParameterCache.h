#ifndef ParameterCache_h
#define ParameterCache_h

#include "Utils.h"

#include <TString.h>

#include <fstream>
#include <map>
#include <vector>

class Combiner;
class MethodAbsScan;
class OptParser;
class RooSlimFitResult;

class ParameterCache {

 public:
  ParameterCache(OptParser* arg);

  void cacheParameters(MethodAbsScan* scanner, TString fileName);
  bool loadPoints(TString fileName);
  void printFitResultToOutStream(std::ofstream& out, RooSlimFitResult* slimFitResult);
  void printPoint();
  int getNPoints();
  void setPoint(Combiner* cmb, int i);
  void setPoint(MethodAbsScan* scanner, int i);
  std::vector<TString> getFixedNames(std::vector<Utils::FixPar> fixPar);
  std::vector<std::map<TString, double>> startingValues;

 private:
  bool m_parametersLoaded = false;
  OptParser* m_arg = nullptr;
};

#endif
