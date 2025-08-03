#ifndef FitResultDump_h
#define FitResultDump_h

#include "MethodAbsScan.h"
#include <fstream>
#include <iostream>

class FitResultDump {

 public:
  FitResultDump();
  ~FitResultDump();

  void dumpResult(std::string ofname, MethodAbsScan* scanner);
  ofstream outf;
};

#endif
