#ifndef FitResultDump_h
#define FitResultDump_h

#include <fstream>
#include <string>

#include "MethodAbsScan.h"

class FitResultDump {

 public:
  void dumpResult(std::string ofname, MethodAbsScan* scanner) const;
  std::ofstream outf;
};

#endif
