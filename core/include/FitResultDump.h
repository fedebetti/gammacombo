#ifndef FitResultDump_h
#define FitResultDump_h

#include <fstream>
#include <string>

class MethodAbsScan;

class FitResultDump {

 public:
  FitResultDump();
  ~FitResultDump();

  void dumpResult(std::string ofname, MethodAbsScan* scanner);
  std::ofstream outf;
};

#endif
