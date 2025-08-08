#ifndef FitResultDump_h
#define FitResultDump_h

#include <fstream>
#include <string>

class MethodAbsScan;

// TODO this class should rather be a utility function
class FitResultDump {

 public:
  void dumpResult(std::string ofname, MethodAbsScan* scanner);
  std::ofstream outf;
};

#endif
