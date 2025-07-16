#ifndef CLIntervalPrinter_h
#define CLIntervalPrinter_h

#include "CLInterval.h"
#include "OptParser.h"

/**
 * Class that prints CL intervals and saves them to disk.
 */
class CLIntervalPrinter {
 public:
  CLIntervalPrinter(const OptParser* arg, TString name, TString var, TString unit, TString method, int CLsType = 0);

  void print() const;
  void savePython() const;
  inline void setDegrees(bool yesno = true) { _convertToDeg = yesno; };
  void addIntervals(const std::vector<CLInterval>& intervals);

 private:
  const OptParser* _arg;                            ///< command line arguments
  TString _name;                                    ///< name of combination
  TString _var;                                     ///< name of scan variable
  TString _unit;                                    ///< unit of scan variable
  TString _method;                                  ///< method name (e.g. Prob)
  bool _convertToDeg = false;                       ///< convert values into degrees
  std::vector<std::vector<CLInterval>> _intervals;  ///< container of intervals
  int _clstype;                                     ///< Type of CLs intervals, 0 means no CLs method
};

#endif
