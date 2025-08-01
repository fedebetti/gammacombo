#ifndef CLIntervalPrinter_h
#define CLIntervalPrinter_h

#include "CLInterval.h"
#include "OptParser.h"

#include <memory>
#include <set>
#include <string>
#include <vector>

#include <TString.h>

/**
 * Class that prints CL intervals and saves them to disk.
 */
class CLIntervalPrinter {
 public:
  CLIntervalPrinter(const OptParser* arg, const TString name, const TString var, const TString unit,
                    const TString method, const int CLsType = 0);

  void addIntervals(const std::vector<CLInterval>& intervals);
  void addIntervals(const std::vector<std::vector<std::unique_ptr<CLInterval>>>& intervals);
  void print() const;
  void savePython() const;
  inline void setDegrees(const bool yesno = true) { _convertToDeg = yesno; };

 private:
  const OptParser* _arg;            ///< Command line arguments
  std::string _name;                ///< Name of combination
  std::string _var;                 ///< Name of scan variable
  std::string _unit;                ///< Unit of scan variable
  std::string _method;              ///< Method name (e.g. Prob)
  bool _convertToDeg = false;       ///< Convert values into degrees
  std::set<CLInterval> _intervals;  ///< Container of intervals, sorted according to CLCompare.
  int _clstype = 0;                 ///< Type of CLs intervals, 0 means no CLs method
};

#endif
