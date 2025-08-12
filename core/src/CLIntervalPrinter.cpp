#include <CLIntervalPrinter.h>

#include <OptParser.h>
#include <Rounder.h>
#include <Utils.h>

#include <TString.h>

#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>

CLIntervalPrinter::CLIntervalPrinter(const OptParser* arg, TString name, TString var, TString unit, TString method,
                                     int CLsType)
    : _name(name), _var(var), _unit(unit), _method(method), _clstype(CLsType) {
  assert(arg);
  _arg = arg;
}

///
/// Set the intervals. If more vectors of intervals are added, each of them will
/// be printed in order.
///
/// \param intervals - vector of confidence intervals, each one corresponding to one solution
///
void CLIntervalPrinter::addIntervals(std::vector<CLInterval>& intervals) { _intervals.push_back(intervals); }

///
/// Helper function to sort the intervals according to their
/// lower boundary.
///
bool CLIntervalPrinter::compareByMin(const CLInterval& a, const CLInterval& b) { return a.min < b.min; }

void CLIntervalPrinter::print() {
  for (int k = 0; k < _intervals.size(); k++) {
    // sort the intervals
    std::vector<CLInterval> sortedIntervals(_intervals[k]);
    std::sort(sortedIntervals.begin(), sortedIntervals.end(), compareByMin);
    for (int j = 0; j < sortedIntervals.size(); j++) {
      CLInterval i = sortedIntervals[j];

      // convert to degrees if necessary
      if (convertToDegrees) {
        using Utils::RadToDeg;
        i.central = RadToDeg(i.central);
        i.min = RadToDeg(i.min);
        i.max = RadToDeg(i.max);
        _unit = "Deg";
      }

      Rounder myRounder(_arg, i.min, i.max, i.central);
      int d = myRounder.getNsubdigits();
      printf("%s = [%7.*f, %7.*f] (%7.*f -%7.*f +%7.*f) @%3.2fCL", _var.Data(), d, myRounder.CLlo(), d,
             myRounder.CLhi(), d, myRounder.central(), d, myRounder.errNeg(), d, myRounder.errPos(), 1. - i.pvalue);
      if (_unit != "") std::cout << ", [" << _unit << "]";
      // \todo remove the following code from quickhack stage once we have switched
      // to the CLIntervalMaker mechanism to get more useful information
      // on the CL intervals
      std::cout << ", " << _method;
      if (_clstype == 1) std::cout << " Simplified CL_s";
      if (_clstype == 2) std::cout << " Standard CL_s";

      if (_arg->isQuickhack(8)) {
        if (_arg->verbose) {
          std::cout << Form(", central: %-7s", i.centralmethod.Data());
          std::cout << Form(", interval: [%-6s, %-6s]", i.minmethod.Data(), i.maxmethod.Data());
          std::cout << ", p(central): " << i.pvalueAtCentral;
        }
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
}

void CLIntervalPrinter::savePython() {
  TString dirname = "plots/cl";
  TString ofname;
  if (_clstype == 0)
    ofname = dirname + "/clintervals_" + _name + "_" + _var + "_" + _method + ".py";
  else
    ofname = dirname + "/clintervals_" + _name + "_" + _var + "_" + _method + "_CLs" + std::to_string(_clstype) + ".py";
  if (_arg->verbose) std::cout << "CLIntervalPrinter::save() : saving " << ofname << std::endl;
  system("mkdir -p " + dirname);
  std::ofstream outf;
  outf.open(ofname);
  outf << "# Confidence Intervals" << std::endl;
  outf << "intervals = {" << std::endl;

  float previousCL = -1.0;

  for (int k = 0; k < _intervals.size(); k++)
    for (int j = 0; j < _intervals[k].size(); j++) {
      CLInterval i = _intervals[k][j];

      // convert to degrees if necessary
      if (convertToDegrees) {
        using Utils::RadToDeg;
        i.central = RadToDeg(i.central);
        i.min = RadToDeg(i.min);
        i.max = RadToDeg(i.max);
        _unit = "Deg";
      }

      Rounder myRounder(_arg, i.min, i.max, i.central);
      int d = myRounder.getNsubdigits();

      float thisCL = 1. - i.pvalue;
      if (previousCL != thisCL) {
        if (previousCL != -1) outf << "  ]," << std::endl;
        outf << Form("  '%.2f' : [", thisCL) << std::endl;
      }

      outf << Form("    {'var':'%s', 'min':'%.*f', 'max':'%.*f', 'central':'%.*f', "
                   "'neg':'%.*f', 'pos':'%.*f', 'cl':'%.2f', 'unit':'%s', 'method':'%s'},\n",
                   _var.Data(), d, myRounder.CLlo(), d, myRounder.CLhi(), d, myRounder.central(), d, myRounder.errNeg(),
                   d, myRounder.errPos(), thisCL, _unit.Data(), _method.Data());
      previousCL = thisCL;
    }
  if (previousCL != -1) { outf << "  ]" << std::endl; }
  outf << "}" << std::endl;
  outf.close();
}
