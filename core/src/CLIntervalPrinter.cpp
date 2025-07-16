#include <CLIntervalPrinter.h>
#include <Rounder.h>
#include <Utils.h>

#include <algorithm>
#include <format>
#include <string>

using namespace std;
using Utils::RadToDeg;

CLIntervalPrinter::CLIntervalPrinter(const OptParser* arg, TString name, TString var, TString unit, TString method,
                                     int CLsType) {
  assert(arg);
  _arg = arg;
  _name = name;
  _var = var;
  _unit = unit;
  _method = method;
  _clstype = CLsType;
}

///
/// Set the intervals. If more vectors of intervals are added, each of them will
/// be printed in order.
///
/// \param intervals - vector of confidence intervals, each one corresponding to one solution
///
void CLIntervalPrinter::addIntervals(const vector<CLInterval>& intervals) { _intervals.push_back(intervals); }

///
/// Helper function to sort the intervals according to their
/// lower boundary.
///

void CLIntervalPrinter::print() const {
  auto compareByMin = [](const CLInterval& a, const CLInterval& b) { return a.min < b.min; };
  for (auto intervals : _intervals) {
    vector<CLInterval> sortedIntervals(intervals);
    std::ranges::sort(sortedIntervals, compareByMin);
    for (auto sInterval : sortedIntervals) {
      CLInterval i = sInterval;

      std::string unit{_unit};
      if (_convertToDeg) {
        i.central = RadToDeg(i.central);
        i.min = RadToDeg(i.min);
        i.max = RadToDeg(i.max);
        unit = "Deg";
      }

      Rounder myRounder(_arg, i.min, i.max, i.central);
      int d = myRounder.getNsubdigits();
      printf("%s = [%7.*f, %7.*f] (%7.*f -%7.*f +%7.*f) @%3.2fCL", _var.Data(), d, myRounder.CLlo(), d,
             myRounder.CLhi(), d, myRounder.central(), d, myRounder.errNeg(), d, myRounder.errPos(), 1. - i.pvalue);
      if (!unit.empty()) cout << std::format(", [{}]", unit);
      // \todo remove the following code from quickhack stage once we have switched
      // to the CLIntervalMaker mechanism to get more useful information
      // on the CL intervals
      cout << ", " << _method;
      if (_clstype == 1) cout << " Simplified CL_s";
      if (_clstype == 2) cout << " Standard CL_s";

      if (_arg->isQuickhack(8)) {
        if (_arg->verbose) {
          cout << Form(", central: %-7s", i.centralmethod.Data());
          cout << Form(", interval: [%-6s, %-6s]", i.minmethod.Data(), i.maxmethod.Data());
          cout << ", p(central): " << i.pvalueAtCentral;
        }
      }
    }
    cout << endl;
  }
}

void CLIntervalPrinter::savePython() const {
  TString dirname = "plots/cl";
  TString ofname;
  if (_clstype == 0)
    ofname = dirname + "/clintervals_" + _name + "_" + _var + "_" + _method + ".py";
  else
    ofname = dirname + "/clintervals_" + _name + "_" + _var + "_" + _method + "_CLs" + std::to_string(_clstype) + ".py";
  if (_arg->verbose) cout << "CLIntervalPrinter::save() : saving " << ofname << endl;
  system("mkdir -p " + dirname);
  ofstream outf;
  outf.open(ofname);
  outf << "# Confidence Intervals" << endl;
  outf << "intervals = {" << endl;

  double previousCL = -1.;

  for (auto intervals : _intervals)
    for (auto interval : intervals) {
      CLInterval i = interval;

      std::string unit{_unit};
      if (_convertToDeg) {
        i.central = RadToDeg(i.central);
        i.min = RadToDeg(i.min);
        i.max = RadToDeg(i.max);
        unit = "Deg";
      }

      Rounder myRounder(_arg, i.min, i.max, i.central);
      int d = myRounder.getNsubdigits();

      double thisCL = 1. - i.pvalue;
      if (previousCL != thisCL) {
        if (previousCL != -1) outf << "  ]," << endl;
        outf << Form("  '%.2f' : [", thisCL) << endl;
      }

      outf << Form("    {'var':'%s', 'min':'%.*f', 'max':'%.*f', 'central':'%.*f', "
                   "'neg':'%.*f', 'pos':'%.*f', 'cl':'%.2f', 'unit':'%s', 'method':'%s'},\n",
                   _var.Data(), d, myRounder.CLlo(), d, myRounder.CLhi(), d, myRounder.central(), d, myRounder.errNeg(),
                   d, myRounder.errPos(), thisCL, unit, _method.Data());
      previousCL = thisCL;
    }
  if (previousCL != -1) { outf << "  ]" << endl; }
  outf << "}" << endl;
  outf.close();
}
