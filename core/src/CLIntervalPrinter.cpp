#include <CLIntervalPrinter.h>
#include <Rounder.h>
#include <Utils.h>

#include <format>
#include <memory>
#include <string>
#include <vector>

CLIntervalPrinter::CLIntervalPrinter(const OptParser* arg, const TString name, const TString var, const TString unit,
                                     const TString method, const int CLsType)
    : _arg(arg), _name(name), _var(var), _unit(unit), _method(method), _clstype(CLsType) {
  assert(_arg);
}

void CLIntervalPrinter::addIntervals(const std::vector<CLInterval>& intervals) {
  for (const auto& i : intervals) _intervals.emplace(i);
}
///
/// Set the intervals. If more vectors of intervals are added, each of them will
/// be printed in order.
///
/// \param intervals - vector of confidence intervals, each one corresponding to one solution
///
void CLIntervalPrinter::addIntervals(const std::vector<std::vector<std::unique_ptr<CLInterval>>>& intervals) {
  for (const auto& ints : intervals) {
    for (const auto& i : ints) {
      if (!i) continue;
      _intervals.emplace(CLInterval(*i));
    }
  }
}

/**
 * Print all valid intervals.
 */
void CLIntervalPrinter::print() const {
  for (auto interval : _intervals) {
    CLInterval i = interval;

    std::string unit{_unit};
    if (_convertToDeg) {
      using Utils::RadToDeg;
      i.central = RadToDeg(i.central);
      i.min = RadToDeg(i.min);
      i.max = RadToDeg(i.max);
      unit = "Deg";
    }

    Rounder rounder(_arg, i.min, i.max, i.central);
    std::string unit_str = unit.empty() ? "" : std::format(", [{}]", unit);
    const int d = rounder.getNsubdigits();
    const auto float_format = std::format("{{:-7.{:d}f}}", d);
    const auto format_str =
        std::format("{{:s}} = [{0}, {0}] ({0} -{0} + {0}) @{{:4.3f}}CL{{:s}}, {{:s}}", float_format);
    std::cout << std::vformat(format_str, std::make_format_args(_var, rounder.CLlo(), rounder.CLhi(), rounder.central(),
                                                                rounder.errNeg(), rounder.errPos(), 1. - i.pvalue,
                                                                unit_str, _method));

    // TODO remove the following code from quickhack stage once we have switched to the CLIntervalMaker mechanism
    // to get more useful information on the CL intervals
    if (_clstype == 1) std::cout << " Simplified CL_s";
    if (_clstype == 2) std::cout << " Standard CL_s";

    if (_arg->isQuickhack(8) && _arg->verbose) {
      std::cout << std::format(", central: {:>7s}, interval: [{:>6s}, {:>6s}], p(central): {:4.3f}",
                               std::string(i.centralmethod), std::string(i.minmethod), std::string(i.maxmethod),
                               i.pvalueAtCentral);
    }
    std::cout << std::endl;
  }
}

void CLIntervalPrinter::savePython() const {
  std::string dirname = "plots/cl";
  std::string clstype_str = (_clstype == 0) ? "" : "_CLs" + std::to_string(_clstype);
  std::string ofname = dirname + "/clintervals_" + _name + "_" + _var + "_" + _method + clstype_str + ".py";
  if (_arg->verbose) std::cout << "CLIntervalPrinter::save() : saving " << ofname << std::endl;
  std::string cmd = "mkdir -p " + dirname;
  system(cmd.c_str());
  std::ofstream outf;
  outf.open(ofname);
  outf << "# Confidence Intervals\n"
       << "intervals = {" << std::endl;

  // TODO
  double previousCL = -1.;

  for (const auto interval : _intervals) {
    CLInterval i = interval;

    std::string unit{_unit};
    if (_convertToDeg) {
      using Utils::RadToDeg;
      i.central = RadToDeg(i.central);
      i.min = RadToDeg(i.min);
      i.max = RadToDeg(i.max);
      unit = "Deg";
    }

    Rounder myRounder(_arg, i.min, i.max, i.central);
    int d = myRounder.getNsubdigits();

    double thisCL = 1. - i.pvalue;
    if (previousCL != thisCL) {
      if (previousCL != -1) outf << "  ]," << std::endl;
      outf << Form("  '%.2f' : [", thisCL) << std::endl;
    }

    outf << Form("    {'var':'%s', 'min':'%.*f', 'max':'%.*f', 'central':'%.*f', "
                 "'neg':'%.*f', 'pos':'%.*f', 'cl':'%.2f', 'unit':'%s', 'method':'%s'},\n",
                 _var, d, myRounder.CLlo(), d, myRounder.CLhi(), d, myRounder.central(), d, myRounder.errNeg(), d,
                 myRounder.errPos(), thisCL, unit, _method);
    previousCL = thisCL;
  }
  if (previousCL != -1) { outf << "  ]" << std::endl; }
  outf << "}" << std::endl;
  outf.close();
}
