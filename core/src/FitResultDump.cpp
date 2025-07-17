#include <FitResultDump.h>
#include <Utils.h>

#include <iostream>

using namespace std;
using Utils::isAngle;
using Utils::RadToDeg;

void FitResultDump::dumpResult(string ofname, MethodAbsScan* scanner) const {

  ofname = "plots/par/" + ofname + ".dat";
  system("mkdir -p plots/par");

  cout << "FitResultDump::dumpResult() : saving " << ofname << endl;

  ofstream outf;
  outf.open(ofname.c_str());
  outf << "# Fit Result Summary" << endl;
  outf << "nSolutions=" << scanner->getSolutions().size() << endl;
  outf << "# pvalue central min max" << endl;

  bool angle = false;
  if (isAngle(scanner->getWorkspace()->var(scanner->getScanVar1Name()))) angle = true;

  for (int cl = 0; cl < 2; ++cl) {
    for (const auto& ci : scanner->clintervals[cl]) {
      auto central = ci->central;
      auto min = ci->min;
      auto max = ci->max;
      if (angle) {
        central = RadToDeg(central);
        min = RadToDeg(min);
        max = RadToDeg(max);
      }
      outf << ci->pvalue << " " << central << " " << min << " " << max << endl;
    }
  }

  outf.close();
}
