#include <FitResultDump.h>

#include <MethodAbsScan.h>
#include <Utils.h>

#include <RooWorkspace.h>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

void FitResultDump::dumpResult(std::string ofname, MethodAbsScan* scanner) {

  ofname = "plots/par/" + ofname + ".dat";
  system("mkdir -p plots/par");

  std::cout << "FitResultDump::dumpResult() : saving " << ofname << std::endl;

  std::ofstream outf;
  outf.open(ofname.c_str());
  outf << "# Fit Result Summary" << std::endl;
  outf << "nSolutions=" << scanner->getSolutions().size() << std::endl;
  outf << "# pvalue central min max" << std::endl;

  bool angle = false;
  if (Utils::isAngle(scanner->getWorkspace()->var(scanner->getScanVar1Name()))) angle = true;

  for (std::vector<CLInterval>::iterator cl = scanner->clintervals1sigma.begin();
       cl != scanner->clintervals1sigma.end(); cl++) {

    auto central = cl->central;
    auto min = cl->min;
    auto max = cl->max;
    if (angle) {
      using Utils::RadToDeg;
      central = RadToDeg(central);
      min = RadToDeg(min);
      max = RadToDeg(max);
    }

    outf << cl->pvalue << " " << central << " " << min << " " << max << std::endl;
  }
  for (std::vector<CLInterval>::iterator cl = scanner->clintervals2sigma.begin();
       cl != scanner->clintervals2sigma.end(); cl++) {

    auto central = cl->central;
    auto min = cl->min;
    auto max = cl->max;
    if (angle) {
      using Utils::RadToDeg;
      central = RadToDeg(central);
      min = RadToDeg(min);
      max = RadToDeg(max);
    }

    outf << cl->pvalue << " " << central << " " << min << " " << max << std::endl;
  }

  outf.close();
}
