#include <Parameter.h>

#include <format>
#include <iostream>

const Parameter::Range Parameter::default_range = {-1e4, 1e4};

void Parameter::Print() const {
  std::cout << std::format(" name       = {:s}\n", name.Data());
  std::cout << std::format(" title      = {:s}\n", title.Data());
  std::cout << std::format(" unit       = {:s}\n", unit.Data());
  std::cout << std::format(" startvalue = {:f}\n", startvalue)
            << std::format(" phys       = {:f} ... {:f}\n", phys.min, phys.max)
            << std::format(" scan       = {:f} ... {:f}\n", scan.min, scan.max)
            << std::format(" force      = {:f} ... {:f}\n", force.min, force.max)
            << std::format(" bboos      = {:f} ... {:f}\n", bboos.min, bboos.max)
            << std::format(" free       = {:f} ... {:f}", free.min, free.max) << std::endl;
}
