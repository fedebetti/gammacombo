#include <Parameter.h>

#include <iostream>

Parameter::Parameter() {
  name = "not initialized";
  title = "not initialized";
  unit = "not initialized";
  startvalue = 0;
  Parameter::Range r = {-1e4, 1e4};
  free = r;
  phys = r;
  scan = r;
  force = r;
  bboos = r;
}

void Parameter::Print() {
  std::cout << " name       = " << name << std::endl;
  std::cout << " title      = " << title << std::endl;
  std::cout << " unit       = " << unit << std::endl;
  std::cout << " startvalue = " << startvalue << std::endl;
  std::cout << " phys       = " << phys.min << " ... " << phys.max << std::endl;
  std::cout << " scan       = " << scan.min << " ... " << scan.max << std::endl;
  std::cout << " force      = " << force.min << " ... " << force.max << std::endl;
  std::cout << " bboos      = " << bboos.min << " ... " << bboos.max << std::endl;
  std::cout << " free       = " << free.min << " ... " << free.max << std::endl;
}
