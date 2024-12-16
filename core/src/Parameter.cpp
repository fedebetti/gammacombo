#include <Parameter.h>

#include <iostream>

using namespace std;

const Parameter::Range Parameter::default_range = {-1e4, 1e4};

void Parameter::Print() const
{
  cout << " name       = " << name << endl;
  cout << " title      = " << title << endl;
  cout << " unit       = " << unit  << endl;
  cout << " startvalue = " << startvalue  << endl;
  cout << " phys       = " << phys.min  << " ... " << phys.max << endl;
  cout << " scan       = " << scan.min  << " ... " << scan.max<< endl;
  cout << " force      = " << force.min << " ... " << force.max << endl;
  cout << " bboos      = " << bboos.min << " ... " << bboos.max << endl;
  cout << " free       = " << free.min  << " ... " << free.max << endl;
}
