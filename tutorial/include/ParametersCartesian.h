#ifndef ParametersCartesian_h
#define ParametersCartesian_h

#include <ParametersAbs.h>

class ParametersCartesian : public ParametersAbs {
 public:
  ParametersCartesian();

 protected:
  void defineParameters() override;
};

#endif
