#ifndef ParametersTutorial_h
#define ParametersTutorial_h

#include <ParametersAbs.h>

class ParametersTutorial : public ParametersAbs {
 public:
  ParametersTutorial();

 private:
  void defineParameters() override;
};

#endif
