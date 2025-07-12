#ifndef ParametersAbs_h
#define ParametersAbs_h

#include "Parameter.h"

#include <RooRealVar.h>

#include <memory>
#include <vector>

/**
 * Abstract class to define the (nuisance) parameters.
 */
class ParametersAbs {
 public:
  Parameter* var(const TString name);
  RooFit::OwningPtr<RooRealVar> get(const TString name);
  Parameter* newParameter(const TString name);
  Parameter::Range range(const double min, const double max) const;

 protected:
  std::vector<std::unique_ptr<Parameter>> m_parameters;
  inline virtual void defineParameters() = 0;

 private:
  Parameter* parameterExists(const TString name);
};

#endif
