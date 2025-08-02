#include <ParametersAbs.h>

#include <RooMsgService.h>
#include <RooRealVar.h>

#include <memory>
#include <stdexcept>
#include <utility>

/**
 * Create a new parameter and add it to the parameters list.
 *
 * Throws a `runtime_error` exception if a parameter with the given name exists already.
 */
Parameter* ParametersAbs::newParameter(const TString name) {
  if (parameterExists(name))
    throw std::runtime_error(std::string("ParametersAbs::newParameter() : parameter " + name + " exist already"));
  auto p = std::make_unique<Parameter>();
  p->name = name;
  m_parameters.push_back(std::move(p));
  return m_parameters.back().get();
}

/**
 * Get raw pointer to a parameter with a given name. Return `nullptr` if name doesn't exist.
 */
Parameter* ParametersAbs::parameterExists(const TString name) {
  for (auto&& par : m_parameters) {
    if (par->name == name) return par.get();
  }
  return nullptr;
}

Parameter::Range ParametersAbs::range(const double min, const double max) const {
  Parameter::Range r = {min, max};
  return r;
}

/**
 * Get raw pointer to a parameter with a given name. Return `nullptr` and raise an error if name doesn't exist.
 */
Parameter* ParametersAbs::var(const TString name) {
  auto par = parameterExists(name);
  if (!par) std::cerr << "ParametersAbs::var() : ERROR : no such parameter '" + name + "'." << std::endl;
  return par;
}

///
/// Ranges:
///  free - here is where Feldman-Cousins-like forbidden regions get implemented
///
RooFit::OwningPtr<RooRealVar> ParametersAbs::get(const TString name) {
  for (auto&& par : m_parameters) {
    if (par->name == name) {
      auto r = new RooRealVar(par->name, par->title, par->startvalue, par->unit);
      RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);  // else we get messages for range creation
      r->setRange("phys", par->phys.min, par->phys.max);
      r->setRange("scan", par->scan.min, par->scan.max);
      r->setRange("force", par->force.min, par->force.max);
      r->setRange("bboos", par->bboos.min, par->bboos.max);
      RooMsgService::instance().setGlobalKillBelow(RooFit::INFO);
      return r;
    }
  }
  std::cout << "ParametersAbs::get() : ERROR : no such parameter '" + name + "'." << std::endl;
  return nullptr;
}
