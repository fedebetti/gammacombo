#include <ParametersAbs.h>

#include <RooMsgService.h>
#include <RooRealVar.h>

using namespace std;
using namespace RooFit;

Parameter* ParametersAbs::newParameter(TString name) {
  auto p = std::make_unique<Parameter>();
  p->name = name;
  m_parameters.push_back(std::move(p));
  return m_parameters.back().get();
}

Parameter::Range ParametersAbs::range(double min, double max) const {
  Parameter::Range r = {min, max};
  return r;
}

Parameter* ParametersAbs::var(TString name) {
  for (auto&& par : m_parameters) {
    if (par->name == name) return par.get();
  }
  cout << "ParametersAbs::var() : ERROR : no such parameter '" + name + "'." << endl;
  return nullptr;
}

///
/// Ranges:
///  free - here is where Feldman-Cousins-like forbidden regions get implemented
///
RooFit::OwningPtr<RooRealVar> ParametersAbs::get(TString name) {
  for (auto&& par : m_parameters) {
    if (par->name == name) {
      auto r = new RooRealVar(par->name, par->title, par->startvalue, par->unit);
      RooMsgService::instance().setGlobalKillBelow(WARNING);  // else we get messages for range creation
      r->setRange("free", par->free.min, par->free.max);
      r->setRange("phys", par->phys.min, par->phys.max);
      r->setRange("scan", par->scan.min, par->scan.max);
      r->setRange("force", par->force.min, par->force.max);
      r->setRange("bboos", par->bboos.min, par->bboos.max);
      RooMsgService::instance().setGlobalKillBelow(INFO);
      return r;
    }
  }
  cout << "ParametersAbs::get() : ERROR : no such parameter '" + name + "'." << endl;
  return nullptr;
}
