#include <PDF_Gaus.h>
#include <ParametersTutorial.h>

#include <Utils.h>

#include <RooArgList.h>
#include <RooArgSet.h>
#include <RooFormulaVar.h>
#include <RooMultiVarGaussian.h>

using namespace std;

PDF_Gaus::PDF_Gaus(TString cObs, TString cErr, TString cCor) : PDF_Abs(1) {
  name = "Gaus";
  initParameters();
  initRelations();
  initObservables();
  setObservables(cObs);
  setUncertainties(cErr);
  setCorrelations(cCor);
  build();
}

void PDF_Gaus::initParameters() {
  ParametersTutorial p;
  parameters = new RooArgList("parameters");
  parameters->add(*(p.get("a_gaus")));
}

void PDF_Gaus::initRelations() {
  theory = new RooArgList("theory");  ///< the order of this list must match that of the COR matrix!
  theory->add(*(Utils::makeTheoryVar("a_gaus_th", "a_gaus_th", "a_gaus", parameters)));
}

void PDF_Gaus::initObservables() {
  observables = new RooArgList("observables");  ///< the order of this list must match that of the COR matrix!
  observables->add(*(new RooRealVar("a_gaus_obs", "a_gaus_obs", 0, -1e4, 1e4)));
}

void PDF_Gaus::setObservables(TString c) {
  if (c.EqualTo("truth")) {
    setObservablesTruth();
  } else if (c.EqualTo("toy")) {
    setObservablesToy();
  } else if (c.EqualTo("year2013")) {
    obsValSource = c;
    setObservable("a_gaus_obs", -0.5);
  } else if (c.EqualTo("year2014")) {
    obsValSource = c;
    setObservable("a_gaus_obs", 1.5);
  } else {
    cout << "PDF_Gaus::setObservables() : ERROR : config " + c + " not found." << endl;
    exit(1);
  }
}

void PDF_Gaus::setUncertainties(TString c) {
  if (c.EqualTo("year2013")) {
    obsErrSource = c;
    StatErr[0] = 1;  // a_gaus
    SystErr[0] = 0;  // a_gaus
  } else if (c.EqualTo("year2014")) {
    obsErrSource = c;
    StatErr[0] = 0.5;   // a_gaus
    SystErr[0] = 0.15;  // a_gaus
  } else {
    cout << "PDF_Gaus::setUncertainties() : ERROR : config " + c + " not found." << endl;
    exit(1);
  }
}

void PDF_Gaus::setCorrelations(TString c) {
  resetCorrelations();
  if (c.EqualTo("year2013") || c.EqualTo("year2014")) {
    corSource = "no correlations for 1 obs";
  } else {
    cout << "PDF_Gaus::setCorrelations() : ERROR : config " + c + " not found." << endl;
    exit(1);
  }
}

void PDF_Gaus::buildPdf() {
  pdf = new RooMultiVarGaussian("pdf_" + name, "pdf_" + name, *(RooArgSet*)observables, *(RooArgSet*)theory, covMatrix);
}
