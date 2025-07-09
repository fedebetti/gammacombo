#ifndef PDF_GausB_h
#define PDF_GausB_h

#include <PDF_Abs.h>

#include <TString.h>

class PDF_GausB : public PDF_Abs {
 public:
  PDF_GausB(TString cObs = "year2014", TString cErr = "year2014", TString cCor = "year2014");
  void buildPdf() override;
  void initObservables() override;
  void initParameters() override;
  void initRelations() override;
  void setCorrelations(TString c) override;
  void setObservables(TString c) override;
  void setUncertainties(TString c) override;
};

#endif
