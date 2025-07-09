#ifndef PDF_Gaus2d_h
#define PDF_Gaus2d_h

#include <PDF_Abs.h>

#include <TString.h>

class PDF_Gaus2d : public PDF_Abs {
 public:
  PDF_Gaus2d(TString cObs = "year2013", TString cErr = "year2013", TString cCor = "year2013");
  void buildPdf() override;
  void initObservables() override;
  void initParameters() override;
  void initRelations() override;
  void setCorrelations(TString c) override;
  void setObservables(TString c) override;
  void setUncertainties(TString c) override;
};

#endif
