/**
 * Gamma Combination
 * Author: Till Moritz Karbach, moritz.karbach@cern.ch
 * Date: Dec 2014
 *
 **/

#ifndef PDF_rb_h
#define PDF_rb_h

#include <PDF_Abs.h>

class PDF_rb : public PDF_Abs {
 public:
  PDF_rb(TString cObs, TString cErr, TString cCor);
  ~PDF_rb();
  void buildPdf() override;
  void initObservables() override;
  void initParameters() override;
  void initRelations() override;
  void setCorrelations(TString c) override;
  void setObservables(TString c) override;
  void setUncertainties(TString c) override;
};

#endif
