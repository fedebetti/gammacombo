/**
 * Gamma Combination
 * Author: Till Moritz Karbach, moritz.karbach@cern.ch
 * Date: November 2014
 *
 **/

#ifndef PDF_Cartesian_h
#define PDF_Cartesian_h

#include <PDF_Abs.h>
#include <Utils.h>

class PDF_Cartesian : public PDF_Abs {
 public:
  PDF_Cartesian(TString cObs, TString cErr, TString cCor);
  void buildPdf() override;
  void initObservables() override;
  void initParameters() override;
  void initRelations() override;
  void setCorrelations(TString c) override;
  void setObservables(TString c) override;
  void setUncertainties(TString c) override;
};

#endif
