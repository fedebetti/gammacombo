#ifndef PDF_CrossCorAbs_h
#define PDF_CrossCorAbs_h

#include <PDF_Abs.h>
#include <UtilsConfig.h>

class PDF_CrossCorAbs : public PDF_Abs {
 public:
  PDF_CrossCorAbs(PDF_Abs* pdf1, PDF_Abs* pdf2);
  void buildPdf() override;
  void initObservables() override;
  void initParameters() override;
  void initRelations() override;
  void setCorrelations(Utils::config c);

 protected:
  void copyMeasurementCovariance();
  PDF_Abs* pdf1;
  PDF_Abs* pdf2;
  int nObs1;
  int nObs2;
};

#endif
