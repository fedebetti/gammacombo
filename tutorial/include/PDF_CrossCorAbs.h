/**
 * Gamma Combination
 * Author: Till Moritz Karbach, moritz.karbach@cern.ch
 * Date: July 2014
 *
 **/

#ifndef PDF_CrossCorAbs_h
#define PDF_CrossCorAbs_h

#include <PDF_Abs.h>
#include <UtilsConfig.h>

class PDF_CrossCorAbs : public PDF_Abs {
 public:
  PDF_CrossCorAbs(PDF_Abs* pdf1, PDF_Abs* pdf2);
  ~PDF_CrossCorAbs();
  void buildPdf();
  void initObservables();
  virtual void initParameters();
  virtual void initRelations();
  virtual void setCorrelations(Utils::config c);

 protected:
  void copyMeasurementCovariance();
  PDF_Abs* pdf1;
  PDF_Abs* pdf2;
  int nObs1;
  int nObs2;
};

#endif
