/**
 * Gamma Combination
 * Author: Till Moritz Karbach, moritz.karbach@cern.ch
 * Date: July 2014
 *
 **/

#ifndef PDF_CrossCorAbs_h
#define PDF_CrossCorAbs_h

#include "RooGenericPdf.h"
#include "RooCrossCorPdf.h"
#include "PDF_Abs.h"
#include "ParametersTutorial.h"

using namespace RooFit;
using namespace std;
using namespace Utils;

class PDF_CrossCorAbs : public PDF_Abs
{
public:
  PDF_CrossCorAbs(PDF_Abs* pdf1, PDF_Abs* pdf2);
  virtual ~PDF_CrossCorAbs();
  void  buildPdf() override;
  void  initObservables() override;
  void  initParameters() override;
  void  initRelations() override;
  void  setCorrelations(config c);

protected:
    void copyMeasurementCovariance();
    PDF_Abs* pdf1;
    PDF_Abs* pdf2;
    int nObs1;
    int nObs2;
};

#endif
