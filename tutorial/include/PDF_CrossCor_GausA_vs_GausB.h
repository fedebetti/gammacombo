#ifndef PDF_CrossCor_GausA_vs_GausB_h
#define PDF_CrossCor_GausA_vs_GausB_h

#include <PDF_CrossCorAbs.h>

class PDF_CrossCor_GausA_vs_GausB : public PDF_CrossCorAbs {
 public:
  PDF_CrossCor_GausA_vs_GausB(PDF_Abs* pdf1, PDF_Abs* pdf2, TString cCor = "lumi1fb");
  ~PDF_CrossCor_GausA_vs_GausB();
  void setCorrelations(TString c) override;
};

#endif
