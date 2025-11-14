#ifndef PDF_DatasetLb2pktaul_h
#define PDF_DatasetLb2pktaul_h

#include <PDF_Datasets.h>

class RooDataSet;
class RooFitResult;
class RooWorkspace;

class PDF_DatasetLb2pktaul : public PDF_Datasets {
 public:
  PDF_DatasetLb2pktaul(RooWorkspace* w);
  RooFitResult* fit(RooDataSet* dataToFit);
  void generateToys(int SeedShift = 0) override;

 private:
  bool drawFitsDebug;  //> for visualizing toys and fit results, only changeable in the code
};

#endif
