#ifndef LatexMaker_h
#define LatexMaker_h

#include <fstream>
#include <iostream>
#include <vector>

#include <TMatrixDSym.h>
#include <TString.h>

#include <RooArgList.h>

#include "PDF_Abs.h"

class LatexMaker {

 public:
  LatexMaker(TString cName, PDF_Abs* _pdf);
  ~LatexMaker();

  void writeFile();
  void writeCorrMatrix(ofstream& file, TMatrixDSym mat, RooArgList* observables, std::vector<TString> labels);

  TString outfname;
  PDF_Abs* pdf;
};

#endif
