#ifndef LatexMaker_h
#define LatexMaker_h

#include <TMatrixDSym.h>
#include <TString.h>

#include <fstream>
#include <vector>

class PDF_Abs;

class RooArgList;

class LatexMaker {

 public:
  LatexMaker(TString cName, PDF_Abs* _pdf);

  void writeFile();
  void writeCorrMatrix(std::ofstream& file, TMatrixDSym mat, RooArgList* observables, std::vector<TString> labels);

  TString outfname;
  PDF_Abs* pdf = nullptr;
};

#endif
