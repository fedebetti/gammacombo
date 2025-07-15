#ifndef RooSlimFitResult_h
#define RooSlimFitResult_h

#include <fstream>

#include <RooArgList.h>
#include <RooFitResult.h>
#include <RooRealVar.h>

/**
 *  Class that mimics the functionality of RooFitResult, but uses less internal memory by not storing the correlation
 *  matrix, if not specifically requested. It is useful to save memory in 2D scans.
 *  It also contains a few extra getters.
 **/
class RooSlimFitResult : public TObject {
 public:
  /// Default constructor (needed for TObject serialisation).
  RooSlimFitResult() = default;
  RooSlimFitResult(const RooFitResult* r, bool storeCorrelation = false);
  RooSlimFitResult(const RooSlimFitResult* other);
  RooSlimFitResult(const RooSlimFitResult& r);

  std::unique_ptr<RooSlimFitResult> Clone();
  RooArgList& constPars() const;
  inline TMatrixDSym correlationMatrix() const { return _correlationMatrix; };
  inline Int_t covQual() const { return _covQual; };
  inline Double_t edm() const { return _edm; };
  RooArgList& floatParsFinal() const;
  double getParVal(TString name) const;
  double getParErr(TString name) const;
  double getConstParVal(TString name) const;
  double getFloatParFinalVal(TString name) const;
  const std::vector<std::string>& getParsNames() const { return _parsNames; };
  double getParVal(const int i) const { return _parsVal[i]; };
  std::vector<double>& getParsVal() { return _parsVal; };
  bool hasParameter(TString name) const;
  bool isAngle(RooRealVar* v) const;
  inline bool isConfirmed() const { return _isConfirmed; };
  inline Double_t minNll() const { return _minNLL; };
  void Print(bool verbose = false, bool printcor = false) const;
  void SaveLatex(std::ofstream& outfile, bool verbose = false, bool printcor = false);
  inline void setConfirmed(bool c) { _isConfirmed = c; };
  inline Int_t status() const { return _status; };

 private:
  template <class FitResult>
  void init(const FitResult* r, bool storeCorrelation = false);

  std::vector<std::string> _parsNames;  // variable names
  std::vector<int> _parsFloatId;  // ID of floating parameter - this corresponds to the correlation matrix position
  std::vector<double> _parsVal;   // values of const parameters, index given by position in _variable names
  std::vector<double> _parsErr;
  std::vector<bool> _parsAngle;  // is it an angle?
  std::vector<bool> _parsConst;  // is it constant?
  double _edm = std::numeric_limits<double>::quiet_NaN();
  double _minNLL = std::numeric_limits<double>::quiet_NaN();
  int _covQual = -9;
  int _status = -9;
  TMatrixDSym _correlationMatrix{0};

  // dummy variables only needed for constPars() and floatParsFinal():
  mutable RooArgList _constParsDummy;  //! <- The exlcamation mark turns off storing in a root file (marks the member
                                       //! transient)
  mutable RooArgList _floatParsFinalDummy;  //! mutables can be changed in const methods

  bool _isConfirmed = false;

 public:
  ClassDef(RooSlimFitResult, 1)  // defines version number, ClassDef is a macro
};

//
// TEMPLATE IMPLEMENTATION
//

template <class FitResult>
void RooSlimFitResult::init(const FitResult* r, bool storeCorrelation) {
  assert(r);
  // copy over const parameters
  int size = r->constPars().getSize();
  for (int i = 0; i < size; i++) {
    RooRealVar* p = (RooRealVar*)r->constPars().at(i);
    _parsNames.push_back(p->GetName());
    _parsVal.push_back(p->getVal());
    _parsErr.push_back(0.);
    _parsAngle.push_back(isAngle(p));
    _parsConst.push_back(true);
    _parsFloatId.push_back(-1);  // floating ID doesn't exist for constant parameters
  }
  // copy over floating parameters
  size = r->floatParsFinal().getSize();
  for (int i = 0; i < size; i++) {
    RooRealVar* p = (RooRealVar*)r->floatParsFinal().at(i);
    _parsNames.push_back(p->GetName());
    _parsVal.push_back(p->getVal());
    _parsErr.push_back(p->getError());
    _parsAngle.push_back(isAngle(p));
    _parsConst.push_back(false);
    _parsFloatId.push_back(i);  // needed to store the parameter's position in the COR matrix (matches
                                // floatParsFinal())
  }
  // copy over numeric values
  _covQual = r->covQual();
  _edm = r->edm();
  _minNLL = r->minNll();
  _status = r->status();
  // store correlation matrix
  if (storeCorrelation) {
    _correlationMatrix.ResizeTo(r->correlationMatrix());
    _correlationMatrix = r->correlationMatrix();
  }
}

#endif
