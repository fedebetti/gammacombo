/**
 * Gamma Combination
 * Author: Till Moritz Karbach, moritz.karbach@cern.ch
 * Date: August 2012
 *
 **/

#ifndef Combiner_h
#define Combiner_h

#include "OptParser.h"
#include "PDF_Abs.h"
#include "Utils.h"

class Combiner {
 public:
  Combiner(const OptParser* arg, TString title);
  Combiner(const OptParser* arg, TString name, TString title);
  ~Combiner();

  /// Add a PDF to the combiner.
  void addPdf(PDF_Abs* p);

  /// Add an arbitrary number of PDFs to the combiner.
  template <typename... Args>
  void addPdf(PDF_Abs* p, Args... args) {
    addPdf(p);
    addPdf(args...);
  }

  void adjustPhysRange(const TString varName, const float min, const float max);
  Combiner* Clone(const TString name, const TString title);
  void combine();
  void fixParameter(const TString var, const float value);
  void fixParameters(TString vars);

  /// Remove a PDF from the combiner.
  void delPdf(PDF_Abs* p);

  /// Remove an arbitrary number of PDFs from the combiner.
  template <typename... Args>
  void delPdf(PDF_Abs* p, Args... args) {
    delPdf(p);
    delPdf(args...);
  }

  inline const OptParser* getArg() { return arg; };
  const RooArgSet* getParameters() const;
  std::vector<std::string>& getParameterNames() const;
  PDF_Abs* getPdfProvidingObservable(const TString obsname);
  const RooArgSet* getObservables() const;
  std::vector<std::string>& getObservableNames() const;
  inline TString getTitle() const { return title; };
  inline TString getName() const { return name; };
  inline TString getPdfName() const { return pdfName; };  ///< Returns name of combined pdf. Call combine() first.
  inline TString getParsName() const {
    return parsName;
  };  ///< Returns name of combined parameter set. Call combine() first.
  inline TString getObsName() const {
    return obsName;
  };  ///< Returns name of combined observables set. Call combine() first.
  RooAbsPdf* getPdf();
  inline std::vector<PDF_Abs*>& getPdfs() { return pdfs; };
  inline RooWorkspace* getWorkspace() { return w; };
  inline bool isCombined() const { return _isCombined; }
  void loadParameterLimits();
  void print() const;
  void replacePdf(PDF_Abs* from, PDF_Abs* to);
  void setName(const TString name);
  void setObservablesToToyValues();
  void setParametersConstant();  // helper function for combine()
  inline void setTitle(const TString title) { this->title = title; };
  inline std::vector<Utils::FixPar> getConstVars() { return constVars; };

 private:
  std::vector<PDF_Abs*> pdfs;         // holds all pdfs to be combined
  TString title;                      // title of the combination, used in plots
  TString name = "";                  // name of the combination, used to refer to it and as part of file names
  TString pdfName = "";               // Name of combined pdf. Call combine() first.
  TString parsName;                   // Name of combined parameter set. Call combine() first.
  TString obsName;                    // Name of combined observables set. Call combine() first.
  const OptParser* arg;               // command line arguments
  RooWorkspace* w;                    // holds all input pdfs, parameters, and observables, as well as the combination
  std::vector<std::string> pdfNames;  // hold all unique names of the pdfs to be combined
  bool _isCombined;                   // make sure we'll only combine once - else all PDFs get double counted!
  std::vector<Utils::FixPar> constVars;  // hold variables that will be set constant (filled by fixParameter())
};

#endif
