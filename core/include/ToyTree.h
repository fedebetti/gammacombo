/*
 * Gamma Combination
 * Author: Till Moritz Karbach, moritz.karbach@cern.ch
 * Date: Dec 2014
 *
 */

#ifndef ToyTree_h
#define ToyTree_h

#include <TChain.h>

#include "OptParser.h"
#include "PDF_Datasets.h"

///
/// Interface class for the root trees that are written
/// by the Plugin method and related functions.
///
class ToyTree {
 public:
  ToyTree(Combiner* c, TChain* t = 0, bool _quiet = false);
  ToyTree(PDF_Datasets* p, const OptParser* opt, TChain* t = 0, bool _quiet = false);
  ~ToyTree(){};

  void activateCoreBranchesOnly();
  void activateAllBranches();
  void activateBranch(const TString& bName);
  void fill();
  void init();
  const OptParser* getArg() const { return arg; };
  Long64_t GetEntries() const;
  void GetEntry(Long64_t i);
  inline TString getName() const { return name; };
  float getScanpointMin();
  float getScanpointMax();
  int getScanpointN();
  float getScanpointyMin();
  float getScanpointyMax();
  int getScanpointyN();
  TTree* getTree() { return t; };
  bool isWsVarAngle(TString var) const;
  void open();
  void setCombiner(Combiner* c);
  void storeParsPll();
  void storeParsFree();
  void storeParsScan();
  void storeParsScan(RooFitResult* values);
  void storeTheory();
  void storeObservables();
  void writeToFile(TString fName);
  void writeToFile();
  void setStoreObs(bool flag) { this->storeObs = flag; };
  void setStoreTh(bool flag) { this->storeTh = flag; };
  void setStoreGlob(bool flag) { this->storeTh = flag; };
  void storeParsGau(const RooArgSet globalConstraintMeans);

  float scanpoint;              ///< the scanpoint for 1D scans, or the x scanpoint for 2D scans
  float scanpointy = 0.;        ///< the y scanpoint for 2D scans
  float chi2min = 0.;           ///< the chi2 of the fit with var fixed to scan point
  float chi2minGlobal = 0.;     ///< the chi2 of the free fit
  float chi2minBkg = 0.;        ///< the chi2 of the fit of the bkg hypothesis (for CLs method)
  float chi2minToy = 0.;        ///< the chi2 of the fit to the toy with var fixed to scan point
  float chi2minGlobalToy = 0.;  ///< the chi2 of the free fit to the toy
  float chi2minBkgToy =
      0.;  ///< the chi2 of the fit of the hypothesis value to the bkg toy distribution (for CLs method)
  float chi2minGlobalBkgToy = 0.;  ///< the chi2 of the free fit to the bkg only toys
  float chi2minBkgBkgToy = 0.;     ///< the chi2 of the bkg fit to the bkg only toys
  float scanbest = 0.;             ///< an alias to the free fit value of the scan variable
  float scanbesty = 0.;            ///< an alias to the free fit value of the scan y variable in 2D scans
  float scanbestBkg = 0.;  ///< an alias to the free fit value of the scan variable on the bkg only toy (for CLs method)
  float scanbestBkgfitBkg =
      0.;  ///< an alias to the fit value of the scan variable of the bkg fit on the bkg only toy (for CLs method)
  float nrun = 0.;    ///< an ID to distinguish different runs, i.e. batch jobs
  float ntoy = 0.;    ///< an ID to distinguish different toys
  float npoint = 0.;  ///< an ID to distinguish different scan point
  float id = 0.;      ///< an ID to distinguish different conditions, e.g. different toys in a coverage test
  float statusFree = -5.;
  float covQualFree = -2.;
  float statusScan = -5.;
  float covQualScan = -2.;
  float statusFreeBkg = -5.;
  float covQualFreeBkg = -2.;
  float statusScanBkg = -5.;
  float covQualScanBkg = -2.;
  float statusBkgBkg = -5.;
  float covQualBkgBkg = -2.;
  float statusScanData = -5.;
  float covQualScanData = -2.;
  int bestIndexScanData = 0;
  float nBergerBoos = 0.;
  float BergerBoos_id = 0.;
  float genericProbPValue = 0.;
  float statusFreePDF = -5.;
  float statusScanPDF = -5.;
  float chi2minToyPDF = 0.;
  float chi2minGlobalToyPDF = 0.;
  float chi2minBkgToyPDF = 0.;
  TTree* t = nullptr;  ///< the tree

 private:
  void computeMinMaxN();
  Combiner* comb = nullptr;        ///< combination bringing in the arg, workspace, and names
  const OptParser* arg = nullptr;  ///< command line arguments
  RooWorkspace* w = nullptr;       ///< holds all input pdfs, parameters, and observables, as well as the combination
  TString name;                    ///< combiner name, ending up in titles and file names
  TString pdfName;                 ///< PDF name in workspace, derived from name
  TString obsName;                 ///< dataset name of observables, derived from name
  TString parsName;                ///< set name of physics parameters, derived from name
  TString thName = "";             ///< set name of theory parameters, derived from name
  TString globName;                ///< set name of explicit set of global observables

  std::map<std::string, float> parametersScan;  ///< fit result of the scan fit
  std::map<std::string, float> parametersFree;  ///< fit result of the free fit
  std::map<std::string, float> parametersPll;   ///< parameters of the profile likelihood curve of the data
  std::map<std::string, float> observables;     ///< values of the observables
  std::map<std::string, float> theory;          ///< theory parameters (=observables at profile likelihood points)
  std::map<TString, float> constraintMeans;     ///< values of global observables

  float scanpointMin = 0.;   ///< minimum of the scanpoint, computed by computeMinMaxN().
  float scanpointMax = 0.;   ///< maximum of the scanpoint, computed by computeMinMaxN().
  int scanpointN = -1.;      ///< number of different values of the scanpoint, computed by computeMinMaxN().
  float scanpointyMin = 0.;  ///< minimum of the scanpointy, computed by computeMinMaxN().
  float scanpointyMax = 0.;  ///< maximum of the scanpointy, computed by computeMinMaxN().
  int scanpointyN = -1.;     ///< number of different values of the scanpointy, computed by computeMinMaxN().

  bool storeObs;   ///< Boolean flag to control storing ToyTree observables, can't store these for DatasetsScans
  bool storeTh;    ///< Boolean flag to control storing ToyTree theory parameters. Not needed in DatasetsScans
  bool storeGlob;  ///< Boolean flag to control storing ToyTree global observables. Extremely handy in DatasetsScans
  bool quiet = false;
};

#endif
