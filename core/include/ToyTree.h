/*
 * Gamma Combination
 * Author: Till Moritz Karbach, moritz.karbach@cern.ch
 * Date: Dec 2014
 *
 */

#ifndef ToyTree_h
#define ToyTree_h

#include <TString.h>

#include <RooArgSet.h>

#include <map>
#include <string>

class Combiner;
class OptParser;
class PDF_Datasets;

class RooFitResult;
class RooWorkspace;

class TChain;
class TTree;

///
/// Interface class for the root trees that are written
/// by the Plugin method and related functions.
///
class ToyTree {
 public:
  ToyTree(Combiner* c, TChain* t = 0, bool _quiet = false);
  ToyTree(PDF_Datasets* p, const OptParser* opt, TChain* t = 0, bool _quiet = false);

  void activateCoreBranchesOnly();
  void activateAllBranches();
  void activateBranch(const TString& bName);
  void fill();
  void init();
  const OptParser* getArg() { return arg; };
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
  bool isWsVarAngle(TString var);
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

  float scanpoint = 0.f;         ///< the scanpoint for 1D scans, or the x scanpoint for 2D scans
  float scanpointy = 0.f;        ///< the y scanpoint for 2D scans
  float chi2min = 0.f;           ///< the chi2 of the fit with var fixed to scan point
  float chi2minGlobal = 0.f;     ///< the chi2 of the free fit
  float chi2minBkg = 0.f;        ///< the chi2 of the fit of the bkg hypothesis (for CLs method)
  float chi2minToy = 0.f;        ///< the chi2 of the fit to the toy with var fixed to scan point
  float chi2minGlobalToy = 0.f;  ///< the chi2 of the free fit to the toy
  /// The chi2 of the fit of the hypothesis value to the bkg toy distribution (for CLs method)
  float chi2minBkgToy = 0.f;
  float chi2minGlobalBkgToy = 0.f;  ///< the chi2 of the free fit to the bkg only toys
  float chi2minBkgBkgToy = 0.f;     ///< the chi2 of the bkg fit to the bkg only toys
  float scanbest = 0.f;             ///< an alias to the free fit value of the scan variable
  float scanbesty = 0.f;            ///< an alias to the free fit value of the scan y variable in 2D scans
  /// An alias to the free fit value of the scan variable on the bkg only toy (for CLs method)
  float scanbestBkg = 0.f;
  /// An alias to the fit value of the scan variable of the bkg fit on the bkg only toy (for CLs method)
  float scanbestBkgfitBkg = 0.f;
  float nrun = 0.f;    ///< an ID to distinguish different runs, i.e. batch jobs
  float ntoy = 0.f;    ///< an ID to distinguish different toys
  float npoint = 0.f;  ///< an ID to distinguish different scan point
  float id = 0.f;      ///< an ID to distinguish different conditions, e.g. different toys in a coverage test
  float statusFree = -5.f;
  float covQualFree = -2.f;
  float statusScan = -5.f;
  float covQualScan = -2.f;
  float statusFreeBkg = -5.f;
  float covQualFreeBkg = -2.f;
  float statusScanBkg = -5.f;
  float covQualScanBkg = -2.f;
  float statusBkgBkg = -5.f;
  float covQualBkgBkg = -2.f;
  float statusScanData = -5.f;
  float covQualScanData = -2.f;
  int bestIndexScanData = 0;
  float nBergerBoos = 0.f;
  float BergerBoos_id = 0.f;
  float genericProbPValue = 0.f;
  float statusFreePDF = -5.f;
  float statusScanPDF = -5.f;
  float chi2minToyPDF = 0.f;
  float chi2minGlobalToyPDF = 0.f;
  float chi2minBkgToyPDF = 0.f;
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
  TString thName;                  ///< set name of theory parameters, derived from name
  TString globName;                ///< set name of explicit set of global observables

  std::map<std::string, float> parametersScan;  ///< fit result of the scan fit
  std::map<std::string, float> parametersFree;  ///< fit result of the free fit
  std::map<std::string, float> parametersPll;   ///< parameters of the profile likelihood curve of the data
  std::map<std::string, float> observables;     ///< values of the observables
  std::map<std::string, float> theory;          ///< theory parameters (=observables at profile likelihood points)
  std::map<TString, float> constraintMeans;     ///< values of global observables

  float scanpointMin = 0.f;   ///< minimum of the scanpoint, computed by computeMinMaxN().
  float scanpointMax = 0.f;   ///< maximum of the scanpoint, computed by computeMinMaxN().
  int scanpointN = -1.f;      ///< number of different values of the scanpoint, computed by computeMinMaxN().
  float scanpointyMin = 0.f;  ///< minimum of the scanpointy, computed by computeMinMaxN().
  float scanpointyMax = 0.f;  ///< maximum of the scanpointy, computed by computeMinMaxN().
  int scanpointyN = -1.f;     ///< number of different values of the scanpointy, computed by computeMinMaxN().

  bool storeObs = true;    ///< Boolean flag to control storing ToyTree observables, can't store these for DatasetsScans
  bool storeTh = true;     ///< Boolean flag to control storing ToyTree theory parameters. Not needed in DatasetsScans
  bool storeGlob = false;  ///< Boolean flag to control storing ToyTree global observables (handy in DatasetsScans)
  bool quiet = false;
};

#endif
