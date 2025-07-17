#ifndef MethodAbsScan_h
#define MethodAbsScan_h

#include <memory>
#include <optional>
#include <type_traits>
#include <utility>
#include <vector>

#include <RooDataSet.h>

#include <TH1.h>
#include <TH2.h>
#include <TRandom3.h>

#include "CLInterval.h"
#include "Combiner.h"
#include "OptParser.h"
#include "PValueCorrection.h"
#include "RooSlimFitResult.h"

class OneMinusClPlotAbs;

/**
 * Abstract class to perform a likelihood scan over the scan range of one or two floating parameters, and save the
 * related results.
 */
class MethodAbsScan {
 public:
  MethodAbsScan() = default;
  MethodAbsScan(Combiner* c);
  MethodAbsScan(const OptParser* opt);

  virtual void calcCLintervals(const int CLsType = 0, const bool calc_expected = false, const bool quiet = false);
  void confirmSolutions();
  void doInitialFit(bool force = false);
  inline const OptParser* getArg() const { return arg; };
  inline const std::vector<std::unique_ptr<RooSlimFitResult>>& getAllResults() const { return allResults; };
  inline const std::vector<RooSlimFitResult*>& getCurveResults() const { return curveResults; };
  inline double getChi2minGlobal() const { return chi2minGlobal; }
  inline double getChi2minBkg() const { return chi2minBkg; }
  double getCL(double val) const;
  const CLInterval* getCLintervalCentral(int sigma = 1, bool quiet = false);
  const CLInterval* getCLinterval(const int iSol = 0, const int index = 0, const bool quiet = false);
  inline Combiner* getCombiner() const { return combiner; };
  int getDrawSolution() const { return drawSolution; }
  inline bool getFilled() const { return drawFilled; };
  inline TH1* getHCL() { return hCL.get(); };
  inline TH1* getHCLs() { return hCLs.get(); };
  inline TH1* getHCLsFreq() { return hCLsFreq.get(); };
  inline TH1* getHCLsExp() { return hCLsExp.get(); };
  inline TH1* getHCLsErr1Up() { return hCLsErr1Up.get(); };
  inline TH1* getHCLsErr1Dn() { return hCLsErr1Dn.get(); };
  inline TH1* getHCLsErr2Up() { return hCLsErr2Up.get(); };
  inline TH1* getHCLsErr2Dn() { return hCLsErr2Dn.get(); };
  inline TH2* getHCL2d() { return hCL2d.get(); };
  inline TH2* getHCLs2d() { return hCLs2d.get(); };
  inline TH1* getHchisq() { return hChi2min.get(); };
  inline TH2* getHchisq2d() { return hChi2min2d.get(); };
  inline int getLineColor() const { return lineColor; };
  inline int getLineStyle() const { return lineStyle; };
  inline int getLineWidth() const { return lineWidth; };
  inline int getFillStyle() const { return fillStyle; };
  inline int getFillColor() const { return fillColor; };
  inline double getFillTransparency() const { return fillTransparency; };
  inline TString getMethodName() const { return methodName; };
  inline TString getName() const { return name; };
  inline int getNObservables() const { return w->set(obsName)->getSize(); }
  inline int getNPoints1d() const { return nPoints1d; }
  inline int getNPoints2dx() const { return nPoints2dx; }
  inline int getNPoints2dy() const { return nPoints2dy; }
  inline const RooArgSet* getObservables() const { return w->set(obsName); }
  inline TString getObsName() const { return obsName; };
  inline TString getParsName() const { return parsName; };
  double getScanVarSolution(const int iVar, const int iSol);
  RooRealVar* getScanVar1();
  TString getScanVar1Name() const { return scanVar1; }
  double getScanVar1Solution(int i = 0);
  RooRealVar* getScanVar2();
  TString getScanVar2Name() const { return scanVar2; }
  double getScanVar2Solution(int i = 0);
  inline const std::vector<std::unique_ptr<RooSlimFitResult>>& getSolutions() { return solutions; };
  inline const int getNSolutions() { return solutions.size(); };
  RooSlimFitResult* getSolution(const int i = 0);
  inline const RooArgSet* getTheory() { return w->set(thName); }
  inline int getTextColor() const { return textColor; };
  inline TString getTitle() const { return title; };
  inline RooWorkspace* getWorkspace() { return w; };
  virtual void initScan();
  void loadParameters(const RooSlimFitResult* r);
  bool loadSolution(const int i = 0);
  virtual bool loadScanner(TString fName = "");
  void plot2d(TString varx, TString vary);
  void plot1d(TString var);
  void plotOn(OneMinusClPlotAbs* plot, int CLsType = 0);
  void plotPulls(int nSolution = 0);
  virtual void print() const;
  void printCLintervals(const int CLsType, const bool calc_expected = false);
  void printLocalMinima() const;
  void saveLocalMinima(TString fName = "") const;
  void saveScanner(TString fName = "") const;
  virtual int scan1d() const;
  virtual int scan2d() const;
  inline void setDrawSolution(int code = 0) { drawSolution = code; };
  inline void setPValueCorrector(PValueCorrection* pvalCor) {
    pvalueCorrector = pvalCor;
    pvalueCorrectorSet = true;
  }
  inline void setScanVar1(TString var) { scanVar1 = var; };
  inline void setScanVar2(TString var) { scanVar2 = var; };
  inline void setNPoints1d(int n) { nPoints1d = n; };
  inline void setNPoints2dx(int n) { nPoints2dx = n; };
  inline void setNPoints2dy(int n) { nPoints2dy = n; };
  inline void setFilled(bool filled) { drawFilled = filled; };
  inline void setLineColor(int c) { lineColor = c; };
  inline void setLineStyle(int c) { lineStyle = c; };
  inline void setLineWidth(int c) { lineWidth = c; };
  inline void setTextColor(int c) { textColor = c; };
  inline void setFillStyle(int c) { fillStyle = c; };
  inline void setFillColor(int c) { fillColor = c; };
  inline void setFillTransparency(double c) { fillTransparency = c; };
  inline void setTitle(TString s) { title = s; };
  void setChi2minGlobal(double x);
  void setSolutions(const std::vector<std::unique_ptr<RooSlimFitResult>>& s);
  inline void setVerbose(bool yesNo = true) { verbose = yesNo; };
  inline void setHCL(TH1* h) { hCL = std::unique_ptr<TH1>(h); };
  inline void setHchisq(TH1* h) { hChi2min = std::unique_ptr<TH1>(h); };
  void setXscanRange(double min, double max);
  void setYscanRange(double min, double max);
  void calcCLintervalsSimple(int CLsType = 0, bool calc_expected = false);
  const std::pair<double, double> getBorders(const TGraph& graph, const double confidence_level,
                                             bool qubic = false) const;
  const std::pair<double, double> getBorders_CLs(const TGraph& graph, const double confidence_level,
                                                 bool qubic = false) const;
  virtual bool checkCLs() const;

  /// All fit results we encounter along the scan.
  std::vector<std::unique_ptr<RooSlimFitResult>> allResults;
  /// All fit results of the the points that make it into the 1-CL curve. Index is the bin number of hCL bins -1.
  std::vector<RooSlimFitResult*> curveResults;
  /// All fit results of the the points that make it into the 1-CL curve. Index is the gobal bin number of hCL2d -1.
  std::vector<std::vector<RooSlimFitResult*>> curveResults2d;

  /// Local minima filled by saveSolutions() and saveSolutions2d().
  /// The names of the CL interval std::vectors might be misleading. They correspond to the default CL intervals.
  /// If the option --CL is given, the 1-3 sigma correspond to the first, second,... given value of the CL.
  std::vector<std::unique_ptr<RooSlimFitResult>> solutions;

  /**
   * All CL intervals found by `calcCLintervals`.
   *
   * The first index corresponds to the CL values set (1, 2, 3 sigma by default).
   * The second index corresponds to the solutions from 0 to n (in case they fall in the scan range, otherwise the value
   * will be nullptr). There may be also two more CL intervals at the end, corresponding e.g. to scans starting from the
   * minimum and maximum values of the scan range.
   */
  std::vector<std::vector<std::unique_ptr<CLInterval>>> clintervals;

  std::unique_ptr<RooFitResult> globalMin;  ///< parameter values at a global minimum

 protected:
  void sortSolutions();

  TString name;                ///< basename, e.g. ggsz
  TString title;               ///< nice string for the legends
  TString methodName = "Abs";  ///< Prob, ...
  TString pdfName;             ///< PDF name in workspace, derived from name
  TString obsName;             ///< dataset name of observables, derived from name
  TString parsName;            ///< set name of physics parameters, derived from name
  TString thName;              ///< set name of theory parameters, derived from name
  TString toysName;            ///< set name of parameters to vary in toys
  TString scanVar1;            ///< scan parameter
  TString scanVar2;            ///< second scan parameter if we're scanning 2d
  int nPoints1d;               ///< number of scan points used by 1d scan
  int nPoints2dx;              ///< number of scan points used by 2d scan, x axis
  int nPoints2dy;              ///< number of scan points used by 2d scan, y axis

  PValueCorrection* pvalueCorrector;  // object which can correct the pvalue for undercoverage if required
  bool pvalueCorrectorSet = false;

  TRandom3 rndm;
  RooWorkspace* w = nullptr;
  std::unique_ptr<RooDataSet> obsDataset;  ///< save the nominal observables so we can restore them after we have fitted
                                           ///< toys
  std::unique_ptr<RooDataSet> startPars;   ///< save the start parameter values before any scan
  std::unique_ptr<TH1> hCL;                ///< 1-CL curve
  std::unique_ptr<TH1> hCLs;               ///< 1-CL curve
  std::unique_ptr<TH1> hCLsFreq;           ///< 1-CL curve
  std::unique_ptr<TH1> hCLsExp;            ///< 1-CL curve
  std::unique_ptr<TH1> hCLsErr1Up;         ///< 1-CL curve
  std::unique_ptr<TH1> hCLsErr1Dn;         ///< 1-CL curve
  std::unique_ptr<TH1> hCLsErr2Up;         ///< 1-CL curve
  std::unique_ptr<TH1> hCLsErr2Dn;         ///< 1-CL curve
  std::unique_ptr<TH2> hCL2d;              ///< 1-CL curve
  std::unique_ptr<TH2> hCLs2d;             ///< 1-CL curve
  std::unique_ptr<TH1> hChi2min;           ///< histogram for the chi2min values before Prob()
  std::unique_ptr<TH2> hChi2min2d;         ///< histogram for the chi2min values before Prob()
  double chi2minGlobal = std::numeric_limits<double>::max();  ///< chi2 value at global minimum
  double chi2minBkg = std::numeric_limits<double>::max();     ///< chi2 value at global minimum
  bool chi2minGlobalFound = false;                            ///< flag to avoid finding minimum twice
  int lineColor = kBlue - 8;
  int textColor = kBlack;  ///< color used for plotted central values
  int lineStyle = 0;
  int lineWidth = 2;
  int fillStyle = 1001;
  int fillColor = kBlue - 8;
  double fillTransparency;
  bool drawFilled = true;  ///< choose if Histogram is drawn filled or not
  int drawSolution = 0;    ///< Configure how to draw solutions on the plots.
  ///< 0=don't plot, 1=plot at central value (1d) or markers (2d)
  ///< Default is taken from arg, unless disabled by setDrawSolution().
  bool verbose;
  int nWarnings = 0;                     ///< number of warnings printed in getScanVarSolution()
  const OptParser* arg;                  ///< command line options
  Combiner* combiner = nullptr;          ///< the combination
  bool m_xrangeset = false;              ///< true if the x range was set manually (setXscanRange())
  bool m_yrangeset = false;              ///< true if the y range was set manually (setYscanRange())
  bool m_initialized = false;            ///< true if initScan() was called
  std::vector<double> ConfidenceLevels;  ///< container of the confidence levels to be computed

 private:
  bool compareSolutions(const RooSlimFitResult* r1, const RooSlimFitResult* r2) const;
  double pq(double p0, double p1, double p2, double y, int whichSol = 0) const;
  void removeDuplicateSolutions();
  std::optional<std::pair<double, double>> interpolate(TH1* h, const int i, const double y, const double central,
                                                       const bool upper) const;
  std::optional<double> interpolateLinear(const TH1* h, const int i, const double y) const;
};

#endif
