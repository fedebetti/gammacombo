#ifndef GammaComboEngine_h
#define GammaComboEngine_h

#include "BatchScriptWriter.h"
#include "Combiner.h"
#include "FileNameBuilder.h"
#include "MethodBergerBoosScan.h"
#include "MethodCoverageScan.h"
#include "MethodPluginScan.h"
#include "MethodProbScan.h"
#include "OneMinusClPlotAbs.h"
#include "OptParser.h"
#include "PDF_Abs.h"
#include "ParameterCache.h"

#include <TApplication.h>
#include <TStopwatch.h>

#include <memory>
#include <vector>

/**
 * Main GammaCombo scanning engine, controlling the application.
 */
class GammaComboEngine {
 public:
  GammaComboEngine(TString name, int argc, char* argv[]);
  GammaComboEngine(TString name, int argc, char* argv[], bool _runOnDataSet);

  void adjustRanges(Combiner* c, const int cId);
  void setupToyVariationSets(Combiner* c, const int cId);
  void addPdf(int id, PDF_Abs* pdf, TString title = "");
  void addSubsetPdf(int id, PDF_Abs* pdf, const std::vector<int>& indices, TString title = "");

  // The following methods are deprecated and maintained only to preserve backward-compatibility.
  void addSubsetPdf(int id, PDF_Abs* pdf, int i1, TString title = "");
  void addSubsetPdf(int id, PDF_Abs* pdf, int i1, int i2, TString title = "");
  void addSubsetPdf(int id, PDF_Abs* pdf, int i1, int i2, int i3, TString title = "");
  void addSubsetPdf(int id, PDF_Abs* pdf, int i1, int i2, int i3, int i4, TString title = "");
  void addSubsetPdf(int id, PDF_Abs* pdf, int i1, int i2, int i3, int i4, int i5, TString title = "");
  void addSubsetPdf(int id, PDF_Abs* pdf, int i1, int i2, int i3, int i4, int i5, int i6, TString title = "");
  void addSubsetPdf(int id, PDF_Abs* pdf, int i1, int i2, int i3, int i4, int i5, int i6, int i7, TString title = "");
  void addSubsetPdf(int id, PDF_Abs* pdf, int i1, int i2, int i3, int i4, int i5, int i6, int i7, int i8,
                    TString title = "");

  void setPdf(PDF_Abs* pdf);
  void addCombiner(int id, std::unique_ptr<Combiner> cmb);
  void cloneCombiner(int newId, int oldId, TString name, TString title);
  Combiner* getCombiner(const int id) const;
  std::vector<int> getCombinersIds() const;
  PDF_Abs* getPdf(int id) const;
  inline OptParser* getArg() const { return arg.get(); };

  void newCombiner(const int id, const TString name, const TString title, const std::vector<int>& pdfIds = {});

  /**
   * This method is deprecated and is maintained only for legacy reasons.
   *
   * Add a new Combiner, consisting of the specified PDFs. The pdf arguments refer to the GammaComboEngine ID of the
   * PDFs that should be combined (add them before using `addPdf()`).
   */
  void newCombiner(const int id, const TString name, const TString title, std::convertible_to<int> auto&&... pdfIds) {
    newCombiner(id, name, title, std::vector<int>{pdfIds...});
  }

  void print() const;
  void printPdfs() const;
  void printCombinations() const;
  void run();
  void runApplication();
  void scanStrategy1d(MethodProbScan* scanner, ParameterCache* pCache);
  void scanStrategy2d(MethodProbScan* scanner, ParameterCache* pCache);
  inline void setRunOnDataSet(bool opt) { runOnDataSet = opt; };
  PDF_Abs* operator[](int idx);

 private:
  void makeAddDelCombinations();
  void checkAsimovArg() const;
  void checkColorArg() const;
  void checkCombinationArg() const;
  void configureAsimovCombinerNames(Combiner* c, const int i);
  bool combinerExists(const int id) const;
  void compareCombinations();
  void customizeCombinerTitles();
  void defineColors();
  void disableSystematics();
  void fixParameters(Combiner* c, const int cId);
  TString getStartParFileName(const int cId) const;
  bool isScanVarObservable(Combiner* c, const TString scanVar) const;
  void loadStartParameters(MethodProbScan* s, ParameterCache* pCache, const int cId);
  void make1dPluginOnlyPlot(std::shared_ptr<MethodPluginScan> sPlugin, const int cId);
  void make1dPluginPlot(std::shared_ptr<MethodPluginScan> sPlugin, std::shared_ptr<MethodProbScan> sProb,
                        const int cId);
  void make1dPluginScan(MethodPluginScan* scannerPlugin, const int cId);
  void make1dProbPlot(std::shared_ptr<MethodProbScan> scanner, const int cId);
  void make1dProbScan(MethodProbScan* scanner, const int cId);
  void make1dCoverageScan(MethodCoverageScan* scanner, const int cId);
  void make1dCoveragePlot(MethodCoverageScan* scanner, const int cId);
  void make1dBergerBoosScan(MethodBergerBoosScan* scanner, const int cId);
  void make2dPluginOnlyPlot(std::shared_ptr<MethodPluginScan> sPlugin, const int cId);
  void make2dPluginPlot(std::shared_ptr<MethodPluginScan> sPlugin, std::shared_ptr<MethodProbScan> sProb,
                        const int cId);
  void make2dPluginScan(std::shared_ptr<MethodPluginScan> scannerPlugin, const int cId);
  void make2dProbPlot(std::shared_ptr<MethodProbScan> scanner, const int cId);
  void make2dProbScan(MethodProbScan* scanner, const int cId);

  void printCombinerStructure(Combiner* c) const;
  void printBanner() const;
  bool pdfExists(const int id) const;
  void savePlot() const;
  void scaleStatErrors();
  void scaleStatAndSystErrors();
  void scaleDownErrors();  // now defunct
  void scan();
  void scanDataSet();
  void setAsimovObservables(Combiner* c);
  void setObservablesFromFile(Combiner* c, const int cId);
  void loadAsimovPoint(Combiner* c, const int cId);
  void setUpPlot();
  void tightenChi2Constraint(Combiner* c, const TString scanVar);
  void usage() const;
  void writebatchscripts();
  void makeLatex(Combiner* c) const;
  void saveWorkspace(Combiner* c, const int i);
  void runToys(Combiner* c);

  std::unique_ptr<OptParser> arg;
  std::vector<std::unique_ptr<Combiner>> cmb;
  std::vector<int> colorsLine;
  std::vector<int> colorsText;
  std::vector<int> fillStyles;
  std::vector<int> fillColors;
  std::vector<double> fillTransparencies;
  std::vector<int> lineColors;
  std::vector<int> lineStyles;
  std::vector<int> lineWidths;
  std::vector<std::shared_ptr<MethodProbScan>> comparisonScanners;
  TString execname;
  std::unique_ptr<FileNameBuilder> m_fnamebuilder;
  std::unique_ptr<BatchScriptWriter> m_batchscriptwriter;
  std::vector<PDF_Abs*> pdf;
  std::unique_ptr<OneMinusClPlotAbs> plot;
  TStopwatch t;
  std::unique_ptr<TApplication> theApp;
  bool runOnDataSet = false;
};

#endif
