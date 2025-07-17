#ifndef OneMinusClPlot_h
#define OneMinusClPlot_h

#include "OneMinusClPlotAbs.h"

#include <TGraph.h>
#include <TH1.h>

#include <memory>

/**
 * Class to plot the distribution of 1 - CL for fit variables.
 */
class OneMinusClPlot : public OneMinusClPlotAbs {
 public:
  OneMinusClPlot(OptParser* arg, TString name = "c1", TString title = "c1");

  void drawSolutions() override;
  void drawCLguideLines();
  TGraph* getGraph(MethodAbsScan* s, const bool first = true, const bool last = false, const bool filled = true,
                   const int CLsType = 0) {
    return scan1dPlot(s, first, last, filled, CLsType);
  };
  inline void setPluginMarkers(bool yesNo = true) { plotPluginMarkers = yesNo; }
  void Draw(const bool beautify) override;

 private:
  void drawCLguideLine(const double pvalue);
  void drawVerticalLine(const double x, const int color, const int style);
  std::unique_ptr<TH1> getHistogram(MethodAbsScan* s, const int CLsType, const bool removeErrs) const;
  TGraph* scan1dPlot(MethodAbsScan* s, const bool first, const bool last, const bool filled, const int CLsType = 0);
  void scan1dPlotSimple(MethodAbsScan* s, const bool first, const int CLsType = 0);
  void scan1dCLsPlot(MethodAbsScan* s, const bool smooth = true, const bool obsError = true);

  bool plotPluginMarkers = true;
};

#endif
