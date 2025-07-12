#ifndef OneMinusClPlot_h
#define OneMinusClPlot_h

#include "OneMinusClPlotAbs.h"

#include <TGraph.h>

/**
 * Class to plot the distribution of 1 - CL for fit variables.
 */
class OneMinusClPlot : public OneMinusClPlotAbs {
 public:
  OneMinusClPlot(OptParser* arg, TString name = "c1", TString title = "c1");

  void drawSolutions() override;
  void drawCLguideLines() const;
  std::unique_ptr<TGraph> getGraph(MethodAbsScan* s, bool first = true, bool last = false, bool filled = true,
                                   int CLsType = 0) {
    return scan1dPlot(s, first, last, filled, CLsType);
  };
  inline void setPluginMarkers(bool yesNo = true) { plotPluginMarkers = yesNo; }
  void Draw() override;

 private:
  void drawCLguideLine(double pvalue) const;
  void drawVerticalLine(double x, int color, int style) const;
  std::unique_ptr<TGraph> scan1dPlot(MethodAbsScan* s, bool first, bool last, bool filled, int CLsType = 0);
  void scan1dPlotSimple(MethodAbsScan* s, bool first, int CLsType = 0);
  void scan1dCLsPlot(MethodAbsScan* s, bool smooth = true, bool obsError = true);

  bool plotPluginMarkers = true;
};

#endif
