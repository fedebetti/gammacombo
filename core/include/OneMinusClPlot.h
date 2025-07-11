/**
 * Gamma Combination
 * Author: Till Moritz Karbach, moritz.karbach@cern.ch
 * Date: August 2012
 *
 **/

#ifndef OneMinusClPlot_h
#define OneMinusClPlot_h

#include "OneMinusClPlotAbs.h"

#include <TGraph.h>

class OneMinusClPlot : public OneMinusClPlotAbs {
 public:
  OneMinusClPlot(OptParser* arg, TString name = "c1", TString title = "c1");

  void drawSolutions();
  void drawCLguideLines() const;
  TGraph* getGraph(MethodAbsScan* s, bool first = true, bool last = false, bool filled = true, int CLsType = 0) {
    return scan1dPlot(s, first, last, filled, CLsType);
  };
  inline TString getName() const { return name; }
  inline void setPluginMarkers(bool yesNo = true) { plotPluginMarkers = yesNo; }
  void Draw() override;

 private:
  void drawCLguideLine(double pvalue) const;
  void drawVerticalLine(double x, int color, int style) const;
  TGraph* scan1dPlot(MethodAbsScan* s, bool first, bool last, bool filled, int CLsType = 0);
  void scan1dPlotSimple(MethodAbsScan* s, bool first, int CLsType = 0);
  void scan1dCLsPlot(MethodAbsScan* s, bool smooth = true, bool obsError = true);

  bool plotPluginMarkers = true;
  TCanvas* m_clsCanvas = nullptr;
};

#endif
