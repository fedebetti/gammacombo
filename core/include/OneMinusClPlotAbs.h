/**
 * Gamma Combination
 * Author: Till Moritz Karbach, moritz.karbach@cern.ch
 * Date: August 2012
 *
 **/

#ifndef OneMinusClPlotAbs_h
#define OneMinusClPlotAbs_h

#include <TCanvas.h>
#include <TString.h>

#include "MethodAbsScan.h"
#include "OptParser.h"

class MethodAbsScan;

class OneMinusClPlotAbs {
 public:
  OneMinusClPlotAbs(OptParser* arg, TString name = "c1", TString title = "c1");
  ~OneMinusClPlotAbs();

  virtual void addScanner(MethodAbsScan* s, int CLsType = 0);
  inline void disableLegend(bool yesNo = false) { plotLegend = yesNo; };
  inline void disableSolution(bool yesNo = false) { plotSolution = yesNo; };
  virtual void drawSolutions() const;
  virtual void drawLabel([[maybe_unused]] float yPos = 0.6) const { std::cout << "nothing yet" << std::endl; }
  virtual void drawGroup(float yPos = 0.6) const;
  inline TString getName() const { return name; };
  void save();
  void setYLogRange(double min = 1.e-3, double max = 1) {
    plotLogYMin = min;
    plotLogYMax = max;
  };
  inline void setFont(int fnum) { font = fnum; };
  inline void setLabelSize(int lnum) { labelsize = lnum; };
  inline void setPlotLabel(TString& lname) { label = lname; };
  inline void Show() const { m_mainCanvas->Show(); };
  virtual void Draw() = 0;

  int font = 133;       ///< font code. The last digit disables scaling with the canvas size.
  int labelsize = 35;   ///< text size of axis labels, numeric solutions, CL guide lines (in pixels)
  int titlesize = 45;   ///< text size of axis titles, group label, "Prliminary" is x0.75 (in pixels)
  int legendsize = 29;  ///< text size of legend entries in 1d and 2d plots (in pixels)

  std::vector<MethodAbsScan*> scanners;
  std::vector<int> do_CLs;   ///< vector, which stores the cls method type to be plotted
  OptParser* arg = nullptr;  ///< command line options
  TCanvas* m_mainCanvas = nullptr;
  TString name;
  TString title;
  TString label;
  bool plotLegend = true;
  bool plotSolution = true;
  double plotLogYMin = 1e-3;
  double plotLogYMax = 1;
};

#endif
