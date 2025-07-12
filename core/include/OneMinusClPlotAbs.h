#ifndef OneMinusClPlotAbs_h
#define OneMinusClPlotAbs_h

#include <TCanvas.h>
#include <TString.h>

#include "MethodAbsScan.h"
#include "OptParser.h"

class MethodAbsScan;

/**
 * Abstract class to plot the CL distributions for fit variables.
 */
class OneMinusClPlotAbs {
 public:
  OneMinusClPlotAbs(OptParser* arg, TString name = "c1", TString title = "c1");

  virtual void addScanner(MethodAbsScan* s, int CLsType = 0);
  inline void disableLegend(bool yesNo = false) { plotLegend = yesNo; };
  virtual void drawSolutions() = 0;
  virtual void drawGroup(double yPos = 0.6) const;
  inline TString getName() const { return name; };
  void save() const;
  inline void setFont(int fnum) { font = fnum; };
  inline void setLabelSize(int lnum) { labelsize = lnum; };
  inline void setPlotLabel(TString& lname) { label = lname; };
  inline void Show() const { m_mainCanvas->Show(); };
  virtual void Draw() = 0;

 protected:
  int font = 133;       ///< font code. The last digit disables scaling with the canvas size.
  int labelsize = 35;   ///< text size of axis labels, numeric solutions, CL guide lines (in pixels)
  int titlesize = 45;   ///< text size of axis titles, group label, "Prliminary" is x0.75 (in pixels)
  int legendsize = 29;  ///< text size of legend entries in 1d and 2d plots (in pixels)

  std::vector<MethodAbsScan*> scanners;
  std::vector<int> do_CLs;   ///< vector, which stores the cls method type to be plotted
  OptParser* arg = nullptr;  ///< command line options
  std::unique_ptr<TCanvas> m_mainCanvas;
  TString name;
  TString title;
  TString label;
  bool plotLegend = true;
};

#endif
