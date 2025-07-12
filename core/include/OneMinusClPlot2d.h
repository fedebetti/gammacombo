#ifndef OneMinusClPlot2d_h
#define OneMinusClPlot2d_h

#include "ConfidenceContours.h"
#include "OneMinusClPlotAbs.h"
#include "Utils.h"

#include <TH2.h>
#include <TLegend.h>
#include <TString.h>

/**
 * Class to plot the 2D contours for the fit variables.
 */
class OneMinusClPlot2d : public OneMinusClPlotAbs {
 public:
  OneMinusClPlot2d(OptParser* arg, TString name = "c1", TString title = "c1");

  void addScanner(MethodAbsScan* s, int CLsType = 0);
  void addFile(TString fName);
  void Draw() override;
  void DrawFull();
  void drawCLcontent(bool isFull = false);
  void drawMarker(double x, double y, int color = 0, int style = 3, double size = 2.0) const;
  void drawGroup() const;
  void drawSolutions() override;
  inline int getNumberOfDefinedColors() const { return linecolor[0].size(); }
  inline void setContoursOnly() { contoursOnly = true; };
  inline void setXaxisTitle(TString s) { xTitle = s; };
  inline void setYaxisTitle(TString s) { yTitle = s; };

 protected:
  std::vector<TH2*> histos;
  TString xTitle;
  TString yTitle;
  bool contoursOnly = false;
  std::vector<std::vector<int>> linecolor;  ///< defines colors of 1 sigma lines and solutions of different scanners
  std::vector<std::vector<int>> fillcolor;  ///< defines colors of 1 sigma areas of different scanners
  std::vector<std::vector<int>> linestyle;  ///< defines the line style of 1 sigma line of different scanners
  std::vector<std::vector<int>> fillstyle;  ///< defines the fill style of
  std::vector<std::vector<int>> linewidth;  ///< defines the line width
  std::vector<std::vector<double>> filltransparency;  ///< defines the fill transparency
  std::vector<int> markerstyle;                       ///< defines marker styles of the solutions of different scanners
  std::vector<double> markersize;

 private:
  void drawLegend();
  bool hasHistoType(Utils::histogramType t) const;
  void makeNewPlotStyle(TString htmlColor, int ROOTColor = -1);
  void makeOneColorPlotStyle(TString htmlColor, int ROOTColor = -1);

  std::vector<Utils::histogramType> histosType;  ///< defines if histogram is interpreted as p-value or chi2
  std::vector<ConfidenceContours*> m_contours;   ///< holds the contours for each scanner
  std::vector<bool> m_contours_computed;  ///< true if the contours were computed for that scanner by computeContours()
  std::unique_ptr<TLegend> m_legend = nullptr;  ///< pointer to the plot legend. Filled by drawLegend().
};

#endif
