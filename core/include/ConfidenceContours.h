/**
 * Gamma Combination
 * Author: Till Moritz Karbach, moritz.karbach@cern.ch
 * Date: Feb 2015
 *
 * Class holding the 1-N sigma confidence contours.
 *
 **/

#ifndef ConfidenceContours_h
#define ConfidenceContours_h

#include "Contour.h"
#include "OptParser.h"
#include "Utils.h"

class ConfidenceContours {
 public:
  ConfidenceContours(const OptParser* arg);
  ~ConfidenceContours();
  void computeContours(TH2F* hist, Utils::histogramType type, int id = 0);
  void Draw();
  void DrawDashedLine();
  void setStyle(std::vector<int>& linecolor, std::vector<int>& linestyle, std::vector<int>& linewidth,
                std::vector<int>& fillcolor, std::vector<int>& fillstyle);
  inline void setTransparency(float percent) { m_transparency = percent; };
  inline void setContoursToPlot(std::vector<int>& contstoplot) { m_contstoplots = contstoplot; };

 private:
  TH2F* addBoundaryBins(TH2F* hist);
  void addFilledPlotArea(TH2F* hist);
  TH2F* transformChi2valleyToHill(TH2F* hist, float offset);
  const OptParser* m_arg;            ///< command line arguments
  std::vector<Contour*> m_contours;  ///< container for the 1,...,N sigma contours
  std::vector<int> m_linecolor;      ///< style for the 1,...,N sigma contours
  std::vector<int> m_linestyle;
  std::vector<int> m_fillcolor;
  std::vector<int> m_fillstyle;
  std::vector<int> m_linewidth;
  float m_transparency;
  std::vector<int> m_contstoplots;  ///< container for which contours to actually draw
  int m_nMaxContours;
};

#endif
