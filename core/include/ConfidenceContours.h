#ifndef ConfidenceContours_h
#define ConfidenceContours_h

#include "Contour.h"
#include "OptParser.h"
#include "Utils.h"

#include <memory>
#include <vector>

/**
 * Class holding the confidence contours from 1 to N sigmas.
 **/
class ConfidenceContours {
 public:
  ConfidenceContours(const OptParser* arg);
  void computeContours(const TH2* hist, const Utils::histogramType type, const int id = 0);
  void Draw();
  void DrawDashedLine();
  void setStyle(const std::vector<int>& linecolor, const std::vector<int>& linestyle, const std::vector<int>& linewidth,
                const std::vector<int>& fillcolor, const std::vector<int>& fillstyle);
  inline void setTransparency(double percent) { m_transparency = percent; };
  inline void setContoursToPlot(const std::vector<int>& contstoplot) { m_contstoplots = contstoplot; };

 private:
  std::unique_ptr<TH2> addBoundaryBins(const TH2* hist);
  void addFilledPlotArea(const TH2* hist);
  std::unique_ptr<TH2> transformChi2valleyToHill(const TH2* hist, const double offset);
  const OptParser* m_arg;                            ///< command line arguments
  std::vector<std::unique_ptr<Contour>> m_contours;  ///< container for the 1,...,N sigma contours
  std::vector<int> m_linecolor;                      ///< style for the 1,...,N sigma contours
  std::vector<int> m_linestyle;
  std::vector<int> m_fillcolor;
  std::vector<int> m_fillstyle;
  std::vector<int> m_linewidth;
  std::vector<int> m_contstoplots;  ///< container for which contours to actually draw
  double m_transparency = 0.;
  const int m_nMaxContours = 9;
};

#endif
