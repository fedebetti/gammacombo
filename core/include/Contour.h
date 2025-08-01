#ifndef Contour_h
#define Contour_h

#include <memory>
#include <utility>
#include <vector>

#include <TH2.h>
#include <TList.h>

#include "OptParser.h"

/**
 * Class holding a 2D contour consisting of several TGraphs.
 */
class Contour {
 public:
  Contour(const OptParser* arg, const TList* listOfGraphs);
  void Draw();
  void DrawFilled();
  void DrawLine() const;
  inline int getSigma() const { return m_sigma; };
  void magneticBoundaries(const TH2* hCL);
  inline void setSigma(int s) { m_sigma = s; };
  void setStyle(const int linecolor, const int linestyle, const int linewidth, const int fillcolor,
                const int fillstyle);
  void setTransparency(const double percent);

 private:
  std::vector<std::unique_ptr<TGraph>> makeHoles(const std::vector<std::unique_ptr<TGraph>>& contours);
  void magneticBoundaries(const std::vector<std::unique_ptr<TGraph>>& contours, const TH2* hCL);

  const OptParser* m_arg;  ///< Command line arguments

  /// Vector of disjoint subcontours.
  std::vector<std::unique_ptr<TGraph>> m_contours;

  /// Vector of contours with holes.
  std::vector<std::unique_ptr<TGraph>> m_contoursHoles;

  std::vector<std::unique_ptr<TGraph>> m_tmpGraphs;

  int m_sigma = -1;
  int m_linecolor = 2;
  int m_linestyle = kSolid;
  int m_fillcolor = 2;
  int m_fillstyle = 1001;
  int m_linewidth = 1;
  double m_alpha = 1.;
};

#endif
