#ifndef Contour_h
#define Contour_h

#include <vector>

#include <TH2.h>

#include "OptParser.h"

/**
 * Class holding a 2D contour consisting of several TGraphs.
 */
class Contour {
 public:
  Contour(const OptParser* arg, TList* listOfGraphs);
  ~Contour();
  void Draw() const;
  void DrawFilled() const;
  void DrawLine() const;
  inline int getSigma() const { return m_sigma; };
  void magneticBoundaries(const TH2* hCL);
  inline void setSigma(int s) { m_sigma = s; };
  void setStyle(const int linecolor, const int linestyle, const int linewidth, const int fillcolor,
                const int fillstyle);
  void setTransparency(const double percent);

 private:
  TGraph* changePointOrder(TGraph* g, int pointId);
  void findClosestPoints(TGraph* g1, TGraph* g2, int& i1, int& i2);
  TGraph* joinIfInside(TGraph* g1, TGraph* g2);
  std::vector<TGraph*> makeHoles(std::vector<TGraph*>& contours);
  void magneticBoundaries(std::vector<TGraph*>& contours, const TH2* hCL);

  const OptParser* m_arg;           ///< Command line arguments
  std::vector<TGraph*> m_contours;  ///< Container for the several disjoint subcontours. Used by DrawLine().

  /// Container for contours with holes (filled by makeHoles(), used by DrawFilled()).
  std::vector<TGraph*> m_contoursHoles;
  int m_sigma = -1;
  int m_linecolor = 2;
  int m_linestyle = kSolid;
  int m_fillcolor = 2;
  int m_fillstyle = 1001;
  int m_linewidth = 1;
  double m_alpha = 1.;
};

#endif
