/**
 * Gamma Combination
 * Author: Till Moritz Karbach, moritz.karbach@cern.ch
 * Date: Feb 2015
 *
 * Class holding a 2D contour consisting of several TGraphs.
 *
 **/

#ifndef Contour_h
#define Contour_h

#include <vector>

#include "OptParser.h"

class Contour {
 public:
  Contour(const OptParser* arg, TList* listOfGraphs);
  ~Contour();
  void Draw() const;
  void DrawFilled() const;
  void DrawLine() const;
  inline int getSigma() const { return m_sigma; };
  void magneticBoundaries(const TH2F* hCL);
  inline void setSigma(int s) { m_sigma = s; };
  void setStyle(int linecolor, int linestyle, int linewidth, int fillcolor, int fillstyle);
  void setTransparency(double percent);

 private:
  TGraph* changePointOrder(TGraph* g, int pointId);
  void findClosestPoints(TGraph* g1, TGraph* g2, int& i1, int& i2);
  TGraph* joinIfInside(TGraph* g1, TGraph* g2);
  std::vector<TGraph*> makeHoles(std::vector<TGraph*>& contours);
  void magneticBoundaries(std::vector<TGraph*>& contours, const TH2F* hCL);

  const OptParser* m_arg;                ///< command line arguments
  std::vector<TGraph*> m_contours;       ///< container for the several disjoint subcontours. Used by DrawLine().
  std::vector<TGraph*> m_contoursHoles;  ///< container for contours with holes. Filled by makeHoles(). Used by
                                         ///< DrawFilled().
  int m_sigma;                           ///< sigma level of the contour
  int m_linecolor;                       ///< style for the contour
  int m_linestyle;
  int m_fillcolor;
  int m_fillstyle;
  int m_linewidth;
  double m_alpha;
};

#endif
