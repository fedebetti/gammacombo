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

#include <TAttLine.h>

class OptParser;

class TH2F;
class TGraph;
class TList;

class Contour {
 public:
  Contour(const OptParser* arg, TList* listOfGraphs);
  ~Contour();
  void Draw();
  void DrawFilled();
  void DrawLine();
  inline int getSigma() { return m_sigma; };
  void magneticBoundaries(const TH2F* hCL);
  inline void setSigma(int s) { m_sigma = s; };
  void setStyle(int linecolor, int linestyle, int linewidth, int fillcolor, int fillstyle);
  void setTransparency(float percent);

 private:
  TGraph* changePointOrder(TGraph* g, int pointId);
  void findClosestPoints(TGraph* g1, TGraph* g2, int& i1, int& i2);
  TGraph* joinIfInside(TGraph* g1, TGraph* g2);
  std::vector<TGraph*> makeHoles(std::vector<TGraph*>& contours);
  void magneticBoundaries(std::vector<TGraph*>& contours, const TH2F* hCL);

  const OptParser* m_arg = nullptr;  ///< command line arguments

  /// Vector of disjoint subcontours.
  std::vector<TGraph*> m_contours;

  /// Vector of contours with holes.
  std::vector<TGraph*> m_contoursHoles;

  int m_sigma = -1;
  int m_linecolor = 2;
  int m_linestyle = kSolid;
  int m_fillcolor = 2;
  int m_fillstyle = 1001;
  int m_linewidth = 1;
  float m_alpha = 1.;
};

#endif
