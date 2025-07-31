#include <Contour.h>
#include <Utils.h>

using namespace std;

///
/// Constructor. The class stores copies of the TGraphs provided in listOfGraphs.
///
/// \param arg - command line options
/// \param listOfGraphs - list of TGraphs that make up the subcontours
///
Contour::Contour(const OptParser* arg, TList* listOfGraphs) {
  assert(arg);
  m_arg = arg;
  for (auto g : *listOfGraphs) { m_contours.push_back(static_cast<TGraph*>(g->Clone())); }

  // compute holes in the contours
  m_contoursHoles = makeHoles(m_contours);
}

Contour::~Contour() {
  for (auto contour : m_contours) { delete contour; }
  for (auto ch : m_contoursHoles) { delete ch; }
}

///
/// Draw the contours into the currently active Canvas.
///
void Contour::Draw() const {
  DrawFilled();
  DrawLine();
}

///
/// Draw filled contours into the currently active Canvas.
/// This plots the contours in m_contoursHoles.
///
void Contour::DrawFilled() const {
  for (auto ch : m_contoursHoles) {
    auto g = dynamic_cast<TGraph*>(ch->Clone());
    g->SetFillStyle(1001);  // solid
    // g->SetFillColor(m_fillcolor);
    g->SetFillColorAlpha(m_fillcolor, m_alpha);  // transparency!
    if (!m_arg->isQuickhack(27) && m_fillstyle != 0) g->Draw("F");
    if (m_fillstyle != 1001) {  // if not solid, add the pattern in the line color
      g = (TGraph*)ch->Clone();
      g->SetFillStyle(m_fillstyle);  // hatched on top
      g->SetFillColor(m_linecolor);
      g->Draw("F");
    }
  }
}

///
/// Draw line contours into the currently active Canvas. This plots the contours in m_contours.
///
void Contour::DrawLine() const {
  for (auto ch : m_contours) {
    auto g = dynamic_cast<TGraph*>(ch->Clone());
    g->SetLineWidth(m_linewidth);
    g->SetLineColor(m_linecolor);
    g->SetLineStyle(m_linestyle);
    g->SetFillStyle(0);  // hollow
    g->Draw("L");
  }
}

///
/// Set the contour style.
///
void Contour::setStyle(const int linecolor, const int linestyle, const int linewidth, const int fillcolor,
                       const int fillstyle) {
  m_linecolor = linecolor;
  m_linestyle = linestyle;
  m_fillcolor = fillcolor;
  m_fillstyle = fillstyle;
  m_linewidth = linewidth;
}

///
/// Make contours that have holes cut away.
/// If a, say, 1sigma contour looks like a ring, it is initially
/// stored as two independent circles in the TList. Here we join
/// them, if one lies inside the other, such that when the resulting
/// contour is plotted, the plot will actually have a hole.
///
/// Fills m_contoursHoles. The reason this is kept separate from
/// the m_contours container, that holds the contours without holes,
/// is the plotting: one can't plot the just-lines version from the
/// contours with holes, else one sees where the contour is artificially
/// closed to be able to fill the inside.
///
/// \param contours - vector containing contours without holes
/// \return - vector with contours with holes
///
vector<TGraph*> Contour::makeHoles(vector<TGraph*>& contours) {
  int n = contours.size();
  bool joined = false;
  int iJoined1;
  int iJoined2;
  TGraph* gJoined = nullptr;
  for (int i1 = 0; i1 < n; i1++) {
    if (joined) break;
    TGraph* g1 = contours[i1];
    for (int i2 = i1 + 1; i2 < n; i2++) {
      TGraph* g2 = contours[i2];
      gJoined = joinIfInside(g1, g2);
      if (gJoined) {
        joined = true;
        iJoined1 = i1;
        iJoined2 = i2;
        break;
      }
    }
  }

  if (joined) {
    vector<TGraph*> newContours;
    newContours.push_back(gJoined);
    for (int i = 0; i < n; i++) {
      if (i != iJoined1 && i != iJoined2) newContours.push_back(contours[i]);
    }
    return makeHoles(newContours);
  }
  return contours;
}

///
/// Helper function for makeHoles().
///
TGraph* Contour::joinIfInside(TGraph* g1, TGraph* g2) {
  // First determine which graph lies inside which, if they include each
  // other at all.
  // Get number of points of g1 that lie inside g2:
  Double_t pointx, pointy;
  int nG1InsideG2 = 0;
  for (int i = 0; i < g1->GetN(); i++) {
    g1->GetPoint(i, pointx, pointy);
    if (g2->IsInside(pointx, pointy)) nG1InsideG2++;
  }
  // reversed: Get number of points of g2 that lie inside g1:
  int nG2InsideG1 = 0;
  for (int i = 0; i < g2->GetN(); i++) {
    g2->GetPoint(i, pointx, pointy);
    if (g1->IsInside(pointx, pointy)) nG2InsideG1++;
  }
  // they don't contain each other
  if (nG1InsideG2 == 0 && nG2InsideG1 == 0) return nullptr;
  // they do contain each other: merge them into g1 so that the line
  // will form a hole! We'll merge at the points that are closest.
  int i1, i2;
  findClosestPoints(g1, g2, i1, i2);
  // cout << "g1 ===========" << endl;
  // g1->Print();
  // cout << "g2 ===========" << endl;
  // g2->Print();
  // change graph order such that it stars and ends with the nearest point
  g1 = changePointOrder(g1, i1);
  g2 = changePointOrder(g2, i2);
  // cout << "g1 ===========" << endl;
  // g1->Print();
  // cout << "g2 ===========" << endl;
  // g2->Print();
  // merge them
  auto gNew = new TGraph(g1->GetN() + g2->GetN());
  for (int i = 0; i < g1->GetN(); i++) {
    g1->GetPoint(i, pointx, pointy);
    gNew->SetPoint(i, pointx, pointy);
  }
  for (int i = 0; i < g2->GetN(); i++) {
    g2->GetPoint(i, pointx, pointy);
    gNew->SetPoint(i + g1->GetN(), pointx, pointy);
  }
  // cout << "gNew ===========" << endl;
  // gNew->Print();
  return gNew;
}

///
/// Helper function for joinIfInside().
///
TGraph* Contour::changePointOrder(TGraph* g, int pointId) {
  Double_t pointx, pointy;
  auto gNew = new TGraph(g->GetN() + 1);
  for (int i = pointId; i < g->GetN() + pointId; i++) {
    g->GetPoint(i < g->GetN() ? i : i - g->GetN(), pointx, pointy);
    gNew->SetPoint(i - pointId, pointx, pointy);
  }
  gNew->GetPoint(0, pointx, pointy);
  gNew->SetPoint(gNew->GetN() - 1, pointx, pointy);
  return gNew;
}

///
/// Helper function for makeHoles().
///
void Contour::findClosestPoints(TGraph* g1, TGraph* g2, int& i1, int& i2) {
  Double_t x1, y1, x2, y2;
  double distance = 1e6;
  for (int ii1 = 0; ii1 < g1->GetN(); ii1++) {
    for (int ii2 = 0; ii2 < g2->GetN(); ii2++) {
      g1->GetPoint(ii1, x1, y1);
      g2->GetPoint(ii2, x2, y2);
      double d = sqrt(Utils::sq(x1 - x2) + Utils::sq(y1 - y2));
      if (d < distance) {
        i1 = ii1;
        i2 = ii2;
        distance = d;
      }
    }
  }
  g1->GetPoint(i1, x1, y1);
  g2->GetPoint(i2, x2, y2);
}

/**
 * Magnetic boundaries. If a contour is closer than half a bin width to a boundary, adjust it to be the boundary.
 *
 * @param contour Input contours.
 * @param hCL     Histogram defining the boundaries.
 */
void Contour::magneticBoundaries(vector<TGraph*>& contours, const TH2* hCL) {
  const double magneticRange = 0.75;
  const double xmin = hCL->GetXaxis()->GetXmin();
  const double xmax = hCL->GetXaxis()->GetXmax();
  const double ymin = hCL->GetYaxis()->GetXmin();
  const double ymax = hCL->GetYaxis()->GetXmax();
  const double xbinwidth = hCL->GetXaxis()->GetBinWidth(1);
  const double ybinwidth = hCL->GetYaxis()->GetBinWidth(1);
  double pointx, pointy;
  for (auto contour : contours) {
    auto g = dynamic_cast<TGraph*>(contour);
    for (int i = 0; i < g->GetN(); i++) {
      g->GetPoint(i, pointx, pointy);
      if (abs(pointx - xmin) < xbinwidth * magneticRange) g->SetPoint(i, xmin, pointy);
      g->GetPoint(i, pointx, pointy);
      if (abs(pointx - xmax) < xbinwidth * magneticRange) g->SetPoint(i, xmax, pointy);
      g->GetPoint(i, pointx, pointy);
      if (abs(pointy - ymin) < ybinwidth * magneticRange) g->SetPoint(i, pointx, ymin);
      g->GetPoint(i, pointx, pointy);
      if (abs(pointy - ymax) < ybinwidth * magneticRange) g->SetPoint(i, pointx, ymax);
    }
  }
}

void Contour::magneticBoundaries(const TH2* hCL) {
  magneticBoundaries(m_contours, hCL);
  magneticBoundaries(m_contoursHoles, hCL);
}

/**
 * Set the transparency of the contour.
 *
 * \param percent 100% means fully transparent, 0% means intransparent.
 */
void Contour::setTransparency(const double percent) {
  if (!(0. <= percent && percent <= 1.)) {
    cout << "Contour::setTransparency() : ERROR : percent not in [0,1]. Skipping." << endl;
    return;
  }
  m_alpha = 1. - percent;
}
