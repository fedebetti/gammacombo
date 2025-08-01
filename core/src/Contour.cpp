#include <Contour.h>
#include <OptParser.h>
#include <Utils.h>

#include <TGraph.h>
#include <TH2.h>
#include <TList.h>

#include <limits>
#include <memory>
#include <utility>
#include <vector>

namespace {

  /// Find the two closest points between two TGraphs
  std::pair<int, int> findClosestPoints(const TGraph* g1, const TGraph* g2) {
    double x1, y1, x2, y2;
    int i1 = -1;
    int i2 = -1;
    auto distance = std::numeric_limits<double>::max();
    for (int ii1 = 0; ii1 < g1->GetN(); ii1++) {
      for (int ii2 = 0; ii2 < g2->GetN(); ii2++) {
        g1->GetPoint(ii1, x1, y1);
        g2->GetPoint(ii2, x2, y2);
        const double d = std::sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
        if (d < distance) {
          i1 = ii1;
          i2 = ii2;
          distance = d;
        }
      }
    }
    return {i1, i2};
  }

  /// Create a new TGraph where the points start from \p pointId.
  std::unique_ptr<TGraph> changePointOrder(const TGraph* g, const int pointId) {
    double pointx, pointy;
    auto gNew = std::make_unique<TGraph>(g->GetN() + 1);
    for (int i = pointId; i < g->GetN() + pointId; i++) {
      g->GetPoint(i < g->GetN() ? i : i - g->GetN(), pointx, pointy);
      gNew->SetPoint(i - pointId, pointx, pointy);
    }
    gNew->GetPoint(0, pointx, pointy);
    gNew->SetPoint(gNew->GetN() - 1, pointx, pointy);
    return gNew;
  }

  /// Join two graphs if one of the two is completely contained in the other.
  std::unique_ptr<TGraph> joinIfInside(const TGraph* g1, const TGraph* g2) {
    // First determine which graph lies inside which, if they include each other at all.
    // Get number of points of g1 that lie inside g2
    double pointx, pointy;
    int nG1InsideG2 = 0;
    for (int i = 0; i < g1->GetN(); i++) {
      g1->GetPoint(i, pointx, pointy);
      if (g2->IsInside(pointx, pointy)) nG1InsideG2++;
    }
    // Get number of points of g2 that lie inside g1
    int nG2InsideG1 = 0;
    for (int i = 0; i < g2->GetN(); i++) {
      g2->GetPoint(i, pointx, pointy);
      if (g1->IsInside(pointx, pointy)) nG2InsideG1++;
    }
    // They don't contain each other
    if (nG1InsideG2 == 0 && nG2InsideG1 == 0) return nullptr;
    // They do contain each other: merge them into g1 so that the line will form a hole!
    // We'll merge at the points that are closest.
    auto [i1, i2] = findClosestPoints(g1, g2);
    // Change graph order such that it starts and ends with the nearest point
    auto g1Reordered = changePointOrder(g1, i1);
    auto g2Reordered = changePointOrder(g2, i2);
    // Merge them
    auto gNew = std::make_unique<TGraph>(g1Reordered->GetN() + g2Reordered->GetN());
    for (int i = 0; i < g1Reordered->GetN(); i++) {
      g1Reordered->GetPoint(i, pointx, pointy);
      gNew->SetPoint(i, pointx, pointy);
    }
    for (int i = 0; i < g2Reordered->GetN(); i++) {
      g2Reordered->GetPoint(i, pointx, pointy);
      gNew->SetPoint(i + g1Reordered->GetN(), pointx, pointy);
    }
    return gNew;
  }

  std::vector<std::unique_ptr<TGraph>> makeHolesHelper(std::vector<std::unique_ptr<TGraph>>& contours) {
    const int n = contours.size();
    int iJoined1;
    int iJoined2;
    std::unique_ptr<TGraph> gJoined;
    for (int i1 = 0; i1 < n; i1++) {
      if (gJoined) break;
      const TGraph* g1 = contours[i1].get();
      for (int i2 = i1 + 1; i2 < n; i2++) {
        const TGraph* g2 = contours[i2].get();
        gJoined = joinIfInside(g1, g2);
        if (gJoined) {
          iJoined1 = i1;
          iJoined2 = i2;
          break;
        }
      }
    }

    if (gJoined) {
      std::vector<std::unique_ptr<TGraph>> newContours;
      contours.erase(contours.begin() + iJoined2);
      contours.erase(contours.begin() + iJoined1);
      contours.push_back(std::move(gJoined));
      return makeHolesHelper(contours);
    }
    return std::move(contours);
  }
}  // namespace

/**
 * Constructor.
 *
 * @param arg          Command line options.
 * @param listOfGraphs List of TGraphs that make up the subcontours.
 */
Contour::Contour(const OptParser* arg, const TList* listOfGraphs) {
  assert(arg);
  m_arg = arg;
  for (auto g : *listOfGraphs) { m_contours.emplace_back(Utils::clone<TGraph>(g)); }

  // compute holes in the contours
  m_contoursHoles = makeHoles(m_contours);
}

///
/// Draw the contours into the currently active Canvas.
///
void Contour::Draw() {
  DrawFilled();
  DrawLine();
}

///
/// Draw filled contours into the currently active Canvas.
/// This plots the contours in m_contoursHoles.
///
void Contour::DrawFilled() {
  for (auto&& ch : m_contoursHoles) {
    m_tmpGraphs.emplace_back(Utils::clone<TGraph>(ch.get()));
    auto g = m_tmpGraphs.back().get();
    g->SetFillStyle(1001);                       // solid
    g->SetFillColorAlpha(m_fillcolor, m_alpha);  // transparency!
    if (!m_arg->isQuickhack(27) && m_fillstyle != 0) g->Draw("F");
    if (m_fillstyle != 1001) {  // if not solid, add the pattern in the line color
      m_tmpGraphs.emplace_back(Utils::clone<TGraph>(ch.get()));
      auto g = m_tmpGraphs.back().get();
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
  for (auto&& ch : m_contours) {
    ch->SetLineWidth(m_linewidth);
    ch->SetLineColor(m_linecolor);
    ch->SetLineStyle(m_linestyle);
    ch->SetFillStyle(0);  // hollow
    ch->Draw("L");
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

/**
 * Make contours that have holes cut away.
 *
 * If a, say, 1sigma contour looks like a ring, it is initially stored as two independent circles in the TList.
 * Here we join them, if one lies inside the other, such that when the resulting contour is plotted, the plot will
 * actually have a hole.
 *
 * Fills m_contoursHoles. The reason this is kept separate from the m_contours container, which holds the contours
 * without holes, is the plotting: one cannot plot the just-lines version from the contours with holes, else one sees
 * where the contour is artificially closed to be able to fill the inside.
 *
 * @param contours Vector containing contours without holes.
 *
 * @return         Vector with contours with holes.
 */
std::vector<std::unique_ptr<TGraph>> Contour::makeHoles(const std::vector<std::unique_ptr<TGraph>>& contours) {
  std::vector<std::unique_ptr<TGraph>> newContours;
  for (auto&& g : contours) { newContours.emplace_back(Utils::clone<TGraph>(g.get())); }
  return makeHolesHelper(newContours);
  // const int n = contours.size();
  // int iJoined1;
  // int iJoined2;
  // std::unique_ptr<TGraph> gJoined;
  // for (int i1 = 0; i1 < n; i1++) {
  //   if (gJoined) break;
  //   const TGraph* g1 = contours[i1].get();
  //   for (int i2 = i1 + 1; i2 < n; i2++) {
  //     const TGraph* g2 = contours[i2].get();
  //     gJoined = joinIfInside(g1, g2);
  //     if (gJoined) {
  //       iJoined1 = i1;
  //       iJoined2 = i2;
  //       break;
  //     }
  //   }
  // }

  // if (gJoined) {
  //   std::vector<std::unique_ptr<TGraph>> newContours;
  //   newContours.push_back(gJoined);
  //   for (int i = 0; i < n; i++) {
  //     if (i != iJoined1 && i != iJoined2) newContours.push_back(contours[i]);
  //   }
  //   return makeHoles(newContours);
  // }
  // return contours;
}

/**
 * Magnetic boundaries. If a contour is closer than half a bin width to a boundary, adjust it to be the boundary.
 *
 * @param contour Input contours.
 * @param hCL     Histogram defining the boundaries.
 */
void Contour::magneticBoundaries(const std::vector<std::unique_ptr<TGraph>>& contours, const TH2* hCL) {
  const double magneticRange = 0.75;
  const double xmin = hCL->GetXaxis()->GetXmin();
  const double xmax = hCL->GetXaxis()->GetXmax();
  const double ymin = hCL->GetYaxis()->GetXmin();
  const double ymax = hCL->GetYaxis()->GetXmax();
  const double xbinwidth = hCL->GetXaxis()->GetBinWidth(1);
  const double ybinwidth = hCL->GetYaxis()->GetBinWidth(1);
  double pointx, pointy;
  for (auto&& g : contours) {
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
    std::cout << "Contour::setTransparency() : ERROR : percent not in [0,1]. Skipping." << std::endl;
    return;
  }
  m_alpha = 1. - percent;
}
