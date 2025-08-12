#include <ConfidenceContours.h>
#include <Utils.h>

#include <TH2.h>
#include <TH2F.h>
#include <TROOT.h>

#include <iomanip>
#include <memory>
#include <vector>

namespace {
  /**
   * Constructs a new 2D histogram that contains 2 more bins in each direction that are set to the minimum,
   * so that the contours will always close.
   *
   * @param hist 2D histogram.
   *
   * @return     New 2D histogram.
   */
  std::unique_ptr<TH2> addBoundaryBins(const TH2* hist) {
    using Utils::getUniqueRootName;
    const double boundary = hist->GetMinimum();
    auto hBoundaries =
        std::make_unique<TH2F>(getUniqueRootName(), getUniqueRootName(), hist->GetNbinsX() + 2,
                               hist->GetXaxis()->GetXmin() - hist->GetXaxis()->GetBinWidth(1),
                               hist->GetXaxis()->GetXmax() + hist->GetXaxis()->GetBinWidth(1), hist->GetNbinsY() + 2,
                               hist->GetYaxis()->GetXmin() - hist->GetYaxis()->GetBinWidth(1),
                               hist->GetYaxis()->GetXmax() + hist->GetYaxis()->GetBinWidth(1));
    /*
        std::make_unique<TH2F>(getUniqueRootName(), getUniqueRootName(), hist->GetNbinsX() + 4,
                 hist->GetXaxis()->GetXmin() - 2 * hist->GetXaxis()->GetBinWidth(1),
                 hist->GetXaxis()->GetXmax() + 2 * hist->GetXaxis()->GetBinWidth(1), hist->GetNbinsY() + 4,
                 hist->GetYaxis()->GetXmin() - 2 * hist->GetYaxis()->GetBinWidth(1),
                 hist->GetYaxis()->GetXmax() + 2 * hist->GetYaxis()->GetBinWidth(1));
    */

    for (int ix = 1; ix <= hBoundaries->GetXaxis()->GetNbins(); ix++) {
      for (int iy = 1; iy <= hBoundaries->GetYaxis()->GetNbins(); iy++) {
        if (ix == 1 || ix == hBoundaries->GetXaxis()->GetNbins() || iy == 1 ||
            iy == hBoundaries->GetYaxis()->GetNbins())
          hBoundaries->SetBinContent(ix, iy, boundary);
        else { hBoundaries->SetBinContent(ix, iy, hist->GetBinContent(ix - 1, iy - 1)); }
      }
    }

    /*
    for (int ix = 1; ix <= hBoundaries->GetXaxis()->GetNbins(); ix++) {
      for (int iy = 1; iy <= hBoundaries->GetYaxis()->GetNbins(); iy++) {
        // Fill outmost extra bins with the boundary value
        if (ix == 1 || ix == hBoundaries->GetXaxis()->GetNbins() || iy == 1 || iy ==
    hBoundaries->GetYaxis()->GetNbins())
        // Fill the outer 4 bins of the inner extra bins with the boundary value
        else if ((ix == 2 && iy == hBoundaries->GetYaxis()->GetNbins() - 1) ||
                 (ix == hBoundaries->GetXaxis()->GetNbins() - 1 && iy == hBoundaries->GetYaxis()->GetNbins() - 1) ||
                 (ix == 2 && iy == 2) || (ix == hBoundaries->GetXaxis()->GetNbins() - 1 && iy == 2)) {
          hBoundaries->SetBinContent(ix, iy, boundary);
        }
        // Fill the rest of the inner bins with the adjacent values of the original histogram
        else if (ix == 2 && (iy >= 3 || iy <= hBoundaries->GetYaxis()->GetNbins() - 2)) {
          hBoundaries->SetBinContent(ix, iy, hist->GetBinContent(ix - 1, iy - 2));
        } else if (ix == hBoundaries->GetXaxis()->GetNbins() - 1 &&
                   (iy >= 3 || iy <= hBoundaries->GetYaxis()->GetNbins() - 2)) {
          hBoundaries->SetBinContent(ix, iy, hist->GetBinContent(ix - 3, iy - 2));
        } else if (iy == 2 && (ix >= 3 || ix <= hBoundaries->GetXaxis()->GetNbins() - 2)) {
          hBoundaries->SetBinContent(ix, iy, hist->GetBinContent(ix - 2, iy - 1));
        } else if (iy == hBoundaries->GetYaxis()->GetNbins() - 1 &&
                   (ix >= 3 || ix <= hBoundaries->GetXaxis()->GetNbins() - 2)) {
          hBoundaries->SetBinContent(ix, iy, hist->GetBinContent(ix - 2, iy - 3));
        }
        // Copy the inner part from the original histogram
        else {
          hBoundaries->SetBinContent(ix, iy, hist->GetBinContent(ix - 2, iy - 2));
        }
      }
    }
    */
    return hBoundaries;
  }

  /**
   * Transform the chi2 valley into a hill to help ROOTs contour mechanism.
   *
   * @param hist   The 2D histogram.
   * @param offset A chi2 offset, usually around 30 units.
   *
   * @return       The transformed 2D histogram.
   */
  std::unique_ptr<TH2> transformChi2valleyToHill(const TH2* hist, const double offset) {
    const double chi2min = hist->GetMinimum();
    auto newHist = Utils::histHardCopy(hist, false, true);
    for (int ix = 1; ix <= hist->GetXaxis()->GetNbins(); ix++) {
      for (int iy = 1; iy <= hist->GetYaxis()->GetNbins(); iy++) {
        newHist->SetBinContent(ix, iy, -hist->GetBinContent(ix, iy) + offset + chi2min);
      }
    }
    return newHist;
  }
}  // namespace

// using namespace std;

ConfidenceContours::ConfidenceContours(const OptParser* arg) {
  assert(arg);
  m_arg = arg;
}

/**
 * Add a new contour that is just the plotted area.
 *
 * @param hist Histogram providing the dimensions of the plotted area.
 */
void ConfidenceContours::addFilledPlotArea(const TH2* hist) {
  const auto xmin = hist->GetXaxis()->GetXmin();
  const auto xmax = hist->GetXaxis()->GetXmax();
  const auto ymin = hist->GetYaxis()->GetXmin();
  const auto ymax = hist->GetYaxis()->GetXmax();
  // make new graph covering the plotted area
  auto g = new TGraph(m_nMaxContours);
  g->SetPoint(0, xmin, ymin);
  g->SetPoint(1, xmin, ymax);
  g->SetPoint(2, xmax, ymax);
  g->SetPoint(3, xmax, ymin);
  for (int i = 4; i < m_nMaxContours; i++) g->SetPoint(i, xmin, ymin);
  // make a new Contour object from it
  auto l = new TList();
  l->Add(g);
  m_contours.push_back(std::make_unique<Contour>(m_arg, l));
}

///
/// Compute the raw N sigma confidence contours from a 2D histogram
/// holding either the chi2 or the p-value curve. The resulting
/// contours are stored into the m_contours member.
///
/// \param hist - the 2D histogram
/// \param type - the type of the 2D histogram, either chi2 or p-value
///
void ConfidenceContours::computeContours(const TH2* inHist, const Utils::histogramType type, const int id) {
  using Utils::histogramType;
  auto hist = Utils::clone<TH2>(inHist);
  if (m_arg->debug)
    std::cout << "ConfidenceContours::computeContours() : making contours of histogram " << hist->GetName() << ", type "
              << (type == histogramType::kChi2 ? "chi2" : "p-value") << std::endl;
  // clean up contours from a previous call
  m_contours.clear();

  // transform chi2 from valley to hill
  const double offset = 100.;
  if (type == histogramType::kChi2) hist = transformChi2valleyToHill(hist.get(), offset);

  // add boundaries
  auto histb = addBoundaryBins(hist.get());

  // make contours
  histb->SetContour(m_nMaxContours);
  if (type == histogramType::kChi2) {
    // chi2 units
    if (m_arg->plot2dcl[id] > 0) {
      for (int i = 0; i < m_nMaxContours; i++) {
        int cLev = m_nMaxContours - 1 - i;
        // hack for >= 9 when ROOT precision fails
        if (i == 8)
          histb->SetContourLevel(cLev, offset - 83.9733);  // 9 sigma
        else if (i == 9)
          histb->SetContourLevel(cLev, offset - 99.2688);  // 10 sigma
        else if (i == 10)
          histb->SetContourLevel(cLev, offset - 114.564);  // 11 sigma
        else
          histb->SetContourLevel(cLev, offset - TMath::ChisquareQuantile(1. - TMath::Prob((i + 1) * (i + 1), 1), 2));
      }
    } else {
      for (int i = 0; i < m_nMaxContours; i++) {
        int cLev = m_nMaxContours - 1 - i;
        histb->SetContourLevel(cLev, offset - (i + 1) * (i + 1));
      }
    }
  } else {
    // p-value units
    if (m_arg->plot2dcl[id] > 0) {
      for (int i = 0; i < m_nMaxContours; i++) {
        int cLev = m_nMaxContours - 1 - i;
        histb->SetContourLevel(cLev, TMath::Prob((i + 1) * (i + 1), 1));
      }
    } else {
      for (int i = 0; i < m_nMaxContours; i++) {
        int cLev = m_nMaxContours - 1 - i;
        histb->SetContourLevel(cLev, TMath::Prob((i + 1) * (i + 1), 2));
      }
    }
  }

  // create and access the contours
  if (m_arg->interactive) gROOT->SetBatch(true);  // don't display the temporary canvas
  const auto previousGPad = gPad;
  auto ctmp = Utils::newNoWarnTCanvas(Utils::getUniqueRootName(), "ctmp");
  histb->Draw("contlist");
  gPad->Update();  // needed to be able to access the contours as TGraphs
  auto contours = dynamic_cast<TObjArray*>(gROOT->GetListOfSpecials()->FindObject("contours"));
  if (m_arg->interactive)
    gROOT->SetBatch(false);  // Else some canvases get screwed up badly resulting in corrupted PDF files.

  // Access contours.
  // They get filled in reverse order, and depend on how many are actually present.
  // If all 5 are filled, index 0 is 5sigma. If only 2 are filled, index 0 is 2 sigma.
  int nEmptyContours = 0;
  for (int ic = m_nMaxContours - 1; ic >= 0; ic--) {
    if (((TList*)contours->At(ic))->IsEmpty()) nEmptyContours++;
  }
  for (int ic = m_nMaxContours - 1; ic >= 0; ic--) {
    if (!(((TList*)contours->At(ic))->IsEmpty())) {
      m_contours.push_back(std::make_unique<Contour>(m_arg, static_cast<TList*>(contours->At(ic))));
      m_contours.back()->setSigma(5 - nEmptyContours - ic);
    }
  }
  // Move gPad to location previous to creation of ctmp to avoid that the caller needs to reset it
  previousGPad->cd();

  // Add the entire plotted area, if one requested contour is empty, i.e. it contains the entire plot range
  if (nEmptyContours > 0) { addFilledPlotArea(hist.get()); }

  // Magnetic boundaries
  if (m_arg->plotmagnetic) {
    for (int ic = m_nMaxContours - 1; ic >= 0; ic--) {
      if (ic >= m_contours.size()) continue;
      m_contours[ic]->magneticBoundaries(hist.get());
    }
  }
}

///
/// Draw the contours into the currently active Canvas.
///
void ConfidenceContours::Draw() {
  if (m_arg->debug) std::cout << "ConfidenceContours::Draw() : Start execution..." << std::endl;
  if (m_contstoplots.empty()) {
    for (int i = m_arg->plotnsigmacont - 1; i >= 0; i--) {
      if (i >= m_contours.size()) continue;
      m_contours[i]->setStyle(m_linecolor[i], m_linestyle[i], m_linewidth[i], m_fillcolor[i], m_fillstyle[i]);
      m_contours[i]->setTransparency(m_transparency);
      m_contours[i]->Draw();
    }
  } else {
    for (int ind = m_contstoplots.size() - 1; ind >= 0; ind--) {
      int i = m_contstoplots[ind] - 1;
      if (i >= m_contours.size()) continue;
      m_contours[i]->setStyle(m_linecolor[i], m_linestyle[i], m_linewidth[i], m_fillcolor[i], m_fillstyle[i]);
      m_contours[i]->setTransparency(m_transparency);
      m_contours[i]->Draw();
    }
  }
}

///
/// Draw the contours into the currently active Canvas.
///
void ConfidenceContours::DrawDashedLine() {
  if (m_arg->debug) std::cout << "ConfidenceContours::DrawDashedLine() : Start execution..." << std::endl;
  if (m_contstoplots.empty()) {
    for (int i = m_arg->plotnsigmacont - 1; i >= 0; i--) {
      if (i >= m_contours.size()) continue;
      m_contours[i]->setStyle(m_linecolor[i], kDashed, m_linewidth[i], 0, 0);
      m_contours[i]->DrawLine();
    }
  } else {
    for (int ind = m_contstoplots.size() - 1; ind >= 0; ind--) {
      int i = m_contstoplots[ind] - 1;
      if (i >= m_contours.size()) continue;
      m_contours[i]->setStyle(m_linecolor[i], m_linestyle[i], m_linewidth[i], m_fillcolor[i], m_fillstyle[i]);
      m_contours[i]->setTransparency(m_transparency);
      m_contours[i]->Draw();
    }
  }
}

///
/// Set the contour style.
///
void ConfidenceContours::setStyle(const std::vector<int>& linecolor, const std::vector<int>& linestyle,
                                  const std::vector<int>& linewidth, const std::vector<int>& fillcolor,
                                  const std::vector<int>& fillstyle) {
  m_linecolor = linecolor;
  m_linestyle = linestyle;
  m_fillcolor = fillcolor;
  m_fillstyle = fillstyle;
  m_linewidth = linewidth;
  if (m_arg->plotnsigmacont > m_linestyle.size()) {
    // Not enough styles were given for the number of contours to be plotted
    std::cout << "ConfidenceContours::setStyle() : ERROR : not enough sigma contour styles defined! ";
    std::cout << "Reusing style of " << m_linestyle.size() << " sigma contour." << std::endl;
    const int laststyle = m_linestyle.size() - 1;
    if (laststyle < 0) {
      std::cout << "ConfidenceContours::setStyle() : ERROR : linestyle is empty. Exit." << std::endl;
      std::exit(1);
    }
    for (int i = m_linestyle.size(); i < m_arg->plotnsigmacont; i++) {
      m_linecolor.push_back(m_linecolor[laststyle]);
      m_linestyle.push_back(m_linestyle[laststyle]);
      m_fillcolor.push_back(m_fillcolor[laststyle]);
      m_fillstyle.push_back(m_fillstyle[laststyle]);
      m_linewidth.push_back(m_linewidth[laststyle]);
    }
  }
}
