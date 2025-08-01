#include <OneMinusClPlot.h>
#include <Rounder.h>
#include <Utils.h>

#include <algorithm>
#include <cmath>
#include <memory>

#include <TColor.h>
#include <TF1.h>
#include <TGaxis.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TGraphSmooth.h>
#include <TH1.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPaveText.h>

#include <cmath>
#include <format>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace {
  auto msgBase = [](const std::string& prefix, const std::string& msg, std::ostream& stream = std::cout) {
    auto msgOut = Utils::replaceAll(msg, "\n", "\n" + std::string(prefix.size(), ' '));
    stream << prefix << msgOut << std::endl;
  };

  auto errBase = [](const std::string& prefix, const std::string& msg) { msgBase(prefix, msg, std::cerr); };
}  // namespace

/// Constructor.
OneMinusClPlot::OneMinusClPlot(OptParser* arg, TString name, TString title) : OneMinusClPlotAbs(arg, name, title) {}

/// Helper function for operation that is repeated several times.
std::unique_ptr<TH1> OneMinusClPlot::getHistogram(MethodAbsScan* s, const int CLsType, const bool removeErrs) const {
  auto error = [](const std::string& msg) {
    errBase("OneMinusClPlot::getHistogram : ERROR : ", msg);
    exit(1);
  };

  if (!s) error("MethodAbsScan is equal to nullptr");
  const TH1* pH = nullptr;
  if (CLsType == 0)
    pH = s->getHCL();
  else if (CLsType == 1)
    pH = s->getHCLs();
  else if (CLsType == 2)
    pH = s->getHCLsFreq();
  if (!pH) error(std::format("Could not retrieve the histogram from MethodAbsScan {:s}", std::string(s->getName())));
  auto h = Utils::clone<TH1>(pH);
  // TODO this should not be needed
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    if (std::isnan(h->GetBinContent(i)) || std::isinf(h->GetBinContent(i))) h->SetBinContent(i, 0.);
  }
  // remove errors the hard way, else root always plots them
  if (removeErrs) h = Utils::histHardCopy(h.get(), true, true);
  return h;
}

/**
 * Make a TGraph out of a 1D histogram holding a 1-CL curve.
 *
 * The strategy is to always convert the 1-CL histogram (hCL) into a TGraph. This way we can add known points
 * (solutions, points at end of scan range) and also have a filled area without line smoothing. This is not possible
 * with histograms due to a Root bug.
 *
 * The function draws the TGraphs, and returns a pointer to the TGraph object that can be used in the TLegend.
 *
 * Markers are plotted if the method name of the scanner is "Plugin" or "BergerBoos" or "DatasetsPlugin".  One can plot
 * a line instead of points even for the Plugin method by using setPluginMarkers().
 *
 * For the angle variables, a new axis is painted that is in Deg.
 *
 * @param s      The scanner to plot.
 * @param first  Flags whether this is the first graph to be plotted.
 * @param last   Flags whether this is the last graph to be plotted.
 * @param filled Flag to steer whether the graph should be plot filled or with a simple line.
 */
TGraph* OneMinusClPlot::scan1dPlot(MethodAbsScan* s, const bool first, const bool last, const bool filled,
                                   const int CLsType) {
  auto debug = [](const std::string& msg) { msgBase("OneMinusClPlot::scan1dPlot : DEBUG : ", msg); };

  if (arg->debug) debug(std::format("Plot {:s} ({:s})", std::string(s->getName()), std::string(s->getMethodName())));
  canvas->cd();
  bool plotPoints = (s->getMethodName() == "Plugin" || s->getMethodName() == "BergerBoos" ||
                     s->getMethodName() == "DatasetsPlugin") &&
                    plotPluginMarkers;

  auto hCL = getHistogram(s, CLsType, !plotPoints);
  // disable any statistics box
  hCL->SetStats(0);

  // Convert the histogram into a TGraph so we can add the solution.
  // Also, the lf2 drawing option is broken in latest root versions.
  using gType = TGraph;
  std::unique_ptr<gType> g;
  if (plotPoints)
    g = std::make_unique<TGraphErrors>(hCL->GetNbinsX());
  else
    g = std::make_unique<TGraph>(hCL->GetNbinsX());
  g->SetName(Utils::getUniqueRootName());
  for (int i = 0; i < hCL->GetNbinsX(); i++) {
    g->SetPoint(i, hCL->GetBinCenter(i + 1), hCL->GetBinContent(i + 1));
    if (plotPoints) dynamic_cast<TGraphErrors*>(g.get())->SetPointError(i, 0., hCL->GetBinError(i + 1));
  }

  // add solution -- this does not make sense for the one-sided test statistic, which only gives non-default values for
  // mu > muhat
  if (arg->teststatistic != 1) {
    if (!s->getSolutions().empty()) {
      g = Utils::addPointToGraphAtFirstMatchingX(g.get(), s->getScanVar1Solution(0), 1.);
    }
  }

  // // set last point to the same p-value as first point by hand
  // // some angle plots sometimes don't manage to do it by themselves...
  // if ( arg->isQuickhack(XX) )
  // {
  //   double pointx0, pointy0, err0;
  //   double pointx1, pointy1, err1;
  //   g->GetPoint(0, pointx0, pointy0);
  //   g->GetPoint(g->GetN()-1, pointx1, pointy1);
  //   g->SetPoint(g->GetN()-1, pointx1, pointy0);
  //   if ( plotPoints ) err0 = ((TGraphErrors*)g)->GetErrorY(0);
  //   if ( plotPoints ) ((TGraphErrors*)g)->SetPointError(g->GetN()-1, 0.0, err0);
  // }

  // add end points of scan range
  if (!plotPoints) {
    double pointx0, pointy0;
    auto gNew = std::make_unique<TGraph>(g->GetN() + 4);
    gNew->SetName(Utils::getUniqueRootName());
    for (int i = 0; i < g->GetN(); i++) {
      g->GetPoint(i, pointx0, pointy0);
      gNew->SetPoint(i + 2, pointx0, pointy0);
    }

    // add origin
    gNew->SetPoint(0, hCL->GetXaxis()->GetXmin(), 0);

    // add a point at first y height but at x=origin.
    g->GetPoint(0, pointx0, pointy0);
    gNew->SetPoint(1, hCL->GetXaxis()->GetXmin(), pointy0);

    // add a point at last y height but at x=xmax.
    g->GetPoint(g->GetN() - 1, pointx0, pointy0);
    gNew->SetPoint(gNew->GetN() - 2, hCL->GetXaxis()->GetXmax(), pointy0);

    // add a point at xmax, 0
    gNew->SetPoint(gNew->GetN() - 1, hCL->GetXaxis()->GetXmax(), 0);
    g = std::move(gNew);
  }

  int color = s->getLineColor();
  if (CLsType > 0 && s->getMethodName().Contains("Plugin") && !arg->plotpluginonly) {
    if (CLsType == 1)
      color = kBlue - 8;
    else if (CLsType == 2)
      color = kBlue - 2;
  } else if (CLsType > 0) {
    if (CLsType == 1) color = s->getLineColor() + 1;
    if (CLsType == 2) color = s->getLineColor() + 1;
  }
  g->SetLineColor(color);

  if (filled) {
    double alpha = arg->isQuickhack(12) ? 0.4 : 1.;
    if (arg->isQuickhack(24)) alpha = 0.;
    if (s->getFillColor() > 0 && CLsType == 0)
      g->SetFillColorAlpha(s->getFillColor(), alpha);
    else
      g->SetFillColorAlpha(color, alpha);
    g->SetFillStyle(s->getFillStyle());
    g->SetLineWidth(s->getLineWidth());
    g->SetLineStyle(s->getLineStyle());
    g->SetLineColor(s->getLineColor());
  } else {
    g->SetLineWidth(s->getLineWidth());
    g->SetLineStyle(s->getLineStyle());
    g->SetLineColor(s->getLineColor());
    g->SetFillStyle(s->getFillStyle());
    if (last && arg->isQuickhack(25)) g->SetLineWidth(3);
  }

  if (CLsType > 0) g->SetLineColor(color);

  if (plotPoints) {
    g->SetLineWidth(1);
    g->SetLineColor(color);
    g->SetMarkerColor(color);
    g->SetMarkerStyle(8);
    g->SetMarkerSize(0.6);
    if (CLsType == 1) {
      g->SetMarkerStyle(33);
      g->SetMarkerSize(1);
    }
    if (CLsType == 2) { g->SetMarkerStyle(25); }
  }

  // build a histogram which holds the axes
  auto min = arg->scanrangeMin == arg->scanrangeMax ? hCL->GetXaxis()->GetXmin() : arg->scanrangeMin;
  auto max = arg->scanrangeMin == arg->scanrangeMax ? hCL->GetXaxis()->GetXmax() : arg->scanrangeMax;
  auto haxes = makeOwnedTObject<TH1F>("haxes" + Utils::getUniqueRootName(), "", 100, min, max);
  haxes->SetStats(0);
  if (arg->xtitle == "")
    haxes->GetXaxis()->SetTitle(s->getScanVar1()->GetTitle());
  else
    haxes->GetXaxis()->SetTitle(arg->xtitle);
  haxes->GetYaxis()->SetTitle("1#minusCL");
  haxes->GetXaxis()->SetLabelFont(font);
  haxes->GetYaxis()->SetLabelFont(font);
  haxes->GetXaxis()->SetTitleFont(font);
  haxes->GetYaxis()->SetTitleFont(font);
  haxes->GetXaxis()->SetTitleOffset(0.9);
  haxes->GetYaxis()->SetTitleOffset(0.95);
  haxes->GetXaxis()->SetLabelSize(labelsize);
  haxes->GetYaxis()->SetLabelSize(labelsize);
  haxes->GetXaxis()->SetTitleSize(titlesize);
  haxes->GetYaxis()->SetTitleSize(titlesize);
  int xndiv = arg->ndiv == -1 ? 407 : abs(arg->ndiv);
  bool optimizeNdiv = arg->ndiv < 0 ? true : false;
  haxes->GetXaxis()->SetNdivisions(xndiv, optimizeNdiv);
  haxes->GetYaxis()->SetNdivisions(407, true);

  // plot y range
  double plotYMax;
  double plotYMin;
  if (plotLegend && !arg->isQuickhack(22)) {
    if (arg->plotlog) {
      plotYMin = 1.e-3;
      plotYMax = 10.;
    } else {
      plotYMin = 0.0;
      plotYMax = 1.3;
    }
  } else {
    if (arg->plotlog) {
      plotYMin = 1.e-3;
      plotYMax = 1.0;
    } else {
      plotYMin = 0.0;
      plotYMax = 1.0;
    }
  }
  // change if passed as option
  plotYMin = arg->plotymin > 0. ? arg->plotymin : plotYMin;
  plotYMax = arg->plotymax > 0. ? arg->plotymax : plotYMax;

  haxes->GetYaxis()->SetRangeUser(plotYMin, plotYMax);
  haxes->Draw("axissame");
  g->SetHistogram(haxes);

  TString drawOption = "";
  if (plotPoints)
    drawOption += " pe";
  else if (filled)
    drawOption += " F";
  else
    drawOption += " L";
  if (first) drawOption += " A";
  auto pG = getTObjectOwnership<gType>(std::move(g));
  pG->Draw(drawOption);
  // if ( drawOption.Contains("F") ) ((TGraph*)g->Clone())->Draw("L");

  gPad->Update();
  const double ymin = gPad->GetUymin();
  const double ymax = gPad->GetUymax();
  const double xmin = gPad->GetUxmin();
  const double xmax = gPad->GetUxmax();

  const double wmin = Utils::isAngle(s->getScanVar1()) ? Utils::RadToDeg(xmin) : xmin;
  const double wmax = Utils::isAngle(s->getScanVar1()) ? Utils::RadToDeg(xmax) : xmax;

  if (last) {
    // Add top axis
    TString chopt = "-U";             // - = downward ticks, U = unlabeled
    if (!optimizeNdiv) chopt += "N";  // n = no bin optimization
    auto axist = makeOwnedTObject<TGaxis>(xmin, 1., xmax, 1., wmin, wmax, xndiv, chopt);
    axist->SetName("axist");
    axist->SetLineWidth(1);
    axist->Draw();
  }

  // for the angles, draw a new axis in units of degrees
  if (Utils::isAngle(s->getScanVar1())) {

    if (arg->xtitle == "")
      haxes->GetXaxis()->SetTitle(s->getScanVar1()->GetTitle() + TString(" [#circ]"));
    else
      haxes->GetXaxis()->SetTitle(arg->xtitle);
    haxes->GetXaxis()->SetNdivisions(0);  // disable old axis
    if (last) {
      using Utils::RadToDeg;

      // new bottom axis
      double axisbMin = RadToDeg(xmin);
      double axisbMax = RadToDeg(xmax);
      if (arg->isQuickhack(3)) {  ///< see documentation of --qh option in OptParser.cpp
        axisbMin += 180.;
        axisbMax += 180.;
      }
      TString chopt = "";               // - = downward ticks, U = unlabeled
      if (!optimizeNdiv) chopt += "N";  // n = no bin optimization
      auto axisb = makeOwnedTObject<TGaxis>(xmin, ymin, xmax, ymin, axisbMin, axisbMax, xndiv, chopt);
      axisb->SetName("axisb");
      axisb->SetLabelFont(font);
      axisb->SetLabelSize(labelsize);
      if (arg->xtitle == "")
        axisb->SetTitle(s->getScanVar1()->GetTitle() + TString(" [#circ]"));
      else
        axisb->SetTitle(arg->xtitle);
      axisb->SetTitleOffset(0.85);
      axisb->SetTitleSize(titlesize);
      axisb->SetTitleFont(font);
      axisb->Draw();
    }
  }

  if (last) {
    // add right axis
    TGaxis* axisr;
    if (arg->plotlog) {
      double f3min = 1e-3;
      double f3max = (plotLegend && !arg->isQuickhack(22)) ? 10. : 1.;
      auto f3 = makeOwnedTObject<TF1>("f3", "log10(x)", f3min, f3max);
      axisr = makeOwnedTObject<TGaxis>(xmax, f3min, xmax, f3max, "f3", 510, "G+");
    } else {
      axisr = makeOwnedTObject<TGaxis>(xmax, ymin, xmax, ymax, 0, (plotLegend && !arg->isQuickhack(22)) ? 1.3 : 1.0,
                                       407, "+");
    }
    axisr->SetLabelSize(0);
    axisr->SetLineWidth(1);
    axisr->SetName("axisr");
    axisr->SetLabelColor(kWhite);
    axisr->SetTitleColor(kWhite);
    axisr->Draw();

    // redraw right axis as well because the 1-CL graph can cover the old one
    haxes->Draw("axissame");
  }

  return pG;
}

/**
 * Make a plot out of a 1D histogram holding a 1-CL curve, without smoothing it.
 *
 * @param s       The scanner to plot.
 * @param first   Set this to true for the first plotted scanner.
 * @param CLsType Steers what histogram to plot (Prob, Plugin, etc.)
 */
void OneMinusClPlot::scan1dPlotSimple(MethodAbsScan* s, const bool first, const int CLsType) {
  if (arg->debug) {
    std::cout << "OneMinusClPlot::scan1dPlotSimple() : plotting ";
    std::cout << s->getName() << " (" << s->getMethodName() << ")" << std::endl;
  }
  canvas->cd();

  auto hCL = getTObjectOwnership<TH1>(getHistogram(s, CLsType, false));

  int color = s->getLineColor();
  if (CLsType == 1) color = color + 2;
  hCL->SetStats(0);
  hCL->SetLineColor(color);
  hCL->SetMarkerColor(color);
  hCL->SetLineWidth(2);
  hCL->SetLineStyle(s->getLineStyle());
  hCL->SetMarkerColor(color);
  hCL->SetMarkerStyle(8);
  hCL->SetMarkerSize(0.6);
  hCL->GetYaxis()->SetNdivisions(407, true);
  if (arg->xtitle == "")
    hCL->GetXaxis()->SetTitle(s->getScanVar1()->GetTitle());
  else
    hCL->GetXaxis()->SetTitle(arg->xtitle);
  hCL->GetYaxis()->SetTitle("1#minusCL");
  hCL->GetXaxis()->SetLabelFont(font);
  hCL->GetYaxis()->SetLabelFont(font);
  hCL->GetXaxis()->SetTitleFont(font);
  hCL->GetYaxis()->SetTitleFont(font);
  hCL->GetXaxis()->SetTitleOffset(0.9);
  hCL->GetYaxis()->SetTitleOffset(0.95);
  hCL->GetXaxis()->SetLabelSize(labelsize);
  hCL->GetYaxis()->SetLabelSize(labelsize);
  hCL->GetXaxis()->SetTitleSize(titlesize);
  hCL->GetYaxis()->SetTitleSize(titlesize);
  if (plotLegend && !arg->isQuickhack(22)) {
    if (arg->plotlog)
      hCL->GetYaxis()->SetRangeUser(1e-3, 10);
    else
      hCL->GetYaxis()->SetRangeUser(0.0, 1.3);
  } else {
    if (arg->plotlog)
      hCL->GetYaxis()->SetRangeUser(1e-3, 1);
    else
      hCL->GetYaxis()->SetRangeUser(0.0, 1.0);
  }
  hCL->Draw(first ? "" : "same");
}

/**
 * Make a plot for the CLs stuff.
 *
 * Convert the hCLExp and hCLErr histogrms into TGraphs and TGraphAsymmErrors
 *
 * @param s The scanner to plot.
 * @param smooth
 */
void OneMinusClPlot::scan1dCLsPlot(MethodAbsScan* s, const bool smooth, const bool obsError) {
  auto info = [](const std::string& msg) { msgBase("OneMinusClPlot::scan1dCLsPlot() : ", msg); };
  if (arg->debug) info(std::format("plotting {:s} ({:s})", std::string(s->getName()), std::string(s->getMethodName())));
  canvas->cd();

  if (!s->checkCLs()) {
    info("Cannot plot");
    return;
  }

  using Utils::getUniqueRootName;
  auto hObs = Utils::clone<TH1>(s->getHCLsFreq());
  auto hExp = Utils::clone<TH1>(s->getHCLsExp());
  auto hErr1Up = Utils::clone<TH1>(s->getHCLsErr1Up());
  auto hErr1Dn = Utils::clone<TH1>(s->getHCLsErr1Dn());
  auto hErr2Up = Utils::clone<TH1>(s->getHCLsErr2Up());
  auto hErr2Dn = Utils::clone<TH1>(s->getHCLsErr2Dn());

  if (!hObs) info("problem - can't find histogram hObs");
  if (!hExp) info("problem - can't find histogram hExp");
  if (!hErr1Up) info("problem - can't find histogram hErr1Up");
  if (!hErr1Dn) info("problem - can't find histogram hErr1Dn");
  if (!hErr2Up) info("problem - can't find histogram hErr2Up");
  if (!hErr2Dn) info("problem - can't find histogram hErr2Dn");

  // convert obs to graph
  using Utils::convertTH1ToTGraph;
  auto gObs = convertTH1ToTGraph(hObs.get(), obsError);
  auto gExpRaw = convertTH1ToTGraph(hExp.get());
  auto gErr1UpRaw = convertTH1ToTGraph(hErr1Up.get());
  auto gErr1DnRaw = convertTH1ToTGraph(hErr1Dn.get());
  auto gErr2UpRaw = convertTH1ToTGraph(hErr2Up.get());
  auto gErr2DnRaw = convertTH1ToTGraph(hErr2Dn.get());

  // smoothing if needed
  std::unique_ptr<TGraph> gExp;
  std::unique_ptr<TGraph> gErr1Up;
  std::unique_ptr<TGraph> gErr1Dn;
  std::unique_ptr<TGraph> gErr2Up;
  std::unique_ptr<TGraph> gErr2Dn;

  const auto xCentral = s->getScanVar1Solution();
  auto smoother = std::make_unique<TGraphSmooth>();
  if (smooth) {
    if (arg->debug) info("smoothing graphs");
    gExp = std::unique_ptr<TGraph>(static_cast<TGraph*>((smoother->SmoothSuper(gExpRaw.get())->Clone("gExp"))));
    // gErr1Up =
    // std::unique_ptr<TGraph>(static_cast<TGraph*>((smoother->SmoothSuper(gErr1UpRaw.get())->Clone("gErr1Up"))));
    gErr1Dn =
        std::unique_ptr<TGraph>(static_cast<TGraph*>((smoother->SmoothSuper(gErr1DnRaw.get())->Clone("gErr1Dn"))));
    // gErr2Up =
    // std::unique_ptr<TGraph>(static_cast<TGraph*>((smoother->SmoothSuper(gErr2UpRaw.get())->Clone("gErr2Up"))));
    gErr2Dn =
        std::unique_ptr<TGraph>(static_cast<TGraph*>((smoother->SmoothSuper(gErr2DnRaw.get())->Clone("gErr2Dn"))));

    // alternative smoothing option, needs more fiddling
    // gExp = std::unique_ptr<TGraph>(static_cast<TGraph*>((smoother->SmoothKern(gExpRaw.get(),
    // "normal",hExp->GetBinWidth(1) * 4)->Clone("gExp"))));
    gErr1Up = std::unique_ptr<TGraph>(static_cast<TGraph*>(
        (smoother->SmoothKern(gErr1UpRaw.get(), "normal", hErr1Up->GetBinWidth(1) * 2, hErr1Up->GetNbinsX())
             ->Clone("gErr1Up"))));
    // gErr1Dn = std::unique_ptr<TGraph>(static_cast<TGraph*>((smoother->SmoothKern(gErr1DnRaw.get(), "normal",
    // hErr1Dn->GetBinWidth(1) * 4)->Clone("gErr1Dn"))));
    gErr2Up = std::unique_ptr<TGraph>(static_cast<TGraph*>(
        (smoother->SmoothKern(gErr2UpRaw.get(), "normal", hErr2Up->GetBinWidth(1) * 2, hErr1Up->GetNbinsX())
             ->Clone("gErr2Up"))));
    // gErr2Dn = std::unique_ptr<TGraph>(static_cast<TGraph*>((smoother->SmoothKern(gErr2DnRaw.get(), "normal",
    // hErr2Dn->GetBinWidth(1) * 4)->Clone("gErr2Dn"))));

    // make sure the CLs=1 points do NOT get smoothed away
    double* xvals = gExp->GetX();
    // double* xvalsRaw = gExpRaw->GetX();
    // double* yvalsRawExp = gExpRaw->GetY();
    double* yvalsRawExpErr1Up = gErr1UpRaw->GetY();
    double* yvalsRawExpErr2Up = gErr2UpRaw->GetY();

    if (arg->teststatistic == 1) {
      for (int i = 0; i < gExp->GetN(); i++) {
        // std::cout << xvalsRaw[i] << "\t" <<xvals[i] << std::endl;
        // if(yvalsRawExp[i]>0.99){
        //  gExp->SetPoint(i,xvals[i], 1.0);
        // }
        if (yvalsRawExpErr1Up[i] > 0.99) { gErr1Up->SetPoint(i, xvals[i], 1.0); }
        if (yvalsRawExpErr2Up[i] > 0.99) { gErr2Up->SetPoint(i, xvals[i], 1.0); }
      }
    }

    // fix point 0 to CLs=1 for all expected curves

    gExp->SetPoint(0, 0., 1.);
    gErr1Up->SetPoint(0, 0., 1.);
    gErr1Dn->SetPoint(0, 0., 1.);
    gErr2Up->SetPoint(0, 0., 1.);
    gErr2Dn->SetPoint(0, 0., 1.);

    if (arg->teststatistic == 1) {
      // remove all observed lines with x<xmeas
      if (arg->debug) info("remove all observed lines with mu<muhat in CLs plot");
      double* xvalsobs = gObs->GetX();
      double* yvalsobs = gObs->GetY();
      double* xerrsobs = gObs->GetEX();
      double* yerrsobs = gObs->GetEY();
      int valabove = gObs->GetN();
      int nentries = gObs->GetN();
      for (int i = 0; i < gObs->GetN(); i++) {
        if (xvalsobs[i] < (xCentral + (hObs->GetBinWidth(1) / 2.))) valabove--;
      }
      // std::cout << "Found entries for obs " << valabove << "\t" << nentries << std::endl;

      auto gObs_new = std::make_unique<TGraphErrors>(valabove);
      int k = 0;
      for (int i = 0; i < nentries; i++) {
        if (xvalsobs[i] < (xCentral + (hObs->GetBinWidth(1) / 2.))) {
          // std::cout << "Ignoring " << xvalsobs[i] << "\t" << yvalsobs[i] << std::endl;
          continue;
        } else {
          // std::cout << "SetPoint " << k << "\t" << xvalsobs[i] << "\t" << yvalsobs[i] << std::endl;
          gObs_new->SetPoint(k, xvalsobs[i], yvalsobs[i]);
          gObs_new->SetPointError(k, xerrsobs[i], yerrsobs[i]);
          k++;
        }
      }
      gObs = std::move(gObs_new);
    }

    if (arg->debug) info("done smoothing graphs");
  } else {
    gExp = std::move(gExpRaw);
    gErr1Up = std::move(gErr1UpRaw);
    gErr1Dn = std::move(gErr1DnRaw);
    gErr2Up = std::move(gErr2UpRaw);
    gErr2Dn = std::move(gErr2DnRaw);
  }

  if (!gObs) info("problem - null graph gObs");
  if (!gExp) info("problem - null graph gExp");
  if (!gErr1Up) info("problem - null graph gErr1Up");
  if (!gErr1Dn) info("problem - null graph gErr1Dn");
  if (!gErr2Up) info("problem - null graph gErr2Up");
  if (!gErr2Dn) info("problem - null graph gErr2Dn");

  gObs->SetName("gObs");
  gExp->SetName("gExp");
  gErr1Up->SetName("gErr1Up");
  gErr1Dn->SetName("gErr1Dn");
  gErr2Up->SetName("gErr2Up");
  gErr2Dn->SetName("gErr2Dn");

  // now make the graphs for the error bands
  auto gErr1 = std::make_unique<TGraphAsymmErrors>(gExp->GetN());
  gErr1->SetName("gErr1");
  auto gErr2 = std::make_unique<TGraphAsymmErrors>(gExp->GetN());
  gErr2->SetName("gErr2");

  double x, y, yerrUp, yerrDn;
  double xerr = (hExp->GetBinCenter(2) - hExp->GetBinCenter(1)) / 2.;

  // protect against smoothing over 1
  for (int i = 0; i < gExp->GetN(); i++) {
    gExp->GetPoint(i, x, y);
    gExp->SetPoint(i, x, TMath::Min(y, 1.));
    gErr1Up->GetPoint(i, x, y);
    gErr1Up->SetPoint(i, x, TMath::Min(y, 1.));
    gErr1Dn->GetPoint(i, x, y);
    gErr1Dn->SetPoint(i, x, TMath::Min(y, 1.));
    gErr2Up->GetPoint(i, x, y);
    gErr2Up->SetPoint(i, x, TMath::Min(y, 1.));
    gErr2Dn->GetPoint(i, x, y);
    gErr2Dn->SetPoint(i, x, TMath::Min(y, 1.));
  }

  for (int i = 0; i < gExp->GetN(); i++) {
    gExp->GetPoint(i, x, y);
    gErr1->SetPoint(i, x, y);
    gErr2->SetPoint(i, x, y);

    gErr1Up->GetPoint(i, x, yerrUp);
    gErr1Dn->GetPoint(i, x, yerrDn);
    gErr1->SetPointError(i, xerr, xerr, y - yerrDn, yerrUp - y);

    gErr2Up->GetPoint(i, x, yerrUp);
    gErr2Dn->GetPoint(i, x, yerrDn);
    gErr2->SetPointError(i, xerr, xerr, y - yerrDn, yerrUp - y);
  }

  gErr2->SetFillColor(TColor::GetColor("#3182bd"));
  gErr2->SetLineColor(TColor::GetColor("#3182bd"));
  gErr1->SetFillColor(TColor::GetColor("#9ecae1"));
  gErr1->SetLineColor(TColor::GetColor("#9ecae1"));
  gExp->SetLineColor(kRed);
  gExp->SetLineWidth(3);
  gObs->SetLineColor(kBlack);
  gObs->SetMarkerColor(kBlack);
  gObs->SetLineWidth(3);
  gObs->SetMarkerSize(1);
  gObs->SetMarkerStyle(20);

  double min = arg->scanrangeMin == arg->scanrangeMax ? hObs->GetXaxis()->GetXmin() : arg->scanrangeMin;
  double max = arg->scanrangeMin == arg->scanrangeMax ? hObs->GetXaxis()->GetXmax() : arg->scanrangeMax;
  auto haxes = std::make_unique<TH1F>("haxes" + getUniqueRootName(), "", 100, min, max);
  haxes->SetDirectory(0);
  haxes->SetStats(0);
  if (arg->xtitle == "")
    haxes->GetXaxis()->SetTitle(s->getScanVar1()->GetTitle());
  else
    haxes->GetXaxis()->SetTitle(arg->xtitle);
  haxes->GetYaxis()->SetTitle("CL_{S}");
  haxes->GetXaxis()->SetLabelFont(font);
  haxes->GetYaxis()->SetLabelFont(font);
  haxes->GetXaxis()->SetTitleFont(font);
  haxes->GetYaxis()->SetTitleFont(font);
  haxes->GetXaxis()->SetTitleOffset(0.9);
  haxes->GetYaxis()->SetTitleOffset(0.95);
  haxes->GetXaxis()->SetLabelSize(labelsize);
  haxes->GetYaxis()->SetLabelSize(labelsize);
  haxes->GetXaxis()->SetTitleSize(titlesize);
  haxes->GetYaxis()->SetTitleSize(titlesize);
  int xndiv = arg->ndiv == -1 ? 407 : abs(arg->ndiv);
  bool optimizeNdiv = arg->ndiv < 0 ? true : false;
  haxes->GetXaxis()->SetNdivisions(xndiv, optimizeNdiv);
  haxes->GetYaxis()->SetNdivisions(407, true);
  haxes->GetYaxis()->SetRangeUser(0., 1.);

  // Legend: make the legend short, the text will extend over the boundary, but the symbol will be shorter
  double legendXmin = 0.68;
  double legendYmin = 0.58;
  double legendXmax = legendXmin + 0.25;
  double legendYmax = legendYmin + 0.22;
  auto leg = makeOwnedTObject<TLegend>(legendXmin, legendYmin, legendXmax, legendYmax);
  leg->SetFillColor(kWhite);
  leg->SetFillStyle(0);
  leg->SetLineColor(kWhite);
  leg->SetBorderSize(0);
  leg->SetTextFont(font);
  leg->SetTextSize(legendsize * 0.75);

  if (obsError)
    leg->AddEntry(gObs.get(), "Observed", "LEP");
  else
    leg->AddEntry(gObs.get(), "Observed", "LP");
  leg->AddEntry(gExp.get(), "Expected", "L");
  leg->AddEntry(gErr1.get(), "#pm 1#sigma", "F");
  leg->AddEntry(gErr2.get(), "#pm 2#sigma", "F");

  haxes->Draw("AXIS+");
  gErr2->Draw("E3same");
  gErr1->Draw("E3same");
  gExp->Draw("Lsame");
  if (obsError)
    gObs->Draw("LEPsame");
  else
    gObs->Draw("LPsame");
  leg->Draw("same");

  if (arg->CL.size() == 0) {
    drawCLguideLine(0.1);
  } else {
    drawCLguideLines();
  }

  // draw the solution
  if (arg->plotsolutions.size() > 0 && arg->plotsolutions[0] != 0) drawVerticalLine(xCentral, kBlack, kDashed);

  double yGroup = 0.83;
  if (arg->plotprelim || arg->plotunoff) yGroup = 0.8;
  drawGroup(yGroup);

  canvas->SetTicks();
  canvas->RedrawAxis();
  canvas->Update();
  canvas->Modified();
  canvas->Show();
  Utils::savePlot(canvas.get(), name + "_expected" + arg->plotext);
  canvas->SetTicks(false);
}

void OneMinusClPlot::drawVerticalLine(const double x, const int color, const int style) {
  canvas->cd();
  auto line = makeOwnedTObject<TLine>(x, 0., x, 1.);
  line->SetLineWidth(1);
  line->SetLineColor(color);
  line->SetLineStyle(style);
  line->Draw();
}

/**
 * Draw a vertical line a the position of the best solution.
 *
 * For the getDrawSolution() code, see OptParser, --ps.
 */
void OneMinusClPlot::drawSolutions() {
  canvas->cd();
  canvas->Update();
  double xmin = gPad->GetUxmin();
  double xmax = gPad->GetUxmax();
  int iDrawn = 0;

  for (int i = 0; i < scanners.size(); i++) {
    if (scanners[i]->getDrawSolution() == 0) continue;
    if (arg->debug)
      std::cout << "OneMinusClPlot::drawSolutions() : adding solution for scanner " << i << " ..." << std::endl;
    double xCentral = scanners[i]->getScanVar1Solution(arg->plotsoln[i]);
    double xCLmin = scanners[i]->getCLinterval(arg->plotsoln[i])->min;
    double xCLmax = scanners[i]->getCLinterval(arg->plotsoln[i])->max;
    int color = scanners[i]->getTextColor();

    // draw vertical lines at central value and upper/lower errors
    if (!arg->isQuickhack(19)) {
      drawVerticalLine(xCentral, color, kSolid);
      if (!arg->isQuickhack(20)) {
        drawVerticalLine(xCLmin, color, kDashed);
        drawVerticalLine(xCLmax, color, kDashed);
      }
    }

    // Draw text box with numerical values after the lines, so that it doesn't get covered -----------------------------

    // compute y position of the printed central value
    double yNumberMin = 0.6 - 0.13 * iDrawn;
    double yNumberMax = yNumberMin + 0.1;
    if (arg->plotlog) {
      double yNumberMinFirst = 0.1;
      if (arg->isQuickhack(1)) yNumberMinFirst = 0.175;
      yNumberMin = yNumberMinFirst / pow(3.0, iDrawn);  // move down by a constant shift on log scale
      yNumberMax = yNumberMin * 2.;
    }
    // if printsoly option then move a bit
    if (arg->printSolY > 0.) {
      yNumberMin += arg->printSolY;
      yNumberMax += arg->printSolY;
    }

    // compute x position of the printed central value
    double xNumberMin, xNumberMax;
    if (scanners[i]->getDrawSolution() == 1) {
      xNumberMin = xCentral + (xmax - xmin) * 0.20;  // draw at central value
      xNumberMax = xCentral + (xmax - xmin) * 0.0;
    } else if (scanners[i]->getDrawSolution() == 2) {
      xNumberMin = xCLmin + (xmax - xmin) * 0.0;  // draw at left CL boundary
      xNumberMax = xCLmin + (xmax - xmin) * 0.25;
    } else if (scanners[i]->getDrawSolution() == 3) {
      xNumberMin = xCLmax + (xmax - xmin) * 0.0;  // draw at right CL boundary
      xNumberMax = xCLmax + (xmax - xmin) * 0.25;
    } else if (scanners[i]->getDrawSolution() == 4) {
      xNumberMin = xCLmin + (xmax - xmin) * -0.20;  // draw a little left of the left CL boundary
      xNumberMax = xCLmin + (xmax - xmin) * 0.0;
    } else {
      std::cout << "OneMinusClPlot::drawSolutions() : ERROR : --ps code ";
      std::cout << scanners[i]->getDrawSolution() << " not found! Use [0,1,2,3]." << std::endl;
      continue;
    }

    // move number a bit to the left so it doesn't cover the right plot border anymore
    if (arg->isQuickhack(4)) {
      xNumberMin -= (xmax - xmin) * 0.225;
      xNumberMax -= (xmax - xmin) * 0.225;
    }

    // If print solution argument is given then overwrite
    if (arg->printSolX > 0.) {
      double diff = xNumberMax - xNumberMin;
      xNumberMin = arg->printSolX;
      xNumberMax = arg->printSolX + diff;
    }

    auto t1 = makeOwnedTObject<TPaveText>(xNumberMin, yNumberMin, xNumberMax, yNumberMax, "BR");
    t1->SetBorderSize(0);
    t1->SetFillStyle(0);
    t1->SetTextAlign(13);
    t1->SetTextFont(font);
    t1->SetTextColor(color);
    t1->SetTextSize(labelsize);
    if (arg->isQuickhack(32)) t1->SetTextSize(1.5 * labelsize);
    if (Utils::isAngle(scanners[i]->getScanVar1())) {
      using Utils::RadToDeg;
      xCentral = RadToDeg(xCentral);
      xCLmin = RadToDeg(xCLmin);
      xCLmax = RadToDeg(xCLmax);
    }
    Rounder myRounder(arg, xCLmin, xCLmax, xCentral);
    int d = myRounder.getNsubdigits();
    double xCentralRd = myRounder.central();
    if (arg->isQuickhack(3)) xCentralRd += 180.;  ///< see documentation of --qh option in OptParser.cpp
    t1->AddText(Form("%.*f^{+%.*f}_{#font[122]{-}%.*f}", d, xCentralRd, d, myRounder.errPos(), d, myRounder.errNeg()));

    if (!arg->isQuickhack(21)) t1->Draw();
    iDrawn += 1;
  }
}

///
/// Draw a horizontal line at given p-value, put a label on top of it stating the corresponding CL.
///
void OneMinusClPlot::drawCLguideLine(const double pvalue) {
  canvas->cd();
  canvas->Update();
  double xmin = gPad->GetUxmin();
  double xmax = gPad->GetUxmax();

  double labelPos = xmin + (xmax - xmin) * 0.10;
  if (arg->isQuickhack(2)) labelPos = xmin + (xmax - xmin) * 0.55;
  if (arg->isQuickhack(23)) labelPos = xmin + (xmax - xmin) * 0.8;
  if (arg->isQuickhack(31)) labelPos = xmin + (xmax - xmin) * 0.01;

  if (arg->CL.size() > 1) {
    std::ranges::sort(arg->CL);
    for (int i = 0; i < arg->CL.size(); i++) {
      if (abs((1 - pvalue) - arg->CL[i] / 100.) < 0.0001 && abs(arg->CL[i] - arg->CL[i - 1]) < 8) {
        if (!arg->isQuickhack(23))
          labelPos = labelPos + (xmax - xmin) * 0.15;
        else
          labelPos = labelPos - (xmax - xmin) * 0.15;
      }
    }
  }

  double labelPosYmin = 0;
  double labelPosYmax = 0;

  if (arg->plotlog) {
    labelPosYmin = pvalue;
    labelPosYmax = labelPosYmin * 2.;
  } else {
    labelPosYmin = pvalue + 0.02;
    labelPosYmax = labelPosYmin + 0.05;
  }

  auto t = makeOwnedTObject<TPaveText>(labelPos, labelPosYmin, labelPos + (xmax - xmin) * 0.5, labelPosYmax, "BR");
  t->SetBorderSize(0);
  t->SetFillStyle(0);
  t->SetTextAlign(12);
  t->SetTextFont(font);
  t->SetTextSize(labelsize);
  t->AddText(Form("%.1f%%", (1. - pvalue) * 100.));
  t->Draw();

  auto l = makeOwnedTObject<TLine>(xmin, pvalue, xmax, pvalue);
  l->SetLineWidth(1);
  l->SetLineColor(kBlack);
  l->SetLineStyle(kDotted);
  l->Draw();
}

///
/// Draw 1, 2, and 3 sigma lines.
///
void OneMinusClPlot::drawCLguideLines() {
  if (arg->CL.size() == 0) {
    drawCLguideLine(0.31731);
    drawCLguideLine(4.550026e-2);
    if (arg->plotlog) {
      drawCLguideLine(2.7e-3);
      if (arg->plotymin < 6.3e-5) { drawCLguideLine(6.3e-5); }
    }
  }
  if (arg->CL.size() > 0) {
    for (auto level : arg->CL) {
      if (level < 99) {
        drawCLguideLine(1. - level / 100.);
      } else if (arg->plotlog) {
        if (arg->plotymin > 6.3e-5 && level < 99.9937) { continue; }
        drawCLguideLine(1. - level / 100.);
      }
    }
  }
}

/**
 * Make the plot for the 1D scan.
 *
 * @param beautify Set to true to use a simpler plot function which directly plots the 1-CL histograms without
 *                 beautification rather than TGraphs.
 */
void OneMinusClPlot::Draw(const bool beautify) {

  // put this in for exponent xaxes
  if (!arg->isQuickhack(30)) canvas->SetRightMargin(0.1);
  if (arg->plotlog) {
    canvas->SetLogy();
    if (!this->name.EndsWith("_log")) this->name = this->name + "_log";
  }
  canvas->cd();

  // Plot the CLs
  for (int i = 0; i < scanners.size(); i++)
    if (do_CLs[i] == 2) scan1dCLsPlot(scanners[i].get(), arg->nsmooth);
  canvas->Update();
  canvas->Show();

  // Legend: make the legend short, the text will extend over the boundary, but the symbol will be shorter
  const double legendXmin = arg->plotlegx != -1. ? arg->plotlegx : 0.19;
  const double legendYmin = arg->plotlegy != -1. ? arg->plotlegy : 0.78;
  const double legendXmax = legendXmin + (arg->plotlegsizex != -1. ? arg->plotlegsizex : 0.31);
  const double legendYmax = legendYmin + (arg->plotlegsizey != -1. ? arg->plotlegsizey : 0.1640559);
  auto leg = makeOwnedTObject<TLegend>(legendXmin, legendYmin, legendXmax, legendYmax);
  leg->SetNColumns(arg->plotlegcols);
  leg->SetFillColorAlpha(kWhite, 1.);
  leg->SetFillStyle(0);
  leg->SetLineColor(kWhite);
  leg->SetBorderSize(0);
  if (arg->isQuickhack(35)) {
    leg->SetFillColorAlpha(0, 0.4);
    leg->SetFillStyle(1001);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetLineColor(kGray + 1);
    leg->SetBorderSize(1);
  }
  leg->SetTextFont(font);
  leg->SetTextSize(legendsize * 0.75);

  std::vector<TString> legTitles;
  for (int i = 0; i < scanners.size(); i++) {
    TString legDrawOption = "f";
    if (plotPluginMarkers &&
        (scanners[i]->getMethodName() == "Plugin" || scanners[i]->getMethodName() == "BergerBoos" ||
         scanners[i]->getMethodName() == "DatasetsPlugin")) {
      legDrawOption = "lep";
    }
    if (arg->plotlegstyle != "default") legDrawOption = arg->plotlegstyle;

    TString legTitle = scanners[i]->getTitle();
    if (legTitle == "default") {
      if (scanners[i]->getMethodName().Contains("Prob")) legTitle = do_CLs[i] ? "Prob CLs" : "Prob";
      if (scanners[i]->getMethodName().Contains("Plugin")) {
        if (do_CLs[i] == 0)
          legTitle = "Plugin";
        else if (do_CLs[i] == 1)
          legTitle = "Simplified CLs";
        else if (do_CLs[i] == 2)
          legTitle = "Standard CLs";
      }
    } else if (!arg->isQuickhack(29)) {
      if (scanners[i]->getMethodName().Contains("Prob")) legTitle += do_CLs[i] ? " (Prob CLs)" : " (Prob)";
      if (scanners[i]->getMethodName().Contains("Plugin")) {
        if (do_CLs[i] == 0)
          legTitle += " (Plugin)";
        else if (do_CLs[i] == 1)
          legTitle += " (Simplified CLs)";
        else if (do_CLs[i] == 2)
          legTitle += " (Standard CLs)";
      }
    }
    legTitles.push_back(legTitle);

    if (!beautify) {
      scan1dPlotSimple(scanners[i].get(), i == 0, do_CLs[i]);
      leg->AddEntry(scanners[i]->getHCL(), legTitle, legDrawOption);
    } else {
      if (scanners[i]->getFillStyle() != 0 || scanners[i]->getFillColor() != 0) {
        const auto g = scan1dPlot(scanners[i].get(), i == 0, false, scanners[i]->getFilled(), do_CLs[i]);
        if (legTitles[i] != "noleg") leg->AddEntry(g, legTitle, legDrawOption);
      }
    }
  }

  // lines only
  if (beautify) {
    for (int i = 0; i < scanners.size(); i++) {
      bool last = i == scanners.size() - 1;
      const auto g = scan1dPlot(scanners[i].get(), false, last, false, do_CLs[i]);
      if (scanners[i]->getFillStyle() == 0 && scanners[i]->getFillColor() == 0 && legTitles[i] != "noleg")
        leg->AddEntry(g, legTitles[i], "L");
    }
  }
  drawSolutions();
  if (plotLegend) leg->Draw();
  if (arg->isQuickhack(22)) leg->Draw();
  canvas->Update();
  if (!arg->isQuickhack(34)) drawCLguideLines();

  // Draw the logo
  double yGroup = 0.6;
  if (plotLegend) {
    yGroup = (arg->plotlog) ? 0.775 : 0.6;
  } else {
    yGroup = (arg->plotlog) ? 0.3 : 0.775;
  }
  drawGroup(yGroup);

  canvas->Update();
  canvas->Show();
}
