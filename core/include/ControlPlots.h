/*
 * GammaCombo
 * Author: Till Moritz Karbach, moritz.karbach@cern.ch
 * Date: Dec 2014
 *
 */

#ifndef ControlPlots_h
#define ControlPlots_h

#include <TCut.h>
#include <TString.h>

#include <vector>

class OptParser;
class MethodProbScan;
class ToyTree;

class TCanvas;
class TTree;
class TVirtualPad;

///
/// Class to make control plots of Plugin toys.
///
class ControlPlots {
 public:
  ControlPlots(ToyTree* tt);

  void ctrlPlotChi2Distribution();
  void ctrlPlotChi2Parabola();
  void ctrlPlotNuisances();
  void ctrlPlotObservables();
  void ctrlPlotChi2();
  void ctrlPlotPvalue();
  void ctrlPlotMore(MethodProbScan* profileLH);
  void saveCtrlPlots();

 private:
  void makePlotsNice(TString htemp = "htemp", TString Graph = "Graph");
  TCanvas* selectNewCanvas(TString title);
  TVirtualPad* selectNewPad();
  void updateCurrentCanvas();

  TString name;                            ///< combiner name, ending up in titles and file names
  ToyTree* tt = nullptr;                   ///< the toy tree
  TTree* t = nullptr;                      ///< the tree
  const OptParser* arg = nullptr;          ///< command line arguments
  std::vector<TCanvas*> ctrlPlotCanvases;  ///< Pointers to the canvases of the control plots, see selectNewCanvas().
  int ctrlPadId = 0;                       ///< ID of currently selected pad, see selectNewPad().
  TCut ctrlPlotCuts;                       ///< Cuts that are applied to all control plots.
};

#endif
