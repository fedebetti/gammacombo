#ifndef ParameterEvolutionPlotter_h
#define ParameterEvolutionPlotter_h

#include <memory>

#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TString.h>
#include <TVirtualPad.h>

#include <RooWorkspace.h>

#include "MethodProbScan.h"

/**
 * Class that plots the parameter evolution of the nuisance parameters from a prob scan.
 **/
class ParameterEvolutionPlotter {
 public:
  ParameterEvolutionPlotter(MethodProbScan* scanner);

  void plotParEvolution();
  void plotObsScanCheck();

 private:
  void getLocalMinPositions();
  void drawLinesAtMinima(TVirtualPad* pad);
  void drawVerticalRedLine(TVirtualPad* pad, const double xpos);
  TGraph* makeChi2Graph(std::vector<RooSlimFitResult*> results);
  TGraph* makeEvolutionGraph(std::vector<RooSlimFitResult*> results, const TString parName);
  TGraphErrors* makeEvolutionGraphErrors(std::vector<RooSlimFitResult*> results, const TString parName);
  void saveEvolutionPlots();
  TCanvas* selectNewCanvas(TString title);
  TVirtualPad* selectNewPad();
  void updateCurrentCanvas();

  const OptParser* arg;                              ///< command line arguments
  std::unique_ptr<RooWorkspace> w;                   ///< a clone of the scanner's workspace
  std::vector<RooSlimFitResult*> allResults;         ///< all results of all scan points
  std::vector<RooSlimFitResult*> curveResults;       ///< only the results of scan points that were accepted into the CL
                                                     ///< curve
  TString title;                                     ///< canvas title
  TString name;                                      ///< scanner name, part of the file name of the plots
  TString parsName;                                  ///< name of parameter set inside the workspace
  TString obsName;                                   ///< name of observables set inside the workspace
  TString scanVar1;                                  ///< name of the can variable
  std::vector<double> m_localMinPositions;           ///< positions of the local minima in scan steps
  std::vector<std::unique_ptr<TCanvas>> m_canvases;  ///< Pointers to the canvases of the plots, see selectNewCanvas().
  int m_padId = 0;                                   ///< ID of currently selected pad, see selectNewPad().
};

#endif
