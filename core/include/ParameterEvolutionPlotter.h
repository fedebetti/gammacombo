/**
 * Gamma Combination
 * Author: Till Moritz Karbach, moritz.karbach@cern.ch
 * Date: Nov 2014
 *
 * Class that plots the parameter evolution of the nuisance parameters
 * from a prob scan.
 *
 **/

#ifndef ParameterEvolutionPlotter_h
#define ParameterEvolutionPlotter_h

#include <TString.h>

#include <vector>

class MethodProbScan;
class OptParser;
class RooSlimFitResult;

class RooWorkspace;

class TCanvas;
class TGraph;
class TGraphErrors;
class TVirtualPad;

class ParameterEvolutionPlotter {
 public:
  ParameterEvolutionPlotter(MethodProbScan* scanner);
  ~ParameterEvolutionPlotter();

  ParameterEvolutionPlotter(ParameterEvolutionPlotter&) = delete;
  ParameterEvolutionPlotter& operator=(const ParameterEvolutionPlotter&) = delete;

  void plotParEvolution();
  void plotObsScanCheck();

 private:
  void getLocalMinPositions();
  void drawLinesAtMinima(TVirtualPad* pad);
  void drawVerticalRedLine(TVirtualPad* pad, double xpos);
  TGraph* makeChi2Graph(const std::vector<RooSlimFitResult*>& results) const;
  TGraph* makeEvolutionGraph(const std::vector<RooSlimFitResult*>& results, TString parName) const;
  TGraphErrors* makeEvolutionGraphErrors(const std::vector<RooSlimFitResult*>& results, TString parName) const;
  void saveEvolutionPlots();
  TCanvas* selectNewCanvas(TString title);
  TVirtualPad* selectNewPad();
  void updateCurrentCanvas();

  const OptParser* arg = nullptr;               ///< command line arguments
  RooWorkspace* w = nullptr;                    ///< a clone of the scanner's workspace
  std::vector<RooSlimFitResult*> allResults;    ///< all results of all scan points
  std::vector<RooSlimFitResult*> curveResults;  ///< only the results of scan points that were accepted into the CL
                                                ///< curve
  TString title;                                ///< canvas title
  TString name;                                 ///< scanner name, part of the file name of the plots
  TString parsName;                             ///< name of parameter set inside the workspace
  TString obsName;                              ///< name of observables set inside the workspace
  TString scanVar1;                             ///< name of the can variable
  std::vector<double> m_localMinPositions;      ///< positions of the local minima in scan steps
  std::vector<TCanvas*> m_canvases;             ///< Pointers to the canvases of the plots, see selectNewCanvas().
  int m_padId = 0;                              ///< ID of currently selected pad, see selectNewPad().
};

#endif
