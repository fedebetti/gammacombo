/**
 * Gamma Combination
 * Author: Till Moritz Karbach, moritz.karbach@cern.ch
 * Date: August 2012
 *
 **/

#ifndef OptParser_h
#define OptParser_h

#include <iostream>
#include <stdlib.h>

#include <TRegexp.h>

#include "Utils.h"
#include "tclap/CmdLine.h"

using namespace TCLAP;

class OptParser {
 public:
  OptParser();
  ~OptParser(){};

  void bookOption(TString opt);
  void bookAllOptions();
  void bookPlottingOptions();
  void bookPluginOptions();
  void bookProbOptions();
  void bookFlowcontrolOptions();
  void parseArguments(int argc, char* argv[]);
  bool isAction(TString s) const;
  bool isAsimovCombiner(int id) const;
  bool isQuickhack(int id) const;

  std::vector<TString> action;
  std::vector<int> asimov;
  std::vector<TString> asimovfile;
  bool cacheStartingValues;
  std::vector<float> CL;
  std::vector<int> cls;
  std::vector<int> color;
  std::vector<int> combid;
  std::vector<std::vector<int>> combmodifications;  // encodes requested modifications to the combiner ID through the -c
                                                    // 26:+12 syntax,format is [cmbid:[+pdf1,-pdf2,...]]
  bool compare = false;
  bool confirmsols = true;
  bool controlplot = false;
  int coverageCorrectionID = 0;
  int coverageCorrectionPoint = 0;
  bool debug = false;
  int digits = -99;
  bool enforcePhysRange = false;
  std::vector<int> fillstyle;
  std::vector<int> fillcolor;
  std::vector<float> filltransparency;
  std::vector<int> linewidth;
  std::vector<int> linecolor;
  std::vector<int> linestyle;
  std::vector<std::string> hexfillcolor;
  std::vector<std::string> hexlinecolor;
  TString filenamechange = "";
  TString filenameaddition;
  std::vector<std::vector<Utils::FixPar>> fixParameters;
  std::vector<std::vector<Utils::StartPar>> startVals;
  std::vector<std::vector<Utils::RangePar>> physRanges;
  std::vector<std::vector<TString>> removeRanges;
  std::vector<std::vector<TString>> randomizeToyVars;
  bool grid = false;
  TString group = "GammaCombo";
  TString groupPos = "";
  TString hfagLabel = "";
  TString hfagLabelPos = "";
  int id = -99;
  bool importance = false;
  bool info = false;
  bool interactive = false;
  std::vector<int> jmax;
  std::vector<int> jmin;
  TString jobdir = ".";
  bool largest = false;
  bool latex = false;
  std::vector<TString> loadParamsFile;
  bool lightfiles = false;
  int batchstartn = 1;
  bool batcheos = false;
  bool batchsubmit = false;
  TString batchout = "";
  TString batchreqs = "";
  int nbatchjobs = -99;
  int nBBpoints;
  int ndiv = 407;
  int ndivy = 407;
  bool nosyst = false;
  int npoints1d = -99;
  int npoints2dx = -99;
  int npoints2dy = -99;
  int npointstoy = -99;
  int ncoveragetoys = -99;
  int nrun = -99;
  int ntoys = -99;
  int nsmooth = 1;
  TString parsavefile;
  bool parevol = false;
  std::vector<int> pevid;
  std::vector<int> plot2dcl;
  TString plotdate = "";
  TString plotext = "";
  int plotid = -99;
  bool plotlog = false;
  bool plotlegend = true;
  float plotlegx = -99;
  float plotlegy = -99;
  float plotlegsizex = -99;
  float plotlegsizey = -99;
  TString plotlegstyle = "default";
  int plotlegcols = 1;
  bool plotlegbox = false;
  float plotlegboxx = -99;
  float plotlegboxy = -99;
  float plotgroupx = -99;
  float plotgroupy = -99;
  Double_t plotHFAGLabelPosX = 0;
  Double_t plotHFAGLabelPosY = 0;
  Double_t plotHFAGLabelScale = 1;
  bool plotmagnetic = false;
  int plotnsigmacont = 2;
  std::map<int, std::vector<int>> contourlabels;
  bool plotpluginonly = false;
  bool plotpulls = false;
  bool plotprelim = false;
  float plotoriginx = -99.;
  float plotoriginy = -99.;
  std::vector<int> plotsolutions;
  std::vector<int> plotsoln;
  bool plotunoff = false;
  float plotymin = -99.;
  float plotymax = -99.;
  bool intprob = false;
  float pluginPlotRangeMin = -100;
  float pluginPlotRangeMax = -100;
  bool probforce = false;
  bool probimprove = false;
  TString probScanResult = "notSet";
  bool printcor = false;
  float printSolX = -999.;
  float printSolY = -999.;
  std::vector<int> qh;
  /*TString         queue;*/
  std::vector<std::vector<TString>> readfromfile;
  std::vector<TString> relation;
  bool runCLs;
  TString save = "";
  bool saveAtMin = false;
  std::vector<float> savenuisances1d;
  std::vector<float> savenuisances2dx;
  std::vector<float> savenuisances2dy;
  bool scanforce = false;
  float scanrangeMin = -101;
  float scanrangeMax = -101;
  float scanrangeyMin = -102;
  float scanrangeyMax = -102;
  float scaleerr = -999.;
  float scalestaterr = -999.;
  bool smooth2d = false;
  bool square = false;
  int teststatistic = 2;
  std::vector<TString> title;
  TString xtitle;
  TString ytitle;
  TString toyFiles = "";
  int updateFreq = 10;
  bool usage = false;
  std::vector<TString> var;
  bool verbose = false;

  CmdLine cmd;

 private:
  int convertToDigitWithCheck(TString parseMe, TString usage) const;
  int convertToIntWithCheck(TString parseMe, TString usage) const;
  void defineOptions();
  void parsePosition(TString parseMe, float& x, float& y, TString usage);
  void parsePositionAndScale(TString parseMe, Double_t& x, Double_t& y, Double_t& scale, TString usage);
  bool parseRange(TString parseMe, float& min, float& max) const;
  bool parseAssignment(TString parseMe, TString& name, TString& value) const;
  bool parseAssignment(TString parseMe, TString& name, float& value) const;
  void parseCombinerString(TString parseMe, int& resultCmbId, std::vector<int>& resultAddDelPdf) const;
  std::vector<TString> availableOptions;
  std::vector<TString> bookedOptions;
};

#endif
