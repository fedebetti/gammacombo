/*
 * Gamma Combination
 * Author: Maximilian Schlupp, maximilian.schlupp@cern.ch
 * Date: January 2013
 */

#ifndef MethodBergerBoosScan_h
#define MethodBergerBoosScan_h

#include "MethodPluginScan.h"

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TTree.h>

#include <memory>

class MethodBergerBoosScan : public MethodPluginScan {
 public:
  MethodBergerBoosScan(MethodProbScan* s, TString d = "XX");
  ~MethodBergerBoosScan();
  std::unique_ptr<TH2> calcPValues(const TH2* better, const TH2* all, const TH2* bg);
  void getBestPValue(TH1* hCL, TH2* pValues);
  int getNBergerBoosPointsPerScanpoint() const { return nBBPoints; };  ///< Return number of BB points per scan point
  void readScan1dTrees(int runMin, int runMax, TString fname) override;
  int scan1d(const int nRun) override;
  inline void setNBergerBoosPointsPerScanpoint(int n) { nBBPoints = n; };  ///< Set number of BB points per scan point
  ///< Set number of Berger Boos points drawn from
  ///< the BergerBoos CL intervals defined in the workspace
  void setNewBergerBoosPoint(int m);  ///< Samples one Berger Boos point
                                      ///< Draws ALL parameters randomly according to their
                                      ///< BergerBoos ranges defined in the workspace
                                      ///< The parameter of interest (scan parameter) has to be set
                                      ///< constant (and to its scan point value) after the function call.

  /**
   * Draw a 2D histogram showing the BB points in varX-varY space.
   *
   * @param save Boolean flag that specifies if a copy of the plot will be saved in the plots folder.
   */
  void drawBBPoints(TString varX, TString varY, int runMin = 1, int runMax = 1, bool save = true) const;

 private:
  int nBBPoints = 1;  ///< number of sampled Berger Boos points per scan point
  std::unique_ptr<TFile> file;
  TTree* BBtree = nullptr;
  TString dir;
};

#endif
