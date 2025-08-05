/**
 * \brief Tools to handle TGraphs
 * \author Till Moritz Karbach, moritz.karbach@cern.ch
 * \date March 2015
 *
 **/

#ifndef TGraphTools_h
#define TGraphTools_h

class TGraph;

class TGraphTools {
 public:
  TGraph* addPointToGraphAtFirstMatchingX(TGraph* g, float xNew, float yNew);
};

#endif
