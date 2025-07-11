#ifndef TGraphTools_h
#define TGraphTools_h

#include <TGraph.h>

#include <memory>

class TGraphTools {
 public:
  std::unique_ptr<TGraph> addPointToGraphAtFirstMatchingX(const TGraph* g, double xNew, double yNew);
};

#endif
