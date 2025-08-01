#ifndef ProgressBar_h
#define ProgressBar_h

#include "OptParser.h"

///
/// Class showing a progress bar.
///
class ProgressBar {
 public:
  ProgressBar(const OptParser* arg, const int n);

  void progress();
  void skipSteps(const int n);

 private:
  void progressBar();
  void progressPercentage();

  const OptParser* _arg;  ///< command line arguments
  int _n;                 ///< maximum number of steps, "100%"
  int _x = 0;             ///< current step, "78%"
  int _width = 50;        ///< width of the progress bar
  int _resolution;        ///< update the display this many times
  bool _batch;            ///< display progress in a log-file compatible way
};

#endif
