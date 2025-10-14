/**
 * Gamma Combination
 * Author: Till Moritz Karbach, moritz.karbach@cern.ch
 * Date: August 2014
 *
 **/

#ifndef ProgressBar_h
#define ProgressBar_h

class OptParser;

///
/// Class showing a progress bar.
///
class ProgressBar {
 public:
  ProgressBar(const OptParser* arg, int n);

  void progress();
  void skipSteps(int n);

 private:
  void progressBar() const;
  void progressPercentage() const;

  int _n = 100;                ///< maximum number of steps, "100%"
  int _x = 0;                  ///< current step, "78%"
  const int _width = 50;       ///< width of the progress bar
  const int _resolution = 50;  ///< update the display this many times
  bool _batch = false;         ///< display progress in a log-file compatible way
};

#endif
