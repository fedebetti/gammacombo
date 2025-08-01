#include <ProgressBar.h>

#include <format>
#include <iomanip>

ProgressBar::ProgressBar(const OptParser* arg, const int n) {
  assert(arg);
  _arg = arg;
  _n = n;
  _resolution = _width;
  _batch = _arg->isAction("pluginbatch") || _arg->isAction("bbbatch");
}

///
/// Call this from inside the loop.
///
void ProgressBar::progress() {
  _x++;
  if ((_x != _n) && (_x % (_n / _resolution + 1) != 0)) return;
  if (_batch)
    progressPercentage();
  else
    progressBar();
}

///
/// Display the progress as a bar:
/// 100% [==================================================]
///
void ProgressBar::progressBar() {
  const auto ratio = static_cast<double>(_x) / _n;
  const auto c = static_cast<int>(ratio * _width);

  std::cout << std::format("{:3d}% [", static_cast<int>(ratio * 100));
  for (int x = 0; x < c; x++) std::cout << "=";
  for (int x = c; x < _width; x++) std::cout << " ";
  std::cout << "]";
  // this code leaves the progress bar in the output once the loop is done
  // if ( _x<_n ) std::cout << "\r" << std::flush;
  // else std::cout << std::endl;
  // and this code removes it
  std::cout << "\r" << std::flush;
  if (_x == _n) {
    int offset = TString("100% []").Sizeof();
    for (int x = 0; x < c + offset; x++) std::cout << " ";
    std::cout << "\r" << std::flush;
  }
}

///
/// Display the progress as percentage without flushing, so that the logfiles of batch jobs look good.
///
void ProgressBar::progressPercentage() {
  const auto ratio = static_cast<double>(_x) / _n;
  std::cout << std::format("{:3d}%", static_cast<int>(ratio * 100)) << std::endl;
}

///
/// Increase the internal count by n steps.
///
void ProgressBar::skipSteps(const int n) { _x += n; }
