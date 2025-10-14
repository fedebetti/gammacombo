#include <ProgressBar.h>

#include <OptParser.h>

#include <cassert>
#include <format>
#include <iomanip>
#include <iostream>

ProgressBar::ProgressBar(const OptParser* arg, int n) : _n(n) {
  assert(arg);
  _batch = arg->isAction("pluginbatch") || arg->isAction("bbbatch");
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
void ProgressBar::progressBar() const {
  const auto ratio = static_cast<double>(_x) / _n;
  const auto c = static_cast<int>(ratio * _width);

  std::cout << std::format("{:3d}% [{:s}{:s}]", static_cast<int>(ratio * 100.), std::string(c, '='),
                           std::string(_width - c, ' '))
            << "\r" << std::flush;
  // This alternative would leave the progress bar in the output once the loop is done:
  // if ( _x<_n ) std::cout << "\r" << std::flush;
  // else std::cout << std::endl;
  if (_x == _n) {
    int offset = TString("100% []").Sizeof();
    std::cout << std::string(c + offset, ' ') << "\r" << std::flush;
  }
}

///
/// Display the progress as percentage without flushing, so that the logfiles of batch jobs look good.
///
void ProgressBar::progressPercentage() const {
  auto ratio = static_cast<double>(_x) / _n;
  std::cout << std::format("{:3d}%", static_cast<int>(ratio * 100.)) << std::endl;
}

///
/// Increase the internal count by n steps.
///
void ProgressBar::skipSteps(int n) { _x += n; }
