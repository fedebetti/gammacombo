#ifndef OneMinusClPlotAbs_h
#define OneMinusClPlotAbs_h

#include "MethodAbsScan.h"
#include "OptParser.h"

#include <concepts>
#include <memory>
#include <utility>
#include <vector>

#include <TCanvas.h>
#include <TObject.h>
#include <TString.h>

/// Concept for objects deriving from TObject.
template <typename T>
concept RootTObject = std::derived_from<T, TObject> == true;

class MethodAbsScan;

/**
 * Abstract class to plot the CL distributions for fit variables.
 */
class OneMinusClPlotAbs {
 public:
  OneMinusClPlotAbs(OptParser* arg, TString name = "c1", TString title = "c1");
  virtual ~OneMinusClPlotAbs() = default;

  virtual void addScanner(MethodAbsScan* s, const int CLsType = 0);
  inline void disableLegend(bool yesNo = false) { plotLegend = yesNo; };
  virtual void Draw(const bool beautify = true) = 0;
  void drawGroup(const double yPos = 0.6);
  virtual void drawSolutions() = 0;
  inline TString getName() const { return name; };
  void save() const;
  inline void setFont(int fnum) { font = fnum; };
  inline void setLabelSize(int lnum) { labelsize = lnum; };
  inline void setPlotLabel(TString& lname) { label = lname; };
  inline void Show() const { canvas->Show(); };

  /**
   * Return a non-owning pointer to a ROOT object of type T.
   *
   * The deallocation of the object will happen at destruction of the OneMinusClPlotAbs object.
   * (This is needed due to the ugly memory model of ROOT).
   */
  template <RootTObject T, typename... Args>
  T* makeOwnedTObject(Args... args) {
    rootObjects.emplace_back(std::unique_ptr<TObject>(std::make_unique<T>(args...)));
    return dynamic_cast<T*>(rootObjects.back().get());
  }

  /// Assign the ownership of a ROOT object to the this isntance of OneMinusClPlotsAbs.
  template <RootTObject T>
  T* getTObjectOwnership(std::unique_ptr<T> obj) {
    rootObjects.emplace_back(std::unique_ptr<TObject>(std::move(obj)));
    return dynamic_cast<T*>(rootObjects.back().get());
  }

 protected:
  int font = 133;       ///< font code. The last digit disables scaling with the canvas size.
  int labelsize = 35;   ///< text size of axis labels, numeric solutions, CL guide lines (in pixels)
  int titlesize = 45;   ///< text size of axis titles, group label, "Prliminary" is x0.75 (in pixels)
  int legendsize = 29;  ///< text size of legend entries in 1d and 2d plots (in pixels)

  std::vector<MethodAbsScan*> scanners;
  std::vector<int> do_CLs;   ///< vector, which stores the cls method type to be plotted
  OptParser* arg = nullptr;  ///< command line options
  std::unique_ptr<TCanvas> canvas;
  TString name;
  TString title;
  TString label;
  bool plotLegend = true;

 private:
  std::vector<std::unique_ptr<TObject>> rootObjects;
};

#endif
