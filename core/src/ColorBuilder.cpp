#include <ColorBuilder.h>

#include <TColor.h>
#include <TObjArray.h>
#include <TROOT.h>

int ColorBuilder::darkcolor(int n) const { return darklightcolor(n, 0.95); }

int ColorBuilder::lightcolor(int n) const { return darklightcolor(n, 1.04); }

///
/// Copied from TColor::GetColorDark(Int_t n), but customized
/// the 'darkness'.
///
int ColorBuilder::darklightcolor(int n, float scale) const {
  if (n < 0) return -1;

  // Get list of all defined colors
  TObjArray* colors = (TObjArray*)gROOT->GetListOfColors();
  Int_t ncolors = colors->GetSize();
  // Get existing color at index n
  TColor* color = nullptr;
  if (n < ncolors) color = (TColor*)colors->At(n);
  if (!color) return -1;

  // Get the rgb of the the new dark color corresponding to color n
  Float_t r, g, b;
  TColor::HLStoRGB(color->GetHue(), scale * color->GetLight(), color->GetSaturation(), r, g, b);

  // Build the dark color (unless the slot nd is already used)
  Int_t nd = scale < 1. ? n + 100 : n + 150;
  TColor* colord = nullptr;
  if (nd < ncolors) colord = (TColor*)colors->At(nd);
  if (colord) return nd;
  colord = new TColor(nd, r, g, b);
  colord->SetName(scale < 1. ? Form("%s_dark", color->GetName()) : Form("%s_light", color->GetName()));
  colors->AddAtAndExpand(colord, nd);
  return nd;
}
