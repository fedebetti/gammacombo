/**
 * Gamma Combination
 * Author: Till Moritz Karbach, moritz.karbach@cern.ch
 * Date: Oct 2014
 *
 **/

#include <GammaComboEngine.h>

#include <PDF_Circle.h>
#include <PDF_CrossCor_GausA_vs_GausB.h>
#include <PDF_Gaus.h>
#include <PDF_Gaus2d.h>
#include <PDF_GausB.h>

int main(int argc, char* argv[]) {
  GammaComboEngine gc("Lb2pktaul", argc, argv);

  ///////////////////////////////////////////////////
  //
  // define PDFs
  //
  ///////////////////////////////////////////////////

  gc.addPdf(1, new PDF_Gaus("year2013", "year2013", "year2013"), "1D Gaussian (a_{obs} = -0.5)");
  gc.addPdf(2, new PDF_Gaus("year2014", "year2014", "year2014"), "1D Gaussian (a_{obs} =  1.5)");
  gc.addPdf(3, new PDF_Gaus2d("year2013", "year2013", "year2013"), "2D Gaussian (a_{obs}, b_{obs})");
  gc.addPdf(4, new PDF_Circle("year2013", "year2013", "year2013"), "circle (a_{obs}, b_{obs})");
  gc.addPdf(5, new PDF_GausB("year2013", "year2013", "year2013"), "1D Gaussian (b_{obs} = 1.0)");
  gc.addPdf(6, new PDF_CrossCor_GausA_vs_GausB(gc[1], gc[5], "year2013"), "Correlation");

  ///////////////////////////////////////////////////
  //
  // Define combinations
  //
  ///////////////////////////////////////////////////

  gc.newCombiner(0, "empty", "empty");
  gc.newCombiner(1, "Lb2pktaul1", "Gaus 1", 1);
  gc.newCombiner(2, "Lb2pktaul2", "Gaus 2", 2);
  gc.newCombiner(3, "Lb2pktaul3", "Gaus 1 & Gaus 2", 1, 2);
  gc.newCombiner(4, "Lb2pktaul4", "2D Gaus", 3);
  gc.newCombiner(5, "Lb2pktaul5", "2D Gaus & Gaus 1", 2, 3);
  gc.newCombiner(6, "Lb2pktaul6", "Circle", 4);
  gc.newCombiner(7, "Lb2pktaul7", "2D Gaus & Circle", 3, 4);
  gc.newCombiner(8, "Lb2pktaul8", "Gaus 1 & Gaus 2 & 2D Gaus", 1, 2, 3);

  ///////////////////////////////////////////////////
  //
  // Run
  //
  ///////////////////////////////////////////////////

  gc.run();
}
