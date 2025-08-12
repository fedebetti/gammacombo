/**
 * Gamma Combination
 * Author: Till Moritz Karbach, moritz.karbach@cern.ch
 * Date: July 2014
 *
 **/

#include <PDF_CrossCorAbs.h>

#include <ParametersTutorial.h>

#include <PDF_Abs.h>
#include <RooCrossCorPdf.h>
#include <UtilsConfig.h>

#include <RooArgList.h>
#include <RooRealVar.h>

PDF_CrossCorAbs::PDF_CrossCorAbs(PDF_Abs* pdf1, PDF_Abs* pdf2) : PDF_Abs(pdf1->getNobs() + pdf2->getNobs()) {
  name = "CrossCorAbs";
  assert(pdf1);
  assert(pdf2);
  this->pdf1 = pdf1;
  this->pdf2 = pdf2;
  title = "cross correlation: [" + pdf1->getTitle() + "] vs [" + pdf2->getTitle() + "]";
  m_isCrossCorPdf = true;
  nObs1 = pdf1->getNobs();
  nObs2 = pdf2->getNobs();
  initParameters();
  initRelations();
  initObservables();
}

void PDF_CrossCorAbs::initParameters() {
  ParametersTutorial p;
  parameters = new RooArgList("parameters");
  // we need the same parameters as both input PDFs:
  // copy over from first PDF
  for (const auto& par : *pdf1->getParameters()) { parameters->add(*(p.get(par->GetName()))); }
  // copy over from second PDF
  for (const auto& par : *pdf2->getParameters()) {
    if (parameters->find(par->GetName())) continue;
    parameters->add(*(p.get(par->GetName())));
  }
}

void PDF_CrossCorAbs::initRelations() {
  theory = new RooArgList("theory");
  // we need the same theory as both input PDFs:
  // copy over from first PDF
  for (const auto& th : *pdf1->getTheory()) { theory->add(*th); }
  // copy over from second PDF
  for (const auto& th : *pdf2->getTheory()) { theory->add(*th); }
}

void PDF_CrossCorAbs::initObservables() {
  observables = new RooArgList("observables");
  // we need the same observables as both input PDFs:
  // copy over from first PDF
  for (const auto& obs : *pdf1->getObservables()) observables->add(*obs);
  // copy over from second PDF
  for (const auto& obs : *pdf2->getObservables()) observables->add(*obs);
}

void PDF_CrossCorAbs::copyMeasurementCovariance() {
  // copy errors
  for (int i = 0; i < nObs; i++) {
    if (i < pdf1->getNobs()) {
      StatErr[i] = pdf1->StatErr[i];
      SystErr[i] = pdf1->SystErr[i];
    } else if (i >= pdf1->getNobs()) {
      int shift = pdf1->getNobs();
      StatErr[i] = pdf2->StatErr[i - shift];
      SystErr[i] = pdf2->SystErr[i - shift];
    }
  }
  // copy correlations
  for (int i = 0; i < nObs; i++)
    for (int j = 0; j < nObs; j++) {
      if (i < pdf1->getNobs() && j < pdf1->getNobs()) {
        corStatMatrix[i][j] = pdf1->corStatMatrix[i][j];
        corSystMatrix[i][j] = pdf1->corSystMatrix[i][j];
      } else if (i >= pdf1->getNobs() && j >= pdf1->getNobs()) {
        int shift = pdf1->getNobs();
        corStatMatrix[i][j] = pdf2->corStatMatrix[i - shift][j - shift];
        corSystMatrix[i][j] = pdf2->corSystMatrix[i - shift][j - shift];
      }
    }
}

void PDF_CrossCorAbs::setCorrelations(Utils::config c) { assert(0); };

void PDF_CrossCorAbs::buildPdf() {
  TMatrixDSym covInverse = covMatrix;
  covInverse.Invert();
  pdf = new RooCrossCorPdf("pdf_" + name, "pdf_" + name, *observables, *theory, covInverse, nObs1);
}
