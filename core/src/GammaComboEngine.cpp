#include <GammaComboEngine.h>

#include <BatchScriptWriter.h>
#include <ColorBuilder.h>
#include <Combiner.h>
#include <FileNameBuilder.h>
#include <Graphviz.h>
#include <LatexMaker.h>
#include <MethodBergerBoosScan.h>
#include <MethodCoverageScan.h>
#include <MethodDatasetsPluginScan.h>
#include <MethodDatasetsProbScan.h>
#include <MethodPluginScan.h>
#include <MethodProbScan.h>
#include <OneMinusClPlot.h>
#include <OneMinusClPlot2d.h>
#include <OptParser.h>
#include <PDF_Abs.h>
#include <PDF_Datasets.h>
#include <PValueCorrection.h>
#include <ParameterCache.h>
#include <ParameterEvolutionPlotter.h>
#include <RooSlimFitResult.h>
#include <Utils.h>

// Needed to define GAMMACOMBO_VERSION. Header created during CMake generation
#include <VersionConfig.h>

#include <RooAbsPdf.h>
#include <RooMsgService.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>

#include <TApplication.h>
#include <TColor.h>
#include <TF1.h>
#include <TFile.h>
#include <TLatex.h>
#include <TLine.h>
#include <TMath.h>
#include <TObjString.h>
#include <TROOT.h>
#include <TString.h>
#include <TTree.h>

#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <string.h>
#include <string>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

GammaComboEngine::GammaComboEngine(TString name, int argc, char* argv[]) : runOnDataSet(false) {
  // time the program
  t.Start();

  // print the copyright banner
  printBanner();

  // parse the command line options
  arg = new OptParser();
  arg->bookAllOptions();
  arg->parseArguments(argc, argv);

  // configure names
  execname = argv[0];
  if (arg->filenameaddition != "") name += "_" + arg->filenameaddition;
  m_fnamebuilder = new FileNameBuilder(arg, name);

  // make batch scripts if appropriate and exit
  m_batchscriptwriter = new BatchScriptWriter(argc, argv);

  // run ROOT in interactive mode, if requested (-i)
  if (arg->interactive) {
    // -a is a reserved option for TApplication
    std::string actionStr = "--action";  // (needed to avoid compiler warnings)
    for (int i = 2; i < argc; ++i)
      if (!strcmp(argv[i], "-a")) { argv[i] = actionStr.data(); }
    theApp = new TApplication("App", &argc, argv);
  } else
    gROOT->SetBatch(false);

  // initialize members
  plot = 0;

  // reconfigure RooFormulaVar output
  RooMsgService::instance().getStream(1).removeTopic(RooFit::InputArguments);
  RooMsgService::instance().getStream(0).addTopic(RooFit::InputArguments);
}

GammaComboEngine::GammaComboEngine(TString name, int argc, char* argv[], bool _runOnDataSet) {
  GammaComboEngine(name, argc, argv);
  runOnDataSet = _runOnDataSet;
}

GammaComboEngine::~GammaComboEngine() {
  delete m_fnamebuilder;
  delete m_batchscriptwriter;
}

///
/// Check if a PDF with a certain ID exits.
///
bool GammaComboEngine::pdfExists(int id) {
  if (id < 0) return false;
  if (id >= this->pdf.size()) return false;
  if (this->pdf[id] == 0) return false;
  return true;
}

///
/// Check if a Combiner with a certain ID exits.
///
bool GammaComboEngine::combinerExists(int id) const {
  if (id < 0) return false;
  if (id >= this->cmb.size()) return false;
  if (this->cmb[id] == 0) return false;
  return true;
}

///
/// Set the PDF (for datasets method) for the GammaComboEngine
///
void GammaComboEngine::setPdf(PDF_Abs* pdf) {
  if (!runOnDataSet) {
    std::cout
        << "It looks like you're trying to set a pdf but you haven't set runOnDataSet=true. I assume this is what you "
           "want so I'm doing it for you"
        << std::endl;
    runOnDataSet = true;
  }
  if (!dynamic_cast<PDF_Datasets*>(pdf)) {
    std::cout << "GammaComboEngine::setPdf() : ERROR : The pdf you are trying to set " << pdf->getName()
              << " cannot be cast to a PDF_Datasets object" << std::endl;
    std::exit(1);
  }
  if (pdf == 0) {
    std::cout << "GammaComboEngine::setPdf() : ERROR : Trying to add zero pointer as the PDF. Exit." << std::endl;
    std::exit(1);
  }
  if (pdfExists(0)) {
    std::cout << "GammaComboEngine::setPdf() : ERROR : You have already set the pdf in GammaComboEngine. Exit."
              << std::endl;
    std::exit(1);
  }
  addPdf(0, pdf);
}

///
/// Add a PDF to the GammaComboEngine object.
///
void GammaComboEngine::addPdf(int id, PDF_Abs* pdf, TString title) {
  if (arg->debug) {
    std::cout << "GammaComboEngine::addPdf() : INFO  : Adding pdf " << id << " = " << title << std::endl;
  }
  if (pdf == 0) {
    std::cout << "GammaComboEngine::addPdf() : ERROR : Trying to add zero pointer as the PDF. Exit." << std::endl;
    std::exit(1);
  }
  // check if requested id exists already
  if (pdfExists(id)) {
    std::cout << "GammaComboEngine::addPdf() : ERROR : Requested PDF id " << id
              << " exists already in GammaComboEngine. Exit." << std::endl;
    std::exit(1);
  }
  // check if storage is large enough, enlarge if necessary
  if (id >= this->pdf.size()) {
    for (int i = this->pdf.size(); i <= id; i++) this->pdf.push_back(0);
  }
  this->pdf[id] = pdf;
  if (title != "") this->pdf[id]->setTitle(title);
  this->pdf[id]->setGcId(id);
}

///
/// Add a pdf with a subset of the observables to the GammaComboEngine
///
void GammaComboEngine::addSubsetPdf(int id, PDF_Abs* pdf, std::vector<int>& indices, TString title) {
  if (indices.size() > pdf->getObservables()->getSize()) {
    std::cout << "GammaComboEngine::addSubsetPdf() : ERROR - the subset size " << indices.size()
              << " is bigger than the observables size " << pdf->getObservables()->getSize() << std::endl;
    std::exit(1);
  }
  for (int i = 0; i < indices.size(); i++) {
    int index = indices[i];
    if (index > pdf->getObservables()->getSize() - 1 || index < 0) {
      std::cout << "GammaComboEngine::addSubsetPdf() : ERROR - one of the subset index values " << index
                << " is larger than the total number of of observables " << pdf->getObservables()->getSize()
                << " or it's less than zero" << std::endl;
      std::exit(1);
    }
  }
  RooArgList* obsToRemove = new RooArgList();
  RooArgList* theoryToRemove = new RooArgList();

  // loop over all observables and remove the ones that aren't in indices
  for (int i = 0; i < pdf->getObservables()->getSize(); i++) {
    if (std::find(indices.begin(), indices.end(), i) == indices.end()) {
      obsToRemove->add(*(pdf->getObservables()->at(i)));
      theoryToRemove->add(*(pdf->getTheory()->at(i)));
    }
  }
  pdf->getObservables()->remove(*obsToRemove);
  pdf->getTheory()->remove(*theoryToRemove);
  delete obsToRemove;
  delete theoryToRemove;

  // now sort out parameters
  RooArgList* paramsToRemove = new RooArgList();
  for (int i = 0; i < pdf->getParameters()->getSize(); i++) {
    bool paramFoundInTheory = false;
    for (int j = 0; j < pdf->getTheory()->getSize(); j++) {
      if (pdf->getTheory()->at(j)->dependsOn(*(pdf->getParameters()->at(i)))) { paramFoundInTheory = true; }
    }
    if (!paramFoundInTheory) paramsToRemove->add(*(pdf->getParameters()->at(i)));
  }
  pdf->getParameters()->remove(*paramsToRemove);
  delete paramsToRemove;

  // now sort out uncertainties
  std::vector<double> oldStatErrs = pdf->StatErr;
  std::vector<double> oldSystErrs = pdf->SystErr;
  pdf->StatErr.clear();
  pdf->SystErr.clear();
  for (int i = 0; i < indices.size(); i++) {
    int index = indices[i];
    pdf->StatErr.push_back(oldStatErrs[index]);
    pdf->SystErr.push_back(oldSystErrs[index]);
  }

  pdf->setNObs(indices.size());
  TMatrixDSym newCorStatMatrix(indices.size());
  TMatrixDSym newCorSystMatrix(indices.size());

  pdf->getSubCorrelationStat(newCorStatMatrix, indices);
  pdf->getSubCorrelationSyst(newCorSystMatrix, indices);

  pdf->corStatMatrix.ResizeTo(indices.size(), indices.size());
  pdf->corSystMatrix.ResizeTo(indices.size(), indices.size());
  pdf->corMatrix.ResizeTo(indices.size(), indices.size());
  pdf->covMatrix.ResizeTo(indices.size(), indices.size());

  pdf->corStatMatrix = newCorStatMatrix;
  pdf->corSystMatrix = newCorSystMatrix;

  pdf->buildCov();
  pdf->buildPdf();

  addPdf(id, pdf, title);
}

void GammaComboEngine::addSubsetPdf(int id, PDF_Abs* pdf, int i1, TString title) {
  std::vector<int> indices;
  indices.push_back(i1);
  addSubsetPdf(id, pdf, indices, title);
}

void GammaComboEngine::addSubsetPdf(int id, PDF_Abs* pdf, int i1, int i2, TString title) {
  std::vector<int> indices;
  indices.push_back(i1);
  indices.push_back(i2);
  addSubsetPdf(id, pdf, indices, title);
}

void GammaComboEngine::addSubsetPdf(int id, PDF_Abs* pdf, int i1, int i2, int i3, TString title) {
  std::vector<int> indices;
  indices.push_back(i1);
  indices.push_back(i2);
  indices.push_back(i3);
  addSubsetPdf(id, pdf, indices, title);
}

void GammaComboEngine::addSubsetPdf(int id, PDF_Abs* pdf, int i1, int i2, int i3, int i4, TString title) {
  std::vector<int> indices;
  indices.push_back(i1);
  indices.push_back(i2);
  indices.push_back(i3);
  indices.push_back(i4);
  addSubsetPdf(id, pdf, indices, title);
}

void GammaComboEngine::addSubsetPdf(int id, PDF_Abs* pdf, int i1, int i2, int i3, int i4, int i5, TString title) {
  std::vector<int> indices;
  indices.push_back(i1);
  indices.push_back(i2);
  indices.push_back(i3);
  indices.push_back(i4);
  indices.push_back(i5);
  addSubsetPdf(id, pdf, indices, title);
}

void GammaComboEngine::addSubsetPdf(int id, PDF_Abs* pdf, int i1, int i2, int i3, int i4, int i5, int i6,
                                    TString title) {
  std::vector<int> indices;
  indices.push_back(i1);
  indices.push_back(i2);
  indices.push_back(i3);
  indices.push_back(i4);
  indices.push_back(i5);
  indices.push_back(i6);
  addSubsetPdf(id, pdf, indices, title);
}

void GammaComboEngine::addSubsetPdf(int id, PDF_Abs* pdf, int i1, int i2, int i3, int i4, int i5, int i6, int i7,
                                    TString title) {
  std::vector<int> indices;
  indices.push_back(i1);
  indices.push_back(i2);
  indices.push_back(i3);
  indices.push_back(i4);
  indices.push_back(i5);
  indices.push_back(i6);
  indices.push_back(i7);
  addSubsetPdf(id, pdf, indices, title);
}

void GammaComboEngine::addSubsetPdf(int id, PDF_Abs* pdf, int i1, int i2, int i3, int i4, int i5, int i6, int i7,
                                    int i8, TString title) {
  std::vector<int> indices;
  indices.push_back(i1);
  indices.push_back(i2);
  indices.push_back(i3);
  indices.push_back(i4);
  indices.push_back(i5);
  indices.push_back(i6);
  indices.push_back(i7);
  indices.push_back(i8);
  addSubsetPdf(id, pdf, indices, title);
}

///
/// Add a Combiner to the GammaComboEngine object.
///
void GammaComboEngine::addCombiner(int id, Combiner* cmb) {
  if (runOnDataSet) {
    std::cout
        << "GammaComboEngine::addCombiner() : ERROR : You're trying to make a combiner but the runOnDataSet flag is "
           "true. You cannot make a combination with this option"
        << std::endl;
    std::exit(1);
  }

  if (cmb == 0) {
    std::cout << "GammaComboEngine::addCombiner() : ERROR : Trying to add zero pointer as the Combiner. Exit."
              << std::endl;
    std::exit(1);
  }
  // check if requested id exists already
  if (combinerExists(id)) {
    std::cout
        << "GammaComboEngine::addCombiner() : ERROR : Requested Combiner id exists already in GammaComboEngine. Exit."
        << std::endl;
    std::exit(1);
  }
  // check if storage is large enough, enlarge if necessary
  if (id >= this->cmb.size()) {
    for (int i = this->cmb.size(); i <= id; i++) this->cmb.push_back(0);
  }
  this->cmb[id] = cmb;
}

///
/// Add a new Combiner to the GammaComboEngine object, that is the
/// clone of an existing one. After that, the clone can be modified
/// using, e.g., getCombiner(newId)->addPdf(...)
///
void GammaComboEngine::cloneCombiner(int newId, int oldId, TString name, TString title) {
  if (runOnDataSet) {
    std::cout
        << "GammaComboEngine::cloneCombiner() : ERROR : You're trying to clone a combiner but the runOnDataSet flag "
           "is true. You can't have combiners when using the dataset option."
        << std::endl;
    std::exit(1);
  }

  if (combinerExists(newId)) {
    std::cout << "GammaComboEngine::cloneCombiner() : ERROR : Requested new Combiner id " << newId
              << " exists already in GammaComboEngine. Exit." << std::endl;
    std::exit(1);
  }
  if (!combinerExists(oldId)) {
    std::cout << "GammaComboEngine::cloneCombiner() : ERROR : Requested old Combiner id " << oldId
              << " doesn't exists in GammaComboEngine. Exit." << std::endl;
    std::exit(1);
  }
  addCombiner(newId, getCombiner(oldId)->Clone(name, title));
}

///
/// Get a combiner.
/// \param id - combiner ID, set when defining the combiner using addCombiner(), cloneCombiner(), or newCombiner()
///
Combiner* GammaComboEngine::getCombiner(int id) const {
  if (!combinerExists(id)) {
    std::cout << "GammaComboEngine::getCombiner() : ERROR : Requested Combiner id " << id
              << " doesn't exist in GammaComboEngine. Exit." << std::endl;
    std::exit(1);
  }
  return cmb[id];
}

///
/// Get a PDF.
/// \param id - PDF ID, set when adding the PDF using addPdf()
///
PDF_Abs* GammaComboEngine::getPdf(int id) {
  if (!pdfExists(id)) {
    std::cout << "GammaComboEngine::getPdf() : ERROR : Requested PDF id doesn't exist in GammaComboEngine. Exit."
              << std::endl;
    std::exit(1);
  }
  return pdf[id];
}

PDF_Abs* GammaComboEngine::operator[](int idx) { return getPdf(idx); }

// const PDF_Abs* GammaComboEngine::operator[](int idx) const
// {
//  const PDF_Abs* p = (const)getPdf(idx)
//  return p;
// }

///
/// Add a new Combiner, consisting of the specified PDFs.
/// The pdf arguments refer to the GammaComboEngine ID of the PDFs
/// that should be combined. Add them before using addPdf().
/// Also an empty combiner is possible, that has no PDFs, so it can be filled
/// later. In that case, leave all pdf arguments at -1.
///
void GammaComboEngine::newCombiner(int id, TString name, TString title, int pdf1, int pdf2, int pdf3, int pdf4,
                                   int pdf5, int pdf6, int pdf7, int pdf8, int pdf9, int pdf10, int pdf11, int pdf12,
                                   int pdf13, int pdf14, int pdf15) {
  if (combinerExists(id)) {
    std::cout
        << "GammaComboEngine::newCombiner() : ERROR : Requested new Combiner id exists already in GammaComboEngine. "
           "Exit."
        << std::endl;
    std::exit(1);
  }
  Combiner* c = new Combiner(arg, name, title);
  if (pdf1 > -1) c->addPdf(getPdf(pdf1));
  if (pdf2 > -1) c->addPdf(getPdf(pdf2));
  if (pdf3 > -1) c->addPdf(getPdf(pdf3));
  if (pdf4 > -1) c->addPdf(getPdf(pdf4));
  if (pdf5 > -1) c->addPdf(getPdf(pdf5));
  if (pdf6 > -1) c->addPdf(getPdf(pdf6));
  if (pdf7 > -1) c->addPdf(getPdf(pdf7));
  if (pdf8 > -1) c->addPdf(getPdf(pdf8));
  if (pdf9 > -1) c->addPdf(getPdf(pdf9));
  if (pdf10 > -1) c->addPdf(getPdf(pdf10));
  if (pdf11 > -1) c->addPdf(getPdf(pdf11));
  if (pdf12 > -1) c->addPdf(getPdf(pdf12));
  if (pdf13 > -1) c->addPdf(getPdf(pdf13));
  if (pdf14 > -1) c->addPdf(getPdf(pdf14));
  if (pdf15 > -1) c->addPdf(getPdf(pdf15));
  addCombiner(id, c);
}

///
/// scale down errors
/// \todo implement a more generic version of this function. As it is,
/// it is only useful for the LHCb gamma combination.
///
void GammaComboEngine::scaleDownErrors() {
  for (int i = 0; i < pdf.size(); i++) {
    // these are the PDFs for the full combination:
    // if (! (i==61 || i==58 || i==60 || i==56 || i==54 || i==40 || i==43) ) continue;
    float scale = 1.;
    if (i == 25) scale = sqrt(5. * 50. / 3.);
    if (i == 60) scale = sqrt(50. / 3.);
    // if ( i==23 ) scale = sqrt(50.);
    if (scale == 1.) continue;
    std::cout << "Configuration: Scaling down LHCb errors by a factor " << scale << ": " << pdf[i]->getTitle()
              << std::endl;
    for (int iObs = 0; iObs < pdf[i]->getNobs(); iObs++) {
      pdf[i]->StatErr[iObs] /= scale;
      pdf[i]->SystErr[iObs] /= scale;
    }
    pdf[i]->setErrorSourceString("Scaled!");
    pdf[i]->buildCov();
    pdf[i]->buildPdf();
  }

  int i = 0;

  // switch off Afav
  i = 25;
  pdf[i]->ScaleError("afav_dk_kpi_obs", 200.);
  pdf[i]->buildCov();
  pdf[i]->buildPdf();

  // scale dD_kpi
  i = 7;
  pdf[i]->ScaleError("dD_kpi_obs", 1. / 20.);
  pdf[i]->buildCov();
  pdf[i]->buildPdf();

  std::cout << std::endl;
}

///
/// scale stat errors
///
void GammaComboEngine::scaleStatErrors() {
  std::cout << "\nConfiguration: Scaling ALL STAT ERRORS by " << arg->scalestaterr << ".\n" << std::endl;
  for (int i = 0; i < pdf.size(); i++) {
    if (pdf[i] == 0) continue;
    for (int iObs = 0; iObs < pdf[i]->getNobs(); iObs++) { pdf[i]->StatErr[iObs] *= arg->scalestaterr; }
    pdf[i]->buildCov();
    pdf[i]->buildPdf();
  }
}

///
/// scale stat+syst errors
///
void GammaComboEngine::scaleStatAndSystErrors() {
  std::cout << "\nConfiguration: Scaling ALL STAT AND SYST ERRORS by " << arg->scaleerr << ".\n" << std::endl;
  for (int i = 0; i < pdf.size(); i++) {
    if (pdf[i] == 0) continue;
    for (int iObs = 0; iObs < pdf[i]->getNobs(); iObs++) {
      pdf[i]->StatErr[iObs] *= arg->scaleerr;
      pdf[i]->SystErr[iObs] *= arg->scaleerr;
    }
    pdf[i]->buildCov();
    pdf[i]->buildPdf();
  }
}

///
/// disable systematics
///
void GammaComboEngine::disableSystematics() {
  std::cout << "\nConfiguration: Setting ALL SYSTEMATICS TO ZERO.\n" << std::endl;
  for (int i = 0; i < pdf.size(); i++) {
    if (pdf[i] == 0) continue;
    for (int iObs = 0; iObs < pdf[i]->getNobs(); iObs++) pdf[i]->SystErr[iObs] = 0.;
    pdf[i]->buildCov();
    pdf[i]->buildPdf();
  }
}

///
/// Make an Asimov toy: set all observables set to truth values.
/// The Asimov point needs to be loaded in the combiner before.
/// \param c - combiner which should be set to an asimov toy
///
void GammaComboEngine::setAsimovObservables(Combiner* c) {
  if (!c->isCombined()) {
    std::cout << "GammaComboEngine::setAsimovObservables() : ERROR : Can't set an Asimov toy before "
                 "the combiner is combined. Call combine() first."
              << std::endl;
    std::exit(1);
  }

  // set observables to asimov values in workspace
  RooWorkspace* w = c->getWorkspace();
  for (const auto& pAbsObs : *c->getObservables()) {
    const auto pObs = static_cast<RooRealVar*>(pAbsObs);
    // get theory name from the observable name
    TString pThName = pObs->GetName();
    pThName.ReplaceAll("obs", "th");
    // get the theory relation
    RooAbsReal* th = w->function(pThName);
    if (th == 0) {
      std::cout << "GammaComboEngine::setAsimovObservables() : ERROR : theory relation not found in workspace: "
                << pThName << std::endl;
      std::exit(1);
    }
    // set the observable to what the theory relation predicts
    pObs->setVal(th->getVal());
  }

  // write back the asimov values to the PDF object so that when
  // the PDF is printed, the asimov values show up
  for (int i = 0; i < c->getPdfs().size(); i++) {
    PDF_Abs* pdf = c->getPdfs()[i];
    pdf->setObservableSourceString("Asimov");
    for (const auto& pObs : *pdf->getObservables()) {
      RooAbsReal* obs = w->var(pObs->GetName());
      if (obs == 0) {
        std::cout << "GammaComboEngine::setAsimovObservables() : ERROR : observable not found in workspace: "
                  << pObs->GetName() << std::endl;
        std::exit(1);
      }
      pdf->setObservable(pObs->GetName(), obs->getVal());
    }
  }
}

///
/// Load start parameters.
///
/// \param s - the scanner
/// \param cId - combiner id
/// \param pCache - parameter cache
///
void GammaComboEngine::loadStartParameters(MethodProbScan* s, ParameterCache* pCache, int cId) {
  std::cout << "Start parameter configuration:\n" << std::endl;
  TString startparfile;
  TString startparfile2 = m_fnamebuilder->getFileNameStartPar(s);
  TString startparfile3;
  if (arg->isAsimovCombiner(cId)) {
    startparfile3 =
        m_fnamebuilder->getFileNameStartPar(s);  // this gets the start parameter file of the Asimov combiner
    startparfile3.ReplaceAll(m_fnamebuilder->getAsimovCombinerNameAddition(arg->asimov[cId]), "");
  }
  bool filefound = false;
  // try the file provided through --parfile
  if (arg->loadParamsFile.size() > cId && !arg->loadParamsFile[cId].EqualTo("default")) {
    startparfile = arg->loadParamsFile[cId];
    filefound = Utils::FileExists(startparfile);
    if (!filefound) {
      std::cout << "\n ERROR: --parfile not found: " << startparfile << std::endl;
      std::cout << "  Will now look for default files.\n" << std::endl;
    }
  }
  // requested file not found, try default
  if (!filefound) {
    startparfile = startparfile2;
    filefound = Utils::FileExists(startparfile);
  }
  // for Asimov combiners, also try non-Asimov start parameters file
  if (!filefound && arg->isAsimovCombiner(cId)) {
    startparfile = startparfile3;
    filefound = Utils::FileExists(startparfile);
  }
  // still not found
  if (!filefound) {
    std::cout << "  No start parameter file was found, will use the default values" << std::endl;
    std::cout << "  configured in ParameterAbs class. Start parameter files are searched" << std::endl;
    std::cout << "  in the following order:" << std::endl;
    std::cout << "  1. --parfile option" << std::endl;
    std::cout << "  2. default start paramter file: " << startparfile2 << std::endl;
    if (arg->isAsimovCombiner(cId)) {
      std::cout << "  3. for Asimov combiners: corresponding non-Asimov start parameter file: " << startparfile3
                << std::endl;
    }
  } else {
    std::cout << "  Loading start parameters from file: " << startparfile << std::endl;
    bool loaded = pCache->loadPoints(startparfile);
    if (!loaded) {
      std::cout << "  Error loading file. Exit." << std::endl;
      std::exit(1);
    }
  }
  std::cout << std::endl;
}

///
/// Configure the names and titles of Asimov combiners.
/// It adds things like "Asimov3" to the combiner name,
/// where 3 means the 3rd Asimov point from the configured
/// Asimov parameter file (typically *_genpoints.dat).
///
void GammaComboEngine::configureAsimovCombinerNames(Combiner* c, int i) {
  if (arg->asimov[i] == 0) {
    std::cout << "\n--asimov 0 : ignoring generator point ID 0" << std::endl;
    return;
  }
  std::cout << "\n--asimov : setting up an ASIMOV TOY based on combination \"" << c->getName() << "\"" << std::endl;
  if (arg->title[i] == TString("default")) c->setTitle(c->getTitle() + " (Asimov)");
  c->setName(c->getName() + m_fnamebuilder->getAsimovCombinerNameAddition(arg->asimov[i]));
  std::cout << "           Asimov combiner name: \"" << c->getName() << "\"" << std::endl;
}

///
/// Make an Asimov toy: set all observables set to truth values.
/// The truth values are loaded from a parameter file.
///
void GammaComboEngine::loadAsimovPoint(Combiner* c, int cId) {
  if (arg->asimov[cId] == 0) return;
  std::cout << "\nAsimov point configuration:\n" << std::endl;
  ParameterCache* pCache = new ParameterCache(arg);
  TString asimovfile;
  TString asimovfile2 = m_fnamebuilder->getFileNameAsimovPar(c);
  TString asimovfile3 =
      m_fnamebuilder->getFileNameStartPar(c);  // this gets the start parameter file of the Asimov combiner
  asimovfile3.ReplaceAll(m_fnamebuilder->getAsimovCombinerNameAddition(arg->asimov[cId]), "");
  TString asimovfile4 =
      m_fnamebuilder->getFileNamePar(c);  // this gets the result parameter file of the Asimov combiner
  asimovfile4.ReplaceAll(m_fnamebuilder->getAsimovCombinerNameAddition(arg->asimov[cId]), "");
  bool filefound = false;
  // try the file provided through --asimovfile
  if (arg->asimovfile.size() > cId && !arg->asimovfile[cId].EqualTo("default")) {
    asimovfile = arg->loadParamsFile[cId];
    filefound = Utils::FileExists(asimovfile);
  }
  // requested file not found, try default
  if (!filefound) {
    asimovfile = asimovfile2;
    filefound = Utils::FileExists(asimovfile);
  }
  // requested file not found, try the start parameter file of the corresponding non-Asimov combiner
  if (!filefound) {
    asimovfile = asimovfile3;
    filefound = Utils::FileExists(asimovfile);
  }
  // requested file not found, try the result parameter file of the corresponding non-Asimov combiner
  if (!filefound) {
    asimovfile = asimovfile4;
    filefound = Utils::FileExists(asimovfile);
  }
  // if no parameter file exits, we use point from the ParameterAbs class
  if (!filefound) {
    std::cout << "  No Asimov point parameter file found. Using start values configured in ParameterAbs class."
              << std::endl;
    std::cout << "  Point files are looked for in the following order:" << std::endl;
    std::cout << "  1. --asimovfile" << std::endl;
    std::cout << "  2. " << asimovfile2 << std::endl;
    std::cout << "  3. " << asimovfile3 << std::endl;
    std::cout << "  4. " << asimovfile4 << std::endl;
  } else {
    std::cout << "  Loading Asimov points from file: " << asimovfile << std::endl;
    bool loaded = pCache->loadPoints(asimovfile);
    if (!loaded) {
      std::cout << "  Error loading file. Exit." << std::endl;
      std::exit(1);
    }
    std::cout << "  Setting point number: " << arg->asimov[cId] << std::endl;
    pCache->setPoint(c, arg->asimov[cId] - 1);
  }
  setAsimovObservables(c);
}

///
/// print usage and exit
///
void GammaComboEngine::usage() {
  if (runOnDataSet) {
    std::cout << "USAGE\n\n"
                 "  # Compute limit on parameter a:\n"
                 "  "
              << execname
              << " --var a --ps 1\n\n"
                 "  # Compute limits on two parameters a and b:\n"
                 "  "
              << execname << " --var a --var b --ps 1\n\n"
              << std::endl;
    std::exit(0);
  }
  std::cout << "USAGE\n\n"
               "  # Compute combination 1, make a 1D Prob scan for variable a_gaus:\n"
               "  "
            << execname
            << " -c 1 -i --var a_gaus --ps 1\n\n"
               "  # Add combinations 1, 2, 3, add to the same plot:\n"
               "  "
            << execname
            << " -c 1 -c 2 -c 3 -i --var a_gaus --ps 1\n\n"
               "  # Add or remove measurements to/from existing cominers:\n"
               "  "
            << execname
            << " -i --var a_gaus -c 3:+1,+2,-3\n\n"
               "  # Make a 2D Prob scan\n"
               "  "
            << execname
            << " -c 4 --var a_gaus --var b_gaus -i\n\n"
               "  # Re-plot a previous 2D Prob scan with 2D CL and best fit points\n"
               "  "
            << execname << " -c 4 --var a_gaus --var b_gaus -i -a plot --2dcl --ps 1\n"
            << std::endl;
  // std::cout << "Available plugin toy control plots:\n"
  //"===================================\n"
  //"\n"
  //"All control plots are produced when option --controlplots is given.\n"
  //"If, in addition, -p <n> is given, only one specific plot is produced.\n"
  //" (1) technical overview\n"
  //" (2) summary\n"
  //" (3) nuisances\n"
  //" (4) observables\n"
  //" (5) chi2 distributions\n"
  //" (6) chi2 parabolas\n" << std::endl;
  print();
  std::exit(0);
}

///
/// Print the available PDFs.
///
void GammaComboEngine::printPdfs() {
  std::cout << "AVAILABLE MEASUREMENTS" << std::endl;
  std::cout << std::endl;
  for (int i = 0; i < pdf.size(); i++) {
    if (pdf[i] == 0) continue;
    if (i < 10)
      printf("   (%i) %s\n", i, pdf[i]->getTitle().Data());
    else if (i < 100)
      printf("  (%2i) %s\n", i, pdf[i]->getTitle().Data());
    else
      printf(" (%3i) %s\n", i, pdf[i]->getTitle().Data());
  }
  std::cout << std::endl;
}

///
/// Print the availabe Combinations.
///
void GammaComboEngine::printCombinations() {
  std::cout << "AVAILABLE COMBINATIONS" << std::endl;
  std::cout << std::endl;
  for (int i = 0; i < cmb.size(); i++) {
    if (cmb[i] == 0) continue;
    if (i < 10)
      printf("   (%i) %s\n", i, cmb[i]->getTitle().Data());
    else if (i < 100)
      printf("  (%2i) %s\n", i, cmb[i]->getTitle().Data());
    else
      printf(" (%3i) %s\n", i, cmb[i]->getTitle().Data());
  }
  std::cout << std::endl;
}

///
/// Print the content of this engine.
///
void GammaComboEngine::print() {
  printPdfs();
  printCombinations();
}

///
/// Check the combination argument (-c), exit if it is bad.
///
void GammaComboEngine::checkCombinationArg() {
  if (runOnDataSet && arg->combid.size() > 0) {
    std::cout << "When running on a dataset do not pass a combination argument (it makes no sense for this use case)"
              << std::endl;
    std::exit(1);
  }
  if (arg->combid.size() == 0 && !runOnDataSet) {
    std::cout << "Please chose a combination ID (-c).\n" << std::endl;
    printCombinations();
    std::exit(1);
  }
  for (int i = 0; i < arg->combid.size(); i++) {
    if (arg->combid[i] >= cmb.size()) {
      std::cout << "Please chose a combination ID (-c) less than " << cmb.size()
                << ".\nUse the -u option to print a list of available combinations." << std::endl;
      std::exit(1);
    }
    if (cmb[arg->combid[i]] == 0) {
      std::cout << "You selected an empty combination.\n"
                << "Use the -u option to print a list of available combinations." << std::endl;
      std::exit(1);
    }
  }
}

///
/// Check Asimov arg: when only one --asimov argument is given
/// with the ID 0, it won't do anything. Print a warning in that
/// case.
///
void GammaComboEngine::checkAsimovArg() {
  if (arg->asimov.size() == 1 && arg->asimov[0] == 0) {
    std::cout << "WARNING : --asimov 0 found, this won't do anything." << std::endl;
    std::cout << "          To run an Asimov toy, the generation point ID" << std::endl;
    std::cout << "          needs to be different from 0." << std::endl;
    std::cout << std::endl;
  }
}

///
/// Check color argument and exit if a non-existing color was requested through
/// the --color argument. Colors for one-dimensional plots are defined in
/// GammaComboEngine::defineColors(). Colors for two-dimensional pltos are
/// defined in OneMinusClPlot2d::OneMinusClPlot2d().
///
void GammaComboEngine::checkColorArg() {
  for (int i = 0; i < arg->color.size(); i++) {
    // colors for one-dimensional plots
    if (arg->var.size() == 1) {
      if (colorsLine.size() <= arg->color[i]) {
        std::cout
            << "Argument error --color: No such color for one-dimensional plots. Please choose a color between 0 and "
            << colorsLine.size() - 1 << std::endl;
        std::exit(1);
      }
    }
    // colors for two-dimensional plots
    else if (arg->var.size() == 2) {
      OneMinusClPlot2d p(arg);
      int nMaxColors = p.getNumberOfDefinedColors();
      if (nMaxColors <= arg->color[i]) {
        std::cout
            << "Argument error --color: No such color for two-dimensional plots. Please choose a color between 0 and "
            << nMaxColors - 1 << std::endl;
        std::exit(1);
      }
    }
  }
}

///
/// Implement the logic needed to create new combinations on the fly
/// using the -c 26:+12 syntax
///
void GammaComboEngine::makeAddDelCombinations() {
  if (runOnDataSet) return;

  // sanity check: the combid and combmodifications vectors should
  // be the same size
  if (arg->combmodifications.size() != arg->combid.size()) {
    std::cout << "GammaComboEngine::makeAddDelCombinations() : ERROR : internal inconsistency. \n"
                 "combid and combmodifications vectors not of same size."
              << std::endl;
    assert(0);
  }
  // loop over list of modifications
  for (int i = 0; i < arg->combmodifications.size(); i++) {
    // get combiner that is to be modified
    Combiner* cOld = cmb[arg->combid[i]];
    // see if anything is to be added or subtracted
    if (arg->combmodifications[i].size() == 0) continue;
    // there are modifications to be done!
    std::cout << "-c : Making a new combination based on combination " << arg->combid[i] << std::endl;
    // compute name and title of new combiner
    TString nameNew = cOld->getName();
    TString titleNew = cOld->getTitle();
    for (int j = 0; j < arg->combmodifications[i].size(); j++) {
      int pdfId = abs(arg->combmodifications[i][j]);
      if (!pdfExists(pdfId)) {
        std::cout << "\nERROR: measurement of given ID does not exist: " << pdfId << std::endl;
        std::cout << "       Here is a list of available measurements:" << std::endl;
        printPdfs();
        std::exit(1);
      }
      if (arg->combmodifications[i][j] > 0) {
        nameNew += Form("+%i", pdfId);
        titleNew += Form(", + Meas.%i", pdfId);
      } else {
        nameNew += Form("-%i", pdfId);
        titleNew += Form(", w/o Meas.%i", pdfId);
      }
    }
    // make the new combiner
    Combiner* cNew = cOld->Clone(nameNew, titleNew);
    // add/delete pdfs
    for (int j = 0; j < arg->combmodifications[i].size(); j++) {
      int pdfId = abs(arg->combmodifications[i][j]);
      if (arg->combmodifications[i][j] > 0) {
        std::cout << "... adding measurement " << pdfId << std::endl;
        cNew->addPdf(pdf[pdfId]);
      } else {
        std::cout << "... deleting measurement " << pdfId << std::endl;
        cNew->delPdf(pdf[pdfId]);
      }
    }
    // add to list of combinations to compute this round
    cmb.push_back(cNew);
    arg->combid[i] = cmb.size() - 1;
  }
}

///
/// print parameter structure of the combinations into
/// .dot file
///
void GammaComboEngine::printCombinerStructure(Combiner* c) {
  Graphviz gviz(arg);
  gviz.printCombiner(c);
  gviz.printCombinerLayer(c);
}

///
/// Override default titles of the combinations, if requested
/// on the command line.
///
void GammaComboEngine::customizeCombinerTitles() {
  for (int i = 0; i < arg->combid.size(); i++) {
    int combinerId = arg->combid[i];
    Combiner* c = cmb[combinerId];
    if (i < arg->title.size()) {
      if (arg->title[i] != TString("default")) c->setTitle(arg->title[i]);
    }
  }
  if (runOnDataSet) { pdf[0]->setTitle(arg->title[0]); }
}

///
/// Set up the plot.
///
void GammaComboEngine::setUpPlot() {
  if (arg->var.size() == 1) {
    plot = new OneMinusClPlot(arg, m_fnamebuilder->getFileNamePlot(cmb), "p-value curves");
  } else {
    plot = new OneMinusClPlot2d(arg, m_fnamebuilder->getFileNamePlot(cmb), "p-value contours");
  }
  plot->disableLegend(arg->plotlegend);
}

///
/// Save the plot to disc.
///
void GammaComboEngine::savePlot() {
  if (arg->hfagLabel != "")
    Utils::HFAGLabel(arg->hfagLabel, arg->plotHFAGLabelPosX, arg->plotHFAGLabelPosY, arg->plotHFAGLabelScale);
  plot->save();
}

///
/// Define the colors of curves and numbers in the 1D plots.
/// It is important to call this only after new combinations were
/// made using the makeAddDelCombinations() mechanism.
///
void GammaComboEngine::defineColors() {
  // no --color option was given on the command line
  // if ( arg->color.size()==0 )
  //{
  // define line colors for 1-CL curves
  colorsLine.push_back(arg->combid.size() <= 1 ? kBlue - 8 : kBlue - 5);  // 0
  colorsLine.push_back(kGreen - 8);                                       // 1
  colorsLine.push_back(kOrange - 8);                                      // 2
  colorsLine.push_back(kMagenta - 6);                                     // 3

  // define text colors for drawn central values
  colorsText.push_back(arg->combid.size() == 1 ? kBlack : TColor::GetColor("#23236b"));
  colorsText.push_back(TColor::GetColor("#234723"));
  colorsText.push_back(kOrange + 3);
  colorsText.push_back(kMagenta - 8);
  //}
  // else
  //{
  colorsLine.push_back(TColor::GetColor("#1b9e77"));  // 4 sea green
  colorsLine.push_back(TColor::GetColor("#d95f02"));  // 5 dark orange
  colorsLine.push_back(TColor::GetColor("#7570b3"));  // 6 medium purple
  colorsLine.push_back(TColor::GetColor("#e7298a"));  // 7 medium violet red
  colorsLine.push_back(TColor::GetColor("#66a61e"));  // 8 forest green
  colorsLine.push_back(TColor::GetColor("#e6ab02"));  // 9 goldenrod
  colorsLine.push_back(TColor::GetColor("#a6761d"));  // 10 chocolate
  colorsLine.push_back(TColor::GetColor("#e31a1c"));  // 11 red
  colorsLine.push_back(TColor::GetColor("#984ea3"));  // 12 darkish purple
  colorsLine.push_back(kBlack);                       // 13 black

  // some color blind safe options
  colorsLine.push_back(TColor::GetColor("#74add1"));  // 14 light blue
  colorsLine.push_back(TColor::GetColor("#f46d43"));  // 15 coral
  colorsLine.push_back(TColor::GetColor("#fdae61"));  // 16 orangey
  colorsLine.push_back(TColor::GetColor("#d73027"));  // 17 red
  colorsLine.push_back(TColor::GetColor("#4575b4"));  // 18 dark blue
  colorsLine.push_back(TColor::GetColor("#fee090"));  // 19 yellow

  // from http://colorbrewer2.org with:
  //   number of data classes: 6
  //   nature of data:         qualitative
  //   second colour scheme

  ColorBuilder cb;

  for (int i = 4; i < colorsLine.size(); i++) {
    // colorsText.push_back(cb.darklightcolor(colorsLine[i], 0.5));
    colorsText.push_back(colorsLine[i]);
  }
  //}

  // default for any additional scanner
  for (int i = colorsLine.size(); i < arg->combid.size(); i++) {
    colorsLine.push_back(kBlue - 8 + i);
    colorsText.push_back(kBlue - 2 + i);
  }

  // sort out the fill style vector
  for (int i = 0; i < arg->combid.size(); i++) {
    if (i >= arg->fillstyle.size())
      fillStyles.push_back(1001);
    else
      fillStyles.push_back(arg->fillstyle[i]);
  }
  // sort out the fill color vector
  for (int i = 0; i < arg->combid.size(); i++) {
    // if --hexfillcolor passed use this
    if (i < arg->hexfillcolor.size()) fillColors.push_back(TColor::GetColor(TString(arg->hexfillcolor[i])));
    // if --fillcolor passed use this
    else if (i < arg->fillcolor.size())
      fillColors.push_back(arg->fillcolor[i]);
    // if --color passed use this
    else if (i < arg->color.size())
      fillColors.push_back(colorsLine[arg->color[i]]);
    // otherwise use default
    else
      fillColors.push_back(colorsLine[i]);
  }
  // sort out the fill transparency vector
  for (int i = 0; i < arg->combid.size(); i++) {
    if (i >= arg->filltransparency.size())
      fillTransparencies.push_back(0.);
    else
      fillTransparencies.push_back(arg->filltransparency[i]);
  }

  // sort out the line width vector
  for (int i = 0; i < arg->combid.size(); i++) {
    if (i >= arg->linewidth.size())
      lineWidths.push_back(2);
    else
      lineWidths.push_back(arg->linewidth[i]);
  }
  // sort out the line style vector
  for (int i = 0; i < arg->combid.size(); i++) {
    if (i >= arg->linestyle.size())
      lineStyles.push_back(1);
    else
      lineStyles.push_back(arg->linestyle[i]);
  }
  // sort out the line color vector
  for (int i = 0; i < arg->combid.size(); i++) {
    // if --hexlinecolor passed use this
    if (i < arg->hexlinecolor.size()) lineColors.push_back(TColor::GetColor(TString(arg->hexlinecolor[i])));
    // if --linecolor passed use this
    else if (i < arg->linecolor.size())
      lineColors.push_back(arg->linecolor[i]);
    // if --color passed use this
    else if (i < arg->color.size())
      lineColors.push_back(colorsLine[arg->color[i]]);
    // otherwise use default
    else
      lineColors.push_back(colorsLine[i]);
  }
  // catch for datasets
  if (arg->combid.size() == 0) {
    // sort out fill style
    if (arg->fillstyle.size() > 0)
      fillStyles.push_back(arg->fillstyle[0]);
    else
      fillStyles.push_back(1001);
    // sort out fill color
    if (arg->fillcolor.size() > 0)
      fillColors.push_back(arg->fillcolor[0]);
    else
      fillColors.push_back(colorsLine[0]);
    // sort out line width
    if (arg->linewidth.size() > 0)
      lineWidths.push_back(arg->linewidth[0]);
    else
      lineWidths.push_back(2);
    // sort out line style
    if (arg->linestyle.size() > 0)
      lineStyles.push_back(arg->linestyle[0]);
    else
      lineStyles.push_back(1);
    // sort out line color
    if (arg->linecolor.size() > 0)
      lineColors.push_back(arg->linecolor[0]);
    else
      lineColors.push_back(colorsLine[0]);
    // sort out fill transparency
    if (arg->filltransparency.size() > 0)
      fillTransparencies.push_back(arg->filltransparency[0]);
    else
      fillTransparencies.push_back(.0);
    // sort out total color if not line and fill color defined before
    if (arg->color.size() > 0 && arg->color[0] < colorsLine.size() && arg->linecolor.size() == 0)
      lineColors[0] = colorsLine[arg->color[0]];
    if (arg->color.size() > 0 && arg->color[0] < colorsLine.size() && arg->fillcolor.size() == 0)
      fillColors[0] = colorsLine[arg->color[0]];
  }
}

///
/// Define scan strategy for a 2D scan.
///
void GammaComboEngine::scanStrategy2d(MethodProbScan* scanner, ParameterCache* pCache) {
  int nStartingPoints = pCache->getNPoints();
  // if no starting values loaded do the default thing
  if (nStartingPoints == 0) {
    std::cout << "\nPerforming default 2D scan:\n"
                 " 1. scan in first variable:  " +
                     scanner->getScanVar1Name() +
                     "\n"
                     " 2. scan in second variable: " +
                     scanner->getScanVar2Name() +
                     "\n"
                     " 3. scan starting from each solution found in 1. and 2."
              << std::endl;

    // setup a scanner for each variable individually
    Combiner* c = scanner->getCombiner();
    MethodProbScan* s1;
    MethodProbScan* s2;

    std::cout << "\n1D scan for X variable, " + scanner->getScanVar1Name() + ":\n" << std::endl;
    if (runOnDataSet) {
      const MethodDatasetsProbScan* temp = dynamic_cast<MethodDatasetsProbScan*>(scanner);
      s1 = new MethodDatasetsProbScan(temp->pdf, arg);
    } else {
      s1 = new MethodProbScan(c);
    }
    s1->setScanVar1(scanner->getScanVar1Name());
    s1->initScan();
    scanStrategy1d(s1, pCache);
    if (arg->verbose) s1->printLocalMinima();

    std::cout << "\n1D scan for Y variable, " + scanner->getScanVar2Name() + ":\n" << std::endl;
    if (runOnDataSet) {
      const MethodDatasetsProbScan* temp = dynamic_cast<MethodDatasetsProbScan*>(scanner);
      s2 = new MethodDatasetsProbScan(temp->pdf, arg);
    } else {
      s2 = new MethodProbScan(c);
    }
    s2->setScanVar1(scanner->getScanVar2Name());
    s2->setXscanRange(arg->scanrangeyMin, arg->scanrangeyMax);
    s2->initScan();
    scanStrategy1d(s2, pCache);
    if (arg->verbose) s2->printLocalMinima();

    // now do the 2D scan from the two starting points
    std::cout << "\n2D scan for " + scanner->getScanVar1Name() + " and " + scanner->getScanVar2Name() + ":\n"
              << std::endl;
    std::vector<RooSlimFitResult*> solutions;
    for (int i = 0; i < s1->getSolutions().size(); i++) solutions.push_back(s1->getSolution(i));
    for (int i = 0; i < s2->getSolutions().size(); i++) solutions.push_back(s2->getSolution(i));
    // \todo remove similar solutions from list
    delete s1;
    delete s2;
    for (int j = 0; j < solutions.size(); j++) {
      std::cout << "2D scan " << j + 1 << " of " << solutions.size() << " ..." << std::endl;
      scanner->loadParameters(solutions[j]);
      scanner->scan2d();
    }
  }
  // otherwise load each starting value found
  else {
    std::cout << "\nPerforming 2D scan from provided starting points." << std::endl;
    std::cout << "Number of scans to run: " << nStartingPoints << std::endl;
    for (int i = 0; i < nStartingPoints; i++) {
      pCache->setPoint(scanner, i);
      scanner->scan2d();
    }
  }
}

///
/// Perform the prob scan.
///
/// \param scanner - the scanner to run the scan with
/// \param cId - the id of this combination on the command line
///
void GammaComboEngine::make1dProbScan(MethodProbScan* scanner, int cId) {
  // load start parameters
  ParameterCache* pCache = new ParameterCache(arg);
  loadStartParameters(scanner, pCache, cId);

  scanner->initScan();
  scanStrategy1d(scanner, pCache);
  std::cout << "\nResults:" << std::endl;
  std::cout << "========\n" << std::endl;
  scanner->printLocalMinima();
  scanner->saveLocalMinima(m_fnamebuilder->getFileNameSolution(scanner));
  scanner->computeCLvalues();
  if (!arg->confirmsols) scanner->calcCLintervals();
  if (arg->cls.size() > 0) scanner->calcCLintervals(1);  // for prob method CLsType>1 doesn't exist
  if (!arg->isAction("pluginbatch") && !arg->plotpluginonly) {
    if (arg->plotpulls) scanner->plotPulls();
    if (arg->parevol) {
      ParameterEvolutionPlotter plotter(scanner);
      plotter.plotParEvolution();
    }
    if (isScanVarObservable(scanner->getCombiner(), arg->var[0])) {
      ParameterEvolutionPlotter plotter(scanner);
      plotter.plotObsScanCheck();
    }
    if (!arg->isAction("plugin")) {
      scanner->saveScanner(m_fnamebuilder->getFileNameScanner(scanner));
      pCache->cacheParameters(scanner, m_fnamebuilder->getFileNamePar(scanner));
    }
  }
}

///
/// Perform the 1D plugin scan. Runs toys in batch mode, and
/// reads them back in.
///
/// \param scannerPlugin - the scanner to run the scan with
/// \param cId - the id of this combination on the command line
///
void GammaComboEngine::make1dPluginScan(MethodPluginScan* scannerPlugin, int cId) {
  scannerPlugin->initScan();
  if (arg->isAction("pluginbatch")) {
    scannerPlugin->scan1d(arg->nrun);
  } else {
    scannerPlugin->readScan1dTrees(arg->jmin[cId], arg->jmax[cId]);
    scannerPlugin->calcCLintervals();
    for (int i = 0; i < arg->cls.size(); i++) {
      scannerPlugin->calcCLintervals(arg->cls[i]);
      if (arg->cls[i] == 2) scannerPlugin->calcCLintervals(arg->cls[i], true);  // calculate expected upper limit
    }
  }
  if (!arg->isAction("pluginbatch")) { scannerPlugin->saveScanner(m_fnamebuilder->getFileNameScanner(scannerPlugin)); }
}

///
/// Perform the 2D plugin scan. Runs toys in batch mode, and
/// reads them back in.
///
/// \param scannerPlugin - the scanner to run the scan with
/// \param cId - the id of this combination on the command line
///
void GammaComboEngine::make2dPluginScan(MethodPluginScan* scannerPlugin, int cId) {
  scannerPlugin->initScan();
  if (arg->isAction("pluginbatch")) {
    scannerPlugin->scan2d(arg->nrun);
  } else {
    scannerPlugin->readScan2dTrees(arg->jmin[cId], arg->jmax[cId]);
    scannerPlugin->saveScanner(m_fnamebuilder->getFileNameScanner(scannerPlugin));
    // plot chi2
    std::cout << "making full chi2 plot ..." << std::endl;
    OneMinusClPlot2d* plotf =
        new OneMinusClPlot2d(arg, plot->getName() + "_plugin_full", "p-value histogram: " + scannerPlugin->getTitle());
    scannerPlugin->plotOn(plotf);
    plotf->DrawFull();
    plotf->save();
  }
}

///
/// Perform the 1D berger boos scan. Runs toys in batch mode, and
/// reads them back in.
///
/// \param scannerPlugin - the scanner to run the scan with
/// \param cId - the id of this combination on the command line
///
void GammaComboEngine::make1dBergerBoosScan(MethodBergerBoosScan* scannerBergerBoos, int cId) {
  scannerBergerBoos->initScan();
  scannerBergerBoos->setNBergerBoosPointsPerScanpoint(arg->nBBpoints);
  if (arg->isAction("bbbatch")) {
    scannerBergerBoos->scan1d(arg->nrun);
  } else {
    scannerBergerBoos->readScan1dTrees(arg->jmin[cId], arg->jmax[cId]);
    scannerBergerBoos->calcCLintervals();
  }
  if (!arg->isAction("bbbatch")) {
    scannerBergerBoos->saveScanner(m_fnamebuilder->getFileNameScanner(scannerBergerBoos));
  }
}

/// Perform the coverage scanner
///
/// \param scanner - the scanner to run the scan with
/// \param cId - the id of this combination on the command line
void GammaComboEngine::make1dCoverageScan(MethodCoverageScan* scanner, int cId) {
  // load coverage point parameters (this can be done automatically)
  ParameterCache* pCache = new ParameterCache(arg);
  if (arg->loadParamsFile.size() != arg->combid.size()) {
    std::cout
        << "\nERROR : For a Coverage scan you must pass a parameter file (--parfile) to throw the toys from. You need "
           "one parfile per combiner"
        << std::endl;
    std::exit(1);
  }
  pCache->loadPoints(arg->loadParamsFile[cId]);

  // do scan
  scanner->initScan();
  scanner->setParameterCache(pCache);  // this can be passed directly to scan
  if (arg->isAction("coveragebatch")) {
    scanner->scan1d(arg->nrun);
  } else {
    scanner->readScan1dTrees(arg->jmin[cId], arg->jmax[cId]);
    scanner->saveScanner(
        m_fnamebuilder->getFileNameScanner(scanner).ReplaceAll("Coverage", Form("id%d_Coverage", arg->id)));
  }
}

///
/// Make a 1D plot of a Prob scanner. The scanner can either be a "fresh"
/// one, as made in make1dProbScan(), or a loaded one.
/// \param scanner - the scanner to plot
/// \param cId - the id of this combination on the command line
///
void GammaComboEngine::make1dProbPlot(MethodProbScan* scanner, int cId) {

  if (!arg->isAction("pluginbatch") && !arg->plotpluginonly) {
    scanner->setDrawSolution(arg->plotsolutions[cId]);
    if (arg->isAction("plugin") || arg->isAction("plot"))
      scanner->computeCLvalues();  // compute new CL values depending on test stat, even if not a rescan is wished
    if (arg->cls.size() > 0) {
      if (runOnDataSet) ((MethodDatasetsProbScan*)scanner)->plotFitRes(m_fnamebuilder->getFileNamePlot(cmb) + "_fit");
      scanner->plotOn(plot, 1);  // for prob ClsType>1 doesn't exist
    }
    scanner->plotOn(plot);
    int colorId = cId;
    if (arg->color.size() > cId) colorId = arg->color[cId];
    // scanner->setLineColor(colorsLine[colorId]);
    scanner->setTextColor(colorsText[colorId]);
    scanner->setLineColor(lineColors[cId]);
    scanner->setLineStyle(lineStyles[cId]);
    scanner->setLineWidth(lineWidths[cId]);
    scanner->setFillColor(fillColors[cId]);
    scanner->setFillStyle(fillStyles[cId]);
    scanner->setFillTransparency(fillTransparencies[cId]);
    plot->Draw();
  }
}

///
/// Perform a 1D scan with the following scan strategy.
/// If no starting values loaded do the default thing:
/// 1. scan once, using the start parameters found in the ParameterAbs-derived parameter class
/// 2. scan again from each solution found in the first step
///
void GammaComboEngine::scanStrategy1d(MethodProbScan* scanner, ParameterCache* pCache) {
  std::cout << "Scan strategy:" << std::endl;
  std::cout << "==============\n" << std::endl;
  int nStartingPoints = pCache->getNPoints();
  if (nStartingPoints == 0) {
    std::cout << "1. perform an initial scan" << std::endl;
    std::cout << "2. perform an additional scan starting from each solution found\n" << std::endl;
    std::cout << "first scan ..." << std::endl;
    scanner->scan1d();
    if (!arg->probforce) {
      std::vector<RooSlimFitResult*> firstScanSolutions = scanner->getSolutions();
      for (int i = 0; i < firstScanSolutions.size(); i++) {
        std::cout << "Scan i: " << i << std::endl;
        // scanner->loadSolution(i);
        scanner->loadParameters(firstScanSolutions[i]);
        scanner->scan1d(true);
      }
    }
  }
  // otherwise load each starting value found
  else {
    std::cout << "Scanning from each point found in start parameter file.\n" << std::endl;
    for (int i = 0; i < nStartingPoints; i++) {
      std::cout << "scan " << i + 1 << " of " << nStartingPoints << " ..." << std::endl;
      pCache->setPoint(scanner, i);
      scanner->scan1d();
    }
  }
}

///
/// Make the default 1D plugin plot:
///  - curve for the prob scan
///  - points with error bars for the plugin scan
///
/// \param sPlugin - the plugin scanner
/// \param sProb - the prob scanner
/// \param cId - the id of this combination on the command line
///
void GammaComboEngine::make1dPluginPlot(MethodPluginScan* sPlugin, MethodProbScan* sProb, int cId) {
  if (arg->isQuickhack(17)) {
    make1dPluginOnlyPlot(sPlugin, cId);
    sProb->setLineColor(kBlack);
    sProb->setDrawSolution(arg->plotsolutions[cId]);
    if (arg->cls.size() > 0) sProb->plotOn(plot, 1);
    sProb->plotOn(plot);
  } else {
    make1dProbPlot(sProb, cId);
    sPlugin->setLineColor(kBlack);
    sPlugin->setDrawSolution(arg->plotsolutions[cId]);
    if (std::count(arg->cls.begin(), arg->cls.end(), 1)) sPlugin->plotOn(plot, 1);
    if (std::count(arg->cls.begin(), arg->cls.end(), 2)) sPlugin->plotOn(plot, 2);
    sPlugin->plotOn(plot);
  }
  plot->Draw();
}

///
/// Make the default 2D plugin plot:
///  - contours for the the prob scan
///  - contours for the the plugin scan
///
/// \param sPlugin - the plugin scanner
/// \param sProb - the prob scanner
/// \param cId - the id of this combination on the command line
///
void GammaComboEngine::make2dPluginPlot(MethodPluginScan* sPlugin, MethodProbScan* sProb, int cId) {
  if (arg->isQuickhack(18)) {
    sProb->setTitle(sProb->getTitle() + "PROB");
    sPlugin->setTitle(sPlugin->getTitle() + "PLUGIN");
  } else {
    sProb->setTitle(sProb->getTitle() + " (Prob)");
    sPlugin->setTitle(sPlugin->getTitle() + " (Plugin)");
  }
  sProb->setDrawSolution(arg->plotsolutions[cId]);
  // sProb->setLineColor(colorsLine[cId]);
  sProb->setLineColor(lineColors[cId]);
  sProb->setLineStyle(lineStyles[cId]);
  sProb->setLineWidth(lineWidths[cId]);
  sProb->setFillColor(fillColors[cId]);
  sProb->setFillStyle(fillStyles[cId]);
  sProb->setFillTransparency(fillTransparencies[cId]);
  sPlugin->setDrawSolution(arg->plotsolutions[cId]);
  if (arg->isQuickhack(17)) {
    sPlugin->plotOn(plot);
    sProb->plotOn(plot);
  } else {
    sProb->plotOn(plot);
    sPlugin->plotOn(plot);
  }
  plot->Draw();
}

///
/// Make the plugin-only plot that only shows the plugin
/// as continuous lines.
///
/// \param sPlugin - the plugin scanner
/// \param cId - the id of this combination on the command line
///
void GammaComboEngine::make1dPluginOnlyPlot(MethodPluginScan* sPlugin, int cId) {
  ((OneMinusClPlot*)plot)->setPluginMarkers(false);
  int colorId = cId;
  if (arg->color.size() > cId) colorId = arg->color[cId];
  sPlugin->setTextColor(colorsText[colorId]);
  sPlugin->setLineColor(lineColors[cId]);
  sPlugin->setLineStyle(lineStyles[cId]);
  sPlugin->setLineWidth(lineWidths[cId]);
  sPlugin->setFillColor(fillColors[cId]);
  sPlugin->setFillStyle(fillStyles[cId]);
  sPlugin->setFillTransparency(fillTransparencies[cId]);
  sPlugin->setDrawSolution(arg->plotsolutions[cId]);
  for (int i = 0; i < arg->cls.size(); i++) sPlugin->plotOn(plot, arg->cls[i]);
  sPlugin->plotOn(plot);
  plot->Draw();
}

///
/// Make the plugin-only 2D plot.
///
/// \param sPlugin - the plugin scanner
/// \param cId - the id of this combination on the command line
///
void GammaComboEngine::make2dPluginOnlyPlot(MethodPluginScan* sPlugin, int cId) {
  sPlugin->setDrawSolution(arg->plotsolutions[cId]);
  sPlugin->plotOn(plot);
  plot->Draw();
}

///
/// Make the 1D coverage plot
///
/// \param scanner - the coverage scanner
/// \param cId - the id of this combination on the command line
///
void GammaComboEngine::make1dCoveragePlot(MethodCoverageScan* scanner, int cId) { scanner->plot(); }

///
/// Make a 2D prob scan.
/// - load start parameters
/// - perform scan
/// - save scanner and parameters
///
/// \param scanner - the scanner
/// \param cId - the id of this combination on the command line
///
void GammaComboEngine::make2dProbScan(MethodProbScan* scanner, int cId) {
  // load start parameters
  ParameterCache* pCache = new ParameterCache(arg);
  loadStartParameters(scanner, pCache, cId);
  // scan
  scanner->initScan();
  scanStrategy2d(scanner, pCache);
  std::cout << std::endl;
  scanner->printLocalMinima();
  // save
  scanner->saveScanner(m_fnamebuilder->getFileNameScanner(scanner));
  pCache->cacheParameters(scanner, m_fnamebuilder->getFileNamePar(scanner));
}

///
/// Make the 2D plot for a prob scanner.
///
void GammaComboEngine::make2dProbPlot(MethodProbScan* scanner, int cId) {
  // plot full
  OneMinusClPlot2d* plotf;
  if (scanner->getMethodName() == "Prob")
    plotf = new OneMinusClPlot2d(arg, m_fnamebuilder->getFileNamePlotSingle(cmb, cId) + "_full",
                                 "p-value histogram: " + scanner->getTitle());
  else if (scanner->getMethodName() == "DatasetsProb")
    plotf = new OneMinusClPlot2d(arg, m_fnamebuilder->getFileNamePlot(cmb) + "_full",
                                 "p-value histogram: " + scanner->getTitle());  // Titus: change to make datasets plot
                                                                                // possible
  else
    std::cout << "The name of the scanner matches neither Prob nor DatasetsProb!" << std::endl;
  scanner->plotOn(plotf, 0);
  plotf->DrawFull();
  plotf->save();
  // plot full CLs
  if (arg->cls.size() > 0) {
    OneMinusClPlot2d* plotfcls;
    if (scanner->getMethodName() == "Prob")
      plotfcls = new OneMinusClPlot2d(arg, m_fnamebuilder->getFileNamePlotSingle(cmb, cId) + "_cls_full",
                                      "p-value histogram: " + scanner->getTitle());
    else if (scanner->getMethodName() == "DatasetsProb")
      plotfcls = new OneMinusClPlot2d(arg, m_fnamebuilder->getFileNamePlot(cmb) + "_cls_full",
                                      "p-value histogram: " + scanner->getTitle());  // Titus: change to make datasets
                                                                                     // plot possible
    else
      std::cout << "The name of the scanner matches neither Prob nor DatasetsProb!" << std::endl;
    scanner->plotOn(plotfcls, 1);
    plotfcls->DrawFull();
    plotfcls->save();
  }
  // contour plot
  scanner->setDrawSolution(arg->plotsolutions[cId]);
  // scanner->setLineColor(colorsLine[cId]);
  scanner->setLineColor(lineColors[cId]);
  scanner->setLineStyle(lineStyles[cId]);
  scanner->setLineWidth(lineWidths[cId]);
  scanner->setFillColor(fillColors[cId]);
  scanner->setFillStyle(fillStyles[cId]);
  scanner->setFillTransparency(fillTransparencies[cId]);
  if (arg->cls.size() > 0) scanner->plotOn(plot, 1);
  scanner->plotOn(plot, 0);
  // only draw the plot once when multiple scanners are plotted,
  // else we end up with too many graphs, and the transparency setting
  // gets screwed up
  // Titus: also draw the plot, if no combiner is set (datasets case)
  if (cId == arg->combid.size() - 1 || arg->combid.empty()) {
    plot->Draw();
    plot->Show();
  }
}

///
/// Helper function for scan(). Fixes parameters, if requested
/// (only possible before combining).
///
void GammaComboEngine::fixParameters(Combiner* c, int cId) {
  if (cId < arg->fixParameters.size()) {
    for (int j = 0; j < arg->fixParameters[cId].size(); j++) {
      c->fixParameter(arg->fixParameters[cId][j].name, arg->fixParameters[cId][j].value);
    }
  }
}

///
/// Helper function for scan(). Adjusts ranges, if requested
/// (only possible before combining).
///
void GammaComboEngine::adjustRanges(Combiner* c, int cId) {
  if (cId < arg->physRanges.size()) {
    for (int j = 0; j < arg->physRanges[cId].size(); j++) {
      c->adjustPhysRange(arg->physRanges[cId][j].name, arg->physRanges[cId][j].min, arg->physRanges[cId][j].max);
    }
  }
  if (cId < arg->removeRanges.size()) {
    for (int j = 0; j < arg->removeRanges[cId].size(); j++) {
      if (arg->removeRanges[cId][j] == "all") {
        for (const auto& par : *c->getParameters()) static_cast<RooRealVar*>(par)->removeRange();
      } else {
        c->adjustPhysRange(arg->removeRanges[cId][j], -999, -999);
      }
    }
  }
}

///
/// Helper function for scan(): Makes named sets for any toy variations that were requested
///
void GammaComboEngine::setupToyVariationSets(Combiner* c, int cId) {
  if (cId < arg->randomizeToyVars.size()) {
    TString toyVarList = "";
    for (int j = 0; j < arg->randomizeToyVars[cId].size(); j++) {
      toyVarList += arg->randomizeToyVars[cId][j];
      if (j < arg->randomizeToyVars[cId].size() - 1) toyVarList += ",";
    }
    c->getWorkspace()->defineSet("toy_" + c->getPdfName(), toyVarList.Data());
  }
}

///
/// Helper function for scan(): Checks if for a given combid (the
/// running index of the -c argument) a start parameter file was
/// configured (-l) argument. If so, it is returned, else the default
/// name is returned.
///
TString GammaComboEngine::getStartParFileName(int cId) {
  if (arg->loadParamsFile.size() <= cId) return m_fnamebuilder->getFileNameStartPar(cmb[cId]);
  if (arg->loadParamsFile[cId].EqualTo("default")) return m_fnamebuilder->getFileNameStartPar(cmb[cId]);
  return arg->loadParamsFile[cId];
}

///
/// Checks if a given variable name is in the list of observables
/// of a combiner. The check is slightly sloppy as it ignores the
/// unique ID.
///
/// \param c        - the combiner
/// \param scanVar  - the scan variable name
/// \return true if included, else false
///
bool GammaComboEngine::isScanVarObservable(Combiner* c, TString scanVar) {
  std::vector<std::string> obs = c->getObservableNames();
  for (int i = 0; i < obs.size(); i++) {
    if (scanVar.Contains(obs[i])) return true;
  }
  return false;
}

///
/// Helper function to set up a scan for an observable, tightens
/// the chi2 constraint.
///
/// \param c        - the combiner
/// \param scanVar  - the scan variable name (must be an observable)
///
void GammaComboEngine::tightenChi2Constraint(Combiner* c, TString scanVar) {
  std::cout << "\n--var " << scanVar << ": Setting up a scan for an observable ..." << std::endl;
  PDF_Abs* pdf = c->getPdfProvidingObservable(scanVar);
  if (pdf == 0) {
    std::cout << "GammaComboEngine::tightenChi2Constraint() : ERROR : no PDF found that contains the observable '"
              << scanVar << "'. Exit." << std::endl;
    std::exit(1);
  }
  float scale = 0.1;
  std::cout << "... observable error is multiplied by a factor " << scale << std::endl;
  pdf->ScaleError(scanVar, scale);
  pdf->buildCov();
  pdf->buildPdf();
}

// FIXME
// WARNING - THIS FUNCTION ALLOWS YOU DO SOME INCREDIBLY STUPID THINGS
//         - SO PLEASE BE CAREFUL WITH IT!
/// Helper function for scan(). Set observables to values from file
/// (only possible before combining).
///
void GammaComboEngine::setObservablesFromFile(Combiner* c, int cId) {

  if (cId >= arg->readfromfile.size()) return;
  if (arg->readfromfile[cId].size() == 0) return;
  if (arg->readfromfile[cId].size() == 1 && arg->readfromfile[cId][0] == TString("default")) return;
  if (arg->readfromfile[cId].size() == 1 && arg->readfromfile[cId][0] == TString("")) return;

  std::vector<PDF_Abs*> pdfs = c->getPdfs();

  if (pdfs.size() != arg->readfromfile[cId].size()) {
    std::cout
        << "ERROR -- I think you are trying to read values from a file but you haven't specified one for each PDF in "
           "the combination"
        << std::endl;
    std::cout << "      -- Pass like this: -c 0:+1,+2,+3 --readfromfile "
                 "default,<path_to_pdf2_vals.dat>,<path_to_pdf3_vals.dat>"
              << std::endl;
    std::exit(1);
  }

  for (int i = 0; i < pdfs.size(); i++) {

    // read from the file for this pdf
    if (arg->readfromfile[cId][i] == TString("default")) continue;

    //
    std::ifstream infile(arg->readfromfile[cId][i].Data());
    if (!infile.is_open()) {
      std::cerr << "No such read file found: " << arg->readfromfile[cId][i] << std::endl;
      std::exit(1);
    }
    std::string line;
    bool pdfExists = false;
    bool pdfFound = false;
    while (getline(infile, line)) {
      if (line.empty()) continue;                   // blank line
      if (boost::starts_with(line, "#")) continue;  // these are comments
      if (boost::starts_with(line, Form("pdf: %s", pdfs[i]->getBaseName().Data()))) {
        pdfFound = true;   // this is the pdf we are looking for
        pdfExists = true;  // keep track of whether it's even there or not
      } else if (boost::starts_with(line, "pdf:"))
        pdfFound = false;  // this is some other pdf after the pdf we are looking for

      if (pdfFound) {
        std::vector<std::string> els;
        boost::split(els, line, boost::is_any_of(" "), boost::token_compress_on);
        TString typ = els[0];
        if (typ == "obs:") {
          TString name = els[1];
          double val = boost::lexical_cast<double>(els[2]);
          pdfs[i]->setObservable(name, val);
          pdfs[i]->obsValSource = "Read from file " + arg->readfromfile[cId][i];
        } else if (typ == "err:") {
          TString name = els[1];
          double stat = boost::lexical_cast<double>(els[2]);
          double syst = boost::lexical_cast<double>(els[3]);
          pdfs[i]->setUncertainty(name, stat, syst);
          pdfs[i]->obsErrSource = "Read from file " + arg->readfromfile[cId][i];
        } else if (typ == "cor:") {
          int mi = boost::lexical_cast<int>(els[1]);
          int mj = boost::lexical_cast<int>(els[2]);
          double corStat = boost::lexical_cast<double>(els[3]);
          double corSyst = boost::lexical_cast<double>(els[4]);
          pdfs[i]->corStatMatrix[mi][mj] = corStat;
          pdfs[i]->corSystMatrix[mi][mj] = corSyst;
          pdfs[i]->corStatMatrix[mj][mi] = corStat;
          pdfs[i]->corSystMatrix[mj][mi] = corSyst;
          pdfs[i]->corSource = "Read from file " + arg->readfromfile[cId][i];
        }
      }
    }
    infile.close();
    if (pdfExists) {
      pdfs[i]->storeErrorsInObs();
      pdfs[i]->buildCov();
      pdfs[i]->buildPdf();
    } else {
      std::cout << "WARNING - did not find any pdf named: " << pdfs[i]->getBaseName()
                << " in file: " << arg->readfromfile[cId][i] << std::endl;
    }
  }
}

///
/// write batch scripts
///
void GammaComboEngine::writebatchscripts() {
  if (runOnDataSet) {
    m_batchscriptwriter->writeScripts_datasets(arg, getPdf(0));
  } else {
    m_batchscriptwriter->writeScripts(arg, &cmb);
  }
  std::exit(0);
}

///
/// make latex
///
void GammaComboEngine::makeLatex(Combiner* c) {
  for (unsigned int p = 0; p < c->getPdfs().size(); p++) {
    PDF_Abs* pdf = c->getPdfs()[p];
    LatexMaker m(c->getName(), pdf);
    m.writeFile();
  }
}

///
/// save workspace
///
void GammaComboEngine::saveWorkspace(Combiner* c, int i) {
  // if --pr then make the ranges
  if (arg->enforcePhysRange) Utils::setLimit(c->getParameters(), "phys");

  // first write a copy THE PDF into the workspace
  const RooAbsPdf* thePdf = c->getPdf();
  const RooArgSet* theParameters = c->getParameters();
  const RooArgSet* theObservables = c->getObservables();

  RooAbsPdf* savePdf = (RooAbsPdf*)thePdf->Clone("ThePdf");
  c->getWorkspace()->import(*savePdf, RooFit::Silence());
  c->getWorkspace()->defineSet("TheParameters", *theParameters);
  c->getWorkspace()->defineSet("TheObservables", *theObservables);

  // set the name of the workspace associated with combiner id
  int cId = arg->combid[i];
  c->getWorkspace()->SetName(Form("w%d", cId));

  // if the first instance then write directly to file
  if (i == 0) {
    c->getWorkspace()->writeToFile(arg->save);
  }
  // otherwise add it to the file
  else {
    TFile* tf = TFile::Open(arg->save, "UPDATE");
    c->getWorkspace()->Write();
    tf->Close();
  }
}

///
/// compare combinations
///
void GammaComboEngine::compareCombinations() {
  for (int i = 0; i < comparisonScanners.size(); i++) {
    for (int j = i + 1; j < comparisonScanners.size(); j++) {
      TH2F* pull_corr = new TH2F(Form("c%d_vs_c%d_corr", i, j),
                                 Form("; Obs pulls for %s [#sigma]; Obs pulls %s [#sigma]",
                                      comparisonScanners[i]->getName().Data(), comparisonScanners[j]->getName().Data()),
                                 10, -5, 5, 10, -5, 5);
      std::vector<double> pullVec1;
      std::vector<double> pullVec2;
      double total_pull = 0.;
      int nmatch = 0;
      comparisonScanners[i]->loadSolution(0);
      comparisonScanners[j]->loadSolution(0);
      const auto sc1obs = comparisonScanners[i]->getObservables();
      const auto sc2obs = comparisonScanners[j]->getObservables();
      for (const auto& pAbsObs1 : *sc1obs) {
        const auto pObs1 = static_cast<RooRealVar*>(pAbsObs1);
        for (const auto& pAbsObs2 : *sc2obs) {
          const auto pObs2 = static_cast<RooRealVar*>(pAbsObs1);

          // look for matches
          TString pTh1Name = pObs1->GetName();
          pTh1Name.ReplaceAll("obs", "th");
          pTh1Name.ReplaceAll("UID", ";");
          pTh1Name = ((TObjString*)pTh1Name.Tokenize(";")->At(0))->GetString();
          TString pTh2Name = pObs2->GetName();
          pTh2Name.ReplaceAll("obs", "th");
          pTh2Name.ReplaceAll("UID", ";");
          pTh2Name = ((TObjString*)pTh2Name.Tokenize(";")->At(0))->GetString();

          if (pTh1Name == pTh2Name) {
            nmatch += 1;
            pTh1Name = pObs1->GetName();
            pTh1Name.ReplaceAll("obs", "th");
            pTh2Name = pObs2->GetName();
            pTh2Name.ReplaceAll("obs", "th");
            RooRealVar* pTh1 = (RooRealVar*)comparisonScanners[i]->getTheory()->find(pTh1Name);
            RooRealVar* pTh2 = (RooRealVar*)comparisonScanners[j]->getTheory()->find(pTh2Name);
            assert(pTh1 && pTh2);
            double pull = (pTh1->getVal() - pTh2->getVal()) / pObs2->getError();
            total_pull += pull * pull;
            double pull1 = (pTh1->getVal() - pObs1->getVal()) / pObs1->getError();
            double pull2 = (pTh2->getVal() - pObs2->getVal()) / pObs2->getError();
            pull_corr->Fill(pull1, pull2);
            pullVec1.push_back(pull1);
            pullVec2.push_back(pull2);
          }
        }
      }
      double chi21 = comparisonScanners[i]->getSolution(0)->minNll();
      double chi22 = comparisonScanners[j]->getSolution(0)->minNll();
      int nObs1 = comparisonScanners[i]->getObservables()->getSize();
      int nObs2 = comparisonScanners[j]->getObservables()->getSize();
      int nPar1 = comparisonScanners[i]->getSolution(0)->floatParsFinal().getSize();
      int nPar2 = comparisonScanners[j]->getSolution(0)->floatParsFinal().getSize();
      CLInterval cl1 = comparisonScanners[i]->getCLinterval(0, 1);
      CLInterval cl2 = comparisonScanners[j]->getCLinterval(0, 1);
      double diff = cl1.central - cl2.central;
      double corr = Utils::getCorrelationFactor(pullVec1, pullVec2);
      double err = 99.;
      if (diff > 0) {  // means value has moved down
        using Utils::sq;
        err = sqrt(sq(cl1.central - cl1.min) + sq(cl2.max - cl2.central) -
                   2. * corr * (cl1.central - cl1.min) * (cl2.max - cl2.central));
      } else {
        using Utils::sq;
        err = sqrt(sq(cl1.max - cl1.central) + sq(cl2.central - cl2.min) -
                   2. * corr * (cl1.max - cl1.central) * (cl2.central - cl2.min));
      }
      std::cout << "Comparison   1): " << Form("%-20s", comparisonScanners[i]->getName().Data())
                << " to 2): " << Form("%-20s", comparisonScanners[j]->getName().Data()) << std::endl;
      std::cout << "        chi2:    " << Form("%-20.3f", chi21) << "        " << Form("%-20.3f", chi22) << std::endl;
      std::cout << "        nObs:    " << Form("%-20d", nObs1) << "        " << Form("%-20d", nObs2) << std::endl;
      std::cout << "        nPar:    " << Form("%-20d", nPar1) << "        " << Form("%-20d", nPar2) << std::endl;
      std::cout << "        pval:    " << Form("%-20.2f", 100. * TMath::Prob(chi21, nObs1 - nPar1)) << "        "
                << Form("%-20.2f", 100. * TMath::Prob(chi22, nObs2 - nPar2)) << std::endl;
      std::cout << "         val:    " << Form("%-6.3f", cl1.central) << " [" << Form("%6.3f", cl1.min) << ","
                << Form("%-6.3f", cl1.max) << "]"
                << "      " << Form("%-6.3f", cl2.central) << "[" << Form("%6.3f", cl2.min) << ","
                << Form("%-6.3f", cl2.max) << "]" << std::endl;
      std::cout << "        dval:    " << diff << " +/- " << err << " (" << TMath::Abs(diff) / err << " sigma)"
                << std::endl;
      std::cout << "PULL PER OBS:  " << total_pull << std::endl;
      std::cout << "CORRELATION:   " << corr << std::endl;
      std::cout << "COMPATIBILITY: " << TMath::Abs(diff) / err << " sigma" << std::endl;

      TCanvas* canv = Utils::newNoWarnTCanvas("pull_corr" + Utils::getUniqueRootName());
      pull_corr->SetMarkerStyle(kMultiply);
      pull_corr->SetMarkerColor(kBlue + 2);
      pull_corr->GetXaxis()->SetTitleSize(0.045);
      pull_corr->GetYaxis()->SetTitleSize(0.045);
      pull_corr->GetXaxis()->SetLabelSize(0.045);
      pull_corr->GetYaxis()->SetLabelSize(0.045);
      pull_corr->Draw("scat");
      TLine* line = new TLine();
      line->DrawLine(-5, 0, 5, 0);
      line->DrawLine(0, -5, 0, 5);
      pull_corr->Draw("scatsame");
      TF1* f1 = new TF1("f1", "[0]*x", -5, 5);
      f1->SetParameter(0, corr);
      f1->Draw("Lsame");
      TLatex* lat = new TLatex();
      lat->DrawLatex(3, 4, Form("#rho = %3.1f", corr));
      lat->DrawLatex(3, 3, Form("#sigma = %3.1f", TMath::Abs(diff) / err));
      Utils::savePlot(canv, Form("pull_corr_%s_%s", comparisonScanners[i]->getName().Data(),
                                 comparisonScanners[j]->getName().Data()));
      total_pull = sqrt(total_pull) / nmatch;
    }
  }
}

///
/// run toys
///
void GammaComboEngine::runToys(Combiner* c) {
  // THIS IS A HACK FOR NOW
  //
  // base scan (overhead here)
  MethodProbScan* probscan = new MethodProbScan(c);
  make1dProbScan(probscan, 0);

  TString toydirname = TString("root/scan1dToys_") + probscan->getName() + TString("_") + probscan->getScanVar1Name();
  TString toyfname = toydirname + TString("/scan1dToys_") + probscan->getName() + TString("_") +
                     probscan->getScanVar1Name() + Form("_run%d.root", arg->nrun);

  TTree* tree = new TTree("toys", "toys");
  std::map<TString, double> vals;
  std::map<TString, double> errs;
  double chi2val = probscan->solutions[0]->minNll();
  int ntoy = -1;
  tree->Branch("ntoy", &ntoy);
  tree->Branch("chi2min", &chi2val);
  RooArgList fitPars = probscan->solutions[0]->floatParsFinal();
  for (const auto& pAbs : fitPars) {
    const auto p = static_cast<RooRealVar*>(pAbs);
    std::cout << "YO:: " << p->GetName() << std::endl;
    vals[p->GetName()] = p->getVal();
    errs[p->GetName()] = p->getError();
    tree->Branch(p->GetName() + TString("_val"), &vals[p->GetName()]);
    tree->Branch(p->GetName() + TString("_err"), &errs[p->GetName()]);
  }
  tree->Fill();
  delete probscan;

  for (int i = 0; i < arg->ntoys; i++) {
    std::cout << "RUNNING TOY " << i << " / " << arg->ntoys << std::endl;
    c->setObservablesToToyValues();
    MethodProbScan* toyscan = new MethodProbScan(c);
    make1dProbScan(toyscan, 0);
    if (toyscan->solutions.size() == 0) continue;
    toyscan->solutions[0]->Print();
    // std::cout << "My stuff I want to save" << std::endl;
    // std::cout << "chi2: " << toyscan->solutions[0]->minNll() << std::endl;

    chi2val = toyscan->solutions[0]->minNll();
    ntoy = i;
    RooArgList toyFitPars = toyscan->solutions[0]->floatParsFinal();
    for (auto const& pAbs : toyFitPars) {
      const auto p = static_cast<RooRealVar*>(pAbs);
      // std::cout << p->GetName() << " " << p->getVal() << " " << p->getError() << std::endl;
      vals[p->GetName()] = p->getVal();
      errs[p->GetName()] = p->getError();
    }
    tree->Fill();
    delete toyscan;
  }
  system("mkdir -p " + toydirname);
  TFile* f = new TFile(toyfname, "recreate");
  tree->Write();
  f->Close();
  delete tree;
  delete f;
}

///
/// scan engine
///
void GammaComboEngine::scan() {
  // if we're running with the dataset option then we go off and do that somewhere else
  if (runOnDataSet) {
    scanDataSet();
    return;
  }

  // combination scanning action happens here
  for (int i = 0; i < arg->combid.size(); i++) {
    int combinerId = arg->combid[i];
    Combiner* c = cmb[combinerId];

    // read observable values, uncertainties and correlations from a file
    setObservablesFromFile(c, i);

    // work with a clone - this way we can easily make plots with the
    // same combination in twice (once with asimov, for example)
    c = c->Clone(c->getName(), c->getTitle());

    // fix parameters according to the command line - only possible before combining
    fixParameters(c, i);

    // configure names to run an Asimov toy - only possible before combining
    if (arg->isAsimovCombiner(i)) configureAsimovCombinerNames(c, i);

    // configure scans for observables - this part is only possible before combining
    if (isScanVarObservable(c, arg->var[0])) { tightenChi2Constraint(c, arg->var[0]); }
    if (arg->var.size() == 2 && isScanVarObservable(c, arg->var[1])) { tightenChi2Constraint(c, arg->var[1]); }

    // combine
    c->combine();
    if (!c->isCombined()) continue;  // error during combining

    // adjust ranges according to the command line - only possible after combining
    adjustRanges(c, i);

    // set up parameter sets for the parameters to vary within the toys (if requested)
    setupToyVariationSets(c, i);

    // make graphviz dot files
    printCombinerStructure(c);

    // set an asimov toy - only possible after combining
    if (arg->isAsimovCombiner(i)) loadAsimovPoint(c, i);

    // configure scans for observables - this part is only possible after combining
    // add the observable(s) to the list of parameters
    if (isScanVarObservable(c, arg->var[0])) { c->getWorkspace()->extendSet(c->getParsName(), arg->var[0]); }
    if (arg->var.size() == 2 && isScanVarObservable(c, arg->var[1])) {
      c->getWorkspace()->extendSet(c->getParsName(), arg->var[1]);
    }

    // printout and latex
    c->print();
    if (arg->debug) c->getWorkspace()->Print("v");
    if (arg->save != "" && !arg->saveAtMin) saveWorkspace(c, i);
    if (arg->latex) makeLatex(c);
    if (arg->info || arg->latex || (arg->save != "" && !arg->saveAtMin)) continue;

    /////////////////////////////////////////////////////
    //
    // PROB
    //
    /////////////////////////////////////////////////////

    if (!arg->isAction("plugin") && !arg->isAction("pluginbatch") && !arg->isAction("coverage") &&
        !arg->isAction("coveragebatch") && !arg->isAction("bb") && !arg->isAction("bbbatch")) {
      MethodProbScan* scannerProb = new MethodProbScan(c);
      // pvalue corrector
      if (arg->coverageCorrectionID > 0) {
        PValueCorrection* pvalueCorrector = new PValueCorrection(arg->coverageCorrectionID, arg->verbose);
        pvalueCorrector->readFiles(m_fnamebuilder->getFileBaseName(c), arg->coverageCorrectionPoint,
                                   false);  // false means for prob
        pvalueCorrector->write("root/pvalueCorrection_prob.root");
        scannerProb->setPValueCorrector(pvalueCorrector);
      }

      // 1D SCANS
      if (arg->var.size() == 1) {
        if (arg->isAction("plot")) {
          scannerProb->loadScanner(m_fnamebuilder->getFileNameScanner(scannerProb));
        } else {
          make1dProbScan(scannerProb, i);
        }
        make1dProbPlot(scannerProb, i);
        if (arg->compare) comparisonScanners.push_back(scannerProb);
      }
      // 2D SCANS
      else if (arg->var.size() == 2) {
        if (arg->isAction("plot")) {
          scannerProb->loadScanner(m_fnamebuilder->getFileNameScanner(scannerProb));
        } else {
          make2dProbScan(scannerProb, i);
        }
        make2dProbPlot(scannerProb, i);
      }
    }

    /////////////////////////////////////////////////////
    //
    // PLUGIN
    //
    /////////////////////////////////////////////////////

    if (arg->isAction("plugin") || arg->isAction("pluginbatch")) {
      // 1D SCANS
      if (arg->var.size() == 1) {
        if (arg->isAction("pluginbatch")) {
          MethodProbScan* scannerProb = new MethodProbScan(c);
          make1dProbScan(scannerProb, i);
          MethodPluginScan* scannerPlugin = new MethodPluginScan(scannerProb);
          make1dPluginScan(scannerPlugin, i);
        }
        // if ( arg->isAction("pluginhybridbatch") ){
        //// Hybrid Plugin: compute a second profile likelihood to define the parameter evolution
        // std::cout << "HYBRID PLUGIN: preparing profile likelihood to be used for parameter evolution:" << std::endl;
        // ParameterCache *pCache = new ParameterCache(arg, m_fnamebuilder->getFileBaseName(cmb[arg->pevid[0]]));
        // pCache->loadPoints();
        // MethodProbScan *scanner3 = new MethodProbScan(cmb[arg->pevid[0]]);
        // scanner3->initScan();
        // scanStrategy1d(scanner3,pCache);
        // scanner3->confirmSolutions();
        // scanner3->printLocalMinima();
        // scanner2->setParevolPLH(scanner3);
        // }
        else if (arg->isAction("plugin")) {
          // Create the Prob scanner: load from disc if it exists, else redo the scan.
          // We don't need the prob scanner for the plugin only plot, if we either just
          // want to replot it.
          MethodProbScan* scannerProb = new MethodProbScan(c);
          if (!arg->plotpluginonly || (arg->plotpluginonly && !arg->isAction("plot"))) {
            if (Utils::FileExists(m_fnamebuilder->getFileNameScanner(scannerProb))) {
              scannerProb->loadScanner(m_fnamebuilder->getFileNameScanner(scannerProb));
            } else {
              std::cout << "\nWARNING : Couldn't load the Prob scanner, will rerun the Prob" << std::endl;
              std::cout << "          scan now. You should have run the Prob scan locally" << std::endl;
              std::cout << "          before running the Plugin scan." << std::endl;
              std::cout << "          missing file: " << m_fnamebuilder->getFileNameScanner(scannerProb) << std::endl;
              std::cout << std::endl;
              make1dProbScan(scannerProb, i);
            }
          }
          // create Plugin scanner
          MethodPluginScan* scannerPlugin = new MethodPluginScan(scannerProb);
          if (arg->isAction("plot")) {
            scannerPlugin->loadScanner(m_fnamebuilder->getFileNameScanner(scannerPlugin));
          } else {
            if (arg->coverageCorrectionID > 0) {
              PValueCorrection* pvalueCorrector = new PValueCorrection(arg->coverageCorrectionID, arg->verbose);
              pvalueCorrector->readFiles(m_fnamebuilder->getFileBaseName(c), arg->coverageCorrectionPoint,
                                         true);  // true means for plugin
              pvalueCorrector->write("root/pvalueCorrection_plugin.root");
              scannerPlugin->setPValueCorrector(pvalueCorrector);
            }
            make1dPluginScan(scannerPlugin, i);
          }
          if (arg->plotpluginonly) {
            make1dPluginOnlyPlot(scannerPlugin, i);
          } else {
            make1dPluginPlot(scannerPlugin, scannerProb, i);
          }
        }

      }
      // 2D SCANS
      else if (arg->var.size() == 2) {
        if (arg->isAction("pluginbatch")) {
          MethodProbScan* scannerProb = new MethodProbScan(c);
          make2dProbScan(scannerProb, i);
          MethodPluginScan* scannerPlugin = new MethodPluginScan(scannerProb);
          make2dPluginScan(scannerPlugin, i);
        } else if (arg->isAction("plugin")) {
          MethodProbScan* scannerProb = new MethodProbScan(c);
          if (!(arg->isAction("plot") && arg->plotpluginonly)) {
            // we don't need the prob scanner if we just want to replot the plugin only
            scannerProb->loadScanner(m_fnamebuilder->getFileNameScanner(scannerProb));
          }
          MethodPluginScan* scannerPlugin = new MethodPluginScan(scannerProb);
          if (arg->isAction("plot")) {
            scannerPlugin->loadScanner(m_fnamebuilder->getFileNameScanner(scannerPlugin));
          } else {
            make2dPluginScan(scannerPlugin, i);
          }
          if (arg->plotpluginonly) {
            make2dPluginOnlyPlot(scannerPlugin, i);
          } else {
            make2dPluginPlot(scannerPlugin, scannerProb, i);
          }
        }
      }
    }

    /////////////////////////////////////////////////////
    //
    // COVERAGE
    //
    /////////////////////////////////////////////////////

    if (arg->isAction("coverage") || arg->isAction("coveragebatch")) {
      if (arg->var.size() != 1) {
        std::cerr << "ERROR -- you can only scan in 1D for a coverage check" << std::endl;
        std::exit(1);
      }
      MethodCoverageScan* coverageScan = new MethodCoverageScan(c);
      if (arg->isAction("coveragebatch")) {
        make1dCoverageScan(coverageScan, i);
      } else if (arg->isAction("coverage")) {
        if (arg->isAction("plot")) {
          if (Utils::FileExists(m_fnamebuilder->getFileNameScanner(coverageScan))) {
            coverageScan->loadScanner(m_fnamebuilder->getFileNameScanner(coverageScan));
          } else {
            std::cout << "\nERROR : Couldn't load the coverage scanner: "
                      << m_fnamebuilder->getFileNameScanner(coverageScan) << std::endl;
            std::exit(1);
          }
        } else {
          make1dCoverageScan(coverageScan, i);
        }
        make1dCoveragePlot(coverageScan, i);
      }
    }

    // RUN TOYS
    if (arg->isAction("runtoys")) runToys(c);
    /////////////////////////////////////////////////////

    // SAVE WORKSPACE
    if (arg->save != "" && arg->saveAtMin) saveWorkspace(c, i);
    /////////////////////////////////////////////////////

    if (i < arg->combid.size() - 1) {
      std::cout << "\n-- now starting -c " << arg->combid[i + 1]
                << " ------------------------------------------------------------------\n"
                << std::endl;
    }
  }
}
//
// special scan engine for datasetss
//
void GammaComboEngine::scanDataSet() {
  if (arg->info || arg->latex) return;

  /////////////////////////////////////////////////////
  //
  // PROB - DATASETS
  //
  /////////////////////////////////////////////////////

  if (!arg->isAction("plugin") && !arg->isAction("pluginbatch") && !arg->isAction("coverage") &&
      !arg->isAction("coveragebatch") && !arg->isAction("bb") && !arg->isAction("bbbatch")) {
    MethodDatasetsProbScan* probScanner = new MethodDatasetsProbScan((PDF_Datasets*)pdf[0], arg);

    // 1D SCANS
    if (arg->var.size() == 1) {
      if (arg->isAction("plot")) {
        probScanner->loadScanner(m_fnamebuilder->getFileNameScanner(probScanner));
      } else {
        make1dProbScan(probScanner, 0);
      }
      make1dProbPlot(probScanner, 0);
    }
    // 2D SCANS
    else if (arg->var.size() == 2) {
      if (arg->isAction("plot")) {
        probScanner->loadScanner(m_fnamebuilder->getFileNameScanner(probScanner));
      } else {
        make2dProbScan(probScanner, 0);
      }
      make2dProbPlot(probScanner, 0);
    }
  }

  /////////////////////////////////////////////////////
  //
  // PLUGIN - DATASETS
  //
  /////////////////////////////////////////////////////

  if (arg->isAction("pluginbatch") || arg->isAction("plugin")) {
    // 1D SCANS
    if (arg->var.size() == 1) {
      if (arg->isAction("pluginbatch")) {
        MethodDatasetsProbScan* scannerProb = new MethodDatasetsProbScan((PDF_Datasets*)pdf[0], arg);
        if (Utils::FileExists(m_fnamebuilder->getFileNameScanner(scannerProb))) {
          scannerProb->initScan();
          scannerProb->loadScanner(m_fnamebuilder->getFileNameScanner(scannerProb));
        } else {
          std::cout << "\nWARNING : Couldn't load the Prob scanner, will rerun the Prob" << std::endl;
          std::cout << "          scan now. You should have run the Prob scan locally" << std::endl;
          std::cout << "          before running the Plugin scan." << std::endl;
          std::cout << "          missing file: " << m_fnamebuilder->getFileNameScanner(scannerProb) << std::endl;
          std::cout << std::endl;
          make1dProbScan(scannerProb, 0);
        }
        MethodDatasetsPluginScan* scannerPlugin = new MethodDatasetsPluginScan(scannerProb, (PDF_Datasets*)pdf[0], arg);
        make1dPluginScan(scannerPlugin, 0);
      } else if (arg->isAction("plugin")) {
        MethodDatasetsProbScan* scannerProb = new MethodDatasetsProbScan((PDF_Datasets*)pdf[0], arg);
        if (!arg->plotpluginonly || (arg->plotpluginonly && !arg->isAction("plot"))) {
          if (Utils::FileExists(m_fnamebuilder->getFileNameScanner(scannerProb))) {
            scannerProb->initScan();
            scannerProb->loadScanner(m_fnamebuilder->getFileNameScanner(scannerProb));
          } else {
            std::cout << "\nWARNING : Couldn't load the Prob scanner, will rerun the Prob" << std::endl;
            std::cout << "          scan now. You should have run the Prob scan locally" << std::endl;
            std::cout << "          before running the Plugin scan." << std::endl;
            std::cout << "          missing file: " << m_fnamebuilder->getFileNameScanner(scannerProb) << std::endl;
            std::cout << std::endl;
            make1dProbScan(scannerProb, 0);
          }
        }
        // create Plugin scanner
        MethodDatasetsPluginScan* scannerPlugin = new MethodDatasetsPluginScan(scannerProb, (PDF_Datasets*)pdf[0], arg);
        if (arg->isAction("plot")) {
          scannerPlugin->loadScanner(m_fnamebuilder->getFileNameScanner(scannerPlugin));
        } else {
          make1dPluginScan(scannerPlugin, 0);
        }
        if (arg->plotpluginonly) {
          make1dPluginOnlyPlot(scannerPlugin, 0);
        } else {
          make1dPluginPlot(scannerPlugin, scannerProb, 0);
        }
      }
    } else if (arg->var.size() == 2) {
      std::cout << "SORRY - 2D plugin scans not yet implemented for datasets - this will probably take a while anyway"
                << std::endl;
      std::exit(1);
    }
  }
  std::cout << "Dataset Scan Done" << std::endl;
}

///
/// run the ROOT application, if the -i flag for interactive
/// mode was set.
///
void GammaComboEngine::runApplication() {
  if (arg->interactive) {
    std::cout << "Exit with Ctrl+c" << std::endl;
    theApp->Run();
  }
}

///
/// print the initial banner
///
void GammaComboEngine::printBanner() {
  std::cout << std::endl
            << "\033[1mGammaCombo " << GAMMACOMBO_VERSION << " \033[0m"
            << "-- All rights reserved under GPLv3, http://www.gnu.org/licenses/gpl.txt" << std::endl
            << std::endl;
}

///
/// run GammaComboEngine, main steering function
///
void GammaComboEngine::run() {
  if (arg->usage) usage();  // print usage and exit
  defineColors();
  checkCombinationArg();
  checkColorArg();
  checkAsimovArg();
  if (arg->scalestaterr > -99) scaleStatErrors();
  if (arg->scaleerr > -99) scaleStatAndSystErrors();
  if (arg->nosyst) disableSystematics();
  makeAddDelCombinations();
  if (arg->nbatchjobs > 0) writebatchscripts();
  customizeCombinerTitles();
  setUpPlot();
  scan();  // most thing gets done here
  if (arg->compare) compareCombinations();
  if (arg->info || arg->latex || (arg->save != "" && !arg->saveAtMin))
    return;  // if only info is requested then we can go home
  if (!arg->isAction("pluginbatch") && !arg->isAction("coveragebatch") && !arg->isAction("coverage")) savePlot();
  std::cout << std::endl;
  t.Stop();
  t.Print();
  runApplication();
}
