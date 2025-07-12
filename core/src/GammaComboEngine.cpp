#include <ColorBuilder.h>
#include <GammaComboEngine.h>
#include <Graphviz.h>
#include <LatexMaker.h>
#include <MethodDatasetsPluginScan.h>
#include <MethodDatasetsProbScan.h>
#include <OneMinusClPlot.h>
#include <OneMinusClPlot2d.h>
#include <PDF_Datasets.h>
#include <ParameterEvolutionPlotter.h>
#include <Utils.h>

#include <RooRealVar.h>

#include <TColor.h>
#include <TDatime.h>
#include <TLatex.h>
#include <TLine.h>
#include <TObjString.h>
#include <TROOT.h>

#include <algorithm>
#include <format>
#include <memory>
#include <vector>

#include <boost/algorithm/string.hpp>
// #include <boost/algorithm/string/predicate.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace RooFit;
using namespace Utils;

GammaComboEngine::GammaComboEngine(TString name, int argc, char* argv[]) {
  // time the program
  t.Start();

  // print the copyright banner
  printBanner();

  // parse the command line options
  arg = make_unique<OptParser>();
  arg->bookAllOptions();
  arg->parseArguments(argc, argv);

  // configure names
  execname = argv[0];
  if (arg->filenameaddition != "") name += "_" + arg->filenameaddition;
  m_fnamebuilder = make_unique<FileNameBuilder>(arg.get(), name);

  // make batch scripts if appropriate and exit
  m_batchscriptwriter = make_unique<BatchScriptWriter>(argc, argv);

  // run ROOT in interactive mode, if requested (-i)
  if (arg->interactive)
    theApp = std::make_unique<TApplication>("App", &argc, argv);
  else
    gROOT->SetBatch(false);

  // reconfigure RooFormulaVar output
  RooMsgService::instance().getStream(1).removeTopic(InputArguments);
  RooMsgService::instance().getStream(0).addTopic(InputArguments);
}

GammaComboEngine::GammaComboEngine(TString name, int argc, char* argv[], bool _runOnDataSet) {
  GammaComboEngine(name, argc, argv);
  runOnDataSet = _runOnDataSet;
}

///
/// Check if a PDF with a certain ID exits.
///
bool GammaComboEngine::pdfExists(int id) const {
  if (id < 0) return false;
  if (id >= this->pdf.size()) return false;
  if (!this->pdf[id]) return false;
  return true;
}

///
/// Check if a Combiner with a certain ID exits.
///
bool GammaComboEngine::combinerExists(int id) const {
  if (id < 0 || id >= this->cmb.size() || !this->cmb[id]) return false;
  return true;
}

///
/// Set the PDF (for datasets method) for the GammaComboEngine
///
void GammaComboEngine::setPdf(PDF_Abs* pdf) {
  if (!runOnDataSet) {
    cout << "It looks like you're trying to set a pdf but you haven't set runOnDataSet=true. I assume this is what you "
            "want so I'm doing it for you"
         << endl;
    runOnDataSet = true;
  }
  if (!dynamic_cast<PDF_Datasets*>(pdf)) {
    cout << "GammaComboEngine::setPdf() : ERROR : The pdf you are trying to set " << pdf->getName()
         << " cannot be cast to a PDF_Datasets object" << endl;
    exit(1);
  }
  if (!pdf) {
    cout << "GammaComboEngine::setPdf() : ERROR : Trying to add zero pointer as the PDF. Exit." << endl;
    exit(1);
  }
  if (pdfExists(0)) {
    cout << "GammaComboEngine::setPdf() : ERROR : You have already set the pdf in GammaComboEngine. Exit." << endl;
    exit(1);
  }
  addPdf(0, pdf);
}

///
/// Add a PDF to the GammaComboEngine object.
///
void GammaComboEngine::addPdf(int id, PDF_Abs* pdf, TString title) {
  if (arg->debug) { cout << "GammaComboEngine::addPdf() : INFO  : Adding pdf " << id << " = " << title << endl; }
  if (!pdf) {
    cout << "GammaComboEngine::addPdf() : ERROR : Trying to add zero pointer as the PDF. Exit." << endl;
    exit(1);
  }
  // check if requested id exists already
  if (pdfExists(id)) {
    cout << "GammaComboEngine::addPdf() : ERROR : Requested PDF id " << id
         << " exists already in GammaComboEngine. Exit." << endl;
    exit(1);
  }
  // check if storage is large enough, enlarge if necessary
  if (id >= this->pdf.size()) {
    for (int i = this->pdf.size(); i <= id; i++) this->pdf.push_back(nullptr);
  }
  this->pdf[id] = pdf;
  if (title != "") this->pdf[id]->setTitle(title);
  this->pdf[id]->setGcId(id);
}

///
/// Add a pdf with a subset of the observables to the GammaComboEngine
///
void GammaComboEngine::addSubsetPdf(const int id, PDF_Abs* pdf, const vector<int>& indices, const TString title) {
  if (indices.size() > pdf->getObservables()->getSize()) {
    cout << "GammaComboEngine::addSubsetPdf() : ERROR - the subset size " << indices.size()
         << " is bigger than the observables size " << pdf->getObservables()->getSize() << endl;
    exit(1);
  }
  for (const auto index : indices) {
    if (index > pdf->getObservables()->getSize() - 1 || index < 0) {
      cout << "GammaComboEngine::addSubsetPdf() : ERROR - one of the subset index values " << index
           << " is larger than the total number of of observables " << pdf->getObservables()->getSize()
           << " or it's less than zero" << endl;
      exit(1);
    }
  }
  RooArgList obsToRemove;
  RooArgList theoryToRemove;

  // loop over all observables and remove the ones that aren't in indices
  for (int i = 0; i < pdf->getObservables()->getSize(); i++) {
    if (std::ranges::find(indices, i) == indices.end()) {
      obsToRemove.add(*(pdf->getObservables()->at(i)));
      theoryToRemove.add(*(pdf->getTheory()->at(i)));
    }
  }
  pdf->getObservables()->remove(obsToRemove);
  pdf->getTheory()->remove(theoryToRemove);

  // now sort out parameters
  RooArgList paramsToRemove;
  for (int i = 0; i < pdf->getParameters()->getSize(); i++) {
    bool paramFoundInTheory = false;
    for (int j = 0; j < pdf->getTheory()->getSize(); j++) {
      if (pdf->getTheory()->at(j)->dependsOn(*(pdf->getParameters()->at(i)))) { paramFoundInTheory = true; }
    }
    if (!paramFoundInTheory) paramsToRemove.add(*(pdf->getParameters()->at(i)));
  }
  pdf->getParameters()->remove(paramsToRemove);

  // now sort out uncertainties
  vector<double> oldStatErrs = pdf->StatErr;
  vector<double> oldSystErrs = pdf->SystErr;
  pdf->StatErr.clear();
  pdf->SystErr.clear();
  for (const auto index : indices) {
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

void GammaComboEngine::addSubsetPdf(int id, PDF_Abs* pdf, int i1, TString title) { addSubsetPdf(id, pdf, {i1}, title); }

void GammaComboEngine::addSubsetPdf(int id, PDF_Abs* pdf, int i1, int i2, TString title) {
  addSubsetPdf(id, pdf, {i1, i2}, title);
}

void GammaComboEngine::addSubsetPdf(int id, PDF_Abs* pdf, int i1, int i2, int i3, TString title) {
  addSubsetPdf(id, pdf, {i1, i2, i3}, title);
}

void GammaComboEngine::addSubsetPdf(int id, PDF_Abs* pdf, int i1, int i2, int i3, int i4, TString title) {
  addSubsetPdf(id, pdf, {i1, i2, i3, i4}, title);
}

void GammaComboEngine::addSubsetPdf(int id, PDF_Abs* pdf, int i1, int i2, int i3, int i4, int i5, TString title) {
  addSubsetPdf(id, pdf, {i1, i2, i3, i4, i5}, title);
}

void GammaComboEngine::addSubsetPdf(int id, PDF_Abs* pdf, int i1, int i2, int i3, int i4, int i5, int i6,
                                    TString title) {
  addSubsetPdf(id, pdf, {i1, i2, i3, i4, i5, i6}, title);
}

void GammaComboEngine::addSubsetPdf(int id, PDF_Abs* pdf, int i1, int i2, int i3, int i4, int i5, int i6, int i7,
                                    TString title) {
  addSubsetPdf(id, pdf, {i1, i2, i3, i4, i5, i6, i7}, title);
}

void GammaComboEngine::addSubsetPdf(int id, PDF_Abs* pdf, int i1, int i2, int i3, int i4, int i5, int i6, int i7,
                                    int i8, TString title) {
  addSubsetPdf(id, pdf, {i1, i2, i3, i4, i5, i6, i7, i8}, title);
}

///
/// Add a Combiner to the GammaComboEngine object.
///
void GammaComboEngine::addCombiner(int id, Combiner* cmb) {
  if (runOnDataSet) {
    cout << "GammaComboEngine::addCombiner() : ERROR : You're trying to make a combiner but the runOnDataSet flag is "
            "true. You cannot make a combination with this option"
         << endl;
    exit(1);
  }
  if (!cmb) {
    cout << "GammaComboEngine::addCombiner() : ERROR : Trying to add nullptr as the Combiner. Exit." << endl;
    exit(1);
  }
  if (combinerExists(id)) {
    cout << "GammaComboEngine::addCombiner() : ERROR : Requested Combiner id exists already in GammaComboEngine. Exit."
         << endl;
    exit(1);
  }
  // check if storage is large enough, enlarge if necessary
  if (id >= this->cmb.size()) {
    for (int i = this->cmb.size(); i <= id; i++) this->cmb.push_back(nullptr);
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
    cout << "GammaComboEngine::cloneCombiner() : ERROR : You're trying to clone a combiner but the runOnDataSet flag "
            "is true. You can't have combiners when using the dataset option."
         << endl;
    exit(1);
  }

  if (combinerExists(newId)) {
    cout << "GammaComboEngine::cloneCombiner() : ERROR : Requested new Combiner id " << newId
         << " exists already in GammaComboEngine. Exit." << endl;
    exit(1);
  }
  if (!combinerExists(oldId)) {
    cout << "GammaComboEngine::cloneCombiner() : ERROR : Requested old Combiner id " << oldId
         << " doesn't exists in GammaComboEngine. Exit." << endl;
    exit(1);
  }
  addCombiner(newId, getCombiner(oldId)->Clone(name, title));
}

///
/// Get a combiner.
/// \param id - combiner ID, set when defining the combiner using addCombiner(), cloneCombiner(), or newCombiner()
///
Combiner* GammaComboEngine::getCombiner(const int id) const {
  if (!combinerExists(id)) {
    cout << "GammaComboEngine::getCombiner() : ERROR : Requested Combiner id " << id
         << " doesn't exist in GammaComboEngine. Exit." << endl;
    exit(1);
  }
  return cmb[id];
}

///
/// Get a vector with the indices of all non-trivial combiners created.
///
vector<int> GammaComboEngine::getCombinersIds() const {
  vector<int> ids;
  for (int id = 0; id < cmb.size(); ++id) {
    if (combinerExists(id)) ids.push_back(id);
  }
  return ids;
}

///
/// Get a PDF.
/// \param id - PDF ID, set when adding the PDF using addPdf()
///
PDF_Abs* GammaComboEngine::getPdf(int id) const {
  if (!pdfExists(id)) {
    cout << "GammaComboEngine::getPdf() : ERROR : Requested PDF id doesn't exist in GammaComboEngine. Exit." << endl;
    exit(1);
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
/// Add a new Combiner, consisting of the PDFs specified by the vector of their GammaComboEngine IDs.
///
/// The PDFs must be added to GammaComboEngine before calling this function.
/// An empty combiner is possible, too, so that it can be filled later.
///
void GammaComboEngine::newCombiner(const int id, const TString name, const TString title, const vector<int>& pdfs) {
  if (combinerExists(id)) {
    cerr << "GammaComboEngine::newCombiner() : ERROR : Requested new Combiner id exists already in GammaComboEngine. "
            "Exit."
         << endl;
    exit(1);
  }
  auto c = new Combiner(arg.get(), name, title);
  for (auto pdf : pdfs) { c->addPdf(getPdf(pdf)); }
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
    double scale = 1.;
    if (i == 25) scale = sqrt(5. * 50. / 3.);
    if (i == 60) scale = sqrt(50. / 3.);
    // if ( i==23 ) scale = sqrt(50.);
    if (scale == 1.) continue;
    cout << "Configuration: Scaling down LHCb errors by a factor " << scale << ": " << pdf[i]->getTitle() << endl;
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

  cout << endl;
}

///
/// scale stat errors
///
void GammaComboEngine::scaleStatErrors() {
  cout << "\nConfiguration: Scaling ALL STAT ERRORS by " << arg->scalestaterr << ".\n" << endl;
  for (auto p : pdf) {
    if (!p) continue;
    for (int iObs = 0; iObs < p->getNobs(); iObs++) { p->StatErr[iObs] *= arg->scalestaterr; }
    p->buildCov();
    p->buildPdf();
  }
}

///
/// scale stat+syst errors
///
void GammaComboEngine::scaleStatAndSystErrors() {
  cout << "\nConfiguration: Scaling ALL STAT AND SYST ERRORS by " << arg->scaleerr << ".\n" << endl;
  for (auto p : pdf) {
    if (!p) continue;
    for (int iObs = 0; iObs < p->getNobs(); iObs++) {
      p->StatErr[iObs] *= arg->scaleerr;
      p->SystErr[iObs] *= arg->scaleerr;
    }
    p->buildCov();
    p->buildPdf();
  }
}

///
/// disable systematics
///
void GammaComboEngine::disableSystematics() {
  cout << "\nConfiguration: Setting ALL SYSTEMATICS TO ZERO.\n" << endl;
  for (auto p : pdf) {
    if (!p) continue;
    for (int iObs = 0; iObs < p->getNobs(); iObs++) p->SystErr[iObs] = 0.;
    p->buildCov();
    p->buildPdf();
  }
}

///
/// Make an Asimov toy: set all observables set to truth values.
/// The Asimov point needs to be loaded in the combiner before.
/// \param c - combiner which should be set to an asimov toy
///
void GammaComboEngine::setAsimovObservables(Combiner* c) {
  if (!c->isCombined()) {
    cout << "GammaComboEngine::setAsimovObservables() : ERROR : Can't set an Asimov toy before "
            "the combiner is combined. Call combine() first."
         << endl;
    exit(1);
  }

  // set observables to asimov values in workspace
  RooWorkspace* w = c->getWorkspace();
  for (const auto pObs : *c->getObservables()) {
    // get theory name from the observable name
    TString pThName = pObs->GetName();
    pThName.ReplaceAll("obs", "th");
    // get the theory relation
    RooAbsReal* th = w->function(pThName);
    if (!th) {
      cout << "GammaComboEngine::setAsimovObservables() : ERROR : theory relation not found in workspace: " << pThName
           << endl;
      exit(1);
    }
    // set the observable to what the theory relation predicts
    static_cast<RooRealVar*>(pObs)->setVal(th->getVal());
  }

  // write back the asimov values to the PDF object so that when
  // the PDF is printed, the asimov values show up
  for (auto pdf : c->getPdfs()) {
    pdf->setObservableSourceString("Asimov");
    for (const auto pObs : *pdf->getObservables()) {
      RooAbsReal* obs = w->var(pObs->GetName());
      if (!obs) {
        cout << "GammaComboEngine::setAsimovObservables() : ERROR : observable not found in workspace: "
             << pObs->GetName() << endl;
        exit(1);
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
  cout << "Start parameter configuration:\n" << endl;
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
    filefound = FileExists(startparfile);
    if (!filefound) {
      cout << "\n ERROR: --parfile not found: " << startparfile << endl;
      cout << "  Will now look for default files.\n" << endl;
    }
  }
  // requested file not found, try default
  if (!filefound) {
    startparfile = startparfile2;
    filefound = FileExists(startparfile);
  }
  // for Asimov combiners, also try non-Asimov start parameters file
  if (!filefound && arg->isAsimovCombiner(cId)) {
    startparfile = startparfile3;
    filefound = FileExists(startparfile);
  }
  // still not found
  if (!filefound) {
    cout << "  No start parameter file was found, will use the default values" << endl;
    cout << "  configured in ParameterAbs class. Start parameter files are searched" << endl;
    cout << "  in the following order:" << endl;
    cout << "  1. --parfile option" << endl;
    cout << "  2. default start paramter file: " << startparfile2 << endl;
    if (arg->isAsimovCombiner(cId)) {
      cout << "  3. for Asimov combiners: corresponding non-Asimov start parameter file: " << startparfile3 << endl;
    }
  } else {
    cout << "  Loading start parameters from file: " << startparfile << endl;
    bool loaded = pCache->loadPoints(startparfile);
    if (!loaded) {
      cout << "  Error loading file. Exit." << endl;
      exit(1);
    }
  }
  cout << endl;
}

///
/// Configure the names and titles of Asimov combiners.
/// It adds things like "Asimov3" to the combiner name,
/// where 3 means the 3rd Asimov point from the configured
/// Asimov parameter file (typically *_genpoints.dat).
///
void GammaComboEngine::configureAsimovCombinerNames(Combiner* c, int i) {
  if (arg->asimov[i] == 0) {
    cout << "\n--asimov 0 : ignoring generator point ID 0" << endl;
    return;
  }
  cout << "\n--asimov : setting up an ASIMOV TOY based on combination \"" << c->getName() << "\"" << endl;
  if (arg->title[i] == TString("default")) c->setTitle(c->getTitle() + " (Asimov)");
  c->setName(c->getName() + m_fnamebuilder->getAsimovCombinerNameAddition(arg->asimov[i]));
  cout << "           Asimov combiner name: \"" << c->getName() << "\"" << endl;
}

///
/// Make an Asimov toy: set all observables set to truth values.
/// The truth values are loaded from a parameter file.
///
void GammaComboEngine::loadAsimovPoint(Combiner* c, int cId) {
  if (arg->asimov[cId] == 0) return;
  cout << "\nAsimov point configuration:\n" << endl;
  auto pCache = new ParameterCache(arg.get());
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
    filefound = FileExists(asimovfile);
  }
  // requested file not found, try default
  if (!filefound) {
    asimovfile = asimovfile2;
    filefound = FileExists(asimovfile);
  }
  // requested file not found, try the start parameter file of the corresponding non-Asimov combiner
  if (!filefound) {
    asimovfile = asimovfile3;
    filefound = FileExists(asimovfile);
  }
  // requested file not found, try the result parameter file of the corresponding non-Asimov combiner
  if (!filefound) {
    asimovfile = asimovfile4;
    filefound = FileExists(asimovfile);
  }
  // if no parameter file exits, we use point from the ParameterAbs class
  if (!filefound) {
    cout << "  No Asimov point parameter file found. Using start values configured in ParameterAbs class." << endl;
    cout << "  Point files are looked for in the following order:" << endl;
    cout << "  1. --asimovfile" << endl;
    cout << "  2. " << asimovfile2 << endl;
    cout << "  3. " << asimovfile3 << endl;
    cout << "  4. " << asimovfile4 << endl;
  } else {
    cout << "  Loading Asimov points from file: " << asimovfile << endl;
    bool loaded = pCache->loadPoints(asimovfile);
    if (!loaded) {
      cout << "  Error loading file. Exit." << endl;
      exit(1);
    }
    cout << "  Setting point number: " << arg->asimov[cId] << endl;
    pCache->setPoint(c, arg->asimov[cId] - 1);
  }
  setAsimovObservables(c);
}

///
/// print usage and exit
///
void GammaComboEngine::usage() const {
  if (runOnDataSet) {
    cout << "USAGE\n\n"
            "  # Compute limit on parameter a:\n"
            "  "
         << execname
         << " --var a --ps 1\n\n"
            "  # Compute limits on two parameters a and b:\n"
            "  "
         << execname << " --var a --var b --ps 1\n\n"
         << endl;
    exit(0);
  }
  cout << "USAGE\n\n"
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
       << endl;
  // cout << "Available plugin toy control plots:\n"
  //"===================================\n"
  //"\n"
  //"All control plots are produced when option --controlplots is given.\n"
  //"If, in addition, -p <n> is given, only one specific plot is produced.\n"
  //" (1) technical overview\n"
  //" (2) summary\n"
  //" (3) nuisances\n"
  //" (4) observables\n"
  //" (5) chi2 distributions\n"
  //" (6) chi2 parabolas\n" << endl;
  print();
  exit(0);
}

///
/// Print the available PDFs.
///
void GammaComboEngine::printPdfs() const {
  cout << "AVAILABLE MEASUREMENTS\n\n";
  for (int i = 0; i < pdf.size(); i++) {
    if (!pdf[i]) continue;
    auto indexStr = format("({:d})", i);
    cout << format("{:>7s} {:s}\n", indexStr, pdf[i]->getTitle().Data());
  }
  cout << endl;
}

///
/// Print the availabe Combinations.
///
void GammaComboEngine::printCombinations() const {
  cout << "AVAILABLE COMBINATIONS\n\n";
  for (int i = 0; i < cmb.size(); i++) {
    if (!cmb[i]) continue;
    auto indexStr = format("({:d})", i);
    cout << format("{:>7s} {:s}\n", indexStr, cmb[i]->getTitle().Data());
  }
  cout << endl;
}

///
/// Print the content of this engine (available PDFs and combinations).
///
void GammaComboEngine::print() const {
  printPdfs();
  printCombinations();
}

///
/// Check the combination argument (-c), exit if it is bad.
///
void GammaComboEngine::checkCombinationArg() const {
  if (runOnDataSet && arg->combid.size() > 0) {
    cout << "When running on a dataset do not pass a combination argument (it makes no sense for this use case)"
         << endl;
    exit(1);
  }
  if (arg->combid.size() == 0 && !runOnDataSet) {
    cout << "Please chose a combination ID (-c).\n" << endl;
    printCombinations();
    exit(1);
  }
  for (auto cid : arg->combid) {
    if (cid >= cmb.size()) {
      cout << "Please chose a combination ID (-c) less than " << cmb.size()
           << ".\nUse the -u option to print a list of available combinations." << endl;
      exit(1);
    }
    if (!cmb[cid]) {
      cout << "You selected an empty combination.\n"
           << "Use the -u option to print a list of available combinations." << endl;
      exit(1);
    }
  }
}

///
/// Check Asimov arg: when only one --asimov argument is given
/// with the ID 0, it won't do anything. Print a warning in that
/// case.
///
void GammaComboEngine::checkAsimovArg() const {
  if (arg->asimov.size() == 1 && arg->asimov[0] == 0) {
    cout << "WARNING : --asimov 0 found, this won't do anything." << endl;
    cout << "          To run an Asimov toy, the generation point ID" << endl;
    cout << "          needs to be different from 0." << endl;
    cout << endl;
  }
}

///
/// Check color argument and exit if a non-existing color was requested through
/// the --color argument. Colors for one-dimensional plots are defined in
/// GammaComboEngine::defineColors(). Colors for two-dimensional pltos are
/// defined in OneMinusClPlot2d::OneMinusClPlot2d().
///
void GammaComboEngine::checkColorArg() const {
  for (int i = 0; i < arg->color.size(); i++) {
    // colors for one-dimensional plots
    if (arg->var.size() == 1) {
      if (colorsLine.size() <= arg->color[i]) {
        cout << "Argument error --color: No such color for one-dimensional plots. Please choose a color between 0 and "
             << colorsLine.size() - 1 << endl;
        exit(1);
      }
    }
    // colors for two-dimensional plots
    else if (arg->var.size() == 2) {
      OneMinusClPlot2d p(arg.get());
      int nMaxColors = p.getNumberOfDefinedColors();
      if (nMaxColors <= arg->color[i]) {
        cout << "Argument error --color: No such color for two-dimensional plots. Please choose a color between 0 and "
             << nMaxColors - 1 << endl;
        exit(1);
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
    cout << "GammaComboEngine::makeAddDelCombinations() : ERROR : internal inconsistency. \n"
            "combid and combmodifications vectors not of same size."
         << endl;
    assert(0);
  }
  // loop over list of modifications
  for (int i = 0; i < arg->combmodifications.size(); i++) {
    // get combiner that is to be modified
    Combiner* cOld = cmb[arg->combid[i]];
    // see if anything is to be added or subtracted
    if (arg->combmodifications[i].size() == 0) continue;
    // there are modifications to be done!
    cout << "-c : Making a new combination based on combination " << arg->combid[i] << endl;
    // compute name and title of new combiner
    TString nameNew = cOld->getName();
    TString titleNew = cOld->getTitle();
    for (auto cm : arg->combmodifications[i]) {
      int pdfId = abs(cm);
      if (!pdfExists(pdfId)) {
        cout << "\nERROR: measurement of given ID does not exist: " << pdfId << endl;
        cout << "       Here is a list of available measurements:" << endl;
        printPdfs();
        exit(1);
      }
      if (cm > 0) {
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
    for (auto cm : arg->combmodifications[i]) {
      int pdfId = abs(cm);
      if (cm > 0) {
        cout << "... adding measurement " << pdfId << endl;
        cNew->addPdf(pdf[pdfId]);
      } else {
        cout << "... deleting measurement " << pdfId << endl;
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
void GammaComboEngine::printCombinerStructure(Combiner* c) const {
  Graphviz gviz(arg.get());
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
    plot = std::make_unique<OneMinusClPlot>(arg.get(), m_fnamebuilder->getFileNamePlot(cmb), "p-value curves");
  } else {
    plot = std::make_unique<OneMinusClPlot2d>(arg.get(), m_fnamebuilder->getFileNamePlot(cmb), "p-value contours");
  }
  plot->disableLegend(arg->plotlegend);
}

///
/// Save the plot to disc.
///
void GammaComboEngine::savePlot() const {
  if (arg->hfagLabel != "")
    HFAGLabel(arg->hfagLabel, arg->plotHFAGLabelPosX, arg->plotHFAGLabelPosY, arg->plotHFAGLabelScale);
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
    cout << "\nPerforming default 2D scan:\n"
            " 1. scan in first variable:  " +
                scanner->getScanVar1Name() +
                "\n"
                " 2. scan in second variable: " +
                scanner->getScanVar2Name() +
                "\n"
                " 3. scan starting from each solution found in 1. and 2."
         << endl;

    // setup a scanner for each variable individually
    Combiner* c = scanner->getCombiner();
    std::unique_ptr<MethodProbScan> s1;
    std::unique_ptr<MethodProbScan> s2;

    cout << "\n1D scan for X variable, " + scanner->getScanVar1Name() + ":\n" << endl;
    if (runOnDataSet) {
      const MethodDatasetsProbScan* temp = dynamic_cast<MethodDatasetsProbScan*>(scanner);
      s1 = std::make_unique<MethodDatasetsProbScan>(temp->pdf, arg.get());
    } else {
      s1 = std::make_unique<MethodProbScan>(c);
    }
    s1->setScanVar1(scanner->getScanVar1Name());
    s1->initScan();
    scanStrategy1d(s1.get(), pCache);
    if (arg->verbose) s1->printLocalMinima();

    cout << "\n1D scan for Y variable, " + scanner->getScanVar2Name() + ":\n" << endl;
    if (runOnDataSet) {
      const MethodDatasetsProbScan* temp = dynamic_cast<MethodDatasetsProbScan*>(scanner);
      s2 = std::make_unique<MethodDatasetsProbScan>(temp->pdf, arg.get());
    } else {
      s2 = std::make_unique<MethodProbScan>(c);
    }
    s2->setScanVar1(scanner->getScanVar2Name());
    s2->setXscanRange(arg->scanrangeyMin, arg->scanrangeyMax);
    s2->initScan();
    scanStrategy1d(s2.get(), pCache);
    if (arg->verbose) s2->printLocalMinima();

    // now do the 2D scan from the two starting points
    cout << "\n2D scan for " + scanner->getScanVar1Name() + " and " + scanner->getScanVar2Name() + ":\n" << endl;
    vector<RooSlimFitResult*> solutions;
    for (int i = 0; i < s1->getSolutions().size(); i++) solutions.push_back(s1->getSolution(i));
    for (int i = 0; i < s2->getSolutions().size(); i++) solutions.push_back(s2->getSolution(i));
    // \todo remove similar solutions from list
    for (int j = 0; j < solutions.size(); j++) {
      cout << "2D scan " << j + 1 << " of " << solutions.size() << " ..." << endl;
      scanner->loadParameters(solutions[j]);
      scanner->scan2d();
    }
  }
  // otherwise load each starting value found
  else {
    cout << "\nPerforming 2D scan from provided starting points." << endl;
    cout << "Number of scans to run: " << nStartingPoints << endl;
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
  auto pCache = new ParameterCache(arg.get());
  loadStartParameters(scanner, pCache, cId);

  scanner->initScan();
  scanStrategy1d(scanner, pCache);
  cout << "\nResults:" << endl;
  cout << "========\n" << endl;
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
    for (const auto cl : arg->cls) {
      scannerPlugin->calcCLintervals(cl);
      if (cl == 2) scannerPlugin->calcCLintervals(cl, true);  // calculate expected upper limit
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
    cout << "making full chi2 plot ..." << endl;
    auto plotf = std::make_unique<OneMinusClPlot2d>(arg.get(), plot->getName() + "_plugin_full",
                                                    "p-value histogram: " + scannerPlugin->getTitle());
    scannerPlugin->plotOn(plotf.get());
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
  auto pCache = new ParameterCache(arg.get());
  if (arg->loadParamsFile.size() != arg->combid.size()) {
    cout << "\nERROR : For a Coverage scan you must pass a parameter file (--parfile) to throw the toys from. You need "
            "one parfile per combiner"
         << endl;
    exit(1);
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
      scanner->plotOn(plot.get(), 1);  // for prob ClsType>1 doesn't exist
    }
    scanner->plotOn(plot.get());
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
  cout << "Scan strategy:" << endl;
  cout << "==============\n" << endl;
  int nStartingPoints = pCache->getNPoints();
  if (nStartingPoints == 0) {
    cout << "1. perform an initial scan" << endl;
    cout << "2. perform an additional scan starting from each solution found\n" << endl;
    cout << "first scan ..." << endl;
    scanner->scan1d();
    if (!arg->probforce) {
      vector<RooSlimFitResult*> firstScanSolutions = scanner->getSolutions();
      for (int i = 0; i < firstScanSolutions.size(); i++) {
        cout << "Scan i: " << i << endl;
        // scanner->loadSolution(i);
        scanner->loadParameters(firstScanSolutions[i]);
        scanner->scan1d(true);
      }
    }
  }
  // otherwise load each starting value found
  else {
    cout << "Scanning from each point found in start parameter file.\n" << endl;
    for (int i = 0; i < nStartingPoints; i++) {
      cout << "scan " << i + 1 << " of " << nStartingPoints << " ..." << endl;
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
    if (arg->cls.size() > 0) sProb->plotOn(plot.get(), 1);
    sProb->plotOn(plot.get());
  } else {
    make1dProbPlot(sProb, cId);
    sPlugin->setLineColor(kBlack);
    sPlugin->setDrawSolution(arg->plotsolutions[cId]);
    if (std::count(arg->cls.begin(), arg->cls.end(), 1)) sPlugin->plotOn(plot.get(), 1);
    if (std::count(arg->cls.begin(), arg->cls.end(), 2)) sPlugin->plotOn(plot.get(), 2);
    sPlugin->plotOn(plot.get());
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
    sPlugin->plotOn(plot.get());
    sProb->plotOn(plot.get());
  } else {
    sProb->plotOn(plot.get());
    sPlugin->plotOn(plot.get());
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
  static_cast<OneMinusClPlot*>(plot.get())->setPluginMarkers(false);
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
  for (const auto cl : arg->cls) sPlugin->plotOn(plot.get(), cl);
  sPlugin->plotOn(plot.get());
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
  sPlugin->plotOn(plot.get());
  plot->Draw();
}

///
/// Make the 1D coverage plot
///
/// \param scanner - the coverage scanner
/// \param cId - the id of this combination on the command line
///
void GammaComboEngine::make1dCoveragePlot(MethodCoverageScan* scanner, [[maybe_unused]] int cId) { scanner->plot(); }

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
  auto pCache = new ParameterCache(arg.get());
  loadStartParameters(scanner, pCache, cId);
  // scan
  scanner->initScan();
  scanStrategy2d(scanner, pCache);
  cout << endl;
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
  std::unique_ptr<OneMinusClPlot2d> plotf;
  if (scanner->getMethodName() == "Prob")
    plotf = std::make_unique<OneMinusClPlot2d>(arg.get(), m_fnamebuilder->getFileNamePlotSingle(cmb, cId) + "_full",
                                               "p-value histogram: " + scanner->getTitle());
  else if (scanner->getMethodName() == "DatasetsProb")
    plotf = std::make_unique<OneMinusClPlot2d>(arg.get(), m_fnamebuilder->getFileNamePlot(cmb) + "_full",
                                               "p-value histogram: " + scanner->getTitle());  // Titus: change to make
                                                                                              // datasets plot possible
  else
    cout << "The name of the scanner matches neither Prob nor DatasetsProb!" << endl;
  scanner->plotOn(plotf.get(), 0);
  plotf->DrawFull();
  plotf->save();
  // plot full CLs
  if (arg->cls.size() > 0) {
    std::unique_ptr<OneMinusClPlot2d> plotfcls;
    if (scanner->getMethodName() == "Prob")
      plotfcls =
          std::make_unique<OneMinusClPlot2d>(arg.get(), m_fnamebuilder->getFileNamePlotSingle(cmb, cId) + "_cls_full",
                                             "p-value histogram: " + scanner->getTitle());
    else if (scanner->getMethodName() == "DatasetsProb")
      plotfcls =
          std::make_unique<OneMinusClPlot2d>(arg.get(), m_fnamebuilder->getFileNamePlot(cmb) + "_cls_full",
                                             "p-value histogram: " + scanner->getTitle());  // Titus: change to make
                                                                                            // datasets plot possible
    else
      cout << "The name of the scanner matches neither Prob nor DatasetsProb!" << endl;
    scanner->plotOn(plotfcls.get(), 1);
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
  if (arg->cls.size() > 0) scanner->plotOn(plot.get(), 1);
  scanner->plotOn(plot.get(), 0);
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
    for (const auto fp : arg->fixParameters[cId]) { c->fixParameter(fp.name, fp.value); }
  }
}

///
/// Helper function for scan(). Adjusts ranges, if requested
/// (only possible before combining).
///
void GammaComboEngine::adjustRanges(Combiner* c, int cId) {
  if (cId < arg->physRanges.size()) {
    for (const auto pr : arg->physRanges[cId]) { c->adjustPhysRange(pr.name, pr.min, pr.max); }
  }
  if (cId < arg->removeRanges.size()) {
    for (const auto rr : arg->removeRanges[cId]) {
      if (rr == "all") {
        const RooArgSet* pars = (RooArgSet*)c->getParameters();
        for (const auto par : *pars) static_cast<RooRealVar*>(par)->removeRange();
      } else {
        c->adjustPhysRange(rr, -999, -999);
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
TString GammaComboEngine::getStartParFileName(int cId) const {
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
bool GammaComboEngine::isScanVarObservable(Combiner* c, TString scanVar) const {
  for (const auto obs : c->getObservableNames()) {
    if (scanVar.Contains(obs)) return true;
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
  cout << "\n--var " << scanVar << ": Setting up a scan for an observable ..." << endl;
  PDF_Abs* pdf = c->getPdfProvidingObservable(scanVar);
  if (!pdf) {
    cout << "GammaComboEngine::tightenChi2Constraint() : ERROR : no PDF found that contains the observable '" << scanVar
         << "'. Exit." << endl;
    exit(1);
  }
  double scale = 0.1;
  cout << "... observable error is multiplied by a factor " << scale << endl;
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

  vector<PDF_Abs*> pdfs = c->getPdfs();

  if (pdfs.size() != arg->readfromfile[cId].size()) {
    cout << "ERROR -- I think you are trying to read values from a file but you haven't specified one for each PDF in "
            "the combination"
         << endl;
    cout << "      -- Pass like this: -c 0:+1,+2,+3 --readfromfile "
            "default,<path_to_pdf2_vals.dat>,<path_to_pdf3_vals.dat>"
         << endl;
    exit(1);
  }

  for (int i = 0; i < pdfs.size(); i++) {

    // read from the file for this pdf
    if (arg->readfromfile[cId][i] == TString("default")) continue;

    //
    ifstream infile(arg->readfromfile[cId][i].Data());
    if (!infile.is_open()) {
      cerr << "No such read file found: " << arg->readfromfile[cId][i] << endl;
      exit(1);
    }
    string line;
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
        vector<string> els;
        boost::split(els, line, boost::is_any_of(" "), boost::token_compress_on);
        TString typ = els[0];
        if (typ == "obs:") {
          TString name = els[1];
          auto val = boost::lexical_cast<double>(els[2]);
          pdfs[i]->setObservable(name, val);
          pdfs[i]->obsValSource = "Read from file " + arg->readfromfile[cId][i];
        } else if (typ == "err:") {
          TString name = els[1];
          auto stat = boost::lexical_cast<double>(els[2]);
          auto syst = boost::lexical_cast<double>(els[3]);
          pdfs[i]->setUncertainty(name, stat, syst);
          pdfs[i]->obsErrSource = "Read from file " + arg->readfromfile[cId][i];
        } else if (typ == "cor:") {
          auto mi = boost::lexical_cast<int>(els[1]);
          auto mj = boost::lexical_cast<int>(els[2]);
          auto corStat = boost::lexical_cast<double>(els[3]);
          auto corSyst = boost::lexical_cast<double>(els[4]);
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
      cout << "WARNING - did not find any pdf named: " << pdfs[i]->getBaseName()
           << " in file: " << arg->readfromfile[cId][i] << endl;
    }
  }
}

///
/// write batch scripts
///
void GammaComboEngine::writebatchscripts() {
  if (runOnDataSet) {
    m_batchscriptwriter->writeScripts_datasets(arg.get(), getPdf(0));
  } else {
    m_batchscriptwriter->writeScripts(arg.get(), &cmb);
  }
  exit(0);
}

///
/// make latex
///
void GammaComboEngine::makeLatex(Combiner* c) const {
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
  if (arg->enforcePhysRange) setLimit(c->getParameters(), "phys");

  // first write a copy THE PDF into the workspace
  const RooAbsPdf* thePdf = c->getPdf();
  const RooArgSet* theParameters = c->getParameters();
  const RooArgSet* theObservables = c->getObservables();

  auto savePdf = (RooAbsPdf*)thePdf->Clone("ThePdf");
  c->getWorkspace()->import(*savePdf, Silence());
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
      auto pull_corr =
          std::make_unique<TH2F>(Form("c%d_vs_c%d_corr", i, j),
                                 Form("; Obs pulls for %s [#sigma]; Obs pulls %s [#sigma]",
                                      comparisonScanners[i]->getName().Data(), comparisonScanners[j]->getName().Data()),
                                 10, -5, 5, 10, -5, 5);
      vector<double> pullVec1;
      vector<double> pullVec2;
      double total_pull = 0.;
      int nmatch = 0;
      comparisonScanners[i]->loadSolution(0);
      comparisonScanners[j]->loadSolution(0);
      const RooArgSet* sc1obs = comparisonScanners[i]->getObservables();
      const RooArgSet* sc2obs = comparisonScanners[j]->getObservables();
      for (const auto pAbsObs1 : *sc1obs) {
        for (const auto pAbsObs2 : *sc2obs) {
          const auto pObs1 = static_cast<RooRealVar*>(pAbsObs1);
          const auto pObs2 = static_cast<RooRealVar*>(pAbsObs2);

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
            auto pTh1 = (RooRealVar*)comparisonScanners[i]->getTheory()->find(pTh1Name);
            auto pTh2 = (RooRealVar*)comparisonScanners[j]->getTheory()->find(pTh2Name);
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
        err = sqrt(sq(cl1.central - cl1.min) + sq(cl2.max - cl2.central) -
                   2. * corr * (cl1.central - cl1.min) * (cl2.max - cl2.central));
      } else {
        err = sqrt(sq(cl1.max - cl1.central) + sq(cl2.central - cl2.min) -
                   2. * corr * (cl1.max - cl1.central) * (cl2.central - cl2.min));
      }
      cout << "Comparison   1): " << Form("%-20s", comparisonScanners[i]->getName().Data())
           << " to 2): " << Form("%-20s", comparisonScanners[j]->getName().Data()) << endl;
      cout << "        chi2:    " << Form("%-20.3f", chi21) << "        " << Form("%-20.3f", chi22) << endl;
      cout << "        nObs:    " << Form("%-20d", nObs1) << "        " << Form("%-20d", nObs2) << endl;
      cout << "        nPar:    " << Form("%-20d", nPar1) << "        " << Form("%-20d", nPar2) << endl;
      cout << "        pval:    " << Form("%-20.2f", 100. * TMath::Prob(chi21, nObs1 - nPar1)) << "        "
           << Form("%-20.2f", 100. * TMath::Prob(chi22, nObs2 - nPar2)) << endl;
      cout << "         val:    " << Form("%-6.3f", cl1.central) << " [" << Form("%6.3f", cl1.min) << ","
           << Form("%-6.3f", cl1.max) << "]"
           << "      " << Form("%-6.3f", cl2.central) << "[" << Form("%6.3f", cl2.min) << "," << Form("%-6.3f", cl2.max)
           << "]" << endl;
      cout << "        dval:    " << diff << " +/- " << err << " (" << TMath::Abs(diff) / err << " sigma)" << endl;
      cout << "PULL PER OBS:  " << total_pull << endl;
      cout << "CORRELATION:   " << corr << endl;
      cout << "COMPATIBILITY: " << TMath::Abs(diff) / err << " sigma" << endl;

      auto canv = newNoWarnTCanvas("pull_corr" + getUniqueRootName());
      pull_corr->SetMarkerStyle(kMultiply);
      pull_corr->SetMarkerColor(kBlue + 2);
      pull_corr->GetXaxis()->SetTitleSize(0.045);
      pull_corr->GetYaxis()->SetTitleSize(0.045);
      pull_corr->GetXaxis()->SetLabelSize(0.045);
      pull_corr->GetYaxis()->SetLabelSize(0.045);
      pull_corr->Draw("scat");
      TLine line;
      line.DrawLine(-5, 0, 5, 0);
      line.DrawLine(0, -5, 0, 5);
      pull_corr->Draw("scatsame");
      TF1 f1("f1", "[0]*x", -5, 5);
      f1.SetParameter(0, corr);
      f1.Draw("Lsame");
      TLatex lat;
      lat.DrawLatex(3, 4, Form("#rho = %3.1f", corr));
      lat.DrawLatex(3, 3, Form("#sigma = %3.1f", TMath::Abs(diff) / err));
      Utils::savePlot(canv.get(), Form("pull_corr_%s_%s", comparisonScanners[i]->getName().Data(),
                                       comparisonScanners[j]->getName().Data()));
      // total_pull = sqrt(total_pull) / nmatch;
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
  auto probscan = std::make_unique<MethodProbScan>(c);
  make1dProbScan(probscan.get(), 0);

  TString toydirname = TString("root/scan1dToys_") + probscan->getName() + TString("_") + probscan->getScanVar1Name();
  TString toyfname = toydirname + TString("/scan1dToys_") + probscan->getName() + TString("_") +
                     probscan->getScanVar1Name() + Form("_run%d.root", arg->nrun);

  auto tree = std::make_unique<TTree>("toys", "toys");
  map<TString, double> vals;
  map<TString, double> errs;
  double chi2val = probscan->solutions[0]->minNll();
  int ntoy = -1;
  tree->Branch("ntoy", &ntoy);
  tree->Branch("chi2min", &chi2val);
  RooArgList fitPars = probscan->solutions[0]->floatParsFinal();
  for (const auto pAbs : fitPars) {
    const auto p = static_cast<RooRealVar*>(pAbs);
    cout << "YO:: " << p->GetName() << endl;
    vals[p->GetName()] = p->getVal();
    errs[p->GetName()] = p->getError();
    tree->Branch(p->GetName() + TString("_val"), &vals[p->GetName()]);
    tree->Branch(p->GetName() + TString("_err"), &errs[p->GetName()]);
  }
  tree->Fill();

  for (int i = 0; i < arg->ntoys; i++) {
    cout << "RUNNING TOY " << i << " / " << arg->ntoys << endl;
    c->setObservablesToToyValues();
    auto toyscan = std::make_unique<MethodProbScan>(c);
    make1dProbScan(toyscan.get(), 0);
    if (toyscan->solutions.size() == 0) continue;
    toyscan->solutions[0]->Print();
    // cout << "My stuff I want to save" << endl;
    // cout << "chi2: " << toyscan->solutions[0]->minNll() << endl;

    chi2val = toyscan->solutions[0]->minNll();
    ntoy = i;
    RooArgList toyFitPars = toyscan->solutions[0]->floatParsFinal();
    for (const auto pAbs : toyFitPars) {
      const auto p = static_cast<RooRealVar*>(pAbs);
      vals[p->GetName()] = p->getVal();
      errs[p->GetName()] = p->getError();
    }
    tree->Fill();
  }
  system("mkdir -p " + toydirname);
  TFile f(toyfname, "recreate");
  tree->Write();
  f.Close();
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
      auto scannerProb = std::make_unique<MethodProbScan>(c);
      // pvalue corrector
      if (arg->coverageCorrectionID > 0) {
        auto pvalueCorrector = new PValueCorrection(arg->coverageCorrectionID, arg->verbose);
        pvalueCorrector->readFiles(m_fnamebuilder->getFileBaseName(c), arg->coverageCorrectionPoint,
                                   false);  // false means for prob
        pvalueCorrector->write("root/pvalueCorrection_prob.root");
        scannerProb->setPValueCorrector(pvalueCorrector);
      }

      // 1D SCANS
      if (arg->var.size() == 1) {
        if (arg->isAction("plot")) {
          scannerProb->loadScanner(m_fnamebuilder->getFileNameScanner(scannerProb.get()));
        } else {
          make1dProbScan(scannerProb.get(), i);
        }
        make1dProbPlot(scannerProb.get(), i);
        if (arg->compare) comparisonScanners.push_back(scannerProb.get());
      }
      // 2D SCANS
      else if (arg->var.size() == 2) {
        if (arg->isAction("plot")) {
          scannerProb->loadScanner(m_fnamebuilder->getFileNameScanner(scannerProb.get()));
        } else {
          make2dProbScan(scannerProb.get(), i);
        }
        make2dProbPlot(scannerProb.get(), i);
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
          auto scannerProb = std::make_unique<MethodProbScan>(c);
          make1dProbScan(scannerProb.get(), i);
          auto scannerPlugin = std::make_unique<MethodPluginScan>(scannerProb.get());
          make1dPluginScan(scannerPlugin.get(), i);
        }
        // if ( arg->isAction("pluginhybridbatch") ){
        //// Hybrid Plugin: compute a second profile likelihood to define the parameter evolution
        // cout << "HYBRID PLUGIN: preparing profile likelihood to be used for parameter evolution:" << endl;
        // ParameterCache *pCache = new ParameterCache(arg.get(), m_fnamebuilder->getFileBaseName(cmb[arg->pevid[0]]));
        // pCache->loadPoints();
        // auto scanner3 = new MethodProbScan(cmb[arg->pevid[0]]);
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
          auto scannerProb = std::make_unique<MethodProbScan>(c);
          if (!arg->plotpluginonly || (arg->plotpluginonly && !arg->isAction("plot"))) {
            if (FileExists(m_fnamebuilder->getFileNameScanner(scannerProb.get()))) {
              scannerProb->loadScanner(m_fnamebuilder->getFileNameScanner(scannerProb.get()));
            } else {
              cout << "\nWARNING : Couldn't load the Prob scanner, will rerun the Prob" << endl;
              cout << "          scan now. You should have run the Prob scan locally" << endl;
              cout << "          before running the Plugin scan." << endl;
              cout << "          missing file: " << m_fnamebuilder->getFileNameScanner(scannerProb.get()) << endl;
              cout << endl;
              make1dProbScan(scannerProb.get(), i);
            }
          }
          // create Plugin scanner
          auto scannerPlugin = std::make_unique<MethodPluginScan>(scannerProb.get());
          if (arg->isAction("plot")) {
            scannerPlugin->loadScanner(m_fnamebuilder->getFileNameScanner(scannerPlugin.get()));
          } else {
            if (arg->coverageCorrectionID > 0) {
              auto pvalueCorrector = new PValueCorrection(arg->coverageCorrectionID, arg->verbose);
              pvalueCorrector->readFiles(m_fnamebuilder->getFileBaseName(c), arg->coverageCorrectionPoint,
                                         true);  // true means for plugin
              pvalueCorrector->write("root/pvalueCorrection_plugin.root");
              scannerPlugin->setPValueCorrector(pvalueCorrector);
            }
            make1dPluginScan(scannerPlugin.get(), i);
          }
          if (arg->plotpluginonly) {
            make1dPluginOnlyPlot(scannerPlugin.get(), i);
          } else {
            make1dPluginPlot(scannerPlugin.get(), scannerProb.get(), i);
          }
        }

      }
      // 2D SCANS
      else if (arg->var.size() == 2) {
        if (arg->isAction("pluginbatch")) {
          auto scannerProb = std::make_unique<MethodProbScan>(c);
          make2dProbScan(scannerProb.get(), i);
          auto scannerPlugin = std::make_unique<MethodPluginScan>(scannerProb.get());
          make2dPluginScan(scannerPlugin.get(), i);
        } else if (arg->isAction("plugin")) {
          auto scannerProb = std::make_unique<MethodProbScan>(c);
          if (!(arg->isAction("plot") && arg->plotpluginonly)) {
            // we don't need the prob scanner if we just want to replot the plugin only
            scannerProb->loadScanner(m_fnamebuilder->getFileNameScanner(scannerProb.get()));
          }
          auto scannerPlugin = std::make_unique<MethodPluginScan>(scannerProb.get());
          if (arg->isAction("plot")) {
            scannerPlugin->loadScanner(m_fnamebuilder->getFileNameScanner(scannerPlugin.get()));
          } else {
            make2dPluginScan(scannerPlugin.get(), i);
          }
          if (arg->plotpluginonly) {
            make2dPluginOnlyPlot(scannerPlugin.get(), i);
          } else {
            make2dPluginPlot(scannerPlugin.get(), scannerProb.get(), i);
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
        cerr << "ERROR -- you can only scan in 1D for a coverage check" << endl;
        exit(1);
      }
      auto coverageScan = std::make_unique<MethodCoverageScan>(c);
      if (arg->isAction("coveragebatch")) {
        make1dCoverageScan(coverageScan.get(), i);
      } else if (arg->isAction("coverage")) {
        if (arg->isAction("plot")) {
          if (FileExists(m_fnamebuilder->getFileNameScanner(coverageScan.get()))) {
            coverageScan->loadScanner(m_fnamebuilder->getFileNameScanner(coverageScan.get()));
          } else {
            cout << "\nERROR : Couldn't load the coverage scanner: "
                 << m_fnamebuilder->getFileNameScanner(coverageScan.get()) << endl;
            exit(1);
          }
        } else {
          make1dCoverageScan(coverageScan.get(), i);
        }
        make1dCoveragePlot(coverageScan.get(), i);
      }
    }

    // RUN TOYS
    if (arg->isAction("runtoys")) runToys(c);
    /////////////////////////////////////////////////////

    // SAVE WORKSPACE
    if (arg->save != "" && arg->saveAtMin) saveWorkspace(c, i);
    /////////////////////////////////////////////////////

    if (i < arg->combid.size() - 1) {
      cout << "\n-- now starting -c " << arg->combid[i + 1]
           << " ------------------------------------------------------------------\n"
           << endl;
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
    auto probScanner = std::make_unique<MethodDatasetsProbScan>((PDF_Datasets*)pdf[0], arg.get());

    // 1D SCANS
    if (arg->var.size() == 1) {
      if (arg->isAction("plot")) {
        probScanner->loadScanner(m_fnamebuilder->getFileNameScanner(probScanner.get()));
      } else {
        make1dProbScan(probScanner.get(), 0);
      }
      make1dProbPlot(probScanner.get(), 0);
    }
    // 2D SCANS
    else if (arg->var.size() == 2) {
      if (arg->isAction("plot")) {
        probScanner->loadScanner(m_fnamebuilder->getFileNameScanner(probScanner.get()));
      } else {
        make2dProbScan(probScanner.get(), 0);
      }
      make2dProbPlot(probScanner.get(), 0);
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
        auto scannerProb = std::make_unique<MethodDatasetsProbScan>((PDF_Datasets*)pdf[0], arg.get());
        if (FileExists(m_fnamebuilder->getFileNameScanner(scannerProb.get()))) {
          scannerProb->initScan();
          scannerProb->loadScanner(m_fnamebuilder->getFileNameScanner(scannerProb.get()));
        } else {
          cout << "\nWARNING : Couldn't load the Prob scanner, will rerun the Prob" << endl;
          cout << "          scan now. You should have run the Prob scan locally" << endl;
          cout << "          before running the Plugin scan." << endl;
          cout << "          missing file: " << m_fnamebuilder->getFileNameScanner(scannerProb.get()) << endl;
          cout << endl;
          make1dProbScan(scannerProb.get(), 0);
        }
        auto scannerPlugin =
            std::make_unique<MethodDatasetsPluginScan>(scannerProb.get(), (PDF_Datasets*)pdf[0], arg.get());
        make1dPluginScan(scannerPlugin.get(), 0);
      } else if (arg->isAction("plugin")) {
        auto scannerProb = std::make_unique<MethodDatasetsProbScan>((PDF_Datasets*)pdf[0], arg.get());
        if (!arg->plotpluginonly || (arg->plotpluginonly && !arg->isAction("plot"))) {
          if (FileExists(m_fnamebuilder->getFileNameScanner(scannerProb.get()))) {
            scannerProb->initScan();
            scannerProb->loadScanner(m_fnamebuilder->getFileNameScanner(scannerProb.get()));
          } else {
            cout << "\nWARNING : Couldn't load the Prob scanner, will rerun the Prob" << endl;
            cout << "          scan now. You should have run the Prob scan locally" << endl;
            cout << "          before running the Plugin scan." << endl;
            cout << "          missing file: " << m_fnamebuilder->getFileNameScanner(scannerProb.get()) << endl;
            cout << endl;
            make1dProbScan(scannerProb.get(), 0);
          }
        }
        // create Plugin scanner
        auto scannerPlugin =
            std::make_unique<MethodDatasetsPluginScan>(scannerProb.get(), (PDF_Datasets*)pdf[0], arg.get());
        if (arg->isAction("plot")) {
          scannerPlugin->loadScanner(m_fnamebuilder->getFileNameScanner(scannerPlugin.get()));
        } else {
          make1dPluginScan(scannerPlugin.get(), 0);
        }
        if (arg->plotpluginonly) {
          make1dPluginOnlyPlot(scannerPlugin.get(), 0);
        } else {
          make1dPluginPlot(scannerPlugin.get(), scannerProb.get(), 0);
        }
      }
    } else if (arg->var.size() == 2) {
      cout << "SORRY - 2D plugin scans not yet implemented for datasets - this will probably take a while anyway"
           << endl;
      exit(1);
    }
  }
  cout << "Dataset Scan Done" << endl;
}

///
/// run the ROOT application, if the -i flag for interactive
/// mode was set.
///
void GammaComboEngine::runApplication() {
  if (arg->interactive) {
    cout << "Exit with Ctrl+c" << endl;
    theApp->Run();
  }
}

///
/// print the initial banner
///
void GammaComboEngine::printBanner() const {
  const char* VTAG = "1.3";
  cout << endl
       << "\033[1mGammaCombo v" << VTAG << " \033[0m"
       << "-- All rights reserved under GPLv3, http://www.gnu.org/licenses/gpl.txt" << endl
       << endl;
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
  cout << endl;
  t.Stop();
  t.Print();
  runApplication();
}
