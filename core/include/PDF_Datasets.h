/**
 * Gamma Combination
 * Author: Maximilian Schlupp, maxschlupp@gmail.com
 * Author: Konstantin Schubert, schubert.konstantin@gmail.com
 * Date: October 2016
 *
 * Abstract class for handling generic PDFs
 * The RooArgLists for Observables and Parameters must be
 * provided by an external workspace. Either there are
 * named sets within the workspace or define the variables
 * for each set manually.
 * The user should inherit from this class.
 * The user has to implement the PDFs methods: fit(), generateToys()
 *
 **/

#ifndef PDF_Datasets_h
#define PDF_Datasets_h

#include "OptParser.h"
#include "PDF_Abs.h"

class PDF_Datasets : public PDF_Abs
{
public:
    PDF_Datasets(RooWorkspace* w, int nObs, const OptParser* opt);
    PDF_Datasets(RooWorkspace* w);
    virtual ~PDF_Datasets();
    void                  deleteNLL() {if (_NLL) {delete _NLL; _NLL = nullptr;}};

    virtual RooFitResult* fit(RooAbsData* dataToFit);
    virtual RooFitResult* fitBkg(RooAbsData* dataToFit, TString signalvar);
    virtual void          generateToys(int SeedShift = 0);
    virtual void          generateToysGlobalObservables(int SeedShift = 0);
    virtual void          generateBkgToys(int SeedShift = 0, TString signalvar="");
    virtual void          generateBkgToysGlobalObservables(int SeedShift = 0, int index = 0);

    void                  initConstraints(const TString& setName);
    void                  initData(const TString& name);
    void                  initObservables(const TString& setName);
    void                  initObservables() override;
    void                  initGlobalObservables(const TString& setName);
    void                  initParameters(const TString& setName);
    void                  initParameters(const std::vector<TString>& parNames);
    void                  initParameters() override;
    void                  initMultipdfCat(const TString& name);
    void                  initPDF(const TString& name);
    void                  initBkgPDF(const TString& name);
    inline  void          addFitObs(TString name) {fitObs.push_back(name);};

    const OptParser*            getArg() const;
    TString               getConstraintName() const {return constraintName;};
    TString               getDataName() const {return dataName;};
    RooAbsData*           getData() {return this->data;};
    inline int            getFitStatus() const {return fitStatus;};
    inline int            getFitStrategy() const {return fitStrategy;};
    inline std::vector<TString>  getFitObs() const {return fitObs;};
    inline std::map<TString,TString> getUnblindRegions() const { return unblindRegions;};
    TString               getGlobalObsName() const {return globalObsName;};
    float                 getMinNll() const {return minNll;};
    float                 getMinNllFree() const {return minNllFree;};
    float                 getMinNllBkg() const {return minNllBkg;};
    float                 getMinNllScan() const {return minNllScan;};
    inline int            getBestIndex() const {return bestIndex;};
    inline int            getBestIndexBkg() const {return bestIndexBkg;};
    inline int            getBestIndexScan() const {return bestIndexScan;};
    TString               getObsName() const {return obsName;};
    TString               getParName() const {return parName;};
    TString               getPdfName() const {return pdfName;};
    TString               getBkgPdfName() const {return pdfBkgName;};
    RooAbsData*           getToyObservables() {return this->toyObservables;};
    RooAbsData*           getBkgToyObservables() {return toyBkgObservables;};
    TString               getMultipdfCatName() const {return multipdfCatName;};
    RooWorkspace*         getWorkspace() {return wspc;};
    // setters
    inline void           setFitStatus(int stat = 0) {fitStatus = stat;};
    inline void           setFitStrategy(int strat = 0) {fitStrategy = strat;};
    inline void           setMinNll(float mnll) {minNll = mnll;};
    inline void           setMinNllFree(float mnll) {minNllFree = mnll;};
    inline void           setMinNllScan(float mnll) {minNllScan = mnll;};
    inline void           setBestIndex(int index) {bestIndex = index;};
    inline void           setBestIndexBkg(int index) {bestIndexBkg = index;};
    inline void           setBestIndexScan(int index) {bestIndexScan = index;};
    void                  setNCPU(int n) {NCPU = n;};
    void                  setVarRange(const TString &varName, const TString &rangeName,
                                      const double &rangeMin, const double &rangeMax);
    void                  setToyData(RooAbsData* ds);
    void                  setBkgToyData(RooAbsData* ds);
    
    void                  setGlobalObsSnapshotBkgToy(TString snapshotname) {globalObsBkgToySnapshotName = snapshotname;};

    void                  unblind(TString var, TString unblindRegs);
    void                  print() const;
    void                  printParameters() const;
    inline  bool          areObservglobalablesSet() const { return areObsSet; };
    inline  bool          areParametersSet() const { return areParsSet; };
    inline  bool          isPdfInitialized() const { return isPdfSet; };
    inline  bool          isMultipdfInitialized() const { return isMultipdfSet; };
    inline  bool          isDataInitialized() const { return isDataSet; };
    inline  bool          isMultipdfCatInitialized() const {return isMultipdfCatSet; };
    inline  bool          notSetupToFit(bool fitToys) const {return (!(isPdfSet && isDataSet) || (fitToys && !(isPdfSet && isToyDataSet))); }; // this comes from a previous if-statement


    int                   NCPU;         //> number of CPU used
    float                 minNll = 0.;

    const TString         globalObsDataSnapshotName = "globalObsDataSnapshotName";
    //> name of a snapshot that stores the values of the global observables in data
    const TString         globalObsToySnapshotName = "globalObsToySnapshotName";
    //> name of a snapshot that stores the latest simulated values for the global observables
   TString         globalObsBkgToySnapshotName = "globalObsBkgToySnapshotName";
    //> name of a snapshot that stores the latest simulated values for the global observables of the bkg-only toy

    //debug counters
    int nbkgfits = 0;
    int nsbfits = 0;

protected:
    void initializeRandomGenerator(int seedShift);
    RooWorkspace* wspc = nullptr;
    RooAbsData*   data = nullptr;
    RooAbsReal*   _NLL = nullptr; // possible pointer to minimization function
    RooAbsPdf*    _constraintPdf = nullptr;
    TString       pdfName = "default_pdf_workspace_name"; //> name of the pdf in the workspace
    TString       pdfBkgName = "default_pdf_bkg_workspace_name"; //> name of the bkg pdf in the workspace
    TString       multipdfCatName; //> name of the multipdf category
    TString       obsName = "default_internal_observables_set_name";
    TString       parName = "default_internal_parameter_set_name";
    TString       dataName = "default_internal_dataset_name";       //> name of the data set in the workspace
    TString       constraintName = "default_internal_constraint_set_name"; //> name of the set with all constraint pdfs
    TString       globalParsName = "default_internal_global_pars_set_name"; //> name of the set of global parameters in the workspace, that is, the parameters that occur (not only) in the  constraints...
    TString       globalObsName = "default_internal_global_obs_set_name";   //> name of the set of global observables in the workspace.
    const OptParser*    arg = nullptr;
    int           fitStrategy = 0;
    int           fitStatus = -10;
    float         minNllFree = 0.;
    float         minNllBkg;
    float         minNllScan = 0.;
    int           bestIndex;
    int           bestIndexBkg;
    int           bestIndexScan;
    bool areObsSet = false;       //> Forces user to set observables
    bool areParsSet = false;      //> Forces user to set parameters
    bool areRangesSet = false;    //> Flag deciding if necessary ranges are set
    bool isPdfSet = false;        //> Flag deciding if PDF is set
    bool isBkgPdfSet = false;     //> Flag deciding if Bkg PDF is set
    bool isMultipdfSet = false;   //> Flag deciding if multipdf is set
    bool isBkgMultipdfSet = false;//> Flag deciding if Bkg multipdf is set
    bool isMultipdfCatSet = false;//> Flag deciding if multipdf category is set
    bool isDataSet = false;       //> Flag deciding if Data is set
    bool isToyDataSet = false;    //> Flag deciding if ToyData is set

    std::vector<TString>    fitObs;
    std::map<TString,TString>   unblindRegions;
};

#endif
