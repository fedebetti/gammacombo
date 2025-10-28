#include <RooAbsArg.h>
#include <RooAddPdf.h>
#include <RooArgList.h>
#include <RooCustomizer.h>
#include <RooDataSet.h>
#include <RooDstD0BG.h>
#include <RooExponential.h>
#include <RooExtendPdf.h>
#include <RooFit.h>
#include <RooFitResult.h>
#include <RooFormulaVar.h>
#include <RooGaussian.h>
#include <RooPlot.h>
#include <RooRealVar.h>
#include <RooWorkspace.h>

#include <TCanvas.h>
#include <TFile.h>
#include <TString.h>

#include <typeinfo>

int main() {
  /////////////////////////////////////////////////////////
  //
  // First, we define all observables and variables and construct the PDF.
  // It will be included in the workspace.
  //
  /////////////////////////////////////////////////////////

  double mL = 4600.0;
  double mR = 10000.0;
  RooRealVar Lb_M_reco("Lb_M_reco", "#it{m}(#it{#Lambda_{b}^{0}}) [MeV / #it{c}^{2}]", mL, mR);

  // First, define the signal peak of the mass model

  // TODO: import signal and bkg pdfs, and bkg yield, from workspaces
  // TODO: make workspace path and name configurable
  TString path_input_workspaces("~/work/Lb_to_taus/vrd-lb2pktaumu/hadronic/fit/");  // path to the input workspaces
  auto file_w_bkg = TFile::Open(path_input_workspaces + "w_bkg_lowBDTCut.root");
  auto w_bkg = (RooWorkspace*)file_w_bkg->Get("w");
  std::cout << "INPUT WORKSPACE FOR BKG:\n";
  w_bkg->Print();
  auto file_w_sig = TFile::Open(path_input_workspaces + "w_sig_lowBDTCut.root");
  auto w_sig = (RooWorkspace*)file_w_sig->Get("w");
  std::cout << "INPUT WORKSPACE FOR SIG:\n";
  w_sig->Print();

  // some gymnastic to map the mass observable to those of sig and bkg functions
  auto m = w_sig->var("m");

  auto sig_tmp = w_sig->pdf("sig");
  auto m_tmp = sig_tmp->getObservables(RooArgSet(*m))->find("m");
  RooCustomizer customizer_sig(*sig_tmp, "sig");
  customizer_sig.replaceArg(*m_tmp, Lb_M_reco);
  RooAbsArg* sig_rebuilt = customizer_sig.build(true);  // true = recycle unchanged nodes
  TObject* clone_obj = sig_rebuilt->clone(Form("%s_clone", sig_rebuilt->GetName()));
  auto* cloned_pdf = dynamic_cast<RooAbsPdf*>(clone_obj);
  std::unique_ptr<RooAbsPdf> signal_model{cloned_pdf};
  signal_model->SetNameTitle("signal_model", "signal_model");

  auto bkg_tmp = w_bkg->pdf("bkg");
  m_tmp = bkg_tmp->getObservables(RooArgSet(*m))->find("m");
  RooCustomizer customizer_bkg(*bkg_tmp, "bkg");
  customizer_bkg.replaceArg(*m_tmp, Lb_M_reco);
  RooAbsArg* bkg_rebuilt = customizer_bkg.build(true);  // true = recycle unchanged nodes
  TObject* clone_obj_bkg = bkg_rebuilt->clone(Form("%s_clone", bkg_rebuilt->GetName()));
  auto* cloned_pdf_bkg = dynamic_cast<RooAbsPdf*>(clone_obj_bkg);
  std::unique_ptr<RooAbsPdf> bkg_only_model{cloned_pdf_bkg};
  bkg_only_model->SetNameTitle("bkg_only_model", "bkg_only_model");

  auto params_signal = signal_model->getParameters(RooArgSet(Lb_M_reco));
  for (auto& param : *params_signal) {
    if (auto* real = dynamic_cast<RooRealVar*>(param)) { real->setConstant(true); }
  }

  RooRealVar n_bkg("Nbkg", "Nbkg", 4981., 0., 10000.);
  RooExtendPdf extended_bkg_model("extended_bkg_model", "extended_bkg_model", *bkg_only_model, n_bkg);

  // The number of signal events is related to the branching ratio via the formula <branching ratio> = <n_sig> *
  // <normalization factor> The normalization factor is not exactly known. Instead, it has to be estimated. The
  // estimator for the normalization factor is a global observable that constrains its value via a Gaussian Constraint.
  double observedValueGlobalObservable =
      2.6e-7;  // estimated value of the normalization constant, taken from slide 30 of //
               // https://indico.cern.ch/event/1527783/contributions/6427746/attachments/3120608/5533689/WG_August20th_VRD_Lbpktaulepton.pdf
  double sigma_observedValueGlobalObservable =
      0.1 * observedValueGlobalObservable;  // assuming 10% relative uncertainty
  RooRealVar norm_constant_obs(
      "norm_constant_glob_obs", "global observable of normalization constant", observedValueGlobalObservable,
      observedValueGlobalObservable - 7 * sigma_observedValueGlobalObservable,
      observedValueGlobalObservable + 7 * sigma_observedValueGlobalObservable);  // this is the observed value, the
                                                                                 // global observable
  RooRealVar norm_constant(
      "norm_constant", "norm_constant", observedValueGlobalObservable - 7 * sigma_observedValueGlobalObservable,
      observedValueGlobalObservable + 7 * sigma_observedValueGlobalObservable);  // the normalization constant is a
                                                                                 // nuisance parameter.
  RooRealVar norm_constant_sigma("norm_constant_sigma", "norm_constant_sigma", sigma_observedValueGlobalObservable);
  RooGaussian norm_constant_constraint("norm_constant_constraint", "norm_constant_constraint", norm_constant_obs,
                                       norm_constant, norm_constant_sigma);

  // Now we can build the mass model by adding the signal and background probability density functions
  RooRealVar branchingRatio("branchingRatio", "branchingRatio", 1e-9, 0,
                            1e-5);  // this is the branching ratio, the parameter of interest
  RooFormulaVar n_sig("Nsig", "branchingRatio/norm_constant", RooArgList(branchingRatio, norm_constant));
  RooExtendPdf extended_sig_model("extended_sig_model", "extended_sig_model", *signal_model, n_sig);

  RooAddPdf mass_model("mass_model", "mass_model", RooArgList(*signal_model, *bkg_only_model),
                       RooArgList(n_sig, n_bkg));

  /////////////////////////////////////////////////////////
  //
  // Open/create the dataset with the data that will be analyzed
  // and set the global observables to the observed values
  //
  /////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////
  //
  // We create the dataset from real data (Wrong Sign for the moment)
  // TODO: set configurable Cut
  TString data_file_name =
      "/eos/lhcb/user/f/fbetti/Lb_to_taus/data/Run2/Lb_pKtaue_3pi_SS_PplusEplus_massConstr_sortPi_vetoes_BDT.root";
  RooRealVar TMVAClassification_BDT_all_noIPCHI2rew("TMVAClassification_BDT_all_noIPCHI2rew",
                                                    "TMVAClassification_BDT_all_noIPCHI2rew", -1.0, 1.0);
  RooDataSet data_tmp("data_tmp", "data_tmp", RooArgSet(Lb_M_reco, TMVAClassification_BDT_all_noIPCHI2rew),
                      RooFit::ImportFromFile(data_file_name, "DecayTree"),
                      RooFit::Cut("TMVAClassification_BDT_all_noIPCHI2rew>-0.05"));  // the name of the dataset MUST be
                                                                                     // "data" in order for the
                                                                                     // framework to identify it.
  RooDataSet data = *(dynamic_cast<RooDataSet*>(data_tmp.reduce(RooFit::SelectVars(RooArgSet(Lb_M_reco)))));
  data.SetNameTitle("data", "data");
  n_bkg.setVal(data.sumEntries());
  n_bkg.setMax(2.0 * data.sumEntries());

  /////////////////////////////////////////////////////////

  // Set the global observables to the observed values and make them constant for the fit.
  norm_constant_obs.setVal(observedValueGlobalObservable);
  norm_constant_obs.setConstant();

  // The workspace must also include the fit result of an initial fit of the model to data.
  RooFitResult& rooFitResult =
      *mass_model.fitTo(data, RooFit::Save(), RooFit::ExternalConstraints(RooArgSet(norm_constant_constraint)));

  /////////////////////////////////////////////////////////
  //
  // Plot it all so we can see what we did
  //
  /////////////////////////////////////////////////////////

  RooPlot* plot = Lb_M_reco.frame();
  data.plotOn(plot);
  extended_bkg_model.plotOn(plot, RooFit::LineColor(kRed));
  mass_model.plotOn(plot);
  TCanvas c("c", "c", 1024, 768);
  plot->Draw();
  c.SaveAs("plots/pdf/data_and_fit_in_workspace.pdf");

  rooFitResult.Print();

  /////////////////////////////////////////////////////////
  //
  // We also must define some RooArgSets:
  //
  /////////////////////////////////////////////////////////

  // One contaninting the constraint PDFs,
  RooArgSet constraint_set(norm_constant_constraint, "constraint_set");

  // one containing the global Observables,
  RooArgSet global_observables_set(norm_constant_obs, "global_observables_set");

  // one containing the normal Observables (the bin variables in the datasets usually) and
  RooArgSet dataset_observables_set(Lb_M_reco, "datasetObservables");

  // one containing the parameters
  auto params = mass_model.getParameters(RooArgSet(Lb_M_reco));
  params->remove(constraint_set);
  params->remove(global_observables_set);
  RooArgSet parameters_set(*params, "parameters");

  /////////////////////////////////////////////////////////
  //
  // Import everything into a workspace and save it
  //
  /////////////////////////////////////////////////////////

  RooWorkspace workspace("dataset_workspace");
  workspace.import(mass_model);
  workspace.import(extended_bkg_model, RooFit::RenameConflictNodes("_internal"));  // TODO check if this is used
  workspace.import(data);
  workspace.import(rooFitResult, "data_fit_result");  // this MUST be called data_fit_result
  workspace.defineSet("constraint_set", constraint_set, true);
  workspace.defineSet("global_observables_set", global_observables_set, true);
  workspace.defineSet("datasetObservables", dataset_observables_set, true);
  workspace.defineSet("parameters", parameters_set, true);

  // Save the workspace to a file
  workspace.SaveAs("workspace.root");

  return 0;
}
