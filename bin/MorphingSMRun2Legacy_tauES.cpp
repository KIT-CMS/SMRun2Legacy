#include "CombineHarvester/CombinePdfs/interface/MorphFunctions.h"
#include "CombineHarvester/CombineTools/interface/Algorithm.h"
#include "CombineHarvester/CombineTools/interface/AutoRebin.h"
#include "CombineHarvester/CombineTools/interface/BinByBin.h"
#include "CombineHarvester/CombineTools/interface/CardWriter.h"
#include "CombineHarvester/CombineTools/interface/CombineHarvester.h"
#include "CombineHarvester/CombineTools/interface/Observation.h"
#include "CombineHarvester/CombineTools/interface/Process.h"
#include "CombineHarvester/CombineTools/interface/Systematics.h"
#include "CombineHarvester/CombineTools/interface/Utilities.h"
#include "CombineHarvester/SMRun2Legacy/interface/BinomialBinByBin.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "TF1.h"
#include "TH2.h"
#include "boost/algorithm/string/predicate.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/program_options.hpp"
#include "boost/regex.hpp"
#include <cstdlib>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include <math.h>

using namespace std;
using boost::starts_with;
using ch::syst::SystMap;
using ch::syst::SystMapAsymm;
using ch::syst::era;
using ch::syst::channel;
using ch::syst::bin_id;
using ch::syst::process;
using ch::JoinStr;
using namespace ch;
namespace po = boost::program_options;

bool debug=false;
void dout() {
    if (debug) std::cout << std::endl;
}
template <typename Head, typename... Tail>
void dout(Head H, Tail... T) {
    if (debug) std::cout << H << ' ';
    dout(T...);
}
template<typename T>
void printVector(const T& t) {
    std::copy(t.cbegin(), t.cend(), std::ostream_iterator<typename T::value_type>(std::cout, ", "));
}
template<typename T>
void printVectorInVector(const T& t) {
    std::for_each(t.cbegin(), t.cend(), printVector<typename T::value_type>);
}
template <class T>
bool contains(std::vector<T> const &v, T const &x) {
    return ! (v.empty() && std::find(v.begin(), v.end(), x) == v.end());
}


int main(int argc, char **argv) {
  typedef vector<string> VString;
  typedef vector<pair<int, string>> Categories;
  using ch::syst::bin_id;
  using ch::JoinStr;

  // Define program options
  string output_folder = "sm_run2";
  string base_path = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/SMRun2Legacy/shapes";
  string input_folder_mt = "Vienna/";
  string postfix = "-ML";
  string midfix = "";
  bool regional_jec = true;
  bool ggh_wg1 = true;
  bool auto_rebin = false;
  bool rebin_categories = true;
  bool manual_rebin_for_yields = false;
  bool real_data = false;
  bool train_ff = true;
  bool train_emb = true;
  bool train_stage0 = false;
  bool classic_bbb = false;
  bool binomial_bbb = false;
  bool verbose = false;
  bool remove_empty_categories = false;
  string stxs_signals = "stxs_stage0"; // "stxs_stage0" or "stxs_stage1p1"
  string categories = "stxs_stage0"; // "stxs_stage0", "stxs_stage1p1" or "gof"
  string gof_category_name = "gof";
  int era = 2016; // 2016 or 2017
  po::variables_map vm;
  po::options_description config("configuration");
  config.add_options()
      ("base_path", po::value<string>(&base_path)->default_value(base_path))
      ("input_folder_mt", po::value<string>(&input_folder_mt)->default_value(input_folder_mt))
      ("postfix", po::value<string>(&postfix)->default_value(postfix))
      ("midfix", po::value<string>(&midfix)->default_value(midfix))
      ("auto_rebin", po::value<bool>(&auto_rebin)->default_value(auto_rebin))
      ("rebin_categories", po::value<bool>(&rebin_categories)->default_value(rebin_categories))
      ("manual_rebin_for_yields", po::value<bool>(&manual_rebin_for_yields)->default_value(manual_rebin_for_yields))
      ("regional_jec", po::value<bool>(&regional_jec)->default_value(regional_jec))
      ("ggh_wg1", po::value<bool>(&ggh_wg1)->default_value(ggh_wg1))
      ("real_data", po::value<bool>(&real_data)->default_value(real_data))
      ("verbose", po::value<bool>(&verbose)->default_value(verbose))
      ("remove_empty_categories", po::value<bool>(&remove_empty_categories)->default_value(remove_empty_categories))
      ("output_folder", po::value<string>(&output_folder)->default_value(output_folder))
      ("stxs_signals", po::value<string>(&stxs_signals)->default_value(stxs_signals))
      ("categories", po::value<string>(&categories)->default_value(categories))
      ("gof_category_name", po::value<string>(&gof_category_name)->default_value(gof_category_name))
      ("train_ff", po::value<bool>(&train_ff)->default_value(train_ff))
      ("train_emb", po::value<bool>(&train_emb)->default_value(train_emb))
      ("train_stage0", po::value<bool>(&train_stage0)->default_value(train_stage0))
      ("classic_bbb", po::value<bool>(&classic_bbb)->default_value(classic_bbb))
      ("binomial_bbb", po::value<bool>(&binomial_bbb)->default_value(binomial_bbb))
      ("era", po::value<int>(&era)->default_value(era));
  po::store(po::command_line_parser(argc, argv).options(config).run(), vm);
  po::notify(vm);
  
  // Prepering the combine run
  RooRealVar taues("taues", "taues", -1.5, -3.0, 0.0);
  vector<string> energy_scales;
  if (era == 2016) {
  energy_scales = {"-2.4","-2.2","-2.0","-1.8","-1.6","-1.4","-1.2", "-1.0", "-1.0", "-0.8", "-0.7","-0.65","-0.6","-0.55","-0.5","-0.45" , "-0.4","-0.35","-0.3","-0.25", "-0.2", "0.0", "0.2","0.4", "0.6","0.8","1.0","1.2","1.4","1.6"};
  }
  else {
  energy_scales = {"-2.4","-2.2","-2.0", "-1.8", "-1.6", "-1.4", "-1.2", "-1.0", "-0.8", "-0.6", "-0.4", "-0.2", "0.0","0.2","0.4","0.6","0.8","1.0","1.2","1.4","1.6"};
  }
  // Define channels
  VString chns;
  chns.push_back("mt");

  // Define background processes
  map<string, VString> bkg_procs;
  VString bkgs;

  bkgs = {"ZL", "TTL", "VVL","jetFakes"};

  std::cout << "[INFO] Considerung the following processes:\n";
  std::cout << "For mt channel : \n";
  for (unsigned int i=0; i < bkgs.size(); i++) std::cout << bkgs[i] << std::endl;
  bkg_procs["mt"] = bkgs;

  // Define categories
  map<string, Categories> cats;
  std::vector<std::string> cats_to_keep; // will be used later for the card writer
  for (auto chn : chns){
      cats[chn]={
        { 1, chn+"_m_vis"},
      };
  }
  for (auto chn : chns){
    for (auto tuple: cats[chn]) cout << tuple.first << ": " << tuple.second << endl;
  }

  // Specify signal processes
  vector<string> sig_procs = {"EMB"};

  // Create combine harverster object
  ch::CombineHarvester cb;

  // Add observations and processes
  std::string era_tag;
  if (era == 2016) era_tag = "Run2016";
  else if (era == 2017) era_tag = "Run2017";
  else if (era == 2018) era_tag = "Run2018";

  else std::runtime_error("Given era is not implemented.");

  for (auto chn : chns) {
    cb.AddObservations({"*"}, {"htt"}, {era_tag}, {chn}, cats[chn]);
    cb.AddProcesses({"*"}, {"htt"}, {era_tag}, {chn}, bkg_procs[chn], cats[chn],
                    false);
    cb.AddProcesses(energy_scales, {"htt"}, {era_tag}, {chn}, sig_procs, cats[chn],
                    true);
  }
  //auto signal = Set2Vec(cb.cp().signals().SetFromProcs(std::mem_fn(&Process::process)));

  // Add systematics

  // MC uncertainties:
  // lumi
  float lumi_unc = 1.025;
  if (era == 2017) lumi_unc = 1.023;
  cb.cp()
      .process({"ZL", "TTL", "VVL"})
      .AddSyst(cb, "lumi", "lnN", SystMap<>::init(lumi_unc));
  cb.cp()
      .process({"EMB"})
      .AddSyst(cb, "embnorm", "lnN", SystMap<>::init(1.04));
  // VV
  cb.cp()
      .process({"VVL"})
      .AddSyst(cb, "vvlXsec", "lnN", SystMap<>::init(1.05));
  // Z
  cb.cp()
      .process({"ZL"})
      .AddSyst(cb, "zlXsec", "lnN", SystMap<>::init(1.04));
  // TT
  cb.cp()
      .process({"TTL"})
      .AddSyst(cb, "ttlXsec", "lnN", SystMap<>::init(1.06));
  // Muon fakes
  cb.cp()
      .channel({"mt"})
      .process({"ZL","TTL"})
      .AddSyst(cb, "CMS_mFakeTau_$ERA", "lnN", SystMap<>::init(1.25));
  // Embedding Uncertainties

  // Embedded Tau ID in 5 pt bins
  std::string tauIDptbins[5] = {"30-35", "35-40", "40-500", "500-1000", "1000-inf"};
  for (auto tauIDbin : tauIDptbins){
    cb.cp()
        .process({"EMB"})
        .AddSyst(cb, "CMS_eff_emb_t_"+tauIDbin+"_$ERA", "shape", SystMap<>::init(1));
  }
  cb.cp()
      .process({"EMB"})
      .AddSyst(cb, "emb_t_$CHANNEL_$ERA", "lnN", SystMap<>::init(1.01));
  // Muon ID
  cb.cp()
      .process({"ZL", "TTL", "VVL","EMB"})
      .AddSyst(cb, "CMS_eff_m", "lnN", SystMap<>::init(1.02));
    // TTbar contamination 
  // cb.cp()
  //   .process({"EMB"})
  //   .AddSyst(cb, "emb_ttbar_$ERA", "shape", SystMap<>::init(1.00));
  // hadronic tau track efficiency 
  cb.cp()
    .process({"EMB"})
    .AddSyst(cb, "CMS_3ProngEff_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
    .process({"EMB"})
    .AddSyst(cb, "CMS_1ProngPi0Eff_$ERA", "shape", SystMap<>::init(1.00));
  // Jet Fake Stuff
// QCD shape stat.
  cb.cp()
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_qcd_njet0_$CHANNEL_stat_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_qcd_njet1_$CHANNEL_stat_$ERA", "shape", SystMap<>::init(1.00));

  // W shape stat.
  cb.cp()
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_w_njet0_$CHANNEL_stat_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_w_njet1_$CHANNEL_stat_$ERA", "shape", SystMap<>::init(1.00));

  // TT shape stat.
  cb.cp()
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_tt_njet1_stat_$ERA", "shape", SystMap<>::init(1.00));

  // Shape syst. of different contributions (QCD/W/tt)
  // uncorrelated between eras
  cb.cp()
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_qcd_mt_syst_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_w_syst_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_tt_syst_$ERA", "shape", SystMap<>::init(1.00));
  // Define the location of the "auxiliaries" directory where we can
  // source the input files containing the datacard shapes
  std::map<string, string> input_dir;
  input_dir["mt"] = base_path + "/" + input_folder_mt + "/";
  // Extract shapes from input ROOT files
  for (string chn : chns) {
    cb.cp().channel({chn}).backgrounds().ExtractShapes(
        input_dir[chn] + std::to_string(era) + midfix + chn + "-synced"+postfix+ ".root",
        "$BIN/$PROCESS", "$BIN/$PROCESS_$SYSTEMATIC");
    cb.cp().channel({chn}).process(sig_procs).ExtractShapes(
        input_dir[chn] + std::to_string(era) + midfix + chn + "-synced"+postfix+ ".root",
        "$BIN/$PROCESS_$MASS", "$BIN/$PROCESS_$MASS_$SYSTEMATIC");
  }

  // Replacing observation with the sum of the backgrounds (Asimov data)
  // useful to be able to check this, so don't do the replacement
  // for these
  if (!real_data) {
    for (auto b : cb.cp().bin_set()) {
      std::cout << "[INFO] Replacing data with asimov in bin " << b << "\n";
      auto background_shape = cb.cp().bin({b}).backgrounds().GetShape();
      auto signal_shape = cb.cp().bin({b}).signals().GetShape();
      auto total_procs_shape = cb.cp().bin({b}).data().GetShape();
      total_procs_shape.Scale(0.0);
      bool no_signal = (signal_shape.GetNbinsX() == 1 && signal_shape.Integral() == 0.0);
      bool no_background = (background_shape.GetNbinsX() == 1 && background_shape.Integral() == 0.0);
      if(no_signal && no_background)
      {
        std::cout << "\t[WARNING] No signal and no background available in bin " << b << std::endl;
      }
      else if(no_background)
      {
        std::cout << "\t[WARNING] No background available in bin " << b << std::endl;
        total_procs_shape = total_procs_shape + signal_shape;
      }
      else if(no_signal)
      {
        std::cout << "\t[WARNING] No signal available in bin " << b << std::endl;
        total_procs_shape = total_procs_shape + background_shape;
      }
      else
      {
        total_procs_shape = total_procs_shape + background_shape + signal_shape;
      }
      cb.cp().bin({b}).ForEachObs([&](ch::Observation *obs) {
        obs->set_shape(total_procs_shape,true);
      });
    }
  }

  // At this point we can fix the negative bins
  std::cout << "[INFO] Fixing negative bins.\n";
  cb.ForEachProc([](ch::Process *p) {
    if (ch::HasNegativeBins(p->shape())) {
      auto newhist = p->ClonedShape();
      ch::ZeroNegativeBins(newhist.get());
      p->set_shape(std::move(newhist), false);
    }
  });

  cb.ForEachSyst([](ch::Systematic *s) {
    if (s->type().find("shape") == std::string::npos)
      return;
    if (ch::HasNegativeBins(s->shape_u()) ||
        ch::HasNegativeBins(s->shape_d())) {
      auto newhist_u = s->ClonedShapeU();
      auto newhist_d = s->ClonedShapeD();
      ch::ZeroNegativeBins(newhist_u.get());
      ch::ZeroNegativeBins(newhist_d.get());
      s->set_shapes(std::move(newhist_u), std::move(newhist_d), nullptr);
    }
  });

  // Merge bins and set bin-by-bin uncertainties if no autoMCStats is used.
  if (classic_bbb) {
    auto bbb = ch::BinByBinFactory()
                   .SetAddThreshold(0.0)
                   .SetMergeThreshold(0.5)
                   .SetFixNorm(false);
    //bbb.MergeBinErrors(cb.cp().backgrounds());
    bbb.AddBinByBin(cb.cp().backgrounds(), cb);
    //bbb.MergeBinErrors(cb.cp().signals());
    bbb.AddBinByBin(cb.cp().signals(), cb);
  }
  if (binomial_bbb) {
    auto bbb = ch::BinomialBinByBinFactory()
                   .SetPattern("CMS_$ANALYSIS_$CHANNEL_$BIN_$ERA_$PROCESS_binomial_bin_$#")
                   .SetBinomialP(0.022)
                   .SetBinomialN(1000.0)
                   .SetFixNorm(false);
    bbb.AddBinomialBinByBin(cb.cp().channel({"mt"}).process({"EMB"}), cb);
  }

  // This function modifies every entry to have a standardised bin name of
  // the form: {analysis}_{channel}_{bin_id}_{era}
  ch::SetStandardBinNames(cb, "$ANALYSIS_$CHANNEL_$BINID_$ERA");
  
  
  dout("First we generate a set of bin names:");
  RooWorkspace ws("htt", "htt");
  string demo_file = "htt_mssm_demo.root";
  TFile demo(demo_file.c_str(), "RECREATE");
  bool do_morphing = true;
  if (do_morphing)
  {
      auto bins = cb.bin_set();
      for (auto b : bins)
      {
          auto procs = cb.cp().bin({b}).signals().process_set();
          for (auto p : procs)
          {
              // ws, cb, std::string const& bin, std::string const& process, RooAbsReal& mass_var, std::string norm_postfix, bool allow_morph, bool verbose, bool force_template_limit, TFile * file) {
              ch::BuildRooMorphing(ws, cb, b, p, taues, "norm", true, true, false, NULL);
              // ch::BuildRooMorphing(ws, cb, b, p, faketaues, "norm", true, true, true, NULL);
          }
      }
  }

  ws.var("taues")->setVal(1.0);
  // ws.var("CMS_th1x_htt_et_2_13TeV")->setVal(0.0);

  demo.Close();
  dout("AddWorkspace:");
  cb.AddWorkspace(ws);
  dout("ExtractPdfs:");
  cb.cp().process({"EMB"}).ExtractPdfs(cb, "htt", "$BIN_$PROCESS_morph");
  dout("cb.PrintAll():");
  cb.PrintAll();

  //string foldercmb = output_dir + "cmb";
  //boost::filesystem::create_directories(foldercmb);

  // Write out datacards. Naming convention important for rest of workflow. We
  // make one directory per chn-cat, one per chn and cmb. In this code we only
  // store the individual datacards for each directory to be combined later.
  ch::CardWriter writer(output_folder + "/$TAG/$MASS/$BIN.txt",
    output_folder +"/$TAG/common/htt_input_" + era_tag + ".root");

  // We're not using mass as an identifier - which we need to tell the
  // CardWriter
  // otherwise it will see "*" as the mass value for every object and skip it
  writer.SetWildcardMasses({});

  // Set verbosity
  if (verbose)
    writer.SetVerbosity(1);

  // Write datacards combined and per channel
  writer.WriteCards("cmb", cb);

  for (auto chn : chns) {
    writer.WriteCards(chn, cb.cp().channel({chn}));
  }

  if (verbose)
    cb.PrintAll();

  cout << "[INFO] Done producing datacards.\n";
}
