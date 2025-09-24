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
#include "CombineHarvester/SMRun2Legacy/interface/HttSystematics_SMRun2.h"
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
namespace po = boost::program_options;

int main(int argc, char **argv) {
  typedef vector<string> VString;
  typedef vector<pair<int, string>> Categories;
  using ch::syst::bin_id;
  using ch::JoinStr;

  // Define program options
  string output_folder = "sm_run2";
  string base_path = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/SMRun2Legacy/shapes";
  string input_folder_et = "Vienna/";
  string input_folder_mt = "Vienna/";
  string input_folder_tt = "Vienna/";
  string chan = "all";
  string postfix = "-ML";
  bool auto_rebin = false;
  bool rebin_categories = true;
  bool manual_rebin_for_yields = false;
  bool jetfakes = false;
  bool embedding = false;
  bool use_automc = true;
  bool verbose = false;
  bool remove_empty_categories = false;
  string gof_category_name = "gof";
  int era = 2016; // 2016 or 2017
  po::variables_map vm;
  po::options_description config("configuration");
  config.add_options()
      ("base_path", po::value<string>(&base_path)->default_value(base_path))
      ("input_folder_et", po::value<string>(&input_folder_et)->default_value(input_folder_et))
      ("input_folder_mt", po::value<string>(&input_folder_mt)->default_value(input_folder_mt))
      ("input_folder_tt", po::value<string>(&input_folder_tt)->default_value(input_folder_tt))
      ("postfix", po::value<string>(&postfix)->default_value(postfix))
      ("channel", po::value<string>(&chan)->default_value(chan))
      ("auto_rebin", po::value<bool>(&auto_rebin)->default_value(auto_rebin))
      ("use_automc", po::value<bool>(&use_automc)->default_value(use_automc))
      ("rebin_categories", po::value<bool>(&rebin_categories)->default_value(rebin_categories))
      ("manual_rebin_for_yields", po::value<bool>(&manual_rebin_for_yields)->default_value(manual_rebin_for_yields))
      ("verbose", po::value<bool>(&verbose)->default_value(verbose))
      ("remove_empty_categories", po::value<bool>(&remove_empty_categories)->default_value(remove_empty_categories))
      ("output_folder", po::value<string>(&output_folder)->default_value(output_folder))
      ("era", po::value<int>(&era)->default_value(era));
  po::store(po::command_line_parser(argc, argv).options(config).run(), vm);
  po::notify(vm);

  // Define channels
  VString chns;
  if (chan.find("mt") != std::string::npos)
    chns.push_back("mt");
  if (chan.find("et") != std::string::npos)
    chns.push_back("et");
  if (chan.find("tt") != std::string::npos)
    chns.push_back("tt");
  if (chan == "all")
    chns = {"mt", "et", "tt"};

  // Define background processes
  map<string, VString> bkg_procs;
  VString bkgs;
  bkgs = {
    "W",
    "EWK",
    "ZTT", "ZL", "ZJ", 
    "TTT", "TTL", "TTJ", 
    "VVJ", "VVT", "VVL", 
    "VVV-VVVJ", "VVV-VVVT", "VVV-VVVL",
    "ST-STT", "ST-STJ", "ST-STL",
    "TTV-TTVJ", "TTV-TTVL", "TTV-TTVT",
    "ggH_htt125", "qqH_htt125", "ttH_htt125", "VH125",
  };

  std::cout << "[INFO] Considerung the following processes:\n";
  if (chan.find("mt") != std::string::npos || chan.find("et") != std::string::npos || chan.find("tt") != std::string::npos) {
    std::cout << "For et,mt,tt channels : \n";
    for (unsigned int i=0; i < bkgs.size(); i++) std::cout << bkgs[i] << std::endl;}
  bkg_procs["et"] = bkgs;
  bkg_procs["mt"] = bkgs;
  bkg_procs["tt"] = bkgs;

  // Define categories
  map<string, Categories> cats;
  std::vector<std::string> cats_to_keep; // will be used later for the card writer
  for (auto chn : chns){
    cats[chn]={
      { 1, chn+"_HH2B2Tau"},
      { 2, chn+"_DY"},
      { 3, chn+"_ST"},
      { 4, chn+"_TT"},
      { 5, chn+"_VV"},
      { 6, chn+"_Other"},
    };

  }
  for (auto chn : chns){
    for (auto tuple: cats[chn]) cout << tuple.first << ": " << tuple.second << endl;
  }

  // Specify signal processes and masses
  vector<string> sig_procs;
  sig_procs = {"HH2B2Tau"};

  vector<string> masses = {"125"};

  // Create combine harverster object
  ch::CombineHarvester cb;

  // Add observations and processes
  std::string era_tag;
  if (era == 2016) era_tag = "2016";
  else if (era == 2017) era_tag = "2017";
  else if (era == 2018) era_tag = "2018";
  else throw std::runtime_error("Given era is not implemented: " + std::to_string(era));

  for (auto chn : chns) {
    cb.AddObservations({"*"}, {"htt"}, {era_tag}, {chn}, cats[chn]);
    cb.AddProcesses({"*"}, {"htt"}, {era_tag}, {chn}, bkg_procs[chn], cats[chn],
                    false);
    // cb.AddProcesses(masses, {"htt"}, {era_tag}, {chn}, sig_procs, cats[chn],
    //                 true);
    cb.AddProcesses({""}, {"htt"}, {era_tag}, {chn}, sig_procs, cats[chn],
                    true);
  }

  // Add systematics
  // ch::AddSMRun2Systematics(cb, jetfakes, embedding, regional_jec, ggh_wg1, qqh_wg1, era);

  // Define the location of the "auxiliaries" directory where we can
  // source the input files containing the datacard shapes
  std::map<string, string> input_dir;
  input_dir["mt"] = base_path + "/" + input_folder_mt + "/";
  input_dir["et"] = base_path + "/" + input_folder_et + "/";
  input_dir["tt"] = base_path + "/" + input_folder_tt + "/";
  // Extract shapes from input ROOT files
  for (string chn : chns) {
    cb.cp().channel({chn}).backgrounds().ExtractShapes(
        input_dir[chn] + "htt_" + chn + ".inputs-sm-Run" + std::to_string(era) + postfix + ".root",
        "$BIN/$PROCESS", "$BIN/$PROCESS_$SYSTEMATIC");
    cb.cp().channel({chn}).process(sig_procs).ExtractShapes(
        input_dir[chn] + "htt_" + chn + ".inputs-sm-Run" + std::to_string(era) + postfix + ".root",
        "$BIN/$PROCESS$MASS", "$BIN/$PROCESS$MASS_$SYSTEMATIC");
  }

  // Rebin categories to predefined binning for binning
  // if (rebin_categories) {
  //   // Rebin background categories
  //   for (auto b : cb.cp().bin_set()) {
  //     TString bstr = b;
  //     std::cout << "[INFO] Rebin background bin " << b << "\n";
  //     auto shape = cb.cp().bin({b}).backgrounds().GetShape();
  //     auto min = shape.GetBinLowEdge(1);
  //     else if(bstr.Contains("et") && bstr.Contains("misc")) cb.cp().bin({b}).VariableRebin({min, 0.4, 0.5, 0.6, 1.0});
  //     else if(bstr.Contains("mt") && bstr.Contains("misc")) cb.cp().bin({b}).VariableRebin({min, 0.4, 0.5, 0.6, 1.0});
  //     else cb.cp().bin({b}).VariableRebin({min, 0.4, 0.5, 0.6, 0.7, 1.0});
  //   }
  // }

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

  // Perform auto-rebinning
  if (auto_rebin) {
    const auto threshold = 10.0;
    for (auto b : cb.cp().bin_set()) {
      std::cout << "[INFO] Rebin bin " << b << "\n";
      // Get shape of this category with sum of backgrounds
      auto shape = cb.cp().bin({b}).backgrounds().GetShape();
      // Push back last bin edge
      vector<double> binning;
      const auto num_bins = shape.GetNbinsX();
      binning.push_back(shape.GetBinLowEdge(num_bins + 1));
      // Now, go backwards through bins (from right to left) and merge a bin if
      // the background yield is below a given threshold.
      auto c = 0.0;
      for(auto i = num_bins; i > 0; i--) {
        auto low_edge = shape.GetBinLowEdge(i);
        c += shape.GetBinContent(i);
        if (c > threshold) { // Set lower edge if the bin content is above the threshold.
          binning.insert(binning.begin(), low_edge);
          c = 0.0;
        }
      }
      if (binning.size() == 1){ // catching case, if the total yield of the histogram is smaller then threshold.
        binning.insert(binning.begin(), shape.GetBinLowEdge(1));
      }
      binning.at(0)=shape.GetBinLowEdge(1); // in case that yield of lowest bin is smaller than threshold merge it with second lowest
      cb.cp().bin({b}).VariableRebin(binning);
    }
    // Remove categories with too little events, if specified
    for (auto b : cb.cp().bin_set()) {
      // Get yield of all backgrounds in this category
      float shape_integral = cb.cp().bin({b}).backgrounds().GetShape().Integral();
      if(shape_integral < threshold && remove_empty_categories){
        std::cout << "[WARNING] Remove category " << b << " due to insufficient population!" << "\n";
      }
      else {
        cats_to_keep.push_back(b);
      }
    }
  }
  else {
    for (auto b : cb.cp().bin_set()) {
        cats_to_keep.push_back(b);
    }
  }
  cb = cb.bin(cats_to_keep);

  if(manual_rebin_for_yields) {
    for(auto b : cb.cp().bin_set()) {
      std::cout << "Rebinning by hand for bin: " << b <<  std::endl;
      cb.cp().bin({b}).VariableRebin({0.0,1.0});
    }
  }

  if (use_automc) {
    std::cout << "[INFO] Adding SetAutoMCStats .\n";
    cb.SetAutoMCStats(cb, 0.);
  }


  // This function modifies every entry to have a standardised bin name of
  // the form: {analysis}_{channel}_{bin_id}_{era}
  ch::SetStandardBinNames(cb, "$ANALYSIS_$CHANNEL_$BINID_$ERA");

  // Write out datacards. Naming convention important for rest of workflow. We
  // make one directory per chn-cat, one per chn and cmb. In this code we only
  // store the individual datacards for each directory to be combined later.
  ch::CardWriter writer(output_folder + "/$TAG/$MASS/$BIN.txt",
    output_folder +"/$TAG/common/htt_input_" + era_tag + ".root");

  // We're not using mass as an identifier - which we need to tell the
  // CardWriter
  // otherwise it will see "*" as the mass value for every object and skip it
  //    writer.SetWildcardMasses({});

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
