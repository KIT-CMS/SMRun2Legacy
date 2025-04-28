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
#include "CombineHarvester/SMRun2Legacy/interface/HttSystematics_NMSSMRun2UL.h"
#include "CombineHarvester/SMRun2Legacy/interface/HttSystematics_NMSSMboostedRun2UL.h"
#include "CombineHarvester/SMRun2Legacy/interface/BinomialBinByBin.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "TF1.h"
#include "TH1.h"
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
  string output_folder = "nmssm_run2";
  string base_path = string(getenv("CMSSW_BASE")) + "/src/CombineHarvester/SMRun2Legacy/shapes";
  string input_folder_em = "Vienna/";
  string input_folder_et = "Vienna/";
  string input_folder_mt = "Vienna/";
  string input_folder_tt = "Vienna/";
  string chan = "all";
  string postfix = "-ML";
  string midfix = "";
  bool regional_jec = true;
  bool auto_rebin = false;
  bool rebin_categories = true;
  bool manual_rebin_for_yields = false;
  bool real_data = false;
  bool jetfakes = true;
  bool train_ff = true;
  bool train_emb = true;
  bool embedding = false;
  bool use_automc = true;
  bool classic_bbb = false;
  bool binomial_bbb = false;
  bool verbose = false;
  bool remove_empty_categories = false;
  bool boosted_tt = false;
  string heavy_mass = "all";
  string light_mass = "all";
  string training_mass = "all";
  string training_batch = "all";
  string categories = "nmssm"; // "stxs_stage0", "stxs_stage1p1" or "gof"
  string gof_category_name = "gof";
  int era = 2018; // 2016 or 2017
  po::variables_map vm;
  po::options_description config("configuration");
  config.add_options()
      ("base_path", po::value<string>(&base_path)->default_value(base_path))
      ("input_folder_et", po::value<string>(&input_folder_et)->default_value(input_folder_et))
      ("input_folder_mt", po::value<string>(&input_folder_mt)->default_value(input_folder_mt))
      ("input_folder_tt", po::value<string>(&input_folder_tt)->default_value(input_folder_tt))
      ("boosted_tt", po::value<bool>(&boosted_tt)->default_value(boosted_tt))
      ("heavy_mass", po::value<string>(&heavy_mass)->default_value(heavy_mass))
      ("light_mass", po::value<string>(&light_mass)->default_value(light_mass))
      ("real_data", po::value<bool>(&real_data)->default_value(real_data))
      ("classic_bbb", po::value<bool>(&classic_bbb)->default_value(classic_bbb))
      ("binomial_bbb", po::value<bool>(&binomial_bbb)->default_value(binomial_bbb))
      ("jetfakes", po::value<bool>(&jetfakes)->default_value(jetfakes))
      ("embedding", po::value<bool>(&embedding)->default_value(embedding))
      ("postfix", po::value<string>(&postfix)->default_value(postfix))
      ("chan", po::value<string>(&chan)->default_value(chan))
      ("auto_rebin", po::value<bool>(&auto_rebin)->default_value(auto_rebin))
      ("era", po::value<int>(&era)->default_value(era))
      ("output_folder", po::value<string>(&output_folder)->default_value(output_folder))
      ("use_automc", po::value<bool>(&use_automc)->default_value(use_automc))
      ("remove_empty_categories", po::value<bool>(&remove_empty_categories)->default_value(remove_empty_categories))
      ("train_ff", po::value<bool>(&train_ff)->default_value(train_ff))
      ("train_emb", po::value<bool>(&train_emb)->default_value(train_emb))
      ("categories", po::value<string>(&categories)->default_value(categories))
      ("verbose", po::value<bool>(&verbose)->default_value(verbose))
      ("gof_category_name", po::value<string>(&gof_category_name)->default_value(gof_category_name));
      // ("midfix", po::value<string>(&midfix)->default_value(midfix))
      // ("rebin_categories", po::value<bool>(&rebin_categories)->default_value(rebin_categories))
      // ("manual_rebin_for_yields", po::value<bool>(&manual_rebin_for_yields)->default_value(manual_rebin_for_yields))
      // ("regional_jec", po::value<bool>(&regional_jec)->default_value(regional_jec))
      // ("training_mass", po::value<string>(&training_mass)->default_value(training_mass))
      // ("training_batch", po::value<string>(&training_batch)->default_value(training_batch))
  po::store(po::command_line_parser(argc, argv).options(config).run(), vm);
  po::notify(vm);

  // output_folder = output_folder + "/" + era + categories + "_" + heavy_mass + "_" + light_mass;

  std::cout << "[INFO] Channels: " << chan << std::endl;
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
  bkgs = {"W", "ZTT_NLO", "QCD", "ZL_NLO", "ZJ_NLO", "TTT", "TTL", "TTJ", "STT", "STL", "STJ", "VVJ", "VVT", "VVL", "ggH_tt", "qqH_tt", "VH_tt", "ggH_bb", "qqH_bb", "VH_bb"};
  // todo added in next iteration: ttH_htt

  if(embedding){
    bkgs.erase(std::remove(bkgs.begin(), bkgs.end(), "ZTT_NLO"), bkgs.end());
    bkgs.erase(std::remove(bkgs.begin(), bkgs.end(), "TTT"), bkgs.end());
    bkgs.erase(std::remove(bkgs.begin(), bkgs.end(), "STT"), bkgs.end());
    bkgs.erase(std::remove(bkgs.begin(), bkgs.end(), "VVT"), bkgs.end());
    bkgs = JoinStr({bkgs,{"EMB"}});
  }
  if(jetfakes){
    bkgs.erase(std::remove(bkgs.begin(), bkgs.end(), "QCD"), bkgs.end());
    bkgs.erase(std::remove(bkgs.begin(), bkgs.end(), "W"), bkgs.end());
    bkgs.erase(std::remove(bkgs.begin(), bkgs.end(), "VVJ"), bkgs.end());
    bkgs.erase(std::remove(bkgs.begin(), bkgs.end(), "TTJ"), bkgs.end());
    bkgs.erase(std::remove(bkgs.begin(), bkgs.end(), "STJ"), bkgs.end());
    bkgs.erase(std::remove(bkgs.begin(), bkgs.end(), "ZJ_NLO"), bkgs.end());
    bkgs = JoinStr({bkgs,{"jetFakes"}});
  }

  std::cout << "[INFO] Considerung the following processes:\n";
  if (chan.find("mt") != std::string::npos || chan.find("et") != std::string::npos || chan.find("tt") != std::string::npos) {
    std::cout << "For et,mt,tt channels : \n";
    for (unsigned int i=0; i < bkgs.size(); i++) std::cout << bkgs[i] << std::endl;
  }
  bkg_procs["et"] = bkgs;
  bkg_procs["mt"] = bkgs;
  bkg_procs["tt"] = bkgs;

  // Specify signal processes and masses
  vector<string> sig_procs;
  sig_procs = {
    "NMSSM_Ytt", "NMSSM_Ybb"
  };
  // sig_procs = {
  //   "NMSSM_"+heavy_mass+"_125_"+light_mass
  // };
  // Define NMSSM model-dependent mass parameters mH, mhprime, mh
  // RooRealVar mX("mX", "mX", 700., 240., 4000.);
  // RooRealVar mY("mY", "mY", 250., 60., 2800.);
  // RooRealVar mH("mH", "mH", 125., 124.9, 125.1);
  // mX.setConstant(true);
  // mY.setConstant(true);
  // mH.setConstant(true);

  // Define categories
  map<string, Categories> cats;
  std::vector<std::string> cats_to_keep; // will be used later for the card writer
  for (auto chn : chns){
    //define mapping for signal categories
    if(categories == "nmssm" && !boosted_tt){
      cats[chn]={
        {0, chn+"_YbbHtt_res"},
        {1, chn+"_YbbHtt_boost"},
        {2, chn+"_YttHbb_res"},
        {3, chn+"_YttHbb_boost"},
        {4, chn+"_genuine_tau"},
        {5, chn+"_tau_fakes"},
        {6, chn+"_ttbar"},
        {7, chn+"_misc"},
      };
    }
    else if(categories == "nmssm" && boosted_tt){
      cats[chn]={
        {0, chn+"_boosted_YbbHtt_res"},
        {1, chn+"_boosted_YbbHtt_boost"},
        {2, chn+"_boosted_YttHbb_res"},
        {3, chn+"_boosted_YttHbb_boost"},
        {4, chn+"_boosted_genuine_tau"},
        {5, chn+"_boosted_tau_fakes"},
        {6, chn+"_boosted_ttbar"},
        {7, chn+"_boosted_misc"},
      };
    }
    else if(categories == "gof") cats[chn] = { {300, chn+"_"+gof_category_name.c_str() }};
    else throw std::runtime_error("Given categorization " + categories + " is not known.");
  }
  for (auto chn : chns){
    for (auto tuple: cats[chn]) cout << tuple.first << ": " << tuple.second << endl;
  }

  // Create combine harverster object
  ch::CombineHarvester cb;
  cb.SetFlag("workspaces-use-clone", true);

  // Add observations and processes
  std::string era_tag;
  if (era == 2016) era_tag = "2016";
  else if (era == 2017) era_tag = "2017";
  else if (era == 2018) era_tag = "2018";

  else std::runtime_error("Given era is not implemented.");

  for (auto chn : chns) {
    cb.AddObservations({light_mass}, {"nmssm"}, {era_tag}, {chn}, cats[chn]);
    cb.AddProcesses({light_mass}, {"nmssm"}, {era_tag}, {chn}, bkg_procs[chn], cats[chn],
                    false);
    cb.AddProcesses({light_mass}, {"nmssm"}, {era_tag}, {chn}, sig_procs, cats[chn],
                    true);
  }

  // Add systematics
  if (boosted_tt) ch::AddRun2BoostedSystematics(cb, jetfakes, embedding, era);
  else ch::AddRun2Systematics(cb, jetfakes, embedding, era);

  // Define the location of the "auxiliaries" directory where we can
  // source the input files containing the datacard shapes
  std::map<string, string> input_dir;
  input_dir["mt"] = base_path + "/" + input_folder_mt; 
  input_dir["et"] = base_path + "/" + input_folder_et;
  input_dir["tt"] = base_path + "/" + input_folder_tt;
  // Extract shapes from input ROOT files
  for (string chn : chns) {
    cb.cp().channel({chn}).backgrounds().ExtractShapes(
        input_dir[chn], 
        "$BIN/$PROCESS", "$BIN/$PROCESS_$SYSTEMATIC");
    cb.cp().channel({chn}).process(sig_procs).ExtractShapes(
        input_dir[chn], 
        "$BIN/$PROCESS", "$BIN/$PROCESS_$SYSTEMATIC");
  }
  
  // Replacing observation with the sum of the backgrounds (Asimov data)
  // useful to be able to check this, so don't do the replacement
  // for these
  if (!real_data) {
    for (auto b : cb.cp().bin_set()) {
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

  // // Rebin categories to predefined binning for binning
  // if (rebin_categories) {
  //   // Rebin background categories
  //   for (auto b : cb.cp().bin_set()) {
  //     TString bstr = b;
  //     if (bstr.Contains("ggh") || bstr.Contains("qqh") || bstr.Contains("vbftopo") || bstr.Contains("xxh")) continue;
  //     std::cout << "[INFO] Rebin background bin " << b << "\n";
  //     auto shape = cb.cp().bin({b}).backgrounds().GetShape();
  //     auto min = shape.GetBinLowEdge(1);
  //     if(bstr.Contains("em") && bstr.Contains("misc")) cb.cp().bin({b}).VariableRebin({min, 0.4, 1.0});
  //     else if(bstr.Contains("em_emb")){
  //       if(categories == "stxs_stage1p1") cb.cp().bin({b}).VariableRebin({min, 0.3, 1.0});
  //       else cb.cp().bin({b}).VariableRebin({min, 0.4, 0.5, 0.6, 1.0});
  //     }
  //     else if(bstr.Contains("et") && bstr.Contains("misc")) cb.cp().bin({b}).VariableRebin({min, 0.4, 0.5, 0.6, 1.0});
  //     else if(bstr.Contains("mt") && bstr.Contains("misc")) cb.cp().bin({b}).VariableRebin({min, 0.4, 0.5, 0.6, 1.0});
  //     else if(bstr.Contains("et") && bstr.Contains("zll") && categories == "stxs_stage1p1") cb.cp().bin({b}).VariableRebin({min, 0.4, 0.5, 1.0});
  //     else if(bstr.Contains("mt") && bstr.Contains("emb")){
  //       if(categories == "stxs_stage1p1") cb.cp().bin({b}).VariableRebin({min, 0.4, 0.45, 0.5, 0.6, 1.0});
  //       else cb.cp().bin({b}).VariableRebin({min, 0.4, 0.5, 0.6, 1.0});
  //     }
  //     else cb.cp().bin({b}).VariableRebin({min, 0.4, 0.5, 0.6, 0.7, 1.0});
  //   }
  //   // Rebin ggh stage 1.1 categories
  //   for (auto b : cb.cp().bin_set()) {
  //     TString bstr = b;
  //     if (bstr.Contains("ggh_10")) {
  //       std::cout << "[INFO] Rebin ggh signal bin " << b << "\n";
  //       auto shape = cb.cp().bin({b}).backgrounds().GetShape();
  //       auto min = shape.GetBinLowEdge(1);
  //       cb.cp().bin({b}).VariableRebin({min, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7, 1.0});
  //     }
  //   }
  //   // Rebin qqh stage 1.1 categories
  //   for (auto b : cb.cp().bin_set()) {
  //     TString bstr = b;
  //     if (bstr.Contains("qqh_20")) {
  //       std::cout << "[INFO] Rebin qqh signal bin " << b << "\n";
  //       auto shape = cb.cp().bin({b}).backgrounds().GetShape();
  //       auto min = shape.GetBinLowEdge(1);
  //       cb.cp().bin({b}).VariableRebin({min, 0.4, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.90, 0.95, 1.0});
  //     }
  //   }
  // }

  // Loop over all shape systematic uncertainties that contain the jes_str in their name to symmetrize them
  // std::string jes_str = "scale_j";

  // cb.ForEachProc([&](ch::Process *p) {
  //   cb.ForEachSyst([&](ch::Systematic *s) {
  //     if (MatchingProcess(*p, *s)) {
  //       // Check if the systematic name contains the target string
  //       if (s->name().find(jes_str) != std::string::npos) {
  //         if (s->type() == "shape") {
  //           // std::cout << "[JES] Uncertainty " << s->name() << " is symmetrized for " << s->process() << "\n";
  //           // Get the nominal histogram from the process
  //           auto nominal = p->ClonedShape();

  //           // Get the up and down histograms
  //           auto newhist_u = s->ClonedShapeU();
  //           auto newhist_d = s->ClonedShapeD();

  //           // Create histograms to hold the symmetrized version
  //           // TH1* nom = (TH1*)nominal->Clone();
  //           TH1* h_sym_up = (TH1*)newhist_u->Clone();
  //           TH1* h_sym_down = (TH1*)newhist_d->Clone();

  //           for (int bin = 1; bin <= nominal->GetNbinsX(); ++bin) {
  //             // double sys_diff = std::abs(h_sym_up->GetBinContent(bin) - h_sym_down->GetBinContent(bin))/ 2.;
  //             double up_diff = h_sym_up->GetBinContent(bin) - nominal->GetBinContent(bin);
  //             double down_diff = nominal->GetBinContent(bin) - h_sym_down->GetBinContent(bin);
  //             double max_diff = std::max(std::abs(up_diff), std::abs(down_diff));
  //             // std::cout << "[JES] Bin " << bin << " up yield before " << h_sym_up->GetBinContent(bin) << "\n";
  //             // std::cout << "[JES] Bin " << bin << " down yield before " << h_sym_down->GetBinContent(bin) << "\n";
  //             // std::cout << "[JES] Bin " << bin << " nominal yield " << nominal->GetBinContent(bin) << "\n";
  //             // Symmetrize by setting the same deviation for both up and down
  //             h_sym_up->SetBinContent(bin, nominal->GetBinContent(bin) + max_diff);
  //             if ((nominal->GetBinContent(bin) - max_diff)>=0) {
  //               h_sym_down->SetBinContent(bin, nominal->GetBinContent(bin) - max_diff);
  //             }
  //             else {
  //               h_sym_down->SetBinContent(bin, 0.);
  //             }
  //             // std::cout << "[JES] Bin " << bin << " up yield after " << h_sym_up->GetBinContent(bin) << "\n";
  //             // std::cout << "[JES] Bin " << bin << " down yield after " << h_sym_down->GetBinContent(bin) << "\n";
  //           }

  //           // const TH1& up = *h_sym_up;
  //           // const TH1& down = *h_sym_down;
  //           // const TH1& nom_hist = *nom;
  //           std::unique_ptr<TH1> up(h_sym_up);
  //           std::unique_ptr<TH1> down(h_sym_down);
  //           // std::unique_ptr<TH1> nom_hist(nom);

  //           // Set the symmetrized shapes back to the systematic
  //           s->set_shapes(std::move(up), std::move(down), nullptr);
  //           // std::cout << ch::Systematic::PrintHeader << *s << "\n";
  //         }
  //       }
  //     }
  //   });
  // });

  // At this point we can fix the negative bins
  // std::cout << "[INFO] Fixing negative bins.\n";
  // cb.ForEachProc([](ch::Process *p) {
  //   // auto newhist = p->ClonedShape();
  //   // const auto num_bins = newhist->GetNbinsX();
  //   // for(auto i = num_bins; i > 0; i--) {
  //   //   if (newhist.get()->GetBinContent(i) < 0. || std::isnan(newhist.get()->GetBinContent(i))) { // Set lower edge if the bin content is above the threshold.
  //   //     std::cout << "[WARNING] Old bin " << newhist.get()->GetBinContent(i) << "\n";
  //   //     newhist.get()->SetBinContent(i, 0.);
  //   //     std::cout << "[WARNING] New bin " << newhist.get()->GetBinContent(i) << "\n";
  //   //   }
  //   // }
  //   // p->set_shape(std::move(newhist), false);
  //   if (ch::HasNegativeBins(p->shape())) {
  //     auto newhist = p->ClonedShape();
  //     ch::ZeroNegativeBins(newhist.get());
  //     p->set_shape(std::move(newhist), false);
  //   }
  // });

  // cb.ForEachSyst([](ch::Systematic *s) {
  //   if (s->type().find("shape") == std::string::npos)
  //     return;
  //   if (ch::HasNegativeBins(s->shape_u()) ||
  //       ch::HasNegativeBins(s->shape_d())) {
  //     auto newhist_u = s->ClonedShapeU();
  //     auto newhist_d = s->ClonedShapeD();
  //     ch::ZeroNegativeBins(newhist_u.get());
  //     ch::ZeroNegativeBins(newhist_d.get());
  //     s->set_shapes(std::move(newhist_u), std::move(newhist_d), nullptr);
  //   }
  // });

  // Fix the negative bins again for JES
  // std::cout << "[INFO] Fixing negative bins again.\n";
  // cb.ForEachProc([&](ch::Process *p) {
  //   cb.ForEachSyst([&](ch::Systematic *s) {
  //     if (MatchingProcess(*p, *s)) {
  //       if (s->type().find("shape") == std::string::npos)
  //         return;
  //       auto newhist_u = s->ClonedShapeU();
  //       TH1* h_up = (TH1*)newhist_u->Clone();
  //       if (ch::HasNegativeBins(s->shape_u())) {
  //         for (int i = 1; i <= newhist_u->GetNbinsX(); ++i) {
  //           if (newhist_u->GetBinContent(i) < 0.) {
  //             h_up->SetBinContent(i, 0.);
  //           }
  //         }
  //       }
  //       auto newhist_d = s->ClonedShapeD();
  //       TH1* h_down = (TH1*)newhist_d->Clone();
  //       if (ch::HasNegativeBins(s->shape_d())) {
  //         for (int i = 1; i <= newhist_d->GetNbinsX(); ++i) {
  //           if (newhist_d->GetBinContent(i) < 0.) {
  //             h_down->SetBinContent(i, 0.);
  //           }
  //         }
  //       }
  //       const TH1& up = *h_up;
  //       const TH1& down = *h_down;
  //       auto nominal = p->ClonedShape();
  //       TH1* nom = (TH1*)nominal->Clone();
  //       const TH1& nom_hist = *nom;
  //       s->set_shapes(up, down, nom_hist);
  //     }
  //   });
  // });

  // Delete processes with 0 yield
  cb.FilterProcs([&](ch::Process *p) {
    bool null_yield = (!(p->rate() > 0.0001) && !(p->signal()));
    if (null_yield) {
      std::cout << "[WARNING] Removing background process with null yield: \n ";
      // std::cout << ch::Process::PrintHeader << *p << "\n";
      cb.FilterSysts([&](ch::Systematic *s) {
        bool remove_syst = (MatchingProcess(*p, *s));
        return remove_syst;
      });
    }
    return null_yield;
  });
  cb.FilterProcs([&](ch::Process *p) {
    bool null_yield = (!(p->rate() > 0.0) && (p->signal()));
    if (null_yield) {
      std::cout << "[WARNING] Removing signal process with null yield: \n ";
      // std::cout << ch::Process::PrintHeader << *p << "\n";
      cb.FilterSysts([&](ch::Systematic *s) {
        bool remove_syst = (MatchingProcess(*p, *s));
        return remove_syst;
      });
    }
    return null_yield;
  });

  // Delete systematics with 0 yield since these result in a bogus norm error in combine
  cb.FilterSysts([&](ch::Systematic *s) {
    if (s->type() == "shape") {
      if (s->shape_u()->Integral() <= 0.001 || std::isnan(s->shape_u()->Integral())) {
        std::cout << "[WARNING] Removing systematic with null yield in up shift:" << std::endl;
        // std::cout << ch::Systematic::PrintHeader << *s << "\n";
        return true;
      }
      if (s->shape_d()->Integral() <= 0.001 || std::isnan(s->shape_d()->Integral())) {
        std::cout << "[WARNING] Removing systematic with null yield in down shift:" << std::endl;
        // std::cout << ch::Systematic::PrintHeader << *s << "\n";
        return true;
      }
    }
    return false;
  });

  // Perform auto-rebinning
  if (auto_rebin) {
    const auto bin_threshold = 0.;
    const auto threshold = 0.001;
    const auto sig_threshold = 2.0;

    // for (auto b : cb.cp().bin_set()) {
    //   std::cout << "[INFO] Prebin bin " << b << "\n";
    //   // Get shape of this category with sum of backgrounds
    //   auto shape = cb.cp().bin({b}).backgrounds().GetShape();
    //   // Push back last bin edge
    //   vector<double> binning;
    //   const auto num_bins = shape.GetNbinsX();
    //   binning.push_back(shape.GetBinLowEdge(num_bins + 1));
    //   // Now, go backwards through bins (from right to left) and merge a bin if
    //   // the background yield is below a given threshold.
    //   auto c = 0.0;
    //   for(auto i = num_bins; i > 0; i--) {
    //     // std::cout << "[INFO] Bin " << i << " in bin " << b << " with yield " << shape.GetBinContent(i) << "\n";

    //     auto low_edge = shape.GetBinLowEdge(i);
    //     c += shape.GetBinContent(i);
    //     if ((i == 6) || (i == 1)) { // Set lower edge if the bin content is above the threshold.
    //       binning.insert(binning.begin(), low_edge);
    //       c = 0.0;
    //     }
    //   }
    //   cb.cp().bin({b}).VariableRebin(binning);
    // }
    for (auto b : cb.cp().bin_set()) {
      std::cout << "[INFO] Rebin bin " << b << "\n";
      // Get shape of this category with sum of backgrounds
      auto shape = cb.cp().bin({b}).backgrounds().GetShape();
      if(std::isnan(shape.Integral())){
        std::cout << "[WARNING]0 category " << b << ", Bkg yield:" << shape.Integral() << "\n";
        
      }
      // Push back last bin edge
      vector<double> binning;
      const auto num_bins = shape.GetNbinsX();
      binning.push_back(shape.GetBinLowEdge(num_bins + 1));
      // Now, go backwards through bins (from right to left) and merge a bin if
      // the background yield is below a given threshold.
      auto c = 0.0;
      for(auto i = num_bins; i > 0; i--) {
        std::cout << "[INFO] Bin " << i << " in bin " << b << " with yield " << shape.GetBinContent(i) << "\n";

        auto low_edge = shape.GetBinLowEdge(i);
        c += shape.GetBinContent(i);
        if (c > bin_threshold) { // Set lower edge if the bin content is above the threshold.
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
      auto shape = cb.cp().bin({b}).backgrounds().GetShape();
      const auto num_bins = shape.GetNbinsX();
      for(auto i = num_bins; i > 0; i--) {
        std::cout << "[INFO] Bin " << i << " in bin " << b << " with yield " << shape.GetBinContent(i) << " rate " << cb.cp().bin({b}).backgrounds().GetRate() << "\n";
      }
      // Get yield of all backgrounds in this category
      auto shape_integral_bkg = cb.cp().bin({b}).backgrounds().GetShape().Integral();
      auto shape_integral_sig = cb.cp().bin({b}).signals().GetShape().Integral();
      auto data_rate = cb.cp().bin({b}).GetObservedRate();
      // std::cout << "[WARNING] Category " << b << " has insufficient population! " << data_rate << "\n";
      if(std::isnan(shape_integral_bkg)){
        std::cout << "[WARNING]1 category " << b << ", Bkg yield: " << shape_integral_bkg << ", signal yield: " << shape_integral_sig << ", data yield: " << data_rate << "\n";

      }
      
      if((shape_integral_bkg < threshold) && remove_empty_categories){
        std::cout << "[WARNING] Remove category " << b << " due to insufficient population!" << "\n";
        std::cout << "[WARNING] Bkg yield: " << shape_integral_bkg << ", signal yield: " << shape_integral_sig << ", data yield: " << data_rate << "\n";
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

  // Merge bins and set bin-by-bin uncertainties if no autoMCStats is used.
  if (classic_bbb) {
    auto bbb = ch::BinByBinFactory()
                   .SetAddThreshold(0.0)
                   .SetMergeThreshold(0.5)
                   .SetFixNorm(false);
    bbb.MergeBinErrors(cb.cp().backgrounds());
    bbb.AddBinByBin(cb.cp().backgrounds(), cb);
  }
  if (binomial_bbb) {
    // Used for statistical fluctuation in embedded weights in em channel
    auto gen_mean = 0.0;
    if (era==2016){
      gen_mean = 0.017;
    }
    else if (era==2017){
      gen_mean = 0.014;
    }
    else if (era==2018){
      gen_mean = 0.019;
    }
    auto bbb = ch::BinomialBinByBinViaAutoMCstatsFactory()
                   .SetPattern("CMS_$ANALYSIS_$CHANNEL_$BIN_$ERA_$PROCESS_binomial_bin_$#")
                   .SetBinomialP(gen_mean)
                   .SetBinomialN(1000.0)
                   .SetFixNorm(false);
    bbb.AddBinomialBinByBin(cb.cp().channel({"em"}).process({"EMB"}), cb);
  }

  if (use_automc) {
    std::cout << "[INFO] Adding SetAutoMCStats .\n";
    cb.SetAutoMCStats(cb, 10.);
  }


  // This function modifies every entry to have a standardised bin name of
  // the form: {analysis}_{channel}_{bin_id}_{era}
  ch::SetStandardBinNames(cb, "$ANALYSIS_$CHANNEL_$BINID_$ERA");

  // Write out datacards. Naming convention important for rest of workflow. We
  // make one directory per chn-cat, one per chn and cmb. In this code we only
  // store the individual datacards for each directory to be combined later.
  string boosted_tag = "";
  if (boosted_tt && categories != "gof")
    boosted_tag = "boosted_";
  
  if (categories == "gof") {
    ch::CardWriter writer(output_folder + "/$TAG/" + heavy_mass + "_" + light_mass + "_" + gof_category_name + "/" + boosted_tag + "$BIN.txt",
      output_folder +"/$TAG/common" + "_" + gof_category_name + "/nmssm_input_" + boosted_tag + era_tag  + "_" + heavy_mass + "_" + light_mass + ".root");
    
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
  } 
  else {
    ch::CardWriter writer(output_folder + "/$TAG/" + heavy_mass + "_" + light_mass + "/" + boosted_tag + "$BIN.txt",
    output_folder +"/$TAG/common/nmssm_input_" + boosted_tag + era_tag + "_" + chan + "_" + heavy_mass + "_" + light_mass + ".root");
  
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
  }

  // if (verbose)
  //   cb.PrintAll();

  cout << "[INFO] Done producing datacards.\n";
}
