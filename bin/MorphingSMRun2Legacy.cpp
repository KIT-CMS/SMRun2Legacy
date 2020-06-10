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
  string input_folder_em = "Vienna/";
  string input_folder_et = "Vienna/";
  string input_folder_mt = "Vienna/";
  string input_folder_tt = "Vienna/";
  string chan = "all";
  string postfix = "-ML";
  string midfix = "";
  bool regional_jec = true;
  bool ggh_wg1 = true;
  bool qqh_wg1 = true;
  bool auto_rebin = false;
  bool rebin_categories = true;
  bool manual_rebin_for_yields = false;
  bool real_data = false;
  bool jetfakes = true;
  bool train_ff = true;
  bool train_emb = true;
  bool train_stage0 = false;
  bool embedding = false;
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
      ("input_folder_em", po::value<string>(&input_folder_em)->default_value(input_folder_em))
      ("input_folder_et", po::value<string>(&input_folder_et)->default_value(input_folder_et))
      ("input_folder_mt", po::value<string>(&input_folder_mt)->default_value(input_folder_mt))
      ("input_folder_tt", po::value<string>(&input_folder_tt)->default_value(input_folder_tt))
      ("postfix", po::value<string>(&postfix)->default_value(postfix))
      ("midfix", po::value<string>(&midfix)->default_value(midfix))
      ("channel", po::value<string>(&chan)->default_value(chan))
      ("auto_rebin", po::value<bool>(&auto_rebin)->default_value(auto_rebin))
      ("rebin_categories", po::value<bool>(&rebin_categories)->default_value(rebin_categories))
      ("manual_rebin_for_yields", po::value<bool>(&manual_rebin_for_yields)->default_value(manual_rebin_for_yields))
      ("regional_jec", po::value<bool>(&regional_jec)->default_value(regional_jec))
      ("ggh_wg1", po::value<bool>(&ggh_wg1)->default_value(ggh_wg1))
      ("qqh_wg1", po::value<bool>(&qqh_wg1)->default_value(qqh_wg1))
      ("real_data", po::value<bool>(&real_data)->default_value(real_data))
      ("verbose", po::value<bool>(&verbose)->default_value(verbose))
      ("remove_empty_categories", po::value<bool>(&remove_empty_categories)->default_value(remove_empty_categories))
      ("output_folder", po::value<string>(&output_folder)->default_value(output_folder))
      ("stxs_signals", po::value<string>(&stxs_signals)->default_value(stxs_signals))
      ("categories", po::value<string>(&categories)->default_value(categories))
      ("gof_category_name", po::value<string>(&gof_category_name)->default_value(gof_category_name))
      ("jetfakes", po::value<bool>(&jetfakes)->default_value(jetfakes))
      ("embedding", po::value<bool>(&embedding)->default_value(embedding))
      ("train_ff", po::value<bool>(&train_ff)->default_value(train_ff))
      ("train_emb", po::value<bool>(&train_emb)->default_value(train_emb))
      ("train_stage0", po::value<bool>(&train_stage0)->default_value(train_stage0))
      ("classic_bbb", po::value<bool>(&classic_bbb)->default_value(classic_bbb))
      ("binomial_bbb", po::value<bool>(&binomial_bbb)->default_value(binomial_bbb))
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
  if (chan.find("em") != std::string::npos)
    chns.push_back("em");
  if (chan == "all")
    chns = {"mt", "et", "tt", "em"};

  // Define background processes
  map<string, VString> bkg_procs;
  VString bkgs, bkgs_em;
  bkgs = {"W", "ZTT", "QCD", "ZL", "ZJ", "TTT", "TTL", "TTJ", "VVJ", "VVT", "VVL", "WH125", "ZH125", "ttH125", "ggHWW125", "qqHWW125", "WHWW125", "ZHWW125"};
  //bkgs = {"W", "ZTT", "QCD", "ZL", "ZJ", "TTT", "TTL", "TTJ", "VVJ", "VVT", "VVL", "WH125", "ZH125", "ggHWW125", "qqHWW125", "WHWW125"};
  bkgs_em = {"W", "ZTT", "TTT","VVT", "QCD", "ZL", "TTL", "VVL", "WH125", "ZH125", "ttH125", "ggHWW125", "qqHWW125", "WHWW125", "ZHWW125"};


  if(embedding){
    bkgs.erase(std::remove(bkgs.begin(), bkgs.end(), "ZTT"), bkgs.end());
    bkgs.erase(std::remove(bkgs.begin(), bkgs.end(), "TTT"), bkgs.end());
    bkgs.erase(std::remove(bkgs.begin(), bkgs.end(), "VVT"), bkgs.end());
    bkgs = JoinStr({bkgs,{"EMB"}});
    bkgs_em.erase(std::remove(bkgs_em.begin(), bkgs_em.end(), "ZTT"), bkgs_em.end());
    bkgs_em.erase(std::remove(bkgs_em.begin(), bkgs_em.end(), "TTT"), bkgs_em.end());
    bkgs_em.erase(std::remove(bkgs_em.begin(), bkgs_em.end(), "VVT"), bkgs_em.end());
    bkgs_em = JoinStr({bkgs_em,{"EMB"}});
  }
  if(jetfakes){
    bkgs.erase(std::remove(bkgs.begin(), bkgs.end(), "QCD"), bkgs.end());
    bkgs.erase(std::remove(bkgs.begin(), bkgs.end(), "W"), bkgs.end());
    bkgs.erase(std::remove(bkgs.begin(), bkgs.end(), "VVJ"), bkgs.end());
    bkgs.erase(std::remove(bkgs.begin(), bkgs.end(), "TTJ"), bkgs.end());
    bkgs.erase(std::remove(bkgs.begin(), bkgs.end(), "ZJ"), bkgs.end());
    bkgs = JoinStr({bkgs,{"jetFakes"}});
  }

  std::cout << "[INFO] Considerung the following processes:\n";
  if (chan.find("em") != std::string::npos) {
    std::cout << "For em channel : \n";
    for (unsigned int i=0; i < bkgs_em.size(); i++) std::cout << bkgs_em[i] << std::endl;
  }
  if (chan.find("mt") != std::string::npos || chan.find("et") != std::string::npos || chan.find("tt") != std::string::npos) {
    std::cout << "For et,mt,tt channels : \n";
    for (unsigned int i=0; i < bkgs.size(); i++) std::cout << bkgs[i] << std::endl;}
  bkg_procs["et"] = bkgs;
  bkg_procs["mt"] = bkgs;
  bkg_procs["tt"] = bkgs;
  bkg_procs["em"] = bkgs_em;

  // Define categories
  map<string, Categories> cats;
  std::vector<std::string> cats_to_keep; // will be used later for the card writer
  for (auto chn : chns){
    //define mapping for signal categories if categories is stage0 /1p1
    if(categories == "stxs_stage0"){
      cats[chn]={
        //{ 1, chn+"_ggh"},
        //{ 2, chn+"_qqh"},
        //use 2D category instead:
        { 1, chn+"_xxh"},
    };
    } else if (categories == "inclusive"){
    cats[chn]={
      { 1, chn+"_xxh"},
    };
    } else if (categories == "stxs_stage1p1") {
        if(train_stage0){
          cats[chn]={
            {100, chn+"_ggh_100"},
            {101, chn+"_ggh_101"},
            {102, chn+"_ggh_102"},
            {103, chn+"_ggh_103"},
            {200, chn+"_qqh_200"},
            {201, chn+"_qqh_201"},
            {202, chn+"_qqh_202"},
            {203, chn+"_qqh_203"},
          };
        }
        else {
          cats[chn]={
            {100, chn+"_ggh_0J_PTH_0_10"},
            {101, chn+"_ggh_0J_PTH_GT10"},
            {102, chn+"_ggh_1J_PTH0to60"},
            {103, chn+"_ggh_1J_PTH60to120"},
            {104, chn+"_ggh_1J_PTH120to200"},
            {105, chn+"_ggh_2J_PTH0to60"},
            {106, chn+"_ggh_2J_PTH60to120"},
            {107, chn+"_ggh_2J_PTH120to200"},
            {108, chn+"_ggh_PTH_200_300"},
            {109, chn+"_ggh_PTHGT300"},
            {110, chn+"_ggh_vbftopo"},
            {200, chn+"_qqh_2J"},
            {201, chn+"_qqh_PTHGT200"},
            //{202, chn+"_vbftopo_highmjj"},
            //{203, chn+"_vbftopo_lowmjj"},
            {202, chn+"_qqh_vbftopo_highmjj"},
            {203, chn+"_qqh_vbftopo_lowmjj"},
          };
        }
    }

    //define mapping for background categories if categories ist stage0 /1p1
    //11 :w , 12: ztt, 13: tt, 14: ss, 15: zll, 16: misc, 17:noniso, 18: st?? 19:db, 20:emb, 21: ff
    if (categories == "stxs_stage0" || categories == "stxs_stage1p1" || categories=="inclusive"){
      cats[chn].push_back({16,chn+"_misc"});

      if(chn=="em") cats[chn].push_back({19, chn+"_db"});
      if(chn!="tt") cats[chn].push_back({13, chn+"_tt"});
      if(chn=="mt" || chn=="et") cats[chn].push_back({15, chn+"_zll"});

      if (train_emb) cats[chn].push_back({20, chn+"_emb"});
      else cats[chn].push_back({12, chn+"_ztt"});

      if (train_ff) {
        if(chn=="em") cats[chn].push_back({14, chn+"_ss"}); // switch to ff does nothing for em so this is added both times
        else cats[chn].push_back({21, chn+"_ff"});
      } else {
        if(chn=="mt" || chn=="et") cats[chn].push_back({11, chn+"_w"});
        if(chn=="tt") cats[chn].push_back({17, chn+"_noniso"});
        else cats[chn].push_back({14, chn+"_ss"});
      }
    }
    else if(categories == "gof") cats[chn] = { {300, gof_category_name.c_str() }};
    else throw std::runtime_error("Given categorization is not known.");
  }
  for (auto chn : chns){
    for (auto tuple: cats[chn]) cout << tuple.first << ": " << tuple.second << endl;
  }

  // Specify signal processes and masses
  vector<string> sig_procs;
  // STXS stage 0: ggH and VBF processes
  if(stxs_signals == "stxs_stage0") sig_procs = {"ggH", "qqH"};
  // STXS stage 1.1: Splits of ggH and VBF processes
  // References:
  // - https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHXSWGFiducialAndSTXS
  // - https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHXSWG2
  // - https://gitlab.cern.ch/LHCHIGGSXS/LHCHXSWG2/STXS
  else if(stxs_signals == "stxs_stage1p1") sig_procs = {
      // ggH
      "ggH_GG2H_FWDH",
      "ggH_GG2H_PTH_200_300",
      "ggH_GG2H_PTH_300_450",
      "ggH_GG2H_PTH_450_650",
      "ggH_GG2H_PTH_GT650",
      "ggH_GG2H_0J_PTH_0_10",
      "ggH_GG2H_0J_PTH_GT10",
      "ggH_GG2H_1J_PTH_0_60",
      "ggH_GG2H_1J_PTH_60_120",
      "ggH_GG2H_1J_PTH_120_200",
      "ggH_GG2H_GE2J_MJJ_0_350_PTH_0_60",
      "ggH_GG2H_GE2J_MJJ_0_350_PTH_60_120",
      "ggH_GG2H_GE2J_MJJ_0_350_PTH_120_200",
      "ggH_GG2H_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_0_25",
      "ggH_GG2H_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_GT25",
      "ggH_GG2H_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_0_25",
      "ggH_GG2H_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_GT25",
      // VBF
      "qqH_QQ2HQQ_FWDH",
      "qqH_QQ2HQQ_0J",
      "qqH_QQ2HQQ_1J",
      "qqH_QQ2HQQ_GE2J_MJJ_0_60",
      "qqH_QQ2HQQ_GE2J_MJJ_60_120",
      "qqH_QQ2HQQ_GE2J_MJJ_120_350",
      "qqH_QQ2HQQ_GE2J_MJJ_GT350_PTH_GT200",
      "qqH_QQ2HQQ_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_0_25",
      "qqH_QQ2HQQ_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_GT25",
      "qqH_QQ2HQQ_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_0_25",
      "qqH_QQ2HQQ_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_GT25"
  };

  else throw std::runtime_error("Given STXS signals are not known.");
  vector<string> masses = {"125"};

  // Create combine harverster object
  ch::CombineHarvester cb;

  // Add observations and processes
  std::string era_tag;
  if (era == 2016) era_tag = "2016";
  else if (era == 2017) era_tag = "2017";
  else if (era == 2018) era_tag = "2018";

  else std::runtime_error("Given era is not implemented.");

  for (auto chn : chns) {
    cb.AddObservations({"*"}, {"htt"}, {era_tag}, {chn}, cats[chn]);
    cb.AddProcesses({"*"}, {"htt"}, {era_tag}, {chn}, bkg_procs[chn], cats[chn],
                    false);
    cb.AddProcesses(masses, {"htt"}, {era_tag}, {chn}, sig_procs, cats[chn],
                    true);
  }

  // Add systematics
  ch::AddSMRun2Systematics(cb, jetfakes, embedding, regional_jec, ggh_wg1, qqh_wg1, era);

  // Define the location of the "auxiliaries" directory where we can
  // source the input files containing the datacard shapes
  std::map<string, string> input_dir;
  input_dir["mt"] = base_path + "/" + input_folder_mt + "/";
  input_dir["et"] = base_path + "/" + input_folder_et + "/";
  input_dir["tt"] = base_path + "/" + input_folder_tt + "/";
  input_dir["em"] = base_path + "/" + input_folder_em + "/";
  // Extract shapes from input ROOT files
  for (string chn : chns) {
    cb.cp().channel({chn}).backgrounds().ExtractShapes(
        input_dir[chn] + std::to_string(era) + midfix + chn + "-synced"+postfix+ ".root",
        "$BIN/$PROCESS", "$BIN/$PROCESS_$SYSTEMATIC");
    cb.cp().channel({chn}).process(sig_procs).ExtractShapes(
        input_dir[chn] + std::to_string(era) + midfix + chn + "-synced"+postfix+ ".root",
        "$BIN/$PROCESS$MASS", "$BIN/$PROCESS$MASS_$SYSTEMATIC");
  }

  // Delete processes with 0 yield
  cb.FilterProcs([&](ch::Process *p) {
    bool null_yield = !(p->rate() > 0.0);
    if (null_yield) {
      std::cout << "[WARNING] Removing process with null yield: \n ";
      std::cout << ch::Process::PrintHeader << *p << "\n";
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
      if (s->shape_u()->Integral() == 0.0) {
        std::cout << "[WARNING] Removing systematic with null yield in up shift:" << std::endl;
        std::cout << ch::Systematic::PrintHeader << *s << "\n";
        return true;
      }
      if (s->shape_d()->Integral() == 0.0) {
        std::cout << "[WARNING] Removing systematic with null yield in down shift:" << std::endl;
        std::cout << ch::Systematic::PrintHeader << *s << "\n";
        return true;
      }
    }
    return false;
  });

  int count_lnN = 0;
  int count_all = 0;
  cb.cp().ForEachSyst([&count_lnN, &count_all](ch::Systematic *s) {
    if (TString(s->name()).Contains("scale")||TString(s->name()).Contains("CMS_htt_boson_reso_met")||TString(s->name()).Contains("res_j")||TString(s->name()).Contains("res_e")){
      count_all++;
      double err_u = 0.0;
      double err_d = 0.0;
      int nbins = s->shape_u()->GetNbinsX();
      double yield_u = s->shape_u()->IntegralAndError(1,nbins,err_u);
      double yield_d = s->shape_d()->IntegralAndError(1,nbins,err_d);
      double value_u = s->value_u();
      double value_d = s->value_d();
      if (std::abs(value_u-1.0)+std::abs(value_d-1.0)<err_u/yield_u+err_d/yield_d){
          count_lnN++;
          std::cout << "[WARNING] Replacing systematic by lnN:" << std::endl;
          std::cout << ch::Systematic::PrintHeader << *s << "\n";
          s->set_type("lnN");
          bool up_is_larger = (value_u>value_d);
          if (value_u < 1.0) value_u = 1.0 / value_u;
          if (value_d < 1.0) value_d = 1.0 / value_d;
          if (up_is_larger){
              value_u = std::sqrt(value_u*value_d);
              value_d = 1.0 / value_u;
          }else{
              value_d = std::sqrt(value_u*value_d);
              value_u = 1.0 / value_d;
          }
          std::cout << "Former relative integral up shift: " << s->value_u() << "; New relative integral up shift: " << value_u << std::endl;
          std::cout << "Former relative integral down shift: " << s->value_d() << "; New relative integral down shift: " << value_d << std::endl;
          s->set_value_u(value_u);
          s->set_value_d(value_d);
      }
    }
  });
  std::cout << "[WARNING] Turned " << count_lnN << " of " << count_all << " checked systematics into lnN:" << std::endl;

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

  // Rebin categories to predefined binning for binning
  if (rebin_categories) {
    // Rebin background categories
    for (auto b : cb.cp().bin_set()) {
      TString bstr = b;
      if (bstr.Contains("ggh") || bstr.Contains("qqh") || bstr.Contains("vbftopo") || bstr.Contains("xxh")) continue;
      std::cout << "[INFO] Rebin background bin " << b << "\n";
      auto shape = cb.cp().bin({b}).backgrounds().GetShape();
      auto min = shape.GetBinLowEdge(1);
      if(bstr.Contains("em") && bstr.Contains("misc")) cb.cp().bin({b}).VariableRebin({min, 0.4, 1.0});
      else if(bstr.Contains("em_emb")){
        if(categories == "stxs_stage1p1") cb.cp().bin({b}).VariableRebin({min, 0.3, 1.0});
        else cb.cp().bin({b}).VariableRebin({min, 0.4, 0.5, 0.6, 1.0});
      }
      else if(bstr.Contains("et") && bstr.Contains("misc")) cb.cp().bin({b}).VariableRebin({min, 0.4, 0.5, 0.6, 1.0});
      else if(bstr.Contains("mt") && bstr.Contains("misc")) cb.cp().bin({b}).VariableRebin({min, 0.4, 0.5, 0.6, 1.0});
      else if(bstr.Contains("et") && bstr.Contains("zll") && categories == "stxs_stage1p1") cb.cp().bin({b}).VariableRebin({min, 0.4, 0.5, 1.0});
      else if(bstr.Contains("mt") && bstr.Contains("emb")){
        if(categories == "stxs_stage1p1") cb.cp().bin({b}).VariableRebin({min, 0.4, 0.45, 0.5, 0.6, 1.0});
        else cb.cp().bin({b}).VariableRebin({min, 0.4, 0.5, 0.6, 1.0});
      }
      else cb.cp().bin({b}).VariableRebin({min, 0.4, 0.5, 0.6, 0.7, 1.0});
    }
    // Rebin ggh stage 1.1 categories
    for (auto b : cb.cp().bin_set()) {
      TString bstr = b;
      if (bstr.Contains("ggh_10")) {
        std::cout << "[INFO] Rebin ggh signal bin " << b << "\n";
        auto shape = cb.cp().bin({b}).backgrounds().GetShape();
        auto min = shape.GetBinLowEdge(1);
        cb.cp().bin({b}).VariableRebin({min, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.7, 1.0});
      }
    }
    // Rebin qqh stage 1.1 categories
    for (auto b : cb.cp().bin_set()) {
      TString bstr = b;
      if (bstr.Contains("qqh_20")) {
        std::cout << "[INFO] Rebin qqh signal bin " << b << "\n";
        auto shape = cb.cp().bin({b}).backgrounds().GetShape();
        auto min = shape.GetBinLowEdge(1);
        cb.cp().bin({b}).VariableRebin({min, 0.4, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.90, 0.95, 1.0});
      }
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
    auto bbb = ch::BinomialBinByBinFactory()
                   .SetPattern("CMS_$ANALYSIS_$CHANNEL_$BIN_$ERA_$PROCESS_binomial_bin_$#")
                   .SetBinomialP(gen_mean)
                   .SetBinomialN(1000.0)
                   .SetFixNorm(false);
    bbb.AddBinomialBinByBin(cb.cp().channel({"em"}).process({"EMB"}), cb);
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
