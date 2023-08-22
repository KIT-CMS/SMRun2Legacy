#include "CombineHarvester/SMRun2Legacy/interface/HttSystematics_SMRun2.h"
#include "CombineHarvester/CombineTools/interface/Process.h"
#include "CombineHarvester/CombineTools/interface/Systematics.h"
#include "CombineHarvester/CombineTools/interface/Utilities.h"
#include <string>
#include <vector>

using namespace std;

namespace ch {

using ch::syst::SystMap;
using ch::syst::SystMapAsymm;
using ch::syst::era;
using ch::syst::channel;
using ch::syst::bin_id;
using ch::syst::process;
using ch::syst::bin;
using ch::JoinStr;

  void AddSMRun2Systematics(CombineHarvester &cb, bool jetfakes, bool embedding, bool regional_jec, bool ggh_wg1, int era) {

  // ##########################################################################
  // Define groups of processes
  // ##########################################################################

  // Signal processes
      // ggH
      // VBF
  std::vector<std::string> signals_ggH = {
      // STXS stage 0
      "ggH_htt",
      // STXS stage 1.1
      "ggH_FWDH_htt",
      "ggH_PTH_200_300_htt",
      "ggH_PTH_300_450_htt",
      "ggH_PTH_450_650_htt",
      "ggH_PTH_GT650_htt",
      "ggH_0J_PTH_0_10_htt",
      "ggH_0J_PTH_GT10_htt",
      "ggH_1J_PTH_0_60_htt",
      "ggH_1J_PTH_60_120_htt",
      "ggH_1J_PTH_120_200_htt",
      "ggH_GE2J_MJJ_0_350_PTH_0_60_htt",
      "ggH_GE2J_MJJ_0_350_PTH_60_120_htt",
      "ggH_GE2J_MJJ_0_350_PTH_120_200_htt",
      "ggH_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_0_25_htt",
      "ggH_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_GT25_htt",
      "ggH_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_0_25_htt",
      "ggH_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_GT25_htt",
      };
  std::vector<std::string> signals_ggZH_had = {
      // STXS stage 0
      "ggZH_had_htt",
      // STXS stage 1.1
      "ggZH_had_FWDH_htt",
      "ggZH_had_PTH_200_300_htt",
      "ggZH_had_PTH_300_450_htt",
      "ggZH_had_PTH_450_650_htt",
      "ggZH_had_PTH_GT650_htt",
      "ggZH_had_0J_PTH_0_10_htt",
      "ggZH_had_0J_PTH_GT10_htt",
      "ggZH_had_1J_PTH_0_60_htt",
      "ggZH_had_1J_PTH_60_120_htt",
      "ggZH_had_1J_PTH_120_200_htt",
      "ggZH_had_GE2J_MJJ_0_350_PTH_0_60_htt",
      "ggZH_had_GE2J_MJJ_0_350_PTH_60_120_htt",
      "ggZH_had_GE2J_MJJ_0_350_PTH_120_200_htt",
      "ggZH_had_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_0_25_htt",
      "ggZH_had_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_GT25_htt",
      "ggZH_had_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_0_25_htt",
      "ggZH_had_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_GT25_htt",
      };
  std::vector<std::string> signals_qqH = {
      // STXS stage 0
      "qqH_htt",
      // STXS stage 1
      "qqH_FWDH_htt",
      "qqH_0J_htt",
      "qqH_1J_htt",
      "qqH_GE2J_MJJ_0_60_htt",
      "qqH_GE2J_MJJ_60_120_htt",
      "qqH_GE2J_MJJ_120_350_htt",
      "qqH_GE2J_MJJ_GT350_PTH_GT200_htt",
      "qqH_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_0_25_htt",
      "qqH_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_GT25_htt",
      "qqH_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_0_25_htt",
      "qqH_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_GT25_htt"
      };
  std::vector<std::string> signals_VH_had = {
      // STXS stage 0
      "WH_had_htt",
      "ZH_had_htt",
      // STXS stage 1
      "WH_had_FWDH_htt",
      "WH_had_0J_htt",
      "WH_had_1J_htt",
      "WH_had_GE2J_MJJ_0_60_htt",
      "WH_had_GE2J_MJJ_60_120_htt",
      "WH_had_GE2J_MJJ_120_350_htt",
      "WH_had_GE2J_MJJ_GT350_PTH_GT200_htt",
      "WH_had_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_0_25_htt",
      "WH_had_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_GT25_htt",
      "WH_had_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_0_25_htt",
      "WH_had_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_GT25_htt",
      "ZH_had_FWDH_htt",
      "ZH_had_0J_htt",
      "ZH_had_1J_htt",
      "ZH_had_GE2J_MJJ_0_60_htt",
      "ZH_had_GE2J_MJJ_60_120_htt",
      "ZH_had_GE2J_MJJ_120_350_htt",
      "ZH_had_GE2J_MJJ_GT350_PTH_GT200_htt",
      "ZH_had_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_0_25_htt",
      "ZH_had_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_GT25_htt",
      "ZH_had_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_0_25_htt",
      "ZH_had_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_GT25_htt"
      };
  std::vector<std::string> signals_VH = {
      // STXS stage 0
      "WH_lep_htt", "ZH_lep_htt", "ggZH_lep_htt", "ttH_htt",
      // STXS stage 1
      "WH_lep_FWDH_htt",
      "WH_lep_PTV_0_75_htt",
      "WH_lep_PTV_75_150_htt",
      "WH_lep_PTV_150_250_0J_htt",
      "WH_lep_PTV_150_250_GE1J_htt",
      "WH_lep_PTV_GT250_htt",
      "ZH_lep_FWDH_htt",
      "ZH_lep_PTV_0_75_htt",
      "ZH_lep_PTV_75_150_htt",
      "ZH_lep_PTV_150_250_0J_htt",
      "ZH_lep_PTV_150_250_GE1J_htt",
      "ZH_lep_PTV_GT250_htt",
      "ggZH_lep_FWDH_htt",
      "ggZH_lep_PTV_0_75_htt",
      "ggZH_lep_PTV_75_150_htt",
      "ggZH_lep_PTV_150_250_0J_htt",
      "ggZH_lep_PTV_150_250_GE1J_htt",
      "ggZH_lep_PTV_GT250_htt"
      };
  std::vector<std::string> signals_ggHToWW = {
     // STXS stage 0
     "ggH_hww"};
  std::vector<std::string> signals_qqHToWW = {
     // STXS stage 0
     "qqH_hww"};
  std::vector<std::string> signals = JoinStr({signals_ggH, signals_ggZH_had, signals_qqH, signals_VH_had, signals_VH}); 
  std::vector<std::string> signals_WH = {"WHplus", "WHminus"};
  // Background processes
  /* // Not used in the function, keep it for documentation purposes.
  std::vector<std::string> backgrounds = {"ZTT",  "W",   "ZL",      "ZJ",
                                          "TTT",  "TTJ", "VVT",     "VVJ",
                                          "EWKZ", "QCD", "jetFakes", "EMB", "TTL"};
  */

  // All processes being taken from simulation
  // FIXME: Adapt for fake factor and embedding
  std::vector<std::string> mc_processes =
      JoinStr({
              signals,
              signals_ggHToWW,
              signals_qqHToWW,
              signals_WH,
              {"WH_hww", "ZH_hww"},
              {"ggZZ", "rem_VH", "WWW", "rem_VV", "WWZ", "ZZZ", "rem_ttbar", "WZ", "Wjets", "ZH", "DY", "ZZ", "TT"}
              });
  // ##########################################################################
  // Uncertainty: Lumi
  // References:
  // - "CMS Luminosity Measurements for the 2016 Data Taking Period"
  //   (PAS, https://cds.cern.ch/record/2257069)
  // - Recommendation twiki
  //    https://twiki.cern.ch/twiki/bin/view/CMS/TWikiLUM#LumiComb  
  // Notes:
  // - FIXME: Adapt for fake factor and embedding
  // ##########################################################################

  float lumi_unc = 1.0;
  float lumi_unc_corr = 1.0;
  float lumi_unc_1718 = 1.0;
  if (era == 2016) {
      lumi_unc = 1.010;
      lumi_unc_corr = 1.006;
  } else if (era == 2017) {
      lumi_unc = 1.020;
      lumi_unc_corr = 1.009;
      lumi_unc_1718 = 1.006;
  } else if (era == 2018) {
      lumi_unc = 1.015;
      lumi_unc_corr = 1.020;
      lumi_unc_1718 = 1.002;
  }
  cb.cp()
      .channel({"emt", "met", "mmt", "mtt", "ett"})
      .process(mc_processes)
      .AddSyst(cb, "lumi_13TeV_$ERA", "lnN", SystMap<>::init(lumi_unc));
  cb.cp()
      .channel({"emt", "met", "mmt", "mtt", "ett"})
      .process(mc_processes)
      .AddSyst(cb, "lumi_13TeV_correlated", "lnN", SystMap<>::init(lumi_unc_corr));
  cb.cp()
      .channel({"emt", "met", "mmt", "mtt", "ett"})
      .process(mc_processes)
      .AddSyst(cb, "lumi_13TeV_1718", "lnN", SystMap<>::init(lumi_unc_1718));

  // ##########################################################################
  // Uncertainty: Prefiring
  // References:
  // - "https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1ECALPrefiringWeightRecipe"
  // Notes:
  // - FIXME: assumed as uncorrelated accross the years for now, what is the recommendation?
  // ##########################################################################
  if (era != 2018) {
  cb.cp()
      .channel({"emt", "met", "mmt", "mtt", "ett"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_prefiring", "shape", SystMap<>::init(1.00));
  }

  // ##########################################################################
  // Uncertainty: Trigger efficiency
  // References:
  // Notes:
  // - FIXME: References?
  // ##########################################################################

  cb.cp()
    .channel({"emt", "met", "ett"})
    .process(mc_processes)
    .AddSyst(cb, "CMS_eff_trigger_et_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"emt", "met", "mtt", "mmt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_eff_trigger_mt_$ERA", "shape", SystMap<>::init(1.00));
          
  // ##########################################################################
  // Uncertainty: Electron, muon and tau ID efficiency
  // References:
  // Notes:
  // - FIXME: Adapt for fake factor and embedding
  // - FIXME: Handling of ZL in fully-hadronic channel?
  // - FIXME: References?
  // ##########################################################################

  // 3% in Tau ID SF with different anti-l fake WP
  cb.cp()
      .channel({"emt", "met", "mmt", "ett", "mtt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_eff_t_wp_$ERA", "lnN", SystMap<>::init(1.03));

  std::string tauIDptbins[4] = {"30-35", "35-40", "40-500", "500-1000"};
  std::string tauIDdmbins[4] = {"0", "1", "10", "11"};

  // Common component acting on MC
  
  // Electron ID
  cb.cp()
      .channel({"emt", "met", "ett"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_eff_e", "lnN", SystMap<>::init(1.02));

  // Muon ID
  cb.cp()
      .channel({"emt", "met", "mmt", "mtt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_eff_m", "lnN", SystMap<>::init(1.02));

  // Tau ID: et and mt with 1 real tau
  for (auto tauIDbin : tauIDptbins){ //first part correlated between channels for IDvsJets
    cb.cp()
        .channel({"emt", "met", "mmt"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_eff_t_"+tauIDbin+"_$ERA", "shape", SystMap<>::init(1.0));
  }
  cb.cp() //second part uncorrelated between channels for IDvsLep
      .channel({"emt", "met", "mmt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_eff_t_$CHANNEL_$ERA", "lnN", SystMap<>::init(1.01));

  // Tau ID: tt with 2 real taus
  for (auto tauIDbin : tauIDdmbins){
    cb.cp()
        .channel({"ett", "mtt"})
        .process(JoinStr({signals, {"ggZZ", "rem_VH", "WWW", "WWZ", "ZZZ", "rem_ttbar", "WZ", "ZZ"}}))
        .AddSyst(cb, "CMS_eff_t_dm"+tauIDbin+"_$ERA", "shape", SystMap<>::init(1.0));
  }
  cb.cp()
      .channel({ "mtt", "ett"})
      .process(JoinStr({signals, {"ggZZ", "rem_VH", "WWW", "WWZ", "ZZZ", "rem_ttbar", "WZ", "ZZ"}}))
      .AddSyst(cb, "CMS_eff_t_$CHANNEL_$ERA", "lnN", SystMap<>::init(1.014));

  // Tau ID: tt with 1 real taus and 1 jet fake
  cb.cp()
      .channel({ "mtt", "ett"})
      .process({"Wjets","DY", "TT", "rem_VV"})
      .AddSyst(cb, "CMS_eff_t_$ERA", "lnN", SystMap<>::init(1.06));


  cb.cp()
      .channel({"ett", "mtt"})
      .process({"Wjets","DY", "TT", "rem_VV"})
      .AddSyst(cb, "CMS_eff_t_$CHANNEL_$ERA", "lnN", SystMap<>::init(1.02));

  // ##########################################################################
  // Uncertainty: b-tag and mistag efficiency
  // References:
  // Notes:
  // - FIXME: References?
  // ##########################################################################

//   cb.cp()
//       .channel({"emt", "met", "mmt", "mtt", "ett"})
//       .process(mc_processes)
//       .AddSyst(cb, "CMS_htt_eff_b_$ERA", "shape", SystMap<>::init(1.00));

//   cb.cp()
//       .channel({"emt", "met", "mmt", "mtt", "ett"})
//       .process(mc_processes)
//       .AddSyst(cb, "CMS_htt_mistag_b_$ERA", "shape", SystMap<>::init(1.00));

  // ##########################################################################
  // Uncertainty: Electron energy scale
  // References:
  // - MC: https://twiki.cern.ch/twiki/bin/view/CMS/EgammaRunIIRecommendations#E_gamma_Energy_Corrections
  // - Embedding: ?
  // Notes:
  // - FIXME: References for embedding missing, need proper correlation accross years for mc, see here: https://twiki.cern.ch/twiki/bin/view/CMS/EgammaRunIIRecommendations#Recommendations_on_Combining_Sys
  // ##########################################################################

  // MC uncorrelated uncertainty

//   cb.cp()
//       .channel({"em", "et"})
//       .process(mc_processes)
//       .AddSyst(cb, "CMS_scale_e", "shape", SystMap<>::init(1.00));
//       //.AddSyst(cb, "CMS_scale_mc_e", "shape", SystMap<>::init(0.71));

//   cb.cp()
//       .channel({"em", "et"})
//       .process(mc_processes)
//       .AddSyst(cb, "CMS_res_e", "shape", SystMap<>::init(1.00));
//       //.AddSyst(cb, "CMS_scale_mc_e", "shape", SystMap<>::init(0.71));
      
  // Embedded uncorrelated uncertainty
      
  // ##########################################################################
  // Uncertainty: Tau energy scale
  // References:
  // Notes:
  // - Tau energy scale is split by decay mode.
  // - FIXME: References?
  // ##########################################################################


  // Common component acting on MC

  cb.cp()
      .channel({ "emt", "met", "mmt", "mtt", "ett"})
      .process(JoinStr({mc_processes}))
      .AddSyst(cb, "CMS_scale_t_1prong_$ERA", "shape", SystMap<>::init(1.0));

  cb.cp()
      .channel({ "emt", "met", "mmt", "mtt", "ett"})
      .process(JoinStr({mc_processes}))
      .AddSyst(cb, "CMS_scale_t_1prong1pizero_$ERA", "shape",
               SystMap<>::init(1.0));

  cb.cp()
      .channel({ "emt", "met", "mmt", "mtt", "ett"})
      .process(JoinStr({mc_processes}))
      .AddSyst(cb, "CMS_scale_t_3prong_$ERA", "shape", SystMap<>::init(1.0));

  cb.cp()
      .channel({ "emt", "met", "mmt", "mtt", "ett"})
      .process(JoinStr({mc_processes}))
      .AddSyst(cb, "CMS_scale_t_3prong1pizero_$ERA", "shape",
               SystMap<>::init(1.0));

  // ##########################################################################
  // Uncertainty: Jet energy scale
  // References:
  // - Talk in CMS Htt meeting by Daniel Winterbottom about regional JES splits:
  //   https://indico.cern.ch/event/740094/contributions/3055870/
  // Notes:
  // ##########################################################################

  if (!regional_jec) {
  cb.cp()
      .channel({"emt", "met", "mmt", "mtt", "ett"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_scale_j_$ERA", "shape", SystMap<>::init(0.71));
  cb.cp()
      .channel({"emt", "met", "mmt", "mtt", "ett"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_scale_j", "shape", SystMap<>::init(0.71));
  }

  // Regional JES
  else {
    // uncorrelated between eras
    cb.cp()
        .channel({"emt", "met", "mmt", "mtt", "ett"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_Absolute_$ERA", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"emt", "met", "mmt", "mtt", "ett"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_BBEC1_$ERA", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"emt", "met", "mmt", "mtt", "ett"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_EC2_$ERA", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"emt", "met", "mmt", "mtt", "ett"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_HF_$ERA", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"emt", "met", "mmt", "mtt", "ett"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_RelativeSample_$ERA", "shape", SystMap<>::init(1.00));
    // correlated between eras
    cb.cp()
        .channel({"emt", "met", "mmt", "mtt", "ett"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_Absolute", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"emt", "met", "mmt", "mtt", "ett"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_BBEC1", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"emt", "met", "mmt", "mtt", "ett"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_EC2", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"emt", "met", "mmt", "mtt", "ett"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_HF", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"emt", "met", "mmt", "mtt", "ett"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_FlavorQCD", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"emt", "met", "mmt", "mtt", "ett"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_RelativeBal", "shape", SystMap<>::init(1.00));
  }

  // JER
  if (era != 2017){
  cb.cp()
      .channel({"emt", "met", "mmt", "mtt", "ett"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_res_j_$ERA", "shape", SystMap<>::init(1.00));
  }else{
  std::vector<std::string> filtered_processes2;
  for (auto element : mc_processes){
      if (element!="ggH_0J_PTH_0_10") filtered_processes2.push_back(element);
  }
  cb.cp()
      .channel({"emt", "met", "mmt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_res_j_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({ "mtt", "ett"})
      .process(filtered_processes2)
      .AddSyst(cb, "CMS_res_j_$ERA", "shape", SystMap<>::init(1.00));
  }

  // ##########################################################################
  // Uncertainty: MET energy scale and Recoil
  // References:
  // Notes:
  // - FIXME: Clustered vs unclustered MET? Inclusion of JES splitting?
  // - FIXME: References?
  // ##########################################################################

  cb.cp()
      .channel({"emt", "met", "mmt", "mtt", "ett"})
      .process(mc_processes)  //Z and W processes are only included due to the EWK fraction. Make sure that there is no contribution to the shift from the DY or Wjets samples.
      .AddSyst(cb, "CMS_scale_met_unclustered_$ERA", "shape", SystMap<>::init(1.00));
    //TODO uncomment for next iteration
  cb.cp()
      .channel({"emt", "met", "mmt", "mtt", "ett"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_htt_boson_scale_met_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"emt", "met", "mmt", "mtt", "ett"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_htt_boson_res_met_$ERA", "shape", SystMap<>::init(1.00));

  // ##########################################################################
  // Uncertainty: Background normalizations
  // References:
  // Notes:
  // - FIXME: Remeasure QCD extrapolation factors for SS and ABCD methods?
  //          Current values are measured by KIT.
  // - FIXME: Adapt for fake factor and embedding
  // - FIXME: W uncertainties: Do we need lnN uncertainties based on the Ersatz
  //          study in Run1 (found in HIG-16043 uncertainty model)
  // - FIXME: References?
  // ##########################################################################

  // VV
  cb.cp()
      .channel({"emt", "met", "mmt", "mtt", "ett"})
      .process({"ggZZ", "WWW", "rem_VV", "WWZ", "ZZZ", "WZ", "ZZ" })
      .AddSyst(cb, "CMS_htt_vvXsec", "lnN", SystMap<>::init(1.05));

  // TT
  cb.cp()
      .channel({"emt", "met", "mmt", "mtt", "ett"})
      .process({"TTT", "TTL", "TTJ", "TT"})
      .AddSyst(cb, "CMS_htt_tjXsec", "lnN", SystMap<>::init(1.06));

  // W
  cb.cp()
      .channel({"emt", "met", "mmt", "mtt", "ett"})
      .process({"Wjets"})
      .AddSyst(cb, "CMS_htt_wjXsec", "lnN", SystMap<>::init(1.04));

  // Z
  cb.cp()
      .channel({"emt", "met", "mmt", "mtt", "ett"})
      .process({"ZTT", "ZL", "ZJ"})
      .AddSyst(cb, "CMS_htt_zjXsec", "lnN", SystMap<>::init(1.02));

  // QCD
  cb.cp()
      .channel({"emt", "met", "ett"})
      .process({"QCD"})
      .AddSyst(cb, "CMS_ExtrapSSOS_$CHANNEL_$ERA", "lnN", SystMap<>::init(1.05));
  cb.cp()
      .channel({"emt", "met", "mtt", "mmt"})
      .process({"QCD"})
      .AddSyst(cb, "CMS_ExtrapSSOS_$CHANNEL_$ERA", "lnN", SystMap<>::init(1.03));
  cb.cp()
      .channel({ "mtt", "ett"})
      .process({"QCD"})
      .AddSyst(cb, "CMS_ExtrapABCD_$CHANNEL_$ERA", "lnN", SystMap<>::init(1.03));

  cb.cp()
      .channel({"em"})
      .process({"QCD"})
      .AddSyst(cb, "CMS_htt_qcd_0jet_rate_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"em"})
      .process({"QCD"})
      .AddSyst(cb, "CMS_htt_qcd_0jet_shape_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"em"})
      .process({"QCD"})
      .AddSyst(cb, "CMS_htt_qcd_0jet_shape2_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"em"})
      .process({"QCD"})
      .AddSyst(cb, "CMS_htt_qcd_1jet_rate_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"em"})
      .process({"QCD"})
      .AddSyst(cb, "CMS_htt_qcd_1jet_shape_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"em"})
      .process({"QCD"})
      .AddSyst(cb, "CMS_htt_qcd_1jet_shape2_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
     .channel({"em"})
     .process({"QCD"})
     .AddSyst(cb, "CMS_htt_qcd_2jet_rate_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
     .channel({"em"})
     .process({"QCD"})
     .AddSyst(cb, "CMS_htt_qcd_2jet_shape_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
     .channel({"em"})
     .process({"QCD"})
     .AddSyst(cb, "CMS_htt_qcd_2jet_shape2_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
     .channel({"em"})
     .process({"QCD"})
     .AddSyst(cb, "CMS_htt_qcd_iso", "shape", SystMap<>::init(1.00));

  // ##########################################################################
  // Uncertainty: Drell-Yan LO->NLO reweighting
  // References:
  // Notes:
  // - FIXME: References?
  // ##########################################################################

  if (era == 2016) {
      cb.cp()
          .channel({"emt", "met", "mmt", "mtt", "ett"})
          .process({"DY"})
          .AddSyst(cb, "CMS_htt_dyShape_$ERA", "shape", SystMap<>::init(0.10));
  } else {
      cb.cp()
          .channel({"emt", "met", "mmt", "mtt", "ett"})
          .process({"DY"})
          .AddSyst(cb, "CMS_htt_dyShape", "shape", SystMap<>::init(0.10));
  }

  // ##########################################################################
  // Uncertainty: TT shape reweighting
  // References:
  // Notes:
  // - FIXME: References?
  // ##########################################################################

  cb.cp()
      .channel({"emt", "met", "mmt", "mtt", "ett"})
      .process({"TT", "rem_ttbar"})
      .AddSyst(cb, "CMS_htt_ttbarShape", "shape", SystMap<>::init(1.00));

  // ##########################################################################
  // Uncertainty: Electron/muon to tau fakes and ZL energy scale
  // References:
  // Notes:
  // - FIXME: References?
  // ##########################################################################

  // ZL energy scale split by decay mode
  cb.cp()
      .channel({"emt"})
      .process({"ZL"})
      .AddSyst(cb, "CMS_ZLShape_$CHANNEL_$ERA", "shape",
               SystMap<>::init(1.00));
  cb.cp()
      .channel({"met"})
      .process({"DY"})
      .AddSyst(cb, "CMS_ZLShape_$CHANNEL_1prong_barrel_$ERA", "shape",
               SystMap<>::init(1.00));

  cb.cp()
      .channel({"met"})
      .process({"DY"})
      .AddSyst(cb, "CMS_ZLShape_$CHANNEL_1prong1pizero_barrel_$ERA", "shape",
               SystMap<>::init(1.00));
  cb.cp()
      .channel({"met"})
      .process({"DY"})
      .AddSyst(cb, "CMS_ZLShape_$CHANNEL_1prong_endcap_$ERA", "shape",
               SystMap<>::init(1.00));

  cb.cp()
      .channel({"met"})
      .process({"DY"})
      .AddSyst(cb, "CMS_ZLShape_$CHANNEL_1prong1pizero_endcap_$ERA", "shape",
               SystMap<>::init(1.00));

  // Electron fakes
  //cb.cp()
  //    .channel({"emt", "met", "ett"})
  //    .process({"ZL"})
  //    .AddSyst(cb, "CMS_fake_e_$ERA", "lnN", SystMap<>::init(1.15));

  cb.cp()
      .channel({"emt", "met", "ett"})
      .process({"DY"})
      .AddSyst(cb, "CMS_fake_e_BA_$ERA", "shape",
               SystMap<>::init(1.00));
  cb.cp()
      .channel({"emt", "met", "ett"})
      .process({"DY"})
      .AddSyst(cb, "CMS_fake_e_EC_$ERA", "shape",
               SystMap<>::init(1.00));

  // Muon fakes
  //cb.cp()
  //    .channel({"emt", "met", "mtt", "mmt"})
  //    .process({"ZL"})
  //    .AddSyst(cb, "CMS_fake_m_$ERA", "lnN", SystMap<>::init(1.25));

  cb.cp()
      .channel({"emt", "met","mmt"})
      .process({"DY"})
      .AddSyst(cb, "CMS_fake_m_WH1_$ERA", "shape",
               SystMap<>::init(1.00));
  cb.cp()
      .channel({"emt", "met","mmt"})
      .process({"DY"})
      .AddSyst(cb, "CMS_fake_m_WH2_$ERA", "shape",
               SystMap<>::init(1.00));
  cb.cp()
      .channel({"emt", "met","mmt"})
      .process({"DY"})
      .AddSyst(cb, "CMS_fake_m_WH3_$ERA", "shape",
               SystMap<>::init(1.00));
  cb.cp()
      .channel({"emt", "met","mmt"})
      .process({"DY"})
      .AddSyst(cb, "CMS_fake_m_WH4_$ERA", "shape",
               SystMap<>::init(1.00));
  cb.cp()
      .channel({"emt", "met","mmt"})
      .process({"DY"})
      .AddSyst(cb, "CMS_fake_m_WH5_$ERA", "shape",
               SystMap<>::init(1.00));

  // ##########################################################################
  // Uncertainty: Jet to tau fakes
  // References:
  // Notes:
  // - FIXME: Adapt for fake factor and embedding
  // - FIXME: References?
  // ##########################################################################

//   cb.cp()
//       .channel({ "emt", "met", "mmt", "mtt", "ett"})
//       .process({"W", "TTJ", "ZJ", "VVJ"})
//       .AddSyst(cb, "CMS_htt_fake_j_$ERA", "shape", SystMap<>::init(1.00));

  // ##########################################################################
  // Uncertainty: Theory uncertainties
  // References:
  // - Gluon-fusion WG1 uncertainty scheme:
  //   https://twiki.cern.ch/twiki/bin/view/CMS/HiggsWG/SignalModelingTools
  // Notes:
  // - FIXME: WG1 scheme currently NOT applied to ggHWW -> on purpose?
  // - FIXME: Add TopMassTreatment from HIG-16043 uncertainty model
  // - FIXME: Compare to HIG-16043 uncertainty model:
  //           - PDF uncertainties split by category?
  //           - QCDUnc uncertainties?
  //           - UEPS uncertainties?
  // - FIXME: Check VH QCD scale uncertainty
  // - FIXME: References?
  // ##########################################################################

  // Uncertainty on branching ratio for HTT at 125 GeV
  cb.cp()
      .channel({"emt", "met", "mmt", "mtt", "ett"})
      .process(signals)
      .AddSyst(cb, "BR_htt_THU", "lnN", SystMap<>::init(1.0117));
  cb.cp()
      .channel({"emt", "met", "mmt", "mtt", "ett"})
      .process(signals)
      .AddSyst(cb, "BR_htt_PU_mq", "lnN", SystMap<>::init(1.0099));
  cb.cp()
      .channel({"emt", "met", "mmt", "mtt", "ett"})
      .process(signals)
      .AddSyst(cb, "BR_htt_PU_alphas", "lnN", SystMap<>::init(1.0061));
  // Uncertainty on branching ratio for HWW at 125 GeV
  cb.cp()
      .channel({"emt", "met", "mmt", "mtt", "ett"})
     .process(JoinStr({signals_ggHToWW,signals_qqHToWW,{"WH_hww", "ZH_hww"}}))
     .AddSyst(cb, "BR_hww_THU", "lnN", SystMap<>::init(1.0098));   
  cb.cp()
      .channel({"emt", "met", "mmt", "mtt", "ett"})
     .process(JoinStr({signals_ggHToWW,signals_qqHToWW,{"WH_hww", "ZH_hww"}}))
     .AddSyst(cb, "BR_hww_PU_mq", "lnN", SystMap<>::init(1.0097));
  cb.cp()
      .channel({"emt", "met", "mmt", "mtt", "ett"})
     .process(JoinStr({signals_ggHToWW,signals_qqHToWW,{"WH_hww", "ZH_hww"}}))
     .AddSyst(cb, "BR_hww_PU_alphas", "lnN", SystMap<>::init(1.0063));
  // QCD scale
  if (!ggh_wg1) {
  cb.cp()
      .channel({"emt", "met", "mmt", "mtt", "ett"})
     .process(JoinStr({signals_ggH,signals_ggHToWW}))
      .AddSyst(cb, "QCDScale_ggH", "lnN", SystMap<>::init(1.039));
  }
  cb.cp()
      .channel({"emt", "met", "mmt", "mtt", "ett"})
      .process({"ZH_had_htt", "ZH_lep_htt", "ggZH_had_htt", "ggZH_lep_htt", "ZH_hww",
      "ZH_lep_FWDH_htt",
      "ZH_lep_PTV_0_75_htt",
      "ZH_lep_PTV_75_150_htt",
      "ZH_lep_PTV_150_250_0J_htt",
      "ZH_lep_PTV_150_250_GE1J_htt",
      "ZH_lep_PTV_GT250_htt",
      "ggZH_lep_FWDH_htt",
      "ggZH_lep_PTV_0_75_htt",
      "ggZH_lep_PTV_75_150_htt",
      "ggZH_lep_PTV_150_250_0J_htt",
      "ggZH_lep_PTV_150_250_GE1J_htt",
      "ggZH_lep_PTV_GT250_htt"})
      .AddSyst(cb, "QCDScale_VH", "lnN", SystMap<>::init(1.009));
  cb.cp()
      .channel({"emt", "met", "mmt", "mtt", "ett"})
      .process(signals_WH)
      .AddSyst(cb, "QCDScale_VH", "lnN", SystMap<>::init(1.008));
  cb.cp()
      .channel({"emt", "met", "mmt", "mtt", "ett"})
      .process({"ttH_htt"})
      .AddSyst(cb, "QCDScale_ttH", "lnN", SystMap<>::init(1.08));

  // PDF
  cb.cp()
      .channel({"emt", "met", "mmt", "mtt", "ett"})
     .process(JoinStr({signals_ggH,signals_ggHToWW}))
      .AddSyst(cb, "pdf_Higgs_gg", "lnN", SystMap<>::init(1.032));
  cb.cp()
      .channel({"emt", "met", "mmt", "mtt", "ett"})
     .process(JoinStr({signals_qqH,signals_qqHToWW}))
      .AddSyst(cb, "pdf_Higgs_qqbar", "lnN", SystMap<>::init(1.021));
  cb.cp()
      .channel({"emt", "met", "mmt", "mtt", "ett"})
      .process({"ZH_had_htt", "ZH_lep_htt", "ggZH_had_htt", "ggZH_lep_htt", "ZH_hww"
      "ZH_lep_FWDH_htt",
      "ZH_lep_PTV_0_75_htt",
      "ZH_lep_PTV_75_150_htt",
      "ZH_lep_PTV_150_250_0J_htt",
      "ZH_lep_PTV_150_250_GE1J_htt",
      "ZH_lep_PTV_GT250_htt",
      "ggZH_lep_FWDH_htt",
      "ggZH_lep_PTV_0_75_htt",
      "ggZH_lep_PTV_75_150_htt",
      "ggZH_lep_PTV_150_250_0J_htt",
      "ggZH_lep_PTV_150_250_GE1J_htt",
      "ggZH_lep_PTV_GT250_htt"})
      .AddSyst(cb, "pdf_Higgs_VH", "lnN", SystMap<>::init(1.013));
  cb.cp()
      .channel({"emt", "met", "mmt", "mtt", "ett"})
      .process(signals_WH)
      .AddSyst(cb, "pdf_Higgs_VH", "lnN", SystMap<>::init(1.018));
  cb.cp()
      .channel({"emt", "met", "mmt", "mtt", "ett"})
      .process({"ttH_htt"})
      .AddSyst(cb, "pdf_Higgs_ttH", "lnN", SystMap<>::init(1.036));

  // Gluon-fusion WG1 uncertainty scheme
  if (ggh_wg1) {
    cb.cp()
      .channel({"emt", "met", "mmt", "mtt", "ett"})
       .process(JoinStr({signals_ggH, signals_ggZH_had}))
      .AddSyst(cb, "THU_ggH_Mig01", "shape", SystMap<>::init(1.00));
    cb.cp()
      .channel({"emt", "met", "mmt", "mtt", "ett"})
       .process(JoinStr({signals_ggH, signals_ggZH_had}))
      .AddSyst(cb, "THU_ggH_Mig12", "shape", SystMap<>::init(1.00));
    cb.cp()
      .channel({"emt", "met", "mmt", "mtt", "ett"})
       .process(JoinStr({signals_ggH, signals_ggZH_had}))
      .AddSyst(cb, "THU_ggH_Mu", "shape", SystMap<>::init(1.00));
    cb.cp()
      .channel({"emt", "met", "mmt", "mtt", "ett"})
       .process(JoinStr({signals_ggH, signals_ggZH_had}))
      .AddSyst(cb, "THU_ggH_PT120", "shape", SystMap<>::init(1.00));
    cb.cp()
      .channel({"emt", "met", "mmt", "mtt", "ett"})
       .process(JoinStr({signals_ggH, signals_ggZH_had}))
      .AddSyst(cb, "THU_ggH_PT60", "shape", SystMap<>::init(1.00));
    cb.cp()
      .channel({"emt", "met", "mmt", "mtt", "ett"})
       .process(JoinStr({signals_ggH, signals_ggZH_had}))
      .AddSyst(cb, "THU_ggH_Res", "shape", SystMap<>::init(1.00));
    cb.cp()
      .channel({"emt", "met", "mmt", "mtt", "ett"})
       .process(JoinStr({signals_ggH, signals_ggZH_had}))
      .AddSyst(cb, "THU_ggH_VBF2j", "shape", SystMap<>::init(1.00));
    cb.cp()
      .channel({"emt", "met", "mmt", "mtt", "ett"})
       .process(JoinStr({signals_ggH, signals_ggZH_had}))
      .AddSyst(cb, "THU_ggH_VBF3j", "shape", SystMap<>::init(1.00));
    cb.cp()
      .channel({"emt", "met", "mmt", "mtt", "ett"})
       .process(JoinStr({signals_ggH, signals_ggZH_had}))
      .AddSyst(cb, "THU_ggH_qmtop", "shape", SystMap<>::init(1.00));
    cb.cp()
      .channel({"emt", "met", "mmt", "mtt", "ett"})
     .process({signals_ggHToWW})
      .AddSyst(cb, "QCDScale_ggHWW", "lnN", SystMap<>::init(1.039));
  }

    cb.cp()
     .channel({"emt", "met", "mmt", "mtt", "ett"})
     .process({signals_qqH})
     .AddSyst(cb, "vbf_scale_0jet", "shape", SystMap<>::init(1.00));
    cb.cp()
     .channel({"emt", "met", "mmt", "mtt", "ett"})
     .process({signals_qqH})
     .AddSyst(cb, "vbf_scale_1jet", "shape", SystMap<>::init(1.00));
    cb.cp()
     .channel({"emt", "met", "mmt", "mtt", "ett"})
     .process({signals_qqH})
     .AddSyst(cb, "vbf_scale_lowmjj", "shape", SystMap<>::init(1.00));
    cb.cp()
     .channel({"emt", "met", "mmt", "mtt", "ett"})
     .process({signals_qqH})
     .AddSyst(cb, "vbf_scale_highmjj_lowpt", "shape", SystMap<>::init(1.00));
    cb.cp()
     .channel({"emt", "met", "mmt", "mtt", "ett"})
     .process({signals_qqH})
     .AddSyst(cb, "vbf_scale_highmjj_highpt", "shape", SystMap<>::init(1.00));
    
    cb.cp()
     .channel({"emt", "met", "mmt", "mtt", "ett"})
     .process({signals_ggH})
     .AddSyst(cb, "ggH_scale_0jet", "shape", SystMap<>::init(1.00));
    cb.cp()
     .channel({"emt", "met", "mmt", "mtt", "ett"})
     .process({signals_ggH})
     .AddSyst(cb, "ggH_scale_1jet_lowpt", "shape", SystMap<>::init(1.00));
    cb.cp()
     .channel({"emt", "met", "mmt", "mtt", "ett"})
     .process({signals_ggH})
     .AddSyst(cb, "ggH_scale_2jet_lowpt", "shape", SystMap<>::init(1.00));
    cb.cp()
     .channel({"emt", "met", "mmt", "mtt", "ett"})
     .process({signals_ggH})
     .AddSyst(cb, "ggH_scale_highpt", "shape", SystMap<>::init(1.00));
    cb.cp()
     .channel({"emt", "met", "mmt", "mtt", "ett"})
     .process({signals_ggH})
     .AddSyst(cb, "ggH_scale_very_highpt", "shape", SystMap<>::init(1.00));
    cb.cp()
     .channel({"emt", "met", "mmt", "mtt", "ett"})
     .process({signals_ggH})
     .AddSyst(cb, "ggH_scale_vbf", "shape", SystMap<>::init(1.00));

  // ##########################################################################
  // Uncertainty: Embedded events
  // References:
  // - https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauTauEmbeddingSamples2016
  // Notes:
  // ##########################################################################

  // Embedded Normalization: No Lumi, Zjxsec information used, instead derived from data using dimuon selection efficiency
  cb.cp()
      .channel({"emt", "met", "mmt", "mtt", "ett"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_htt_doublemutrg_$ERA", "lnN", SystMap<>::init(1.04));

  // TTbar contamination in embedded events: 10% shape uncertainty of assumed ttbar->tautau event shape
  cb.cp()
    .channel({"emt", "met", "mmt", "mtt", "ett"})
    .process({"EMB"})
    .AddSyst(cb, "CMS_htt_emb_ttbar_$ERA", "shape", SystMap<>::init(1.00));

  // Uncertainty of hadronic tau track efficiency correction
  // uncorrelated between eras
  cb.cp()
    .channel({ "emt", "met", "mmt", "mtt", "ett"})
    .process({"EMB"})
    .AddSyst(cb, "CMS_3ProngEff_$ERA", "shape", SystMap<>::init(0.71));

  cb.cp()
    .channel({ "emt", "met", "mmt", "mtt", "ett"})
    .process({"EMB"})
    .AddSyst(cb, "CMS_1ProngPi0Eff_$ERA", "shape", SystMap<>::init(0.71));
  // correlated between eras
  cb.cp()
    .channel({ "emt", "met", "mmt", "mtt", "ett"})
    .process({"EMB"})
    .AddSyst(cb, "CMS_3ProngEff", "shape", SystMap<>::init(0.71));

  cb.cp()
    .channel({ "emt", "met", "mmt", "mtt", "ett"})
    .process({"EMB"})
    .AddSyst(cb, "CMS_1ProngPi0Eff", "shape", SystMap<>::init(0.71));

  // ##########################################################################
  // Uncertainty: Jet fakes
  // ##########################################################################

  cb.cp()
      .channel({"emt", "met", "mmt", "mtt", "ett"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_stat_$ERA", "shape", SystMap<>::init(1.00));  
  cb.cp()
      .channel({"emt", "met", "mmt", "mtt", "ett"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_syst_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"mtt", "ett"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_njet", "lnN", SystMap<>::init(1.10));  
  cb.cp()
      .channel({"emt", "met", "mmt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_njet_bkgcomp_$CHANNEL", "lnN", SystMap<>::init(1.15));

  // ##########################################################################
  // Uncertainty: Jet fakes
  // References:
  // - https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauJet2TauFakes
  // Notes:
  // - FIXME: add 2017 norm uncertainties, and properly correlate across years
  // ##########################################################################
//   cb.cp()
//       .channel({"mtt", "ett"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_njet", "lnN", SystMap<>::init(1.10));
//   cb.cp()
//       .channel({"emt", "met", "mmt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_njet_bkgcomp_$CHANNEL", "lnN", SystMap<>::init(1.15));

  // QCD shape stat.
//   cb.cp()
//       .channel({ "emt", "met", "mmt", "mtt", "ett"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_qcd_njet0_morphed_stat_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
//   cb.cp()
//       .channel({ "emt", "met", "mmt", "mtt", "ett"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_qcd_njet1_morphed_stat_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
//   cb.cp()
//       .channel({ "emt", "met", "mmt", "mtt", "ett"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_qcd_njet2_morphed_stat_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));


//   // W shape stat.
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_w_lowdR_njet0_morphed_stat_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_w_lowdR_njet1_morphed_stat_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_w_lowdR_njet2_morphed_stat_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_w_highdR_njet0_morphed_stat_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_w_highdR_njet1_morphed_stat_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_w_highdR_njet2_morphed_stat_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));



//   // TT shape stat.
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_tt_njet0_morphed_stat_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));

//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_tt_njet1_morphed_stat_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
  
//   // MC subtraction uncertainty
//   // uncorrelated between eras
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_w_mc_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({ "emt", "met", "mmt", "mtt", "ett"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_qcd_mc_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_frac_w_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));

      
//   // Shape syst. of different contributions (QCD/W/tt)
//   // uncorrelated between eras
//   cb.cp()
//       .channel({ "emt", "met", "mmt", "mtt", "ett"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_qcd_mvis_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({ "emt", "met", "mmt", "mtt", "ett"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_qcd_mvis_osss_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({ "emt", "met", "mmt", "mtt", "ett"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_corr_qcd_mvis_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({ "emt", "met", "mmt", "mtt", "ett"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_corr_qcd_mvis_osss_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_qcd_muiso_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_corr_qcd_muiso_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({ "mtt", "ett"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_qcd_tau2_pt_0jet_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({ "mtt", "ett"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_corr_qcd_tau2_pt_0jet_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({ "mtt", "ett"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_qcd_tau2_pt_1jet_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({ "mtt", "ett"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_corr_qcd_tau2_pt_1jet_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({ "mtt", "ett"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_w_syst_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({ "mtt", "ett"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_tt_syst_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_tt_morphed_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_corr_tt_syst_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_tt_sf_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_w_lepPt_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_corr_w_lepPt_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_w_mt_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_corr_w_mt_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
            

//   //below: jetFakes norm uncertainties. Current values are for 2016, which are probably a good approx. for 2017. To be updated.


//   // Stat. norm (uncorrelated across years)
//   cb.cp()
//       .channel({ "emt", "met", "mmt", "mtt", "ett"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_norm_stat_$CHANNEL_$BIN_$ERA", "lnN", SystMap<channel, bin_id>::init
// 	       ({"emt", "met", "mtt", "mmt"}, {11},  1.04) //w
// 	       ({"emt", "met", "mtt", "mmt"}, {12},  1.052) //ztt
// 	       ({"emt", "met", "mtt", "mmt"}, {13},  1.051) //tt
// 	       ({"emt", "met", "mtt", "mmt"}, {14},  1.047) //ss
// 	       ({"emt", "met", "mtt", "mmt"}, {15},  1.04) //zll
// 	       ({"emt", "met", "mtt", "mmt"}, {16},  1.059) //misc
// 	       ({"emt", "met", "mtt", "mmt"}, {20},  1.052) //emb
// 	       ({"emt", "met", "mtt", "mmt"}, {21},  1.047) //ff
// 	       ({"emt", "met", "mtt", "mmt"}, {300}, 1.037) //incl
// 	       ({"emt", "met", "ett"}, {11},  1.066) //w
// 	       ({"emt", "met", "ett"}, {12},  1.095) //ztt
// 	       ({"emt", "met", "ett"}, {13},  1.083) //tt
// 	       ({"emt", "met", "ett"}, {14},  1.054) //ss
// 	       ({"emt", "met", "ett"}, {15},  1.095) //zll
// 	       ({"emt", "met", "ett"}, {16},  1.107) //misc
// 	       ({"emt", "met", "ett"}, {20},  1.095) //emb
// 	       ({"emt", "met", "ett"}, {21},  1.066) //ff
// 	       ({"emt", "met", "ett"}, {300}, 1.065) //incl
// 	       ({ "mtt", "ett"}, {12},  1.049) //ztt
// 	       ({ "mtt", "ett"}, {16},  1.028) //misc
// 	       ({ "mtt", "ett"}, {17},  1.041) //noniso
// 	       ({ "mtt", "ett"}, {20},  1.049) //emb
// 	       ({ "mtt", "ett"}, {21},  1.041) //ff
// 	       ({ "mtt", "ett"}, {300}, 1.041) //incl
// 	       );
//   // ggH and qqH categories
//   cb.cp()
//       .channel({ "emt", "met", "mmt", "mtt", "ett"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_norm_stat_$CHANNEL_ggH_$ERA", "lnN", SystMap<channel, bin_id>::init
// 	       ({"emt", "met", "mtt", "mmt"}, {1, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110},  1.049)
// 	       ({"emt", "met", "ett"}, {1, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110},  1.074)
// 	       ({ "mtt", "ett"}, {1, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110},  1.041)
// 	       );

//   cb.cp()
//       .channel({ "emt", "met", "mmt", "mtt", "ett"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_norm_stat_$CHANNEL_qqH_$ERA", "lnN", SystMap<channel, bin_id>::init
// 	       ({"emt", "met", "mtt", "mmt"}, {2, 200, 201, 202, 203},  1.068)
// 	       ({"emt", "met", "ett"}, {2, 200, 201, 202, 203},  1.112)
// 	       ({ "mtt", "ett"}, {2, 200, 201, 202, 203},  1.052)
// 	       );
    
//   // Syst. norm: Bin-correlated
//   // uncorrelated between eras

//   cb.cp()
//       .channel({ "emt", "met", "mmt", "mtt", "ett"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_jetbinned_stat_0jet_norm_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({ "emt", "met", "mmt", "mtt", "ett"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_jetbinned_stat_1jet_norm_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({ "emt", "met", "mmt", "mtt", "ett"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_jetbinned_stat_2jet_norm_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));

//   /*
//   cb.cp()
//       .channel({ "emt", "met", "mmt", "mtt", "ett"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_norm_syst_$CHANNEL_$ERA", "lnN", SystMap<channel, bin_id>::init
// 	       ({"emt", "met", "mtt", "mmt"}, {1},     1.069) //ggh
// 	       ({"emt", "met", "mtt", "mmt"}, {100},   1.069) //ggh
// 	       ({"emt", "met", "mtt", "mmt"}, {101},   1.069) //ggh
// 	       ({"emt", "met", "mtt", "mmt"}, {102},   1.069) //ggh
// 	       ({"emt", "met", "mtt", "mmt"}, {103},   1.069) //ggh
// 	       ({"emt", "met", "mtt", "mmt"}, {104},   1.069) //ggh
//                ({"emt", "met", "mtt", "mmt"}, {105},   1.069) //ggh
//                ({"emt", "met", "mtt", "mmt"}, {106},   1.069) //ggh
//                ({"emt", "met", "mtt", "mmt"}, {107},   1.069) //ggh
//                ({"emt", "met", "mtt", "mmt"}, {108},   1.069) //ggh
//                ({"emt", "met", "mtt", "mmt"}, {109},   1.069) //ggh
//                ({"emt", "met", "mtt", "mmt"}, {110},   1.069) //ggh
// 	       ({"emt", "met", "mtt", "mmt"}, {2},     1.058) //qqh
// 	       ({"emt", "met", "mtt", "mmt"}, {200},   1.058) //qqh
// 	       ({"emt", "met", "mtt", "mmt"}, {201},   1.058) //qqh
// 	       ({"emt", "met", "mtt", "mmt"}, {202},   1.058) //qqh
// 	       ({"emt", "met", "mtt", "mmt"}, {203},   1.058) //qqh
// 	       ({"emt", "met", "mtt", "mmt"}, {11},  1.054) //w
// 	       ({"emt", "met", "mtt", "mmt"}, {12},  1.098) //ztt
// 	       ({"emt", "met", "mtt", "mmt"}, {13},  1.052) //tt
// 	       ({"emt", "met", "mtt", "mmt"}, {14},  1.091) //ss
// 	       ({"emt", "met", "mtt", "mmt"}, {15},  1.068) //zll
// 	       ({"emt", "met", "mtt", "mmt"}, {16},  1.091) //misc
// 	       ({"emt", "met", "mtt", "mmt"}, {20},  1.098) //emb
// 	       ({"emt", "met", "mtt", "mmt"}, {21},  1.064) //ff
// 	       ({"emt", "met", "mtt", "mmt"}, {300}, 1.059) //incl
// 	       ({"emt", "met", "ett"}, {1},     1.059) //ggh
// 	       ({"emt", "met", "ett"}, {100},   1.059) //ggh
// 	       ({"emt", "met", "ett"}, {101},   1.059) //ggh
// 	       ({"emt", "met", "ett"}, {102},   1.059) //ggh
// 	       ({"emt", "met", "ett"}, {103},   1.059) //ggh
// 	       ({"emt", "met", "ett"}, {104},   1.059) //ggh
//                ({"emt", "met", "ett"}, {105},   1.059) //ggh
//                ({"emt", "met", "ett"}, {106},   1.059) //ggh
//                ({"emt", "met", "ett"}, {107},   1.059) //ggh
//                ({"emt", "met", "ett"}, {108},   1.059) //ggh
//                ({"emt", "met", "ett"}, {109},   1.059) //ggh
//                ({"emt", "met", "ett"}, {110},   1.059) //ggh
// 	       ({"emt", "met", "ett"}, {2},     1.057) //qqh
// 	       ({"emt", "met", "ett"}, {200},   1.057) //qqh
// 	       ({"emt", "met", "ett"}, {201},   1.057) //qqh
// 	       ({"emt", "met", "ett"}, {202},   1.057) //qqh
// 	       ({"emt", "met", "ett"}, {203},   1.057) //qqh
// 	       ({"emt", "met", "ett"}, {11},  1.052) //w
// 	       ({"emt", "met", "ett"}, {12},  1.088) //ztt
// 	       ({"emt", "met", "ett"}, {13},  1.057) //tt
// 	       ({"emt", "met", "ett"}, {14},  1.064) //ss
// 	       ({"emt", "met", "ett"}, {15},  1.072) //zll
// 	       ({"emt", "met", "ett"}, {16},  1.058) //misc
// 	       ({"emt", "met", "ett"}, {20},  1.088) //ztt
// 	       ({"emt", "met", "ett"}, {21},  1.057) //ff
// 	       ({"emt", "met", "ett"}, {300}, 1.059) //incl
// 	       ({ "mtt", "ett"}, {1},     1.096) //ggh
// 	       ({ "mtt", "ett"}, {100},   1.096) //ggh
// 	       ({ "mtt", "ett"}, {101},   1.096) //ggh
// 	       ({ "mtt", "ett"}, {102},   1.096) //ggh
// 	       ({ "mtt", "ett"}, {103},   1.096) //ggh
// 	       ({ "mtt", "ett"}, {104},   1.096) //ggh
//                ({ "mtt", "ett"}, {105},   1.096) //ggh
//                ({ "mtt", "ett"}, {106},   1.096) //ggh
//                ({ "mtt", "ett"}, {107},   1.096) //ggh
//                ({ "mtt", "ett"}, {108},   1.096) //ggh
//                ({ "mtt", "ett"}, {109},   1.096) //ggh
//                ({ "mtt", "ett"}, {110},   1.096) //ggh
// 	       ({ "mtt", "ett"}, {2},     1.095) //qqh
// 	       ({ "mtt", "ett"}, {200},   1.095) //qqh
// 	       ({ "mtt", "ett"}, {201},   1.095) //qqh
// 	       ({ "mtt", "ett"}, {202},   1.095) //qqh
// 	       ({ "mtt", "ett"}, {203},   1.095) //qqh
// 	       ({ "mtt", "ett"}, {12},  1.095) //ztt
// 	       ({ "mtt", "ett"}, {16},  1.11) //misc
// 	       ({ "mtt", "ett"}, {17},  1.099) //noniso
// 	       ({ "mtt", "ett"}, {20},  1.095) //emb
// 	       ({ "mtt", "ett"}, {21},  1.099) //ff
// 	       ({ "mtt", "ett"}, {300}, 1.095) //incl
// 	       );
//     */
//   // Syst. norm: Bin-dependent, correlated across years
//   // uncorrelated between eras
//   cb.cp()
//       .channel({ "emt", "met", "mmt", "mtt", "ett"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_sub_syst_$CHANNEL_$BIN_$ERA", "lnN", SystMap<channel, bin_id>::init
// 	       ({"emt", "met", "mtt", "mmt"}, {11},  1.025) //w
// 	       ({"emt", "met", "mtt", "mmt"}, {12},  1.045) //ztt
// 	       ({"emt", "met", "mtt", "mmt"}, {13},  1.03) //tt
// 	       ({"emt", "met", "mtt", "mmt"}, {14},  1.02) //ss
// 	       ({"emt", "met", "mtt", "mmt"}, {15},  1.04) //zll
// 	       ({"emt", "met", "mtt", "mmt"}, {16},  1.035) //misc
// 	       ({"emt", "met", "mtt", "mmt"}, {20},  1.045) //emb
// 	       ({"emt", "met", "mtt", "mmt"}, {21},  1.024) //ss
// 	       ({"emt", "met", "mtt", "mmt"}, {300}, 1.035) //incl
// 	       ({"emt", "met", "ett"}, {11},  1.02) //w
// 	       ({"emt", "met", "ett"}, {12},  1.04) //ztt
// 	       ({"emt", "met", "ett"}, {13},  1.03) //tt
// 	       ({"emt", "met", "ett"}, {14},  1.02) //ss
// 	       ({"emt", "met", "ett"}, {15},  1.04) //zll
// 	       ({"emt", "met", "ett"}, {16},  1.035) //misc
// 	       ({"emt", "met", "ett"}, {20},  1.04) //emb
// 	       ({"emt", "met", "ett"}, {21},  1.023) //ff
// 	       ({"emt", "met", "ett"}, {300}, 1.035) //incl
// 	       ({ "mtt", "ett"}, {12},  1.035) //ztt
// 	       ({ "mtt", "ett"}, {16},  1.03) //misc
// 	       ({ "mtt", "ett"}, {17},  1.02) //noniso
// 	       ({ "mtt", "ett"}, {20},  1.035) //emb
// 	       ({ "mtt", "ett"}, {21},  1.02) //ff
// 	       ({ "mtt", "ett"}, {300}, 1.03) //incl
// 	       );

//   // ggH and qqH categories
//   cb.cp()
//       .channel({ "emt", "met", "mmt", "mtt", "ett"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_sub_syst_$CHANNEL_ggH_$ERA", "lnN", SystMap<channel, bin_id>::init
// 	       ({"emt", "met", "mtt", "mmt"}, {1, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110},  1.04)
// 	       ({"emt", "met", "ett"}, {1, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110},  1.04)
// 	       ({ "mtt", "ett"}, {1, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110},  1.03)
// 	       );

//   cb.cp()
//       .channel({ "emt", "met", "mmt", "mtt", "ett"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_sub_syst_$CHANNEL_qqH_$ERA", "lnN", SystMap<channel, bin_id>::init
// 	       ({"emt", "met", "mtt", "mmt"}, {2, 200, 201, 202, 203},  1.04)
// 	       ({"emt", "met", "ett"}, {2, 200, 201, 202, 203},  1.035)
// 	       ({ "mtt", "ett"}, {2, 200, 201, 202, 203},  1.03)
// 	       );


}
} // namespace ch
