#include "CombineHarvester/SMRun2Legacy/interface/HttSystematics_NMSSMRun2UL.h"
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

  void AddRun2Systematics(CombineHarvester &cb, bool jetfakes, bool embedding, int era) {

  // ##########################################################################
  // Define groups of processes
  // ##########################################################################

  // Signal processes
      // NMSSM
    
  std::vector<std::string> signals = {
    "NMSSM_Ytt", "NMSSM_Ybb"
  };
//   std::vector<std::string> signals_ggH = {
//       // STXS stage 0
//       "ggH_htt",
//       // STXS stage 1.1
//       "ggH_FWDH_htt",
//       "ggH_PTH_200_300_htt",
//       "ggH_PTH_300_450_htt",
//       "ggH_PTH_450_650_htt",
//       "ggH_PTH_GT650_htt",
//       "ggH_0J_PTH_0_10_htt",
//       "ggH_0J_PTH_GT10_htt",
//       "ggH_1J_PTH_0_60_htt",
//       "ggH_1J_PTH_60_120_htt",
//       "ggH_1J_PTH_120_200_htt",
//       "ggH_GE2J_MJJ_0_350_PTH_0_60_htt",
//       "ggH_GE2J_MJJ_0_350_PTH_60_120_htt",
//       "ggH_GE2J_MJJ_0_350_PTH_120_200_htt",
//       "ggH_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_0_25_htt",
//       "ggH_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_GT25_htt",
//       "ggH_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_0_25_htt",
//       "ggH_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_GT25_htt",
//       };
//   std::vector<std::string> signals_ggZH_had = {
//       // STXS stage 0
//       "ggZH_had_htt",
//       // STXS stage 1.1
//       "ggZH_had_FWDH_htt",
//       "ggZH_had_PTH_200_300_htt",
//       "ggZH_had_PTH_300_450_htt",
//       "ggZH_had_PTH_450_650_htt",
//       "ggZH_had_PTH_GT650_htt",
//       "ggZH_had_0J_PTH_0_10_htt",
//       "ggZH_had_0J_PTH_GT10_htt",
//       "ggZH_had_1J_PTH_0_60_htt",
//       "ggZH_had_1J_PTH_60_120_htt",
//       "ggZH_had_1J_PTH_120_200_htt",
//       "ggZH_had_GE2J_MJJ_0_350_PTH_0_60_htt",
//       "ggZH_had_GE2J_MJJ_0_350_PTH_60_120_htt",
//       "ggZH_had_GE2J_MJJ_0_350_PTH_120_200_htt",
//       "ggZH_had_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_0_25_htt",
//       "ggZH_had_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_GT25_htt",
//       "ggZH_had_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_0_25_htt",
//       "ggZH_had_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_GT25_htt",
//       };
//   std::vector<std::string> signals_qqH = {
//       // STXS stage 0
//       "qqH_htt",
//       // STXS stage 1
//       "qqH_FWDH_htt",
//       "qqH_0J_htt",
//       "qqH_1J_htt",
//       "qqH_GE2J_MJJ_0_60_htt",
//       "qqH_GE2J_MJJ_60_120_htt",
//       "qqH_GE2J_MJJ_120_350_htt",
//       "qqH_GE2J_MJJ_GT350_PTH_GT200_htt",
//       "qqH_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_0_25_htt",
//       "qqH_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_GT25_htt",
//       "qqH_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_0_25_htt",
//       "qqH_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_GT25_htt"
//       };
//   std::vector<std::string> signals_VH_had = {
//       // STXS stage 0
//       "WH_had_htt",
//       "ZH_had_htt",
//       // STXS stage 1
//       "WH_had_FWDH_htt",
//       "WH_had_0J_htt",
//       "WH_had_1J_htt",
//       "WH_had_GE2J_MJJ_0_60_htt",
//       "WH_had_GE2J_MJJ_60_120_htt",
//       "WH_had_GE2J_MJJ_120_350_htt",
//       "WH_had_GE2J_MJJ_GT350_PTH_GT200_htt",
//       "WH_had_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_0_25_htt",
//       "WH_had_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_GT25_htt",
//       "WH_had_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_0_25_htt",
//       "WH_had_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_GT25_htt",
//       "ZH_had_FWDH_htt",
//       "ZH_had_0J_htt",
//       "ZH_had_1J_htt",
//       "ZH_had_GE2J_MJJ_0_60_htt",
//       "ZH_had_GE2J_MJJ_60_120_htt",
//       "ZH_had_GE2J_MJJ_120_350_htt",
//       "ZH_had_GE2J_MJJ_GT350_PTH_GT200_htt",
//       "ZH_had_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_0_25_htt",
//       "ZH_had_GE2J_MJJ_350_700_PTH_0_200_PTHJJ_GT25_htt",
//       "ZH_had_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_0_25_htt",
//       "ZH_had_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_GT25_htt"
//       };
//   std::vector<std::string> signals_VH = {
//       // STXS stage 0
//       "WH_lep_htt", "ZH_lep_htt", "ggZH_lep_htt", "ttH_htt",
//       "WH_htt", "ZH_htt",
//       // STXS stage 1
//       "WH_lep_FWDH_htt",
//       "WH_lep_PTV_0_75_htt",
//       "WH_lep_PTV_75_150_htt",
//       "WH_lep_PTV_150_250_0J_htt",
//       "WH_lep_PTV_150_250_GE1J_htt",
//       "WH_lep_PTV_GT250_htt",
//       "ZH_lep_FWDH_htt",
//       "ZH_lep_PTV_0_75_htt",
//       "ZH_lep_PTV_75_150_htt",
//       "ZH_lep_PTV_150_250_0J_htt",
//       "ZH_lep_PTV_150_250_GE1J_htt",
//       "ZH_lep_PTV_GT250_htt",
//       "ggZH_lep_FWDH_htt",
//       "ggZH_lep_PTV_0_75_htt",
//       "ggZH_lep_PTV_75_150_htt",
//       "ggZH_lep_PTV_150_250_0J_htt",
//       "ggZH_lep_PTV_150_250_GE1J_htt",
//       "ggZH_lep_PTV_GT250_htt"
//       };
//   std::vector<std::string> signals_ggHToWW = {
//      // STXS stage 0
//      "ggH_hww"};
//   std::vector<std::string> signals_qqHToWW = {
//      // STXS stage 0
//      "qqH_hww"};
//   std::vector<std::string> signals = JoinStr({signals_ggH, signals_ggZH_had, signals_qqH, signals_VH_had, signals_VH});

  // Background processes
  /* // Not used in the function, keep it for documentation purposes.
  std::vector<std::string> backgrounds = {"W", "QCD", "ZL_NLO",
                                          "TTL", "STL", "VVL", "jetFakes", "EMB", 
                                          "ggH125", "qqH125", "VH125"};

  std::vector<std::string> backgrounds = {"ZTT",  "W",   "ZL",      "ZJ",
                                          "TTT",  "TTJ", "VVT",     "VVJ",
                                          "EWKZ", "QCD", "jetFakes", "EMB", "TTL"};
  */

  // All processes being taken from simulation
  // FIXME: Adapt for fake factor and embedding
  std::vector<std::string> mc_processes =
      JoinStr({
              signals,
            //   signals_ggHToWW,
            //   signals_qqHToWW,
            //   {"WH_hww", "ZH_hww"},
              {"ZTT_NLO", "ZJ_NLO", "ZL_NLO", "TTT", "TTL", "TTJ", "VVT", "VVL", "VVJ", "STT", "STL", "STJ", "W", "ggH_tt", "qqH_tt", "VH_tt", "ggH_bb", "qqH_bb", "VH_bb"}
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

//   float lumi_unc = 1.0;
//   float lumi_unc_corr = 1.0;
//   float lumi_unc_1718 = 1.0;
//   if (era == 2016) {
//       lumi_unc = 1.010;
//       lumi_unc_corr = 1.006;
//   } else if (era == 2017) {
//       lumi_unc = 1.020;
//       lumi_unc_corr = 1.009;
//       lumi_unc_1718 = 1.006;
//   } else if (era == 2018) {
//       lumi_unc = 1.015;
//       lumi_unc_corr = 1.020;
//       lumi_unc_1718 = 1.002;
//   }
//   cb.cp()
//       .channel({"et", "mt", "tt"})
//       .process(mc_processes)
//       .AddSyst(cb, "lumi_13TeV_$ERA", "lnN", SystMap<>::init(lumi_unc));
//   cb.cp()
//       .channel({"et", "mt", "tt"})
//       .process(mc_processes)
//       .AddSyst(cb, "lumi_13TeV_correlated", "lnN", SystMap<>::init(lumi_unc_corr));
//   cb.cp()
//       .channel({"et", "mt", "tt"})
//       .process(mc_processes)
//       .AddSyst(cb, "lumi_13TeV_1718", "lnN", SystMap<>::init(lumi_unc_1718));
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process(mc_processes)
      .AddSyst(cb, "lumi_13TeV_$ERA", "lnN", SystMap<>::init(1.025));

  // ##########################################################################
  // Uncertainty: Pileup
  // References:
  // Notes:
  // - FIXME: References?
  // ##########################################################################
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_PileUp", "shape", SystMap<>::init(1.00));

  // ##########################################################################
  // Uncertainty: Prefiring
  // References:
  // - "https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1ECALPrefiringWeightRecipe"
  // Notes:
  // - FIXME: assumed as uncorrelated accross the years for now, what is the recommendation?
  // ##########################################################################
  if (era != 2018) {
  cb.cp()
      .channel({"et", "mt", "tt"})
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
    .channel({"et"})
    .process(mc_processes)
    .AddSyst(cb, "CMS_eff_trigger_et_$ERA", "shape", SystMap<>::init(1.00));
// TODO add xtrigger for et
//   cb.cp()
//     .channel({"et"})
//     .process(mc_processes)
//     .AddSyst(cb, "CMS_eff_xtrigger_l_et_$ERA", "shape", SystMap<>::init(1.00));
  // 100% uncorrelated for embedded
  cb.cp()
    .channel({"et"})
    .process({"EMB"})
    .AddSyst(cb, "CMS_eff_trigger_emb_et_$ERA", "shape", SystMap<>::init(1.00));
// TODO add xtrigger for et
//   cb.cp()
//     .channel({"et"})
//     .process({"EMB"})
//     .AddSyst(cb, "CMS_eff_xtrigger_l_emb_et_$ERA", "shape", SystMap<>::init(1.00));

  cb.cp()
      .channel({"mt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_eff_trigger_mt_$ERA", "shape", SystMap<>::init(1.00));
// TODO add xtrigger for mt
//   cb.cp()
//       .channel({"mt"})
//       .process(mc_processes)
//       .AddSyst(cb, "CMS_eff_xtrigger_l_mt_$ERA", "shape", SystMap<>::init(1.00));
  // 100% uncorrelated for embedded
  cb.cp()
      .channel({"mt"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_eff_trigger_emb_mt_$ERA", "shape", SystMap<>::init(1.00));
// TODO add xtrigger for mt
//   cb.cp()
//       .channel({"mt"})
//       .process({"EMB"})
//       .AddSyst(cb, "CMS_eff_xtrigger_l_emb_mt_$ERA", "shape", SystMap<>::init(1.00));

  cb.cp()
      .channel({"tt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_eff_trigger_tt_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"tt"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_eff_trigger_emb_tt_$ERA", "shape", SystMap<>::init(0.866));
  // 50% correlation for MC and EMB
  cb.cp()
      .channel({"tt"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_eff_trigger_tt_$ERA", "shape", SystMap<>::init(0.5));

  // Tau trigger efficiencies implemented as shape uncertainties in all channels.
//   std::string tauTriggerdmbins[4] = {"0", "1", "10", "11"};
//   for (auto tauTriggerbin: tauTriggerdmbins)
//   {
//       // mt cross trigger
//     // TODO add xtrigger for et and mt
//     //   cb.cp()
//     //       .channel({"mt", "et"})
//     //       .process(mc_processes)
//     //       .AddSyst(cb, "CMS_eff_xtrigger_t_$CHANNEL_dm"+tauTriggerbin+"_$ERA", "shape", SystMap<>::init(1.00));

//     //   cb.cp()
//     //       .channel({"mt", "et"})
//     //       .process({"EMB"})
//     //       .AddSyst(cb, "CMS_eff_xtrigger_t_emb_$CHANNEL_dm"+tauTriggerbin+"_$ERA", "shape", SystMap<>::init(0.866));

//     //   // Correlated component acting on Embedded
//     //   cb.cp()
//     //       .channel({"mt", "et"})
//     //       .process({"EMB"})
//     //       .AddSyst(cb, "CMS_eff_xtrigger_t_$CHANNEL_dm"+tauTriggerbin+"_$ERA", "shape", SystMap<>::init(0.5));

//       // di-tau trigger
//       cb.cp()
//           .channel({"tt"})
//           .process(mc_processes)
//           .AddSyst(cb, "CMS_eff_trigger_tt_dm"+tauTriggerbin+"_$ERA", "shape", SystMap<>::init(1.00));
          
//       cb.cp()
//           .channel({"tt"})
//           .process({"EMB"})
//           .AddSyst(cb, "CMS_eff_trigger_emb_tt_dm"+tauTriggerbin+"_$ERA", "shape", SystMap<>::init(0.866));

//       // Correlated component acting on Embedded
//       cb.cp()
//           .channel({"tt"})
//           .process({"EMB"})
//           .AddSyst(cb, "CMS_eff_trigger_tt_dm"+tauTriggerbin+"_$ERA", "shape", SystMap<>::init(0.5));
//   }

  // ##########################################################################
  // Uncertainty: Electron, muon ID/Iso and tau ID efficiency
  // References:
  // Notes:
  // - FIXME: Adapt for fake factor and embedding
  // - FIXME: Handling of ZL in fully-hadronic channel?
  // - FIXME: References?
  // ##########################################################################

  // Electron ID
  cb.cp()
      .channel({"et"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_eff_e", "lnN", SystMap<>::init(1.02));
  cb.cp()
      .channel({"et"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_eff_e_emb", "lnN", SystMap<>::init(1.017));
  // 50% correlated between MC and EMB
  cb.cp()
      .channel({"et"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_eff_e", "lnN", SystMap<>::init(1.01));
  // Electron Iso
  cb.cp()
      .channel({"et"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_iso_e", "lnN", SystMap<>::init(1.02));
  cb.cp()
      .channel({"et"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_iso_e_emb", "lnN", SystMap<>::init(1.017));
  // 50% correlated between MC and EMB
  cb.cp()
      .channel({"et"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_iso_e", "lnN", SystMap<>::init(1.01));

  // Muon ID
  cb.cp()
      .channel({"mt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_eff_m", "lnN", SystMap<>::init(1.02));
  cb.cp()
      .channel({"mt"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_eff_m_emb", "lnN", SystMap<>::init(1.017));
  // 50% correlated between MC and EMB
  cb.cp()
      .channel({"mt"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_eff_m", "lnN", SystMap<>::init(1.01));
  // Muon Iso
  cb.cp()
      .channel({"mt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_iso_m", "lnN", SystMap<>::init(1.02));
  cb.cp()
      .channel({"mt"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_iso_m_emb", "lnN", SystMap<>::init(1.017));
  // 50% correlated between MC and EMB
  cb.cp()
      .channel({"mt"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_iso_m", "lnN", SystMap<>::init(1.01));


  std::string tauIDptbins[5] = {"30-35", "35-40", "40-500", "500-1000", "1000-Inf"};
  std::string tauIDptbins_emb[5] = {"20-25", "25-30", "30-35", "35-40", "40-Inf"};
  std::string tauIDptbins_emb_corr[3] = {"30-35", "35-40", "40-500"};
  std::string tauIDdmbins[4] = {"0", "1", "10", "11"};

  // Common component acting on MC

  // 3% in Tau ID SF with different anti-l fake WP
  cb.cp()
      .channel({"mt", "tt"})
      .process(JoinStr({signals, {"EMB", "ZTT_NLO", "TTT", "TTL", "VVT", "VVL", "STT", "STL", "ggH_tt", "qqH_tt", "VH_tt", "ggH_bb", "qqH_bb", "VH_bb"}}))
      .AddSyst(cb, "CMS_eff_t_wp_$ERA", "lnN", SystMap<>::init(1.03));
  
  // Tau ID: et and mt with 1 real tau
      
  for (auto tauIDbin : tauIDptbins){ //first part correlated between channels for IDvsJets
    cb.cp()
        .channel({"et", "mt"})
        .process(JoinStr({signals, {"ZTT_NLO", "TTT", "TTL", "VVT", "VVL", "STT", "STL", "ggH_tt", "qqH_tt", "VH_tt", "ggH_bb", "qqH_bb", "VH_bb"}}))
        .AddSyst(cb, "CMS_eff_t_"+tauIDbin+"_$ERA", "shape", SystMap<>::init(1.0));
  }
  cb.cp() //second part uncorrelated between channels for IDvsLep
      .channel({"et", "mt"})
      .process(JoinStr({signals, {"ZTT_NLO", "TTT", "TTL", "VVT", "VVL", "STT", "STL", "ggH_tt", "qqH_tt", "VH_tt", "ggH_bb", "qqH_bb", "VH_bb"}}))
      .AddSyst(cb, "CMS_eff_t_$CHANNEL_$ERA", "lnN", SystMap<>::init(1.01));

  // Tau ID: tt with 2 real taus
  for (auto tauIDbin : tauIDdmbins){
    cb.cp()
        .channel({"tt"})
        .process(JoinStr({signals, {"ZTT_NLO", "TTT", "TTL", "VVT", "VVL", "STT", "STL", "ggH_tt", "qqH_tt", "VH_tt", "ggH_bb", "qqH_bb", "VH_bb"}}))
        .AddSyst(cb, "CMS_eff_t_dm"+tauIDbin+"_$ERA", "shape", SystMap<>::init(1.0));
  }
  cb.cp()
      .channel({"tt"})
      .process(JoinStr({signals, {"ZTT_NLO", "TTT", "TTL", "VVT", "VVL", "STT", "STL", "ggH_tt", "qqH_tt", "VH_tt", "ggH_bb", "qqH_bb", "VH_bb"}}))
      .AddSyst(cb, "CMS_eff_t_$CHANNEL_$ERA", "lnN", SystMap<>::init(1.014));

  // Component for EMB only

  // Tau ID: et and mt with 1 real tau
  for (auto tauIDbin : tauIDptbins_emb){
    cb.cp()
        .channel({"et", "mt"})
        .process({"EMB"})
        .AddSyst(cb, "CMS_eff_t_emb_"+tauIDbin+"_$ERA", "shape", SystMap<>::init(0.866));
  }
  cb.cp()
      .channel({"et", "mt"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_eff_t_emb_$CHANNEL_$ERA", "lnN", SystMap<>::init(1.0087));

  // Tau ID: tt with 2 real taus
  for (auto tauIDbin : tauIDdmbins){
    cb.cp()
        .channel({"tt"})
        .process({"EMB"})
        .AddSyst(cb, "CMS_eff_t_emb_dm"+tauIDbin+"_$ERA", "shape", SystMap<>::init(0.866));
  }
  cb.cp()
      .channel({"tt"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_eff_t_emb_$CHANNEL_$ERA", "lnN", SystMap<>::init(1.012));


  // Common NP acting on EMB
  
  // Tau ID: et and mt with 1 real tau
  for (auto tauIDbin : tauIDptbins_emb_corr){
    cb.cp()
        .channel({"et", "mt"})
        .process({"EMB"})
        .AddSyst(cb, "CMS_eff_t_"+tauIDbin+"_$ERA", "shape", SystMap<>::init(0.5));
  }
  cb.cp()
      .channel({"et", "mt"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_eff_t_$CHANNEL_$ERA", "lnN", SystMap<>::init(1.005));

  // Tau ID: tt with 2 real taus
  for (auto tauIDbin : tauIDdmbins){
    cb.cp()
        .channel({"tt"})
        .process({"EMB"})
        .AddSyst(cb, "CMS_eff_t_dm"+tauIDbin+"_$ERA", "shape", SystMap<>::init(0.5));
  }
  cb.cp()
      .channel({"tt"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_eff_t_$CHANNEL_$ERA", "lnN", SystMap<>::init(1.007));

  // Tau ID: tt with 1 real taus and 1 jet fake 
  // NEEDED?
  cb.cp()
      .channel({"tt"})
      .process({"W", "ZJ_NLO", "TTJ", "VVJ", "STJ"})
      .AddSyst(cb, "CMS_eff_t_$ERA", "lnN", SystMap<>::init(1.06));

  cb.cp()
      .channel({"tt"})
      .process({"W", "ZJ_NLO", "TTJ", "VVJ", "STJ"})
      .AddSyst(cb, "CMS_eff_t_$CHANNEL_$ERA", "lnN", SystMap<>::init(1.02));


  // repeat tt channel for correlated part between 2016 and 2017
  
  // MC uncorrelated uncertainty
  // Tau ID: et and mt with 1 real tau
  /*cb.cp()
      .channel({"et", "mt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_eff_mc_t", "lnN", SystMap<>::init(tauID_corr));

  cb.cp()
      .channel({"et", "mt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_eff_mc_t_$CHANNEL", "lnN", SystMap<>::init(tauID_uncorr));

  // Tau ID: tt with 2 real taus
  cb.cp()
      .channel({"tt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_eff_mc_t", "lnN", SystMap<>::init(ditauID_corr));

  cb.cp()
      .channel({"tt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_eff_mc_t_$CHANNEL", "lnN", SystMap<>::init(ditauID_uncorr));

  // Embedded uncorrelated uncertainty
  // Tau ID: et and mt with 1 real tau
  cb.cp()
      .channel({"et", "mt"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_eff_emb_t", "lnN", SystMap<>::init(tauID_corr));

  cb.cp()
      .channel({"et", "mt"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_eff_emb_t_$CHANNEL", "lnN", SystMap<>::init(tauID_uncorr));

  // Tau ID: tt with 2 real taus
  cb.cp()
      .channel({"tt"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_eff_emb_t", "lnN", SystMap<>::init(ditauID_corr));

  cb.cp()
      .channel({"tt"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_eff_emb_t_$CHANNEL", "lnN", SystMap<>::init(ditauID_uncorr));

  // MC + embedded correlated uncertainty
  // Tau ID: et and mt with 1 real tau
  cb.cp()
      .channel({"et", "mt"})
      .process(JoinStr({mc_processes, {"EMB"}}))
      .AddSyst(cb, "CMS_eff_t", "lnN", SystMap<>::init(tauID_corr));

  cb.cp()
      .channel({"et", "mt"})
      .process(JoinStr({mc_processes, {"EMB"}}))
      .AddSyst(cb, "CMS_eff_t_$CHANNEL", "lnN", SystMap<>::init(tauID_uncorr));

  // Tau ID: tt with 2 real taus
  cb.cp()
      .channel({"tt"})
      .process(JoinStr({mc_processes, {"EMB"}}))
      .AddSyst(cb, "CMS_eff_t", "lnN", SystMap<>::init(ditauID_corr));

  cb.cp()
      .channel({"tt"})
      .process(JoinStr({mc_processes, {"EMB"}}))
      .AddSyst(cb, "CMS_eff_t_$CHANNEL", "lnN", SystMap<>::init(ditauID_uncorr));

  // Tau ID: tt with 1 real taus and 1 jet fake
  cb.cp()
      .channel({"tt"})
      .process({"W", "ZJ", "TTJ", "VVJ"})
      .AddSyst(cb, "CMS_eff_t", "lnN", SystMap<>::init(1.06));

  cb.cp()
      .channel({"tt"})
      .process({"W", "ZJ", "TTJ", "VVJ"})
      .AddSyst(cb, "CMS_eff_t_$CHANNEL", "lnN", SystMap<>::init(1.02));*/

  // ##########################################################################
  // Uncertainty: b-tag and particleNet tag efficiency
  // References:
  // Notes:
  // - FIXME: References?
  // ##########################################################################

// TODO what about those ?
  
  std::string btag_uncs[8] = {"b_HF", "b_HFstats1_$ERA", "b_HFstats2_$ERA", "j_LF", "j_LFstats1_$ERA", "j_LFstats2_$ERA", "c_CFerr1", "c_CFerr2"};
  
  for (auto btag_unc : btag_uncs){
    cb.cp()
        .channel({"et", "mt", "tt"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_btag_"+btag_unc, "shape", SystMap<>::init(1.00));
  }
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_XbbTag_fj_$ERA", "shape", SystMap<>::init(1.00));
//   cb.cp()
//       .channel({"et", "mt", "tt"})
//       .process(mc_processes)
//       .AddSyst(cb, "CMS_XbbTag_fj_$ERA", "lnN", SystMap<>::init(1.1));

  // ##########################################################################
  // Uncertainty: Electron energy scale
  // References:
  // - MC: https://twiki.cern.ch/twiki/bin/view/CMS/EgammaRunIIRecommendations#E_gamma_Energy_Corrections
  // - Embedding: ?
  // Notes:
  // - FIXME: References for embedding missing, need proper correlation accross years for mc, see here: https://twiki.cern.ch/twiki/bin/view/CMS/EgammaRunIIRecommendations#Recommendations_on_Combining_Sys
  // ##########################################################################

  // MC uncorrelated uncertainty
  cb.cp()
      .channel({"et"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_scale_e", "shape", SystMap<>::init(1.00));
      //.AddSyst(cb, "CMS_scale_mc_e", "shape", SystMap<>::init(0.71));
  cb.cp()
      .channel({"et"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_res_e", "shape", SystMap<>::init(1.00));
      //.AddSyst(cb, "CMS_scale_mc_e", "shape", SystMap<>::init(0.71));
      
  // Embedded uncorrelated uncertainty   
  cb.cp()
      .channel({"et"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_scale_e_barrel_emb", "shape", SystMap<>::init(1.00));
      //.AddSyst(cb, "CMS_scale_emb_e", "shape", SystMap<>::init(0.71));
  cb.cp()
      .channel({"et"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_scale_e_endcap_emb", "shape", SystMap<>::init(1.00));

  // MC + embedded correlated uncertainty

  //cb.cp()
  //    .channel({"em", "et"})
  //    .process(JoinStr({mc_processes, {"EMB"}}))
  //    .AddSyst(cb, "CMS_scale_e", "shape", SystMap<>::init(0.71));


  // ##########################################################################
  // Uncertainty: Tau energy scale
  // References:
  // Notes:
  // - Tau energy scale is split by decay mode.
  // - FIXME: References?
  // ##########################################################################


  // Common component acting on MC

  cb.cp()
      .channel({"et", "mt", "tt"})
      .process(JoinStr({signals, {"ZTT", "TTT", "TTL", "VVT", "VVL", "STT", "STL"}}))
      .AddSyst(cb, "CMS_scale_t_1prong_$ERA", "shape", SystMap<>::init(1.0));

  cb.cp()
      .channel({"et", "mt", "tt"})
      .process(JoinStr({signals, {"ZTT", "TTT", "TTL", "VVT", "VVL", "STT", "STL"}}))
      .AddSyst(cb, "CMS_scale_t_1prong1pizero_$ERA", "shape",
               SystMap<>::init(1.0));

  cb.cp()
      .channel({"et", "mt", "tt"})
      .process(JoinStr({signals, {"ZTT", "TTT", "TTL", "VVT", "VVL", "STT", "STL"}}))
      .AddSyst(cb, "CMS_scale_t_3prong_$ERA", "shape", SystMap<>::init(1.0));

  cb.cp()
      .channel({"et", "mt", "tt"})
      .process(JoinStr({signals, {"ZTT", "TTT", "TTL", "VVT", "VVL", "STT", "STL"}}))
      .AddSyst(cb, "CMS_scale_t_3prong1pizero_$ERA", "shape",
               SystMap<>::init(1.0));

  // Component for EMB only

  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_scale_t_emb_1prong_$ERA", "shape", SystMap<>::init(0.866));

  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_scale_t_emb_1prong1pizero_$ERA", "shape", SystMap<>::init(0.866));

  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_scale_t_emb_3prong_$ERA", "shape", SystMap<>::init(0.866));

  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_scale_t_emb_3prong1pizero_$ERA", "shape", SystMap<>::init(0.866));

  // Common component acting on EMB

  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_scale_t_1prong_$ERA", "shape", SystMap<>::init(0.5));

  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_scale_t_1prong1pizero_$ERA", "shape", SystMap<>::init(0.5));

  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_scale_t_3prong_$ERA", "shape", SystMap<>::init(0.5));

  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_scale_t_3prong1pizero_$ERA", "shape", SystMap<>::init(0.5));

  // ##########################################################################
  // Uncertainty: Jet energy scale
  // References:
  // - Talk in CMS Htt meeting by Daniel Winterbottom about regional JES splits:
  //   https://indico.cern.ch/event/740094/contributions/3055870/
  // Notes:
  // ##########################################################################

  // uncorrelated between eras
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_scale_j_Absolute_$ERA", "shape", SystMap<>::init(1.00));

  cb.cp()
      .channel({"et", "mt", "tt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_scale_j_BBEC1_$ERA", "shape", SystMap<>::init(1.00));

  cb.cp()
      .channel({"et", "mt", "tt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_scale_j_EC2_$ERA", "shape", SystMap<>::init(1.00));

  cb.cp()
      .channel({"et", "mt", "tt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_scale_j_HF_$ERA", "shape", SystMap<>::init(1.00));

  cb.cp()
      .channel({"et", "mt", "tt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_scale_j_RelativeSample_$ERA", "shape", SystMap<>::init(1.00));
  // correlated between eras
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_scale_j_Absolute", "shape", SystMap<>::init(1.00));

  cb.cp()
      .channel({"et", "mt", "tt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_scale_j_BBEC1", "shape", SystMap<>::init(1.00));

  cb.cp()
      .channel({"et", "mt", "tt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_scale_j_EC2", "shape", SystMap<>::init(1.00));

  cb.cp()
      .channel({"et", "mt", "tt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_scale_j_HF", "shape", SystMap<>::init(1.00));

  cb.cp()
      .channel({"et", "mt", "tt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_scale_j_FlavorQCD", "shape", SystMap<>::init(1.00));

  cb.cp()
      .channel({"et", "mt", "tt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_scale_j_RelativeBal", "shape", SystMap<>::init(1.00));

  if (era == 2018){
    cb.cp()
        .channel({"et", "mt", "tt"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_HEMIssue_$ERA", "shape", SystMap<>::init(1.00));
  }
  // JER
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_res_j_$ERA", "shape", SystMap<>::init(1.00));


  // ##########################################################################
  // Uncertainty: MET energy scale and Recoil
  // References:
  // Notes:
  // - FIXME: Clustered vs unclustered MET? Inclusion of JES splitting?
  // - FIXME: References?
  // ##########################################################################

  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"TTT", "TTL", "TTJ", "VV", "VVT", "VVL", "VVJ", "STT", "STL", "STJ"})  // change to mc_processes for next iteration
      .AddSyst(cb, "CMS_scale_met_unclustered_$ERA", "shape", SystMap<>::init(1.00));
//   cb.cp()
//       .channel({"et", "mt", "tt"})
//       .process(JoinStr({signals, {"ZTT_NLO", "ZL_NLO", "ZJ_NLO", "W", "ggH_tt", "qqH_tt", "VH_tt", "ggH_bb", "qqH_bb", "VH_bb"}}))
//       .AddSyst(cb, "CMS_htt_boson_scale_met_$ERA", "shape", SystMap<>::init(1.00));
//   cb.cp()
//       .channel({"et", "mt", "tt"})
//       .process(JoinStr({signals, {"ZTT_NLO", "ZL_NLO", "ZJ_NLO", "W", "ggH_tt", "qqH_tt", "VH_tt", "ggH_bb", "qqH_bb", "VH_bb"}}))
//       .AddSyst(cb, "CMS_htt_boson_res_met_$ERA", "shape", SystMap<>::init(1.00));
// TODO do we need this still ?
//   cb.cp()
//       .channel({"et", "mt", "tt"})
//       .process({"EMB"})
//       .AddSyst(cb, "CMS_scale_met_emb", "shape", SystMap<>::init(1.00));
//   cb.cp()
//       .channel({"em"})
//       .process({"EMB"})
//       .AddSyst(cb, "CMS_scale_met_emb_em", "shape", SystMap<>::init(1.00));

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
      .channel({"et", "mt", "tt"})
      .process({"VVT", "VVJ", "VVL"})
      .AddSyst(cb, "CMS_VV_xsec", "lnN", SystMap<>::init(1.056));

  // ST
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"STT", "STL", "STJ"})
      .AddSyst(cb, "CMS_ST_xsec", "lnN", SystMap<>::init(1.027));

  // TT
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"TTT", "TTL", "TTJ"})
      .AddSyst(cb, "CMS_ttbar_xsec", "lnN", SystMap<>::init(1.044));

  // W
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"W"})
      .AddSyst(cb, "CMS_Wj_xsec", "lnN", SystMap<>::init(1.008));

  // Z
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"ZTT", "ZL", "ZJ"})
      .AddSyst(cb, "CMS_Zj_xsec", "lnN", SystMap<>::init(1.02));

  // QCD
//   cb.cp()
//       .channel({"et"})
//       .process({"QCD"})
//       .AddSyst(cb, "CMS_ExtrapSSOS_$CHANNEL_$ERA", "lnN", SystMap<>::init(1.05));
//   cb.cp()
//       .channel({"mt"})
//       .process({"QCD"})
//       .AddSyst(cb, "CMS_ExtrapSSOS_$CHANNEL_$ERA", "lnN", SystMap<>::init(1.03));
//   cb.cp()
//       .channel({"tt"})
//       .process({"QCD"})
//       .AddSyst(cb, "CMS_ExtrapABCD_$CHANNEL_$ERA", "lnN", SystMap<>::init(1.03));

//   cb.cp()
//       .channel({"em"})
//       .process({"QCD"})
//       .AddSyst(cb, "CMS_htt_qcd_0jet_rate_$ERA", "shape", SystMap<>::init(1.00));
//   cb.cp()
//       .channel({"em"})
//       .process({"QCD"})
//       .AddSyst(cb, "CMS_htt_qcd_0jet_shape_$ERA", "shape", SystMap<>::init(1.00));
//   cb.cp()
//       .channel({"em"})
//       .process({"QCD"})
//       .AddSyst(cb, "CMS_htt_qcd_0jet_shape2_$ERA", "shape", SystMap<>::init(1.00));
//   cb.cp()
//       .channel({"em"})
//       .process({"QCD"})
//       .AddSyst(cb, "CMS_htt_qcd_1jet_rate_$ERA", "shape", SystMap<>::init(1.00));
//   cb.cp()
//       .channel({"em"})
//       .process({"QCD"})
//       .AddSyst(cb, "CMS_htt_qcd_1jet_shape_$ERA", "shape", SystMap<>::init(1.00));
//   cb.cp()
//       .channel({"em"})
//       .process({"QCD"})
//       .AddSyst(cb, "CMS_htt_qcd_1jet_shape2_$ERA", "shape", SystMap<>::init(1.00));
//   cb.cp()
//      .channel({"em"})
//      .process({"QCD"})
//      .AddSyst(cb, "CMS_htt_qcd_2jet_rate_$ERA", "shape", SystMap<>::init(1.00));
//   cb.cp()
//      .channel({"em"})
//      .process({"QCD"})
//      .AddSyst(cb, "CMS_htt_qcd_2jet_shape_$ERA", "shape", SystMap<>::init(1.00));
//   cb.cp()
//      .channel({"em"})
//      .process({"QCD"})
//      .AddSyst(cb, "CMS_htt_qcd_2jet_shape2_$ERA", "shape", SystMap<>::init(1.00));
//   cb.cp()
//      .channel({"em"})
//      .process({"QCD"})
//      .AddSyst(cb, "CMS_htt_qcd_iso", "shape", SystMap<>::init(1.00));

  // ##########################################################################
  // Uncertainty: Drell-Yan LO->NLO reweighting
  // References:
  // Notes:
  // - FIXME: References?
  // ##########################################################################

//   if (era == 2016) {
//       cb.cp()
//           .channel({"et", "mt", "tt", "em"})
//           .process({"ZTT", "ZL", "ZJ"})
//           .AddSyst(cb, "CMS_htt_dyShape_$ERA", "shape", SystMap<>::init(0.10));
//   } else {
//       cb.cp()
//           .channel({"et", "mt", "tt", "em"})
//           .process({"ZTT", "ZL", "ZJ"})
//           .AddSyst(cb, "CMS_htt_dyShape", "shape", SystMap<>::init(0.10));
//   }

  // ##########################################################################
  // Uncertainty: TT shape reweighting (top pT)
  // References:
  // Notes:
  // - FIXME: References?
  // ##########################################################################

  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"TTT", "TTL", "TTJ"})
      .AddSyst(cb, "CMS_topPt_Shape", "shape", SystMap<>::init(1.00));

  // ##########################################################################
  // Uncertainty: Electron/muon to tau fakes and ZL energy scale
  // References:
  // Notes:
  // - FIXME: References?
  // ##########################################################################

  // ZL energy scale split by decay mode, for mt, no split by decay mode is avialable
  cb.cp()
      .channel({"mt"})
      .process({"ZL_NLO"})
      .AddSyst(cb, "CMS_ZLShape_$CHANNEL_$ERA", "shape",
               SystMap<>::init(1.00));

  cb.cp()
      .channel({"et"})
      .process({"ZL_NLO"})
      .AddSyst(cb, "CMS_ZLShape_$CHANNEL_1prong_barrel_$ERA", "shape",
               SystMap<>::init(1.00));

  cb.cp()
      .channel({"et"})
      .process({"ZL_NLO"})
      .AddSyst(cb, "CMS_ZLShape_$CHANNEL_1prong1pizero_barrel_$ERA", "shape",
               SystMap<>::init(1.00));
  cb.cp()
      .channel({"et"})
      .process({"ZL_NLO"})
      .AddSyst(cb, "CMS_ZLShape_$CHANNEL_1prong_endcap_$ERA", "shape",
               SystMap<>::init(1.00));

  cb.cp()
      .channel({"et"})
      .process({"ZL_NLO"})
      .AddSyst(cb, "CMS_ZLShape_$CHANNEL_1prong1pizero_endcap_$ERA", "shape",
               SystMap<>::init(1.00));

  // Electron fakes ID

  cb.cp()
      .channel({"et"})
      .process({"ZL_NLO"})
      .AddSyst(cb, "CMS_fake_e_Barrel_$ERA", "shape",
               SystMap<>::init(1.00));
  cb.cp()
      .channel({"et"})
      .process({"ZL_NLO"})
      .AddSyst(cb, "CMS_fake_e_Endcap_$ERA", "shape",
               SystMap<>::init(1.00));

  // Muon fakes ID

  cb.cp()
      .channel({"mt"})
      .process({"ZL_NLO"})
      .AddSyst(cb, "CMS_fake_m_Wheel1_$ERA", "shape",
               SystMap<>::init(1.00));
  cb.cp()
      .channel({"mt"})
      .process({"ZL_NLO"})
      .AddSyst(cb, "CMS_fake_m_Wheel2_$ERA", "shape",
               SystMap<>::init(1.00));
  cb.cp()
      .channel({"mt"})
      .process({"ZL_NLO"})
      .AddSyst(cb, "CMS_fake_m_Wheel3_$ERA", "shape",
               SystMap<>::init(1.00));
  cb.cp()
      .channel({"mt"})
      .process({"ZL_NLO"})
      .AddSyst(cb, "CMS_fake_m_Wheel4_$ERA", "shape",
               SystMap<>::init(1.00));
  cb.cp()
      .channel({"mt"})
      .process({"ZL_NLO"})
      .AddSyst(cb, "CMS_fake_m_Wheel5_$ERA", "shape",
               SystMap<>::init(1.00));

  // ##########################################################################
  // Uncertainty: Jet to tau fakes
  // References:
  // Notes:
  // - FIXME: Adapt for fake factor and embedding
  // - FIXME: References?
  // ##########################################################################

  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"W", "TTJ", "ZJ_NLO", "VVJ", "STJ"})
      .AddSyst(cb, "CMS_htt_fake_j_$ERA", "shape", SystMap<>::init(1.00));

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
  
  // NMSSM signal
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"NMSSM_Ybb"})
      .AddSyst(cb, "LHE_scale_norm", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"NMSSM_Ytt"})
      .AddSyst(cb, "LHE_scale_norm", "shape", SystMap<>::init(1.00));
//   cb.cp()
//       .channel({"et", "mt", "tt"})
//       .process(signals)
//       .AddSyst(cb, "PDF_scale", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"NMSSM_Ybb"})
      .AddSyst(cb, "PDF_scale_Ybb", "lnN", SystMap<>::init(1.18));
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"NMSSM_Ytt"})
      .AddSyst(cb, "PDF_scale_Ytt", "lnN", SystMap<>::init(1.18));

  // Uncertainty on branching ratio for HTT at 125 GeV
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"NMSSM_Ybb", "ggH_tt", "qqH_tt", "VH_tt"})
      .AddSyst(cb, "BR_htt_THU", "lnN", SystMap<>::init(1.0117));
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"NMSSM_Ybb", "ggH_tt", "qqH_tt", "VH_tt"})
      .AddSyst(cb, "BR_htt_PU_mq", "lnN", SystMap<>::init(1.0098));
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"NMSSM_Ybb", "ggH_tt", "qqH_tt", "VH_tt"})
      .AddSyst(cb, "BR_htt_PU_alphaS", "lnN", SystMap<>::init(1.0062));

  // Uncertainty on branching ratio for HBB at 125 GeV
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"NMSSM_Ytt", "ggH_bb", "qqH_bb", "VH_bb"})
      .AddSyst(cb, "BR_hbb_THU", "lnN", SystMap<>::init(1.0065));
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"NMSSM_Ytt", "ggH_bb", "qqH_bb", "VH_bb"})
      .AddSyst(cb, "BR_hbb_PU_mq", "lnN", SystMap<>::init(1.0074));
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"NMSSM_Ytt", "ggH_bb", "qqH_bb", "VH_bb"})
      .AddSyst(cb, "BR_hbb_PU_alphaS", "lnN", SystMap<>::init(1.0079));
  
  // QCD scale
  cb.cp()
      .channel({"et", "mt", "tt"})
     .process({"ggH_tt", "ggH_bb"})
      .AddSyst(cb, "QCDScale_ggH", "lnN", SystMap<>::init(1.039));
  cb.cp()
      .channel({"et", "mt", "tt"})
     .process({"qqH_tt", "qqH_bb"})
      .AddSyst(cb, "QCDScale_qqH", "lnN", SystMap<>::init(1.005));
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"VH_tt", "VH_bb"})
      .AddSyst(cb, "QCDScale_VH", "lnN", SystMap<>::init(1.01));
//   cb.cp()
//       .channel({"et", "mt", "tt"})
//       .process({"ttH125"})
//       .AddSyst(cb, "QCDScale_ttH", "lnN", SystMap<>::init(1.08));

  // PDF
  cb.cp()
      .channel({"et", "mt", "tt"})
     .process({"ggH_tt", "ggH_bb"})
      .AddSyst(cb, "pdf_Higgs_gg", "lnN", SystMap<>::init(1.032));
  cb.cp()
      .channel({"et", "mt", "tt"})
     .process({"qqH_tt", "qqH_bb"})
      .AddSyst(cb, "pdf_Higgs_qqbar", "lnN", SystMap<>::init(1.021));
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"VH_tt", "VH_bb"})
      .AddSyst(cb, "pdf_Higgs_VH", "lnN", SystMap<>::init(1.018));
//   cb.cp()
//       .channel({"et", "mt", "tt"})
//       .process({"ttH125"})
//       .AddSyst(cb, "pdf_Higgs_ttH", "lnN", SystMap<>::init(1.036));

  // ##########################################################################
  // Uncertainty: Embedded events
  // References:
  // - https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauTauEmbeddingSamples2016
  // Notes:
  // ##########################################################################

  // Embedded Normalization: No Lumi, Zjxsec information used, instead derived from data using dimuon selection efficiency
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_emb_doublemutrg_$ERA", "lnN", SystMap<>::init(1.04));

  // TTbar contamination in embedded events: 10% shape uncertainty of assumed ttbar->tautau event shape
  cb.cp()
    .channel({"et", "mt", "tt"})
    .process({"EMB"})
    .AddSyst(cb, "CMS_emb_ttbar_$ERA", "shape", SystMap<>::init(1.00));


  // ##########################################################################
  // Uncertainty: Jet fakes
  // References:
  // - https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauJet2TauFakes
  // Notes:
  // - FIXME: add 2017 norm uncertainties, and properly correlate across years
  // ##########################################################################

  // QCD shape stat.
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_QCDFFslopeUnc_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_QCDFFnormUnc_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_QCDFFmcSubUnc_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_QCDSubleadingFFslopeUnc_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_QCDSubleadingFFnormUnc_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_QCDSubleadingFFmcSubUnc_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));


  // W shape stat.
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_WjetsFFslopeUnc_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_WjetsFFnormUnc_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_WjetsFFmcSubUnc_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_w_highdR_njet0_morphed_stat_", "shape", SystMap<>::init(1.00));
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_w_highdR_njet1_morphed_stat_", "shape", SystMap<>::init(1.00));
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_w_highdR_njet2_morphed_stat_", "shape", SystMap<>::init(1.00));

//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_tt_lowdR_njet0_morphed_stat_", "shape", SystMap<>::init(1.00));

//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_tt_lowdR_njet1_morphed_stat_", "shape", SystMap<>::init(1.00));

  // ttbar shape stat.
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ttbarFFslopeUnc_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ttbarFFnormUnc_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ttbarSubleadingFFslopeUnc_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ttbarSubleadingFFnormUnc_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_frac_w_", "shape", SystMap<>::init(1.0));
  // Process fraction variations
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_fracQCDUnc_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_fracWjetsUnc_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_fracTTbarUnc_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_fracQCDSubleadingUnc_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_fracWjetsSubleadingUnc_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_fracTTbarSubleadingUnc_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));

  // Shape syst. of different contributions (QCD/W/tt)
  // uncorrelated between eras
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_QCDClosureLeadingLepPtCorr_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_QCDClosureSubleadingLepPtCorr_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_QCDSubleadingClosureLeadingLepPtCorr_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_QCDClosureSubleadingTauMassCorr_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_QCDClosureLeadingTauMassCorr_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_QCDSubleadingClosureSubleadingTauMassCorr_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_QCDDRtoSRCorr_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_QCDSubleadingDRtoSRCorr_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));    
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_WjetsClosureLeadingLepPtCorr_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_WjetsClosureSubleadingTauMassCorr_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_WjetsDRtoSRCorr_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ttbarClosureLeadingLepPtCorr_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ttbarClosureSubleadingTauMassCorr_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ttbarClosureSubleadingLepPtCorr_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ttbarClosureLeadingTauMassCorr_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ttbarSubleadingClosureLeadingLepPtCorr_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ttbarSubleadingClosureSubleadingTauMassCorr_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"tt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_qcd_tau2_pt_0jet_", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"tt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_corr_qcd_tau2_pt_0jet_", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"tt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_qcd_tau2_pt_1jet_", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"tt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_corr_qcd_tau2_pt_1jet_", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"tt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_w_syst_", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"tt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_tt_syst_", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_tt_morphed_", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_corr_tt_syst_", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_tt_sf_", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_w_lepPt_", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_corr_w_lepPt_", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_w_mt_", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_corr_w_mt_", "shape", SystMap<>::init(1.0));
// TODO in next iteration add channel and era for unorrelated unc
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



// //   TT shape stat.
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_tt_lowdR_njet0_morphed_stat_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));

//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_tt_lowdR_njet1_morphed_stat_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));

//   // MC subtraction uncertainty
//   // uncorrelated between eras
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_w_mc_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"et", "mt", "tt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_qcd_mc_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_frac_w_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));


//   // Shape syst. of different contributions (QCD/W/tt)
//   // uncorrelated between eras
//   cb.cp()
//       .channel({"et", "mt", "tt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_qcd_mvis_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"et", "mt", "tt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_qcd_mvis_osss_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"et", "mt", "tt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_corr_qcd_mvis_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"et", "mt", "tt"})
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
//       .channel({"tt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_qcd_tau2_pt_0jet_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"tt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_corr_qcd_tau2_pt_0jet_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"tt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_qcd_tau2_pt_1jet_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"tt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_corr_qcd_tau2_pt_1jet_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"tt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_w_syst_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"tt"})
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
// TODO fix in next next iteration
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_$CHANNEL_$ERA_ff_w_lowdR_njet0_morphed_stat_", "shape", SystMap<>::init(1.00));
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_$CHANNEL_$ERA_ff_w_lowdR_njet1_morphed_stat_", "shape", SystMap<>::init(1.00));
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_$CHANNEL_$ERA_ff_w_lowdR_njet2_morphed_stat_", "shape", SystMap<>::init(1.00));
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_$CHANNEL_$ERA_ff_w_highdR_njet0_morphed_stat_", "shape", SystMap<>::init(1.00));
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_$CHANNEL_$ERA_ff_w_highdR_njet1_morphed_stat_", "shape", SystMap<>::init(1.00));
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_$CHANNEL_$ERA_ff_w_highdR_njet2_morphed_stat_", "shape", SystMap<>::init(1.00));



// //   TT shape stat.
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_$CHANNEL_$ERA_ff_tt_lowdR_njet0_morphed_stat_", "shape", SystMap<>::init(1.00));

//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_$CHANNEL_$ERA_ff_tt_lowdR_njet1_morphed_stat_", "shape", SystMap<>::init(1.00));

//   // MC subtraction uncertainty
//   // uncorrelated between eras
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_$CHANNEL_$ERA_ff_w_mc_", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"et", "mt", "tt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_$CHANNEL_$ERA_ff_qcd_mc_", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_$CHANNEL_$ERA_ff_frac_w_", "shape", SystMap<>::init(1.0));


//   // Shape syst. of different contributions (QCD/W/tt)
//   // uncorrelated between eras
//   cb.cp()
//       .channel({"et", "mt", "tt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_$CHANNEL_$ERA_ff_qcd_mvis_", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"et", "mt", "tt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_$CHANNEL_$ERA_ff_qcd_mvis_osss_", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"et", "mt", "tt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_$CHANNEL_$ERA_ff_corr_qcd_mvis_", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"et", "mt", "tt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_$CHANNEL_$ERA_ff_corr_qcd_mvis_osss_", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_$CHANNEL_$ERA_ff_qcd_muiso_", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_$CHANNEL_$ERA_ff_corr_qcd_muiso_", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"tt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_$CHANNEL_$ERA_ff_qcd_tau2_pt_0jet_", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"tt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_$CHANNEL_$ERA_ff_corr_qcd_tau2_pt_0jet_", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"tt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_$CHANNEL_$ERA_ff_qcd_tau2_pt_1jet_", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"tt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_$CHANNEL_$ERA_ff_corr_qcd_tau2_pt_1jet_", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"tt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_$CHANNEL_$ERA_ff_w_syst_", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"tt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_$CHANNEL_$ERA_ff_tt_syst_", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_$CHANNEL_$ERA_ff_tt_morphed_", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_$CHANNEL_$ERA_ff_corr_tt_syst_", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_$CHANNEL_$ERA_ff_tt_sf_", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_$CHANNEL_$ERA_ff_w_lepPt_", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_$CHANNEL_$ERA_ff_corr_w_lepPt_", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_$CHANNEL_$ERA_ff_w_mt_", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_$CHANNEL_$ERA_ff_corr_w_mt_", "shape", SystMap<>::init(1.0));


//   //below: jetFakes norm uncertainties. Current values are for 2016, which are probably a good approx. for 2017. To be updated.


//   // Stat. norm (uncorrelated across years)
//   cb.cp()
//       .channel({"et", "mt", "tt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_norm_stat_$CHANNEL_$BIN_$ERA", "lnN", SystMap<channel, bin_id>::init
// 	       ({"mt"}, {11},  1.04) //w
// 	       ({"mt"}, {12},  1.052) //ztt
// 	       ({"mt"}, {13},  1.051) //tt
// 	       ({"mt"}, {14},  1.047) //ss
// 	       ({"mt"}, {15},  1.04) //zll
// 	       ({"mt"}, {16},  1.059) //misc
// 	       ({"mt"}, {20},  1.052) //emb
// 	       ({"mt"}, {21},  1.047) //ff
// 	       ({"mt"}, {300}, 1.037) //incl
// 	       ({"et"}, {11},  1.066) //w
// 	       ({"et"}, {12},  1.095) //ztt
// 	       ({"et"}, {13},  1.083) //tt
// 	       ({"et"}, {14},  1.054) //ss
// 	       ({"et"}, {15},  1.095) //zll
// 	       ({"et"}, {16},  1.107) //misc
// 	       ({"et"}, {20},  1.095) //emb
// 	       ({"et"}, {21},  1.066) //ff
// 	       ({"et"}, {300}, 1.065) //incl
// 	       ({"tt"}, {12},  1.049) //ztt
// 	       ({"tt"}, {16},  1.028) //misc
// 	       ({"tt"}, {17},  1.041) //noniso
// 	       ({"tt"}, {20},  1.049) //emb
// 	       ({"tt"}, {21},  1.041) //ff
// 	       ({"tt"}, {300}, 1.041) //incl
// 	       );
//   // ggH and qqH categories
//   cb.cp()
//       .channel({"et", "mt", "tt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_norm_stat_$CHANNEL_ggH_$ERA", "lnN", SystMap<channel, bin_id>::init
// 	       ({"mt"}, {1, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110},  1.049)
// 	       ({"et"}, {1, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110},  1.074)
// 	       ({"tt"}, {1, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110},  1.041)
// 	       );

//   cb.cp()
//       .channel({"et", "mt", "tt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_norm_stat_$CHANNEL_qqH_$ERA", "lnN", SystMap<channel, bin_id>::init
// 	       ({"mt"}, {2, 200, 201, 202, 203},  1.068)
// 	       ({"et"}, {2, 200, 201, 202, 203},  1.112)
// 	       ({"tt"}, {2, 200, 201, 202, 203},  1.052)
// 	       );
    
  // Syst. norm: Bin-correlated
  // uncorrelated between eras
// TODO needed ?
//   cb.cp()
//       .channel({"et", "mt", "tt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_jetbinned_stat_0jet_norm_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"et", "mt", "tt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_jetbinned_stat_1jet_norm_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"et", "mt", "tt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_jetbinned_stat_2jet_norm_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));

  /*
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_ff_norm_syst_$CHANNEL_$ERA", "lnN", SystMap<channel, bin_id>::init
	       ({"mt"}, {1},     1.069) //ggh
	       ({"mt"}, {100},   1.069) //ggh
	       ({"mt"}, {101},   1.069) //ggh
	       ({"mt"}, {102},   1.069) //ggh
	       ({"mt"}, {103},   1.069) //ggh
	       ({"mt"}, {104},   1.069) //ggh
               ({"mt"}, {105},   1.069) //ggh
               ({"mt"}, {106},   1.069) //ggh
               ({"mt"}, {107},   1.069) //ggh
               ({"mt"}, {108},   1.069) //ggh
               ({"mt"}, {109},   1.069) //ggh
               ({"mt"}, {110},   1.069) //ggh
	       ({"mt"}, {2},     1.058) //qqh
	       ({"mt"}, {200},   1.058) //qqh
	       ({"mt"}, {201},   1.058) //qqh
	       ({"mt"}, {202},   1.058) //qqh
	       ({"mt"}, {203},   1.058) //qqh
	       ({"mt"}, {11},  1.054) //w
	       ({"mt"}, {12},  1.098) //ztt
	       ({"mt"}, {13},  1.052) //tt
	       ({"mt"}, {14},  1.091) //ss
	       ({"mt"}, {15},  1.068) //zll
	       ({"mt"}, {16},  1.091) //misc
	       ({"mt"}, {20},  1.098) //emb
	       ({"mt"}, {21},  1.064) //ff
	       ({"mt"}, {300}, 1.059) //incl
	       ({"et"}, {1},     1.059) //ggh
	       ({"et"}, {100},   1.059) //ggh
	       ({"et"}, {101},   1.059) //ggh
	       ({"et"}, {102},   1.059) //ggh
	       ({"et"}, {103},   1.059) //ggh
	       ({"et"}, {104},   1.059) //ggh
               ({"et"}, {105},   1.059) //ggh
               ({"et"}, {106},   1.059) //ggh
               ({"et"}, {107},   1.059) //ggh
               ({"et"}, {108},   1.059) //ggh
               ({"et"}, {109},   1.059) //ggh
               ({"et"}, {110},   1.059) //ggh
	       ({"et"}, {2},     1.057) //qqh
	       ({"et"}, {200},   1.057) //qqh
	       ({"et"}, {201},   1.057) //qqh
	       ({"et"}, {202},   1.057) //qqh
	       ({"et"}, {203},   1.057) //qqh
	       ({"et"}, {11},  1.052) //w
	       ({"et"}, {12},  1.088) //ztt
	       ({"et"}, {13},  1.057) //tt
	       ({"et"}, {14},  1.064) //ss
	       ({"et"}, {15},  1.072) //zll
	       ({"et"}, {16},  1.058) //misc
	       ({"et"}, {20},  1.088) //ztt
	       ({"et"}, {21},  1.057) //ff
	       ({"et"}, {300}, 1.059) //incl
	       ({"tt"}, {1},     1.096) //ggh
	       ({"tt"}, {100},   1.096) //ggh
	       ({"tt"}, {101},   1.096) //ggh
	       ({"tt"}, {102},   1.096) //ggh
	       ({"tt"}, {103},   1.096) //ggh
	       ({"tt"}, {104},   1.096) //ggh
               ({"tt"}, {105},   1.096) //ggh
               ({"tt"}, {106},   1.096) //ggh
               ({"tt"}, {107},   1.096) //ggh
               ({"tt"}, {108},   1.096) //ggh
               ({"tt"}, {109},   1.096) //ggh
               ({"tt"}, {110},   1.096) //ggh
	       ({"tt"}, {2},     1.095) //qqh
	       ({"tt"}, {200},   1.095) //qqh
	       ({"tt"}, {201},   1.095) //qqh
	       ({"tt"}, {202},   1.095) //qqh
	       ({"tt"}, {203},   1.095) //qqh
	       ({"tt"}, {12},  1.095) //ztt
	       ({"tt"}, {16},  1.11) //misc
	       ({"tt"}, {17},  1.099) //noniso
	       ({"tt"}, {20},  1.095) //emb
	       ({"tt"}, {21},  1.099) //ff
	       ({"tt"}, {300}, 1.095) //incl
	       );
    */
  // Syst. norm: Bin-dependent, correlated across years
  // uncorrelated between eras
//   cb.cp()
//       .channel({"et", "mt", "tt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_sub_syst_$CHANNEL_$BIN_$ERA", "lnN", SystMap<channel, bin_id>::init
// 	       ({"mt"}, {11},  1.025) //w
// 	       ({"mt"}, {12},  1.045) //ztt
// 	       ({"mt"}, {13},  1.03) //tt
// 	       ({"mt"}, {14},  1.02) //ss
// 	       ({"mt"}, {15},  1.04) //zll
// 	       ({"mt"}, {16},  1.035) //misc
// 	       ({"mt"}, {20},  1.045) //emb
// 	       ({"mt"}, {21},  1.024) //ss
// 	       ({"mt"}, {300}, 1.035) //incl
// 	       ({"et"}, {11},  1.02) //w
// 	       ({"et"}, {12},  1.04) //ztt
// 	       ({"et"}, {13},  1.03) //tt
// 	       ({"et"}, {14},  1.02) //ss
// 	       ({"et"}, {15},  1.04) //zll
// 	       ({"et"}, {16},  1.035) //misc
// 	       ({"et"}, {20},  1.04) //emb
// 	       ({"et"}, {21},  1.023) //ff
// 	       ({"et"}, {300}, 1.035) //incl
// 	       ({"tt"}, {12},  1.035) //ztt
// 	       ({"tt"}, {16},  1.03) //misc
// 	       ({"tt"}, {17},  1.02) //noniso
// 	       ({"tt"}, {20},  1.035) //emb
// 	       ({"tt"}, {21},  1.02) //ff
// 	       ({"tt"}, {300}, 1.03) //incl
// 	       );

//   // ggH and qqH categories
//   cb.cp()
//       .channel({"et", "mt", "tt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_sub_syst_$CHANNEL_ggH_$ERA", "lnN", SystMap<channel, bin_id>::init
// 	       ({"mt"}, {1, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110},  1.04)
// 	       ({"et"}, {1, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110},  1.04)
// 	       ({"tt"}, {1, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110},  1.03)
// 	       );

//   cb.cp()
//       .channel({"et", "mt", "tt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_sub_syst_$CHANNEL_qqH_$ERA", "lnN", SystMap<channel, bin_id>::init
// 	       ({"mt"}, {2, 200, 201, 202, 203},  1.04)
// 	       ({"et"}, {2, 200, 201, 202, 203},  1.035)
// 	       ({"tt"}, {2, 200, 201, 202, 203},  1.03)
// 	       );


}
} // namespace ch