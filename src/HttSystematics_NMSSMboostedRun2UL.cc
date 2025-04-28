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

  void AddRun2BoostedSystematics(CombineHarvester &cb, bool jetfakes, bool embedding, int era) {

  // ##########################################################################
  // Define groups of processes
  // ##########################################################################

  // Signal processes
      // NMSSM
    
  std::vector<std::string> signals = {
    "NMSSM_Ytt", "NMSSM_Ybb"
  };

  // All processes being taken from simulation
  // FIXME: Adapt for fake factor and embedding
  std::vector<std::string> mc_processes =
      JoinStr({
              signals,
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

//   cb.cp()
//     .channel({"et"})
//     .process(mc_processes)
//     .AddSyst(cb, "CMS_eff_trigger_boosted_et_$ERA", "shape", SystMap<>::init(1.00));
//   cb.cp()
//       .channel({"mt"})
//       .process(mc_processes)
//       .AddSyst(cb, "CMS_eff_trigger_boosted_mt_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"tt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_eff_trigger_boosted_tt_$ERA", "shape", SystMap<>::init(1.00));

  cb.cp()
    .channel({"et"})
    .process(mc_processes)
    .AddSyst(cb, "CMS_eff_trigger_boosted_et_$ERA", "lnN", SystMap<>::init(1.02));
  cb.cp()
      .channel({"mt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_eff_trigger_boosted_mt_$ERA", "lnN", SystMap<>::init(1.02));


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

  // Electron Iso
  cb.cp()
      .channel({"et"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_iso_e", "lnN", SystMap<>::init(1.02));

  // Muon ID
  cb.cp()
      .channel({"mt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_eff_m", "lnN", SystMap<>::init(1.02));

  // Muon Iso
  cb.cp()
      .channel({"mt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_iso_m", "lnN", SystMap<>::init(1.02));


  std::string tauIDptbins[3] = {"40-500", "500-1000", "1000-Inf"};
  std::string tauIDdmbins[3] = {"0", "1", "10"};

  // Common component acting on MC

  // 3% in Tau ID SF with different anti-l fake WP
  cb.cp()
      .channel({"mt", "tt"})
      .process(JoinStr({signals, {"EMB", "ZTT_NLO", "TTT", "TTL", "VVT", "VVL", "STT", "STL", "ggH_tt", "qqH_tt", "VH_tt", "ggH_bb", "qqH_bb", "VH_bb"}}))
      .AddSyst(cb, "CMS_eff_boosted_t_wp_$ERA", "lnN", SystMap<>::init(1.03));
  
  // Tau ID: et and mt with 1 real tau
      
  for (auto tauIDbin : tauIDptbins){ //first part correlated between channels for IDvsJets
    cb.cp()
        .channel({"et", "mt"})
        .process(JoinStr({signals, {"ZTT_NLO", "TTT", "TTL", "VVT", "VVL", "STT", "STL", "ggH_tt", "qqH_tt", "VH_tt", "ggH_bb", "qqH_bb", "VH_bb"}}))
        .AddSyst(cb, "CMS_eff_boosted_t_"+tauIDbin+"_$ERA", "shape", SystMap<>::init(1.0));
  }
  cb.cp() //second part uncorrelated between channels for IDvsLep
      .channel({"et", "mt"})
      .process(JoinStr({signals, {"ZTT_NLO", "TTT", "TTL", "VVT", "VVL", "STT", "STL", "ggH_tt", "qqH_tt", "VH_tt", "ggH_bb", "qqH_bb", "VH_bb"}}))
      .AddSyst(cb, "CMS_eff_boosted_t_$CHANNEL_$ERA", "lnN", SystMap<>::init(1.01));

  // Tau ID: tt with 2 real taus
  for (auto tauIDbin : tauIDdmbins){
    cb.cp()
        .channel({"tt"})
        .process(JoinStr({signals, {"ZTT_NLO", "TTT", "TTL", "VVT", "VVL", "STT", "STL", "ggH_tt", "qqH_tt", "VH_tt", "ggH_bb", "qqH_bb", "VH_bb"}}))
        .AddSyst(cb, "CMS_eff_boosted_t_dm"+tauIDbin+"_$ERA", "shape", SystMap<>::init(1.0));
  }
  cb.cp()
      .channel({"tt"})
      .process(JoinStr({signals, {"ZTT_NLO", "TTT", "TTL", "VVT", "VVL", "STT", "STL", "ggH_tt", "qqH_tt", "VH_tt", "ggH_bb", "qqH_bb", "VH_bb"}}))
      .AddSyst(cb, "CMS_eff_boosted_t_$CHANNEL_$ERA", "lnN", SystMap<>::init(1.014));

  // Tau ID: tt with 1 real taus and 1 jet fake 
  // NEEDED?
  cb.cp()
      .channel({"tt"})
      .process({"W", "ZJ_NLO", "TTJ", "VVJ", "STJ"})
      .AddSyst(cb, "CMS_eff_boosted_t_$ERA", "lnN", SystMap<>::init(1.06));

  cb.cp()
      .channel({"tt"})
      .process({"W", "ZJ_NLO", "TTJ", "VVJ", "STJ"})
      .AddSyst(cb, "CMS_eff_boosted_t_$CHANNEL_$ERA", "lnN", SystMap<>::init(1.02));

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
  cb.cp()
      .channel({"et"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_res_e", "shape", SystMap<>::init(1.00));

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
      .AddSyst(cb, "CMS_scale_boosted_t_1prong_$ERA", "shape", SystMap<>::init(1.0));

  cb.cp()
      .channel({"et", "mt", "tt"})
      .process(JoinStr({signals, {"ZTT", "TTT", "TTL", "VVT", "VVL", "STT", "STL"}}))
      .AddSyst(cb, "CMS_scale_boosted_t_1prong1pizero_$ERA", "shape",
               SystMap<>::init(1.0));

  cb.cp()
      .channel({"et", "mt", "tt"})
      .process(JoinStr({signals, {"ZTT", "TTT", "TTL", "VVT", "VVL", "STT", "STL"}}))
      .AddSyst(cb, "CMS_scale_boosted_t_3prong_$ERA", "shape", SystMap<>::init(1.0));


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
//   cb.cp()
//       .channel({"mt"})
//       .process({"ZL_NLO"})
//       .AddSyst(cb, "CMS_ZLShape_$CHANNEL_$ERA", "shape",
//                SystMap<>::init(1.00));

//   cb.cp()
//       .channel({"et"})
//       .process({"ZL_NLO"})
//       .AddSyst(cb, "CMS_ZLShape_$CHANNEL_1prong_barrel_$ERA", "shape",
//                SystMap<>::init(1.00));

//   cb.cp()
//       .channel({"et"})
//       .process({"ZL_NLO"})
//       .AddSyst(cb, "CMS_ZLShape_$CHANNEL_1prong1pizero_barrel_$ERA", "shape",
//                SystMap<>::init(1.00));
//   cb.cp()
//       .channel({"et"})
//       .process({"ZL_NLO"})
//       .AddSyst(cb, "CMS_ZLShape_$CHANNEL_1prong_endcap_$ERA", "shape",
//                SystMap<>::init(1.00));

//   cb.cp()
//       .channel({"et"})
//       .process({"ZL_NLO"})
//       .AddSyst(cb, "CMS_ZLShape_$CHANNEL_1prong1pizero_endcap_$ERA", "shape",
//                SystMap<>::init(1.00));

  // Electron fakes ID

  cb.cp()
      .channel({"et"})
      .process({"ZL_NLO"})
      .AddSyst(cb, "CMS_fake_boosted_e_Barrel_$ERA", "shape",
               SystMap<>::init(1.00));
  cb.cp()
      .channel({"et"})
      .process({"ZL_NLO"})
      .AddSyst(cb, "CMS_fake_boosted_e_Endcap_$ERA", "shape",
               SystMap<>::init(1.00));

  // Muon fakes ID

  cb.cp()
      .channel({"mt"})
      .process({"ZL_NLO"})
      .AddSyst(cb, "CMS_fake_boosted_m_Wheel1_$ERA", "shape",
               SystMap<>::init(1.00));
  cb.cp()
      .channel({"mt"})
      .process({"ZL_NLO"})
      .AddSyst(cb, "CMS_fake_boosted_m_Wheel2_$ERA", "shape",
               SystMap<>::init(1.00));
  cb.cp()
      .channel({"mt"})
      .process({"ZL_NLO"})
      .AddSyst(cb, "CMS_fake_boosted_m_Wheel3_$ERA", "shape",
               SystMap<>::init(1.00));
  cb.cp()
      .channel({"mt"})
      .process({"ZL_NLO"})
      .AddSyst(cb, "CMS_fake_boosted_m_Wheel4_$ERA", "shape",
               SystMap<>::init(1.00));
  cb.cp()
      .channel({"mt"})
      .process({"ZL_NLO"})
      .AddSyst(cb, "CMS_fake_boosted_m_Wheel5_$ERA", "shape",
               SystMap<>::init(1.00));

  // ##########################################################################
  // Uncertainty: Jet to tau fakes
  // References:
  // Notes:
  // - FIXME: Adapt for fake factor and embedding
  // - FIXME: References?
  // ##########################################################################

//   cb.cp()
//       .channel({"et", "mt", "tt"})
//       .process({"W", "TTJ", "ZJ_NLO", "VVJ", "STJ"})
//       .AddSyst(cb, "CMS_htt_boosted_fake_j_$ERA", "shape", SystMap<>::init(1.00));

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
      .AddSyst(cb, "CMS_boostedtau_QCDFFslopeUnc_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_boostedtau_QCDFFnormUnc_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_boostedtau_QCDFFmcSubUnc_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_boostedtau_QCDSubleadingFFslopeUnc_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_boostedtau_QCDSubleadingFFnormUnc_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_boostedtau_QCDSubleadingFFmcSubUnc_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));


  // W shape stat.
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_boostedtau_WjetsFFslopeUnc_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_boostedtau_WjetsFFnormUnc_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_boostedtau_WjetsFFmcSubUnc_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));

  // ttbar shape stat.
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_boostedtau_ttbarFFslopeUnc_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_boostedtau_ttbarFFnormUnc_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_boostedtau_ttbarSubleadingFFslopeUnc_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_boostedtau_ttbarSubleadingFFnormUnc_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));


  // Process fraction variations
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_boostedtau_fracQCDUnc_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_boostedtau_fracWjetsUnc_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_boostedtau_fracTTbarUnc_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_boostedtau_fracQCDSubleadingUnc_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_boostedtau_fracWjetsSubleadingUnc_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_boostedtau_fracTTbarSubleadingUnc_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));

  // Shape syst. of different contributions (QCD/W/tt)
  // uncorrelated between eras
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_boostedtau_QCDClosureLeadingLepPtCorr_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_boostedtau_QCDClosureSubleadingLepPtCorr_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_boostedtau_QCDSubleadingClosureLeadingLepPtCorr_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_boostedtau_QCDClosureSubleadingTauMassCorr_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_boostedtau_QCDClosureLeadingTauMassCorr_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_boostedtau_QCDSubleadingClosureSubleadingTauMassCorr_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"et", "mt", "tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_boostedtau_QCDDRtoSRCorr_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_boostedtau_QCDSubleadingDRtoSRCorr_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));    
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_boostedtau_WjetsClosureLeadingLepPtCorr_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_boostedtau_WjetsClosureSubleadingTauMassCorr_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_boostedtau_WjetsDRtoSRCorr_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_boostedtau_ttbarClosureLeadingLepPtCorr_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"et", "mt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_boostedtau_ttbarClosureSubleadingTauMassCorr_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_boostedtau_ttbarClosureSubleadingLepPtCorr_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_boostedtau_ttbarClosureLeadingTauMassCorr_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_boostedtau_ttbarSubleadingClosureLeadingLepPtCorr_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
  cb.cp()
      .channel({"tt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_boostedtau_ttbarSubleadingClosureSubleadingTauMassCorr_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));


}
} // namespace ch