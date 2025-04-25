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

  void AddSMRun2Systematics(CombineHarvester &cb, bool jetfakes, bool embedding, bool regional_jec, bool ggh_wg1, string era) {

  // ##########################################################################
  // Define groups of processes
  // ##########################################################################
  std::vector<std::string> signals_WH = {"WH_htt_plus", "WH_htt_minus", "WH_hww_plus", "WH_hww_minus"};
  std::vector<std::string> mc_processes =
      JoinStr({
              signals_WH,
              {"WH_hww", "ZH_hww"},
              {"ggZZ", "ggH", "qqH", "ttH", "ggZH", "ZH", "VVV", "rem_VV", "rem_ttbar", "WZ", "Wjets", "DY", "ZZ", "TT"}
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
  if (era == "2016preVFP") {
      lumi_unc = 1.010;
      lumi_unc_corr = 1.006;
  } 
  else if (era == "2016postVFP") {
      lumi_unc = 1.010;
      lumi_unc_corr = 1.006;
  } 
  else if (era == "2017") {
      lumi_unc = 1.020;
      lumi_unc_corr = 1.009;
      lumi_unc_1718 = 1.006;
  } else if (era == "2018") {
      lumi_unc = 1.015;
      lumi_unc_corr = 1.020;
      lumi_unc_1718 = 1.002;
  }
  cb.cp()
      .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
      .process(mc_processes)
      .AddSyst(cb, "lumi_$ERA", "lnN", SystMap<>::init(lumi_unc));
  cb.cp()
      .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
      .process(mc_processes)
      .AddSyst(cb, "lumi_13TeV", "lnN", SystMap<>::init(lumi_unc_corr));
  cb.cp()
      .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
      .process(mc_processes)
      .AddSyst(cb, "lumi_13TeV_1718", "lnN", SystMap<>::init(lumi_unc_1718));

  // ##########################################################################
  // Uncertainty: Prefiring
  // References:
  // - "https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1ECALPrefiringWeightRecipe"
  // Notes:
  // - FIXME: assumed as uncorrelated accross the years for now, what is the recommendation?
  // ##########################################################################
  if (era != "2018") {
  cb.cp()
      .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_prefiring", "shape", SystMap<>::init(1.00));
  }
// ##########################################################################
  // Uncertainty: pileup
  // ##########################################################################
    cb.cp()
        .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_pileup_$ERA", "shape", SystMap<>::init(1.00));
  // ##########################################################################
  // Uncertainty: Trigger efficiency
  // References:
  // Notes:
  // - FIXME: References?
  // ##########################################################################

  cb.cp()
    .channel({"emt", "llt", "met", "ett", "ltt"})
    .process(mc_processes)
    .AddSyst(cb, "CMS_eff_e_trigger_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"emt", "llt", "met", "mtt", "mmt", "ltt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_eff_m_trigger_$ERA", "shape", SystMap<>::init(1.00));
          
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
      .channel({"emt", "llt", "met", "mmt", "ett", "mtt", "ltt", "ltt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_eff_t_wp_$ERA", "lnN", SystMap<>::init(1.03));

  std::string tauIDptbins[4] = {"30-35", "35-40", "40-500", "500-1000"};
  std::string tauIDdmbins[4] = {"0", "1", "10", "11"};

  // Common component acting on MC
  
  // Electron ID
//   cb.cp()
//       .channel({"emt", "llt", "met", "ett", "ltt"})
//       .process(mc_processes)
//       .AddSyst(cb, "CMS_eff_e", "lnN", SystMap<>::init(1.02));
    cb.cp()
    .channel({"emt", "llt", "met", "ett", "ltt"})
    .process(mc_processes)
    .AddSyst(cb, "CMS_eff_e", "shape", SystMap<>::init(1.0));
// Electron Iso
  cb.cp()
      .channel({"emt", "llt", "met", "ett", "ltt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_eff_e_iso", "lnN", SystMap<>::init(1.005));
// Electron Reco
cb.cp()
.channel({"emt", "llt", "met", "ett", "ltt"})
.process(mc_processes)
.AddSyst(cb, "CMS_eff_e_reco", "shape", SystMap<>::init(1.0));
  // Muon ID
  cb.cp()
      .channel({"emt", "llt", "met", "mmt", "mtt", "ltt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_eff_m", "lnN", SystMap<>::init(1.02));

  // Tau ID: et and mt with 1 real tau
  for (auto tauIDbin : tauIDptbins){ //first part correlated between channels for IDvsJets
    cb.cp()
        .channel({"emt", "llt", "met", "mmt"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_eff_t_"+tauIDbin+"_$ERA", "shape", SystMap<>::init(1.0));
  }
//   cb.cp() //second part uncorrelated between channels for IDvsLep
//       .channel({"emt", "llt", "met", "mmt"})
//       .process(mc_processes)
//       .AddSyst(cb, "CMS_eff_t_$CHANNEL_$ERA", "lnN", SystMap<>::init(1.01));

  // Tau ID: tt with 2 real taus
  for (auto tauIDbin : tauIDdmbins){
    cb.cp()
        .channel({"ett", "mtt", "ltt"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_eff_t_DM"+tauIDbin+"_$ERA", "shape", SystMap<>::init(1.0));
  }
//   cb.cp()
//       .channel({ "mtt", "ett", "ltt"})
//       .process(mc_processes)
//       .AddSyst(cb, "CMS_eff_t_$CHANNEL_$ERA", "lnN", SystMap<>::init(1.014));

  // Tau ID: tt with 1 real taus and 1 jet fake
  cb.cp()
      .channel({ "mtt", "ett", "ltt"})
      .process({"Wjets","DY", "TT", "rem_VV"})
      .AddSyst(cb, "CMS_eff_t_$ERA", "lnN", SystMap<>::init(1.06));


  cb.cp()
      .channel({"ett", "mtt", "ltt"})
      .process({"Wjets","DY", "TT", "rem_VV"})
      .AddSyst(cb, "CMS_eff_t_$CHANNEL_$ERA", "lnN", SystMap<>::init(1.02));

  // ##########################################################################
  // Uncertainty: b-tag and mistag efficiency
  // References:
  // Notes:
  // - FIXME: btag jes
  // ##########################################################################

  cb.cp()
      .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_btag_shape_hf", "shape", SystMap<>::init(1.00));

  cb.cp()
      .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_btag_shape_hfstats1_$ERA", "shape", SystMap<>::init(1.00));

  cb.cp()
      .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_btag_shape_hfstats2_$ERA", "shape", SystMap<>::init(1.00));

  cb.cp()
      .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_btag_shape_lf", "shape", SystMap<>::init(1.00));

  cb.cp()
      .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_btag_shape_lfstats1_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_btag_shape_cferr1", "shape", SystMap<>::init(1.00));

  cb.cp()
      .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_btag_shape_cferr2", "shape", SystMap<>::init(1.00));

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
      .channel({"emt", "met", "ett", "ltt", "llt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_scale_e", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"emt", "met", "ett", "ltt", "llt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_res_e", "shape", SystMap<>::init(1.00));

//   cb.cp()
//       .channel({"emt", "met", "ett", "ltt", "llt"})
//       .process(mc_processes)
//       .AddSyst(cb, "CMS_res_e", "shape", SystMap<>::init(1.00));
      //.AddSyst(cb, "CMS_scale_mc_e", "shape", SystMap<>::init(0.71));
      
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
      .channel({ "emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
      .process(JoinStr({mc_processes}))
      .AddSyst(cb, "CMS_scale_t_DM0_$ERA", "shape", SystMap<>::init(1.0));

  cb.cp()
      .channel({ "emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
      .process(JoinStr({mc_processes}))
      .AddSyst(cb, "CMS_scale_t_DM1_$ERA", "shape",
               SystMap<>::init(1.0));

  cb.cp()
      .channel({ "emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
      .process(JoinStr({mc_processes}))
      .AddSyst(cb, "CMS_scale_t_DM10_$ERA", "shape", SystMap<>::init(1.0));

  cb.cp()
      .channel({ "emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
      .process(JoinStr({mc_processes}))
      .AddSyst(cb, "CMS_scale_t_DM11_$ERA", "shape",
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
      .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_scale_j_$ERA", "shape", SystMap<>::init(0.71));
  cb.cp()
      .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_scale_j", "shape", SystMap<>::init(0.71));
  }

  // Regional JES
  else {
    // uncorrelated between eras
    cb.cp()
        .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_Absolute_$ERA", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_BBEC1_$ERA", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_EC2_$ERA", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_HF_$ERA", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_RelativeSample_$ERA", "shape", SystMap<>::init(1.00));
    // correlated between eras
    cb.cp()
        .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_Absolute", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_BBEC1", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_EC2", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_HF", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_FlavorQCD", "shape", SystMap<>::init(1.00));

    cb.cp()
        .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
        .process(mc_processes)
        .AddSyst(cb, "CMS_scale_j_RelativeBal", "shape", SystMap<>::init(1.00));
  }

  // JER
  cb.cp()
      .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
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
      .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
      .process(mc_processes)  //Z and W processes are only included due to the EWK fraction. Make sure that there is no contribution to the shift from the DY or Wjets samples.
      .AddSyst(cb, "CMS_scale_met_$ERA", "shape", SystMap<>::init(1.00));
//   cb.cp()
    //   .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
    //   .process(mc_processes)
    //   .AddSyst(cb, "CMS_htt_boson_scale_met_$ERA", "shape", SystMap<>::init(1.00));
  cb.cp()
      .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
      .process(mc_processes)
      .AddSyst(cb, "CMS_res_met_$ERA", "shape", SystMap<>::init(1.00));

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
      .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
      .process({"ggZZ", "rem_VV", "WZ", "ZZ" })
      .AddSyst(cb, "cross_section_VV", "lnN", SystMap<>::init(1.075));

// VVV
  cb.cp()
      .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
      .process({"VVV" })
      .AddSyst(cb, "cross_section_VVV", "lnN", SystMap<>::init(1.10));

  // TT
  cb.cp()
      .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
      .process({"TTT", "TTL", "TTJ", "TT", "rem_ttbar"})
      .AddSyst(cb, "cross_section_TTV", "lnN", SystMap<>::init(1.06));

  // W
  cb.cp()
      .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
      .process({"Wjets"})
      .AddSyst(cb, "CMS_htt_wjXsec", "lnN", SystMap<>::init(1.04));

  // Z
  cb.cp()
      .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
      .process({"ZTT", "ZL", "ZJ"})
      .AddSyst(cb, "CMS_htt_zjXsec", "lnN", SystMap<>::init(1.02));

  // QCD
  cb.cp()
      .channel({"emt", "llt", "met", "ett", "ltt"})
      .process({"QCD"})
      .AddSyst(cb, "CMS_ExtrapSSOS_$CHANNEL_$ERA", "lnN", SystMap<>::init(1.05));
  cb.cp()
      .channel({"emt", "llt", "met", "mtt", "mmt", "ltt"})
      .process({"QCD"})
      .AddSyst(cb, "CMS_ExtrapSSOS_$CHANNEL_$ERA", "lnN", SystMap<>::init(1.03));
  cb.cp()
      .channel({ "mtt", "ett", "ltt"})
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

  if (era == "2016preVFP" or era == "2016postVFP") {
      cb.cp()
          .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
          .process({"DY"})
          .AddSyst(cb, "CMS_htt_dyShape_$ERA", "shape", SystMap<>::init(0.10));
  } else {
      cb.cp()
          .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
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
      .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
      .process({"TT", "rem_ttbar"})
      .AddSyst(cb, "top_pt_reweighting", "shape", SystMap<>::init(1.00));

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
  //    .channel({"emt", "llt", "met", "ett"})
  //    .process({"ZL"})
  //    .AddSyst(cb, "CMS_fake_e_$ERA", "lnN", SystMap<>::init(1.15));

  cb.cp()
      .channel({"emt", "llt", "met", "ett", "ltt"})
      .process({"DY"})
      .AddSyst(cb, "CMS_fake_e_BA_$ERA", "shape",
               SystMap<>::init(1.00));
  cb.cp()
      .channel({"emt", "llt", "met", "ett", "ltt"})
      .process({"DY"})
      .AddSyst(cb, "CMS_fake_e_EC_$ERA", "shape",
               SystMap<>::init(1.00));

  // Muon fakes
  //cb.cp()
  //    .channel({"emt", "llt", "met", "mtt", "mmt"})
  //    .process({"ZL"})
  //    .AddSyst(cb, "CMS_fake_m_$ERA", "lnN", SystMap<>::init(1.25));

  cb.cp()
      .channel({"emt", "llt", "met","mmt"})
      .process({"DY"})
      .AddSyst(cb, "CMS_fake_m_WH1_$ERA", "shape",
               SystMap<>::init(1.00));
  cb.cp()
      .channel({"emt", "llt", "met","mmt"})
      .process({"DY"})
      .AddSyst(cb, "CMS_fake_m_WH2_$ERA", "shape",
               SystMap<>::init(1.00));
  cb.cp()
      .channel({"emt", "llt", "met","mmt"})
      .process({"DY"})
      .AddSyst(cb, "CMS_fake_m_WH3_$ERA", "shape",
               SystMap<>::init(1.00));
  cb.cp()
      .channel({"emt", "llt", "met","mmt"})
      .process({"DY"})
      .AddSyst(cb, "CMS_fake_m_WH4_$ERA", "shape",
               SystMap<>::init(1.00));
  cb.cp()
      .channel({"emt", "llt", "met","mmt"})
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
//       .channel({ "emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
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
  // Uncertainty due to the missing ggZH(tautau)-UL sample
//   cb.cp()
//       .channel({"emt", "llt", "met", "mmt"})
//       .process({"ZH"})
//       .AddSyst(cb, "ggZH_htt_yield", "lnN", SystMap<>::init(1.2));
//   cb.cp()
//       .channel({"ett", "mtt", "ltt"})
//       .process({"ZH"})
//       .AddSyst(cb, "ggZH_htt_yield", "lnN", SystMap<>::init(1.24));
    // Uncertainty on branching ratio for HTT at 125 GeV
  cb.cp()
      .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
      .process({"ggH", "qqH", "ttH", "ZH", "WH_htt_plus", "WH_htt_minus"})
      .AddSyst(cb, "BR_htt", "lnN", SystMap<>::init(1.0117));
  cb.cp()
      .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
      .process({"ggH", "qqH", "ttH", "ZH", "WH_htt_plus", "WH_htt_minus"})
      .AddSyst(cb, "BR_htt_mq", "lnN", SystMap<>::init(1.0099));
  cb.cp()
      .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
      .process({"ggH", "qqH", "ttH", "ZH", "WH_htt_plus", "WH_htt_minus"})
      .AddSyst(cb, "BR_htt_alphas", "lnN", SystMap<>::init(1.0061));
  // QCD scale

cb.cp()
    .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
    .process({"WH_htt_plus", "WH_htt_minus","WH_hww_plus", "WH_hww_minus", "ZH"})
    .AddSyst(cb, "QCD_ren_scale_VH", "shape", SystMap<>::init(1.0));
cb.cp()
    .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
    .process({"WH_htt_plus", "WH_htt_minus","WH_hww_plus", "WH_hww_minus", "ZH"})
    .AddSyst(cb, "QCD_fac_scale_VH", "shape", SystMap<>::init(1.0));
cb.cp()
    .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
    .process({"ggZH"})
    .AddSyst(cb, "QCD_ren_scale_ggZH", "shape", SystMap<>::init(1.0));
cb.cp()
    .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
    .process({"ggZH"})
    .AddSyst(cb, "QCD_fac_scale_ggZH", "shape", SystMap<>::init(1.0));
cb.cp()
    .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
    .process({"WZ", "ZZ"})
    .AddSyst(cb, "QCD_ren_scale_VV", "shape", SystMap<>::init(1.0));
cb.cp()
    .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
    .process({"WZ", "ZZ"})
    .AddSyst(cb, "QCD_fac_scale_VV", "shape", SystMap<>::init(1.0));
cb.cp()
    .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
    .process({"rem_ttbar"})
    .AddSyst(cb, "QCD_ren_scale_TTV", "shape", SystMap<>::init(1.0));
cb.cp()
    .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
    .process({"rem_ttbar"})
    .AddSyst(cb, "QCD_fac_scale_TTV", "shape", SystMap<>::init(1.0));
cb.cp()
    .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
    .process({"VVV"})
    .AddSyst(cb, "QCD_ren_scale_VVV", "shape", SystMap<>::init(1.0));
cb.cp()
    .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
    .process({"VVV"})
    .AddSyst(cb, "QCD_fac_scale_VVV", "shape", SystMap<>::init(1.0));
  // PDF
cb.cp()
      .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
     .process({"WH_htt_plus", "WH_htt_minus","WH_hww_plus", "WH_hww_minus"})
      .AddSyst(cb, "pdf_WH", "shape", SystMap<>::init(1.0));
// xsec normalization unc.
cb.cp()
    .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
    .process({"WH_htt_plus", "WH_htt_minus","WH_hww_plus", "WH_hww_minus"})
    .AddSyst(cb, "cross_section_WH", "lnN", SystMap<>::init(1.0065));
//   cb.cp()
//       .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
//      .process({"ggH"})
//       .AddSyst(cb, "pdf_Higgs_gg", "lnN", SystMap<>::init(1.032));
//   cb.cp()
//       .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
//      .process({"qqH"})
//       .AddSyst(cb, "pdf_Higgs_qqbar", "lnN", SystMap<>::init(1.021));
//   cb.cp()
//       .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
//       .process(JoinStr({{"ggZH","ZH"}}))
//       .AddSyst(cb, "pdf_Higgs_VH", "lnN", SystMap<>::init(1.013));
//   cb.cp()
//       .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
//       .process(signals_WH)
//       .AddSyst(cb, "pdf_Higgs_VH", "lnN", SystMap<>::init(1.018));
//   cb.cp()
//       .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
//       .process({"ttH"})
//       .AddSyst(cb, "pdf_Higgs_ttH", "lnN", SystMap<>::init(1.036));

  // ##########################################################################
  // Uncertainty: Embedded events
  // References:
  // - https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauTauEmbeddingSamples2016
  // Notes:
  // ##########################################################################

  // Embedded Normalization: No Lumi, Zjxsec information used, instead derived from data using dimuon selection efficiency
  cb.cp()
      .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
      .process({"EMB"})
      .AddSyst(cb, "CMS_htt_doublemutrg_$ERA", "lnN", SystMap<>::init(1.04));

  // TTbar contamination in embedded events: 10% shape uncertainty of assumed ttbar->tautau event shape
  cb.cp()
    .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
    .process({"EMB"})
    .AddSyst(cb, "CMS_htt_emb_ttbar_$ERA", "shape", SystMap<>::init(1.00));

  // Uncertainty of hadronic tau track efficiency correction
  // uncorrelated between eras
  cb.cp()
    .channel({ "emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
    .process({"EMB"})
    .AddSyst(cb, "CMS_3ProngEff_$ERA", "shape", SystMap<>::init(0.71));

  cb.cp()
    .channel({ "emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
    .process({"EMB"})
    .AddSyst(cb, "CMS_1ProngPi0Eff_$ERA", "shape", SystMap<>::init(0.71));
  // correlated between eras
  cb.cp()
    .channel({ "emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
    .process({"EMB"})
    .AddSyst(cb, "CMS_3ProngEff", "shape", SystMap<>::init(0.71));

  cb.cp()
    .channel({ "emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
    .process({"EMB"})
    .AddSyst(cb, "CMS_1ProngPi0Eff", "shape", SystMap<>::init(0.71));

  // ##########################################################################
  // Uncertainty: Jet fakes
  // ##########################################################################
  std::cout<<"oleeeeee"<<std::endl;
    cb.cp()
    .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
    .process({"jetFakes"})
    .AddSyst(cb, "CMS_fake_stat_$ERA", "shape", SystMap<>::init(1.0));  
    cb.cp()
    .channel({"emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
    .process({"jetFakes"})
    .AddSyst(cb, "CMS_fake_irredbkg_$ERA", "shape", SystMap<>::init(1.0));
    cb.cp()
    .channel({"ltt", "ett", "mtt"})
    .process({"jetFakes"})
    .AddSyst(cb, "CMS_fake_metltt_$ERA", "shape", SystMap<>::init(1.0));
    cb.cp()
    .channel({"llt", "emt", "met", "mmt"})
    .process({"jetFakes"})
    .AddSyst(cb, "CMS_fake_metllt_$ERA", "shape", SystMap<>::init(1.0));
    cb.cp()
    .channel({"ltt", "ett", "mtt"})
    .process({"jetFakes"})
    .AddSyst(cb, "CMS_fake_pt1ltt_$ERA", "shape", SystMap<>::init(1.0));
    cb.cp()
    .channel({"llt", "emt", "met", "mmt"})
    .process({"jetFakes"})
    .AddSyst(cb, "CMS_fake_pt1llt_$ERA", "shape", SystMap<>::init(1.0));
    //////////////
    //  cb.cp()
    // if (era == "2017") {
    // cb.cp()
    // .channel({"ltt", "ett", "mtt"})
    // .process({"jetFakes"})
    // .AddSyst(cb, "CMS_fake_nnscore_$ERA", "shape", SystMap<>::init(1.0));
    // }
    
    //   .channel({"llt", "ltt"})
    //   .process({"jetFakes"})
    //   .AddSyst(cb, "CMS_fake_systdet_$ERA", "shape", SystMap<>::init(1.0));
    //////////////
  cb.cp()
      .channel({"mtt", "ett", "ltt","emt", "llt", "met", "mmt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_fake_comp_$CHANNEL", "lnN", SystMap<>::init(1.20));  
  cb.cp()
      .channel({"emt", "llt", "met", "mmt"})
      .process({"jetFakes"})
      .AddSyst(cb, "CMS_fake_mtcut", "lnN", SystMap<>::init(1.03));  
  // ##########################################################################
  // Uncertainty: Jet fakes
  // References:
  // - https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauJet2TauFakes
  // Notes:
  // - FIXME: add 2017 norm uncertainties, and properly correlate across years
  // ##########################################################################
//   cb.cp()
//       .channel({"mtt", "ett", "ltt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_njet", "lnN", SystMap<>::init(1.10));
//   cb.cp()
//       .channel({"emt", "llt", "met", "mmt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_njet_bkgcomp_$CHANNEL", "lnN", SystMap<>::init(1.15));

  // QCD shape stat.
//   cb.cp()
//       .channel({ "emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_qcd_njet0_morphed_stat_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
//   cb.cp()
//       .channel({ "emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_qcd_njet1_morphed_stat_$CHANNEL_$ERA", "shape", SystMap<>::init(1.00));
//   cb.cp()
//       .channel({ "emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
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
//       .channel({ "emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_qcd_mc_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({"et", "mt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_frac_w_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));

      
//   // Shape syst. of different contributions (QCD/W/tt)
//   // uncorrelated between eras
//   cb.cp()
//       .channel({ "emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_qcd_mvis_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({ "emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_qcd_mvis_osss_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({ "emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_corr_qcd_mvis_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({ "emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
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
//       .channel({ "mtt", "ett", "ltt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_qcd_tau2_pt_0jet_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({ "mtt", "ett", "ltt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_corr_qcd_tau2_pt_0jet_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({ "mtt", "ett", "ltt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_qcd_tau2_pt_1jet_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({ "mtt", "ett", "ltt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_corr_qcd_tau2_pt_1jet_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({ "mtt", "ett", "ltt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_w_syst_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({ "mtt", "ett", "ltt"})
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
//       .channel({ "emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_norm_stat_$CHANNEL_$BIN_$ERA", "lnN", SystMap<channel, bin_id>::init
// 	       ({"emt", "llt", "met", "mtt", "mmt"}, {11},  1.04) //w
// 	       ({"emt", "llt", "met", "mtt", "mmt"}, {12},  1.052) //ztt
// 	       ({"emt", "llt", "met", "mtt", "mmt"}, {13},  1.051) //tt
// 	       ({"emt", "llt", "met", "mtt", "mmt"}, {14},  1.047) //ss
// 	       ({"emt", "llt", "met", "mtt", "mmt"}, {15},  1.04) //zll
// 	       ({"emt", "llt", "met", "mtt", "mmt"}, {16},  1.059) //misc
// 	       ({"emt", "llt", "met", "mtt", "mmt"}, {20},  1.052) //emb
// 	       ({"emt", "llt", "met", "mtt", "mmt"}, {21},  1.047) //ff
// 	       ({"emt", "llt", "met", "mtt", "mmt"}, {300}, 1.037) //incl
// 	       ({"emt", "llt", "met", "ett"}, {11},  1.066) //w
// 	       ({"emt", "llt", "met", "ett"}, {12},  1.095) //ztt
// 	       ({"emt", "llt", "met", "ett"}, {13},  1.083) //tt
// 	       ({"emt", "llt", "met", "ett"}, {14},  1.054) //ss
// 	       ({"emt", "llt", "met", "ett"}, {15},  1.095) //zll
// 	       ({"emt", "llt", "met", "ett"}, {16},  1.107) //misc
// 	       ({"emt", "llt", "met", "ett"}, {20},  1.095) //emb
// 	       ({"emt", "llt", "met", "ett"}, {21},  1.066) //ff
// 	       ({"emt", "llt", "met", "ett"}, {300}, 1.065) //incl
// 	       ({ "mtt", "ett", "ltt"}, {12},  1.049) //ztt
// 	       ({ "mtt", "ett", "ltt"}, {16},  1.028) //misc
// 	       ({ "mtt", "ett", "ltt"}, {17},  1.041) //noniso
// 	       ({ "mtt", "ett", "ltt"}, {20},  1.049) //emb
// 	       ({ "mtt", "ett", "ltt"}, {21},  1.041) //ff
// 	       ({ "mtt", "ett", "ltt"}, {300}, 1.041) //incl
// 	       );
//   // ggH and qqH categories
//   cb.cp()
//       .channel({ "emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_norm_stat_$CHANNEL_ggH_$ERA", "lnN", SystMap<channel, bin_id>::init
// 	       ({"emt", "llt", "met", "mtt", "mmt"}, {1, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110},  1.049)
// 	       ({"emt", "llt", "met", "ett"}, {1, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110},  1.074)
// 	       ({ "mtt", "ett", "ltt"}, {1, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110},  1.041)
// 	       );

//   cb.cp()
//       .channel({ "emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_norm_stat_$CHANNEL_qqH_$ERA", "lnN", SystMap<channel, bin_id>::init
// 	       ({"emt", "llt", "met", "mtt", "mmt"}, {2, 200, 201, 202, 203},  1.068)
// 	       ({"emt", "llt", "met", "ett"}, {2, 200, 201, 202, 203},  1.112)
// 	       ({ "mtt", "ett", "ltt"}, {2, 200, 201, 202, 203},  1.052)
// 	       );
    
//   // Syst. norm: Bin-correlated
//   // uncorrelated between eras

//   cb.cp()
//       .channel({ "emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_jetbinned_stat_0jet_norm_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({ "emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_jetbinned_stat_1jet_norm_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));
//   cb.cp()
//       .channel({ "emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_jetbinned_stat_2jet_norm_$CHANNEL_$ERA", "shape", SystMap<>::init(1.0));

//   /*
//   cb.cp()
//       .channel({ "emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_norm_syst_$CHANNEL_$ERA", "lnN", SystMap<channel, bin_id>::init
// 	       ({"emt", "llt", "met", "mtt", "mmt"}, {1},     1.069) //ggh
// 	       ({"emt", "llt", "met", "mtt", "mmt"}, {100},   1.069) //ggh
// 	       ({"emt", "llt", "met", "mtt", "mmt"}, {101},   1.069) //ggh
// 	       ({"emt", "llt", "met", "mtt", "mmt"}, {102},   1.069) //ggh
// 	       ({"emt", "llt", "met", "mtt", "mmt"}, {103},   1.069) //ggh
// 	       ({"emt", "llt", "met", "mtt", "mmt"}, {104},   1.069) //ggh
//                ({"emt", "llt", "met", "mtt", "mmt"}, {105},   1.069) //ggh
//                ({"emt", "llt", "met", "mtt", "mmt"}, {106},   1.069) //ggh
//                ({"emt", "llt", "met", "mtt", "mmt"}, {107},   1.069) //ggh
//                ({"emt", "llt", "met", "mtt", "mmt"}, {108},   1.069) //ggh
//                ({"emt", "llt", "met", "mtt", "mmt"}, {109},   1.069) //ggh
//                ({"emt", "llt", "met", "mtt", "mmt"}, {110},   1.069) //ggh
// 	       ({"emt", "llt", "met", "mtt", "mmt"}, {2},     1.058) //qqh
// 	       ({"emt", "llt", "met", "mtt", "mmt"}, {200},   1.058) //qqh
// 	       ({"emt", "llt", "met", "mtt", "mmt"}, {201},   1.058) //qqh
// 	       ({"emt", "llt", "met", "mtt", "mmt"}, {202},   1.058) //qqh
// 	       ({"emt", "llt", "met", "mtt", "mmt"}, {203},   1.058) //qqh
// 	       ({"emt", "llt", "met", "mtt", "mmt"}, {11},  1.054) //w
// 	       ({"emt", "llt", "met", "mtt", "mmt"}, {12},  1.098) //ztt
// 	       ({"emt", "llt", "met", "mtt", "mmt"}, {13},  1.052) //tt
// 	       ({"emt", "llt", "met", "mtt", "mmt"}, {14},  1.091) //ss
// 	       ({"emt", "llt", "met", "mtt", "mmt"}, {15},  1.068) //zll
// 	       ({"emt", "llt", "met", "mtt", "mmt"}, {16},  1.091) //misc
// 	       ({"emt", "llt", "met", "mtt", "mmt"}, {20},  1.098) //emb
// 	       ({"emt", "llt", "met", "mtt", "mmt"}, {21},  1.064) //ff
// 	       ({"emt", "llt", "met", "mtt", "mmt"}, {300}, 1.059) //incl
// 	       ({"emt", "llt", "met", "ett"}, {1},     1.059) //ggh
// 	       ({"emt", "llt", "met", "ett"}, {100},   1.059) //ggh
// 	       ({"emt", "llt", "met", "ett"}, {101},   1.059) //ggh
// 	       ({"emt", "llt", "met", "ett"}, {102},   1.059) //ggh
// 	       ({"emt", "llt", "met", "ett"}, {103},   1.059) //ggh
// 	       ({"emt", "llt", "met", "ett"}, {104},   1.059) //ggh
//                ({"emt", "llt", "met", "ett"}, {105},   1.059) //ggh
//                ({"emt", "llt", "met", "ett"}, {106},   1.059) //ggh
//                ({"emt", "llt", "met", "ett"}, {107},   1.059) //ggh
//                ({"emt", "llt", "met", "ett"}, {108},   1.059) //ggh
//                ({"emt", "llt", "met", "ett"}, {109},   1.059) //ggh
//                ({"emt", "llt", "met", "ett"}, {110},   1.059) //ggh
// 	       ({"emt", "llt", "met", "ett"}, {2},     1.057) //qqh
// 	       ({"emt", "llt", "met", "ett"}, {200},   1.057) //qqh
// 	       ({"emt", "llt", "met", "ett"}, {201},   1.057) //qqh
// 	       ({"emt", "llt", "met", "ett"}, {202},   1.057) //qqh
// 	       ({"emt", "llt", "met", "ett"}, {203},   1.057) //qqh
// 	       ({"emt", "llt", "met", "ett"}, {11},  1.052) //w
// 	       ({"emt", "llt", "met", "ett"}, {12},  1.088) //ztt
// 	       ({"emt", "llt", "met", "ett"}, {13},  1.057) //tt
// 	       ({"emt", "llt", "met", "ett"}, {14},  1.064) //ss
// 	       ({"emt", "llt", "met", "ett"}, {15},  1.072) //zll
// 	       ({"emt", "llt", "met", "ett"}, {16},  1.058) //misc
// 	       ({"emt", "llt", "met", "ett"}, {20},  1.088) //ztt
// 	       ({"emt", "llt", "met", "ett"}, {21},  1.057) //ff
// 	       ({"emt", "llt", "met", "ett"}, {300}, 1.059) //incl
// 	       ({ "mtt", "ett", "ltt"}, {1},     1.096) //ggh
// 	       ({ "mtt", "ett", "ltt"}, {100},   1.096) //ggh
// 	       ({ "mtt", "ett", "ltt"}, {101},   1.096) //ggh
// 	       ({ "mtt", "ett", "ltt"}, {102},   1.096) //ggh
// 	       ({ "mtt", "ett", "ltt"}, {103},   1.096) //ggh
// 	       ({ "mtt", "ett", "ltt"}, {104},   1.096) //ggh
//                ({ "mtt", "ett", "ltt"}, {105},   1.096) //ggh
//                ({ "mtt", "ett", "ltt"}, {106},   1.096) //ggh
//                ({ "mtt", "ett", "ltt"}, {107},   1.096) //ggh
//                ({ "mtt", "ett", "ltt"}, {108},   1.096) //ggh
//                ({ "mtt", "ett", "ltt"}, {109},   1.096) //ggh
//                ({ "mtt", "ett", "ltt"}, {110},   1.096) //ggh
// 	       ({ "mtt", "ett", "ltt"}, {2},     1.095) //qqh
// 	       ({ "mtt", "ett", "ltt"}, {200},   1.095) //qqh
// 	       ({ "mtt", "ett", "ltt"}, {201},   1.095) //qqh
// 	       ({ "mtt", "ett", "ltt"}, {202},   1.095) //qqh
// 	       ({ "mtt", "ett", "ltt"}, {203},   1.095) //qqh
// 	       ({ "mtt", "ett", "ltt"}, {12},  1.095) //ztt
// 	       ({ "mtt", "ett", "ltt"}, {16},  1.11) //misc
// 	       ({ "mtt", "ett", "ltt"}, {17},  1.099) //noniso
// 	       ({ "mtt", "ett", "ltt"}, {20},  1.095) //emb
// 	       ({ "mtt", "ett", "ltt"}, {21},  1.099) //ff
// 	       ({ "mtt", "ett", "ltt"}, {300}, 1.095) //incl
// 	       );
//     */
//   // Syst. norm: Bin-dependent, correlated across years
//   // uncorrelated between eras
//   cb.cp()
//       .channel({ "emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_sub_syst_$CHANNEL_$BIN_$ERA", "lnN", SystMap<channel, bin_id>::init
// 	       ({"emt", "llt", "met", "mtt", "mmt"}, {11},  1.025) //w
// 	       ({"emt", "llt", "met", "mtt", "mmt"}, {12},  1.045) //ztt
// 	       ({"emt", "llt", "met", "mtt", "mmt"}, {13},  1.03) //tt
// 	       ({"emt", "llt", "met", "mtt", "mmt"}, {14},  1.02) //ss
// 	       ({"emt", "llt", "met", "mtt", "mmt"}, {15},  1.04) //zll
// 	       ({"emt", "llt", "met", "mtt", "mmt"}, {16},  1.035) //misc
// 	       ({"emt", "llt", "met", "mtt", "mmt"}, {20},  1.045) //emb
// 	       ({"emt", "llt", "met", "mtt", "mmt"}, {21},  1.024) //ss
// 	       ({"emt", "llt", "met", "mtt", "mmt"}, {300}, 1.035) //incl
// 	       ({"emt", "llt", "met", "ett"}, {11},  1.02) //w
// 	       ({"emt", "llt", "met", "ett"}, {12},  1.04) //ztt
// 	       ({"emt", "llt", "met", "ett"}, {13},  1.03) //tt
// 	       ({"emt", "llt", "met", "ett"}, {14},  1.02) //ss
// 	       ({"emt", "llt", "met", "ett"}, {15},  1.04) //zll
// 	       ({"emt", "llt", "met", "ett"}, {16},  1.035) //misc
// 	       ({"emt", "llt", "met", "ett"}, {20},  1.04) //emb
// 	       ({"emt", "llt", "met", "ett"}, {21},  1.023) //ff
// 	       ({"emt", "llt", "met", "ett"}, {300}, 1.035) //incl
// 	       ({ "mtt", "ett", "ltt"}, {12},  1.035) //ztt
// 	       ({ "mtt", "ett", "ltt"}, {16},  1.03) //misc
// 	       ({ "mtt", "ett", "ltt"}, {17},  1.02) //noniso
// 	       ({ "mtt", "ett", "ltt"}, {20},  1.035) //emb
// 	       ({ "mtt", "ett", "ltt"}, {21},  1.02) //ff
// 	       ({ "mtt", "ett", "ltt"}, {300}, 1.03) //incl
// 	       );

//   // ggH and qqH categories
//   cb.cp()
//       .channel({ "emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_sub_syst_$CHANNEL_ggH_$ERA", "lnN", SystMap<channel, bin_id>::init
// 	       ({"emt", "llt", "met", "mtt", "mmt"}, {1, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110},  1.04)
// 	       ({"emt", "llt", "met", "ett"}, {1, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110},  1.04)
// 	       ({ "mtt", "ett", "ltt"}, {1, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110},  1.03)
// 	       );

//   cb.cp()
//       .channel({ "emt", "llt", "met", "mmt", "mtt", "ett", "ltt"})
//       .process({"jetFakes"})
//       .AddSyst(cb, "CMS_ff_sub_syst_$CHANNEL_qqH_$ERA", "lnN", SystMap<channel, bin_id>::init
// 	       ({"emt", "llt", "met", "mtt", "mmt"}, {2, 200, 201, 202, 203},  1.04)
// 	       ({"emt", "llt", "met", "ett"}, {2, 200, 201, 202, 203},  1.035)
// 	       ({ "mtt", "ett", "ltt"}, {2, 200, 201, 202, 203},  1.03)
// 	       );


}
} // namespace ch
