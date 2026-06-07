#include "CombineHarvester/SMRun2Legacy/interface/HttSystematics_SMRun2.h"
#include "CombineHarvester/CombineTools/interface/Process.h"
#include "CombineHarvester/CombineTools/interface/Systematics.h"
#include "CombineHarvester/CombineTools/interface/Utilities.h"
#include <string>
#include <vector>
#include <functional>

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

namespace {

    class SystematicBuilder {
        public:
            SystematicBuilder(CombineHarvester& cb) : cb_(cb) {}

            // Add a systematic with a simple value (for lnN or shape).
            void AddSyst(const std::string& name, const std::string& type,
                        const std::vector<std::string>& processes,
                        const std::vector<std::string>& channels,
                        double value = 1.0) {
                cb_.cp().process(processes).channel(channels).AddSyst(cb_, name, type, SystMap<>::init(value));
            }

            // Add a systematic with a complex SystMap (e.g., binned by channel and bin_id).
            template <typename... Args>
            void AddSyst(const std::string& name, const std::string& type,
                        const std::vector<std::string>& processes,
                        const std::vector<std::string>& channels,
                        const SystMap<Args...>& map) {
                cb_.cp().process(processes).channel(channels).AddSyst(cb_, name, type, map);
            }

        private:
            CombineHarvester& cb_;
        };
    
    namespace channels {
        const std::vector<std::string> all    = {"et", "mt", "tt", "em"};
        const std::vector<std::string> lt_tt  = {"et", "mt", "tt"};
        const std::vector<std::string> mt_tt  = {"mt", "tt"};
        const std::vector<std::string> et_tt  = {"et", "tt"};
        const std::vector<std::string> em_et  = {"em", "et"};
        const std::vector<std::string> em_mt  = {"em", "mt"};
        const std::vector<std::string> em_lt  = {"em", "mt", "et"};
        const std::vector<std::string> lt     = {"et", "mt"};
        const std::vector<std::string> et     = {"et"};
        const std::vector<std::string> mt     = {"mt"};
        const std::vector<std::string> tt     = {"tt"};
        const std::vector<std::string> em     = {"em"};
    };

    namespace processes {
        const std::vector<std::string> ggH = {
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
        const std::vector<std::string> ggZH_had = {
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
        const std::vector<std::string> qqH = {
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
        const std::vector<std::string> VH_had = {
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
        const std::vector<std::string> VH = {
            // STXS stage 0
            "WH_lep_htt", "ZH_lep_htt", "ggZH_lep_htt", "ttH_htt",
            "WH_htt", "ZH_htt",
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
        const std::vector<std::string> ggHToWW = {"ggH_hww"};
        const std::vector<std::string> qqHToWW = {"qqH_hww"};
        const std::vector<std::string> signals = JoinStr({ggH, ggZH_had, qqH, VH_had, VH});

        const std::vector<std::string> htt = JoinStr({ggH, ggZH_had, qqH, VH_had, VH});
        const std::vector<std::string> hww = JoinStr({ggHToWW, qqHToWW, {"WH_hww", "ZH_hww"}});

        const std::vector<std::string> z = {"ZTT", "ZL", "ZJ"};
        const std::vector<std::string> ttbar = {"TTT", "TTL", "TTJ", "TT"};
        const std::vector<std::string> vv = {"VVT", "VVJ", "VVL", "VV", "ST"};
        const std::vector<std::string> w = {"W"};
        const std::vector<std::string> real_tau_bkg_mc = {"ZTT", "TTT", "TTL", "VVT", "VVL"};
        const std::vector<std::string> jetFakes_bkg_mc = {"W", "TTJ", "ZJ", "VVJ"};
        const std::vector<std::string> qcd = {"QCD"};
        const std::vector<std::string> emb = {"EMB"};
        const std::vector<std::string> jetFakes = {"jetFakes"};

        const std::vector<std::string> mc_bkg = JoinStr({z, ttbar, vv, w});
        const std::vector<std::string> mc = JoinStr({htt, hww, mc_bkg});
        const std::vector<std::string> real_taus_mc = JoinStr({htt, real_tau_bkg_mc});
        const std::vector<std::string> tes_affected_processes = JoinStr({signals, real_tau_bkg_mc, jetFakes});

        const std::vector<std::string> mc_gte1j = [] {
            std::vector<std::string> procs;
            for (const auto& proc : mc) {
                if (proc.find("_0J") == std::string::npos) { // std::string::npos means not found
                    procs.push_back(proc);
                }
            }
            return procs;
        }();
    };
}

// private:
//     CombineHarvester& cb_;
// };

  void AddSMRun2Systematics(CombineHarvester &cb, bool jetfakes, bool embedding, bool regional_jec, bool ggh_wg1, bool qqh_wg1, int era) {

    bool split_tau_id_and_es_by_pt = true;
    bool use_ml_ff_scheme = true;
    bool correlate_emb = true;

    using namespace std::string_literals;
    const std::vector<std::string> tau_decaymodes = {"1prong0pizero", "1prong1pizero", "3prong0pizero", "3prong1pizero"};

    SystematicBuilder builder(cb);

    // Uncertainty Lumi ; References:
    // - "CMS Luminosity Measurements for the 2016 Data Taking Period" (PAS, https://cds.cern.ch/record/2257069)
    // - Recommendation twiki https://twiki.cern.ch/twiki/bin/view/CMS/TWikiLUM#LumiComb  
    
    float lumi_unc = (era == 2016) ? 1.012 : (era == 2017) ? 1.0082 : (era == 2018) ? 1.0084 : 1.0;
    float lumi_unc_corr = (era == 2016) ? 1.006 : (era == 2017) ? 1.009 : (era == 2018) ? 1.020 : 1.0;
    float lumi_unc_1718 = (era == 2017) ? 1.006 : (era == 2018) ? 1.002 : 1.0;

    builder.AddSyst("lumi_13TeV_Run$ERA", "lnN", processes::mc, channels::all, lumi_unc);
    builder.AddSyst("lumi_13TeV_correlated", "lnN", processes::mc, channels::all, lumi_unc_corr);
    builder.AddSyst("lumi_13TeV_1718", "lnN", processes::mc, channels::all, lumi_unc_1718);

    // Uncertainty: Prefiring; References:
    // - "https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1ECALPrefiringWeightRecipe", Note: assumed uncorrelated accross years

    if (era != 2018){
        builder.AddSyst("CMS_prefiring", "shape", processes::mc, channels::all);
    }

    // Trigger efficiencies, TODO: References

    builder.AddSyst("CMS_eff_trigger_et_Run$ERA", "shape", processes::mc, channels::et);
    builder.AddSyst("CMS_eff_trigger_mt_Run$ERA", "shape", processes::mc, channels::mt);
    builder.AddSyst("CMS_eff_trigger_em_Run$ERA", "lnN", processes::mc, channels::em, 1.02);

    if (embedding){
        builder.AddSyst("CMS_eff_trigger_emb_et_Run$ERA", "shape", processes::emb, channels::et);
        builder.AddSyst("CMS_eff_trigger_emb_mt_Run$ERA", "shape", processes::emb, channels::mt);
        builder.AddSyst("CMS_eff_trigger_emb_em_Run$ERA", "lnN", processes::emb, channels::em, 1.02);
    }

    // TODO: Adjust those to the new naming scheme.
    for (const auto& _dm : tau_decaymodes) {
        builder.AddSyst("CMS_eff_trigger_tt_dm"s + _dm + "_Run$ERA", "shape", processes::mc, channels::tt);
        builder.AddSyst("CMS_eff_trigger_tt_dm"s + _dm + "_Run$ERA", "shape", processes::emb, channels::tt, 0.5);
        builder.AddSyst("CMS_eff_trigger_emb_tt_dm"s + _dm + "_Run$ERA", "shape", processes::emb, channels::tt, 0.866);
    }

    // 3% in Tau ID SF with different anti-l fake WP

    if (correlate_emb) {
        builder.AddSyst("CMS_eff_t_wp_Run$ERA", "lnN", JoinStr({processes::htt, processes::emb, processes::real_tau_bkg_mc}), channels::mt_tt, 1.03);
    } else {
        // Decorrelated: MC gets the standard name, EMB gets an emb-specific name
        builder.AddSyst("CMS_eff_t_wp_Run$ERA", "lnN", JoinStr({processes::htt, processes::real_tau_bkg_mc}), channels::mt_tt, 1.03);
        builder.AddSyst("CMS_eff_t_emb_wp_Run$ERA", "lnN", processes::emb, channels::mt_tt, 1.03);
    }

    // Lepton ID

    builder.AddSyst("CMS_eff_e", "lnN", processes::mc, channels::em_et, 1.02);
    builder.AddSyst("CMS_eff_m", "lnN", processes::mc, channels::em_mt, 1.02);
    if (embedding) {
        if (correlate_emb) {
            builder.AddSyst("CMS_eff_e", "lnN", processes::emb, channels::em_et, 1.01);
            builder.AddSyst("CMS_eff_m", "lnN", processes::emb, channels::em_mt, 1.01);
            builder.AddSyst("CMS_eff_e_emb", "lnN", processes::emb, channels::em_et, 1.017);
            builder.AddSyst("CMS_eff_m_emb", "lnN", processes::emb, channels::em_mt, 1.017);
        } else {
            builder.AddSyst("CMS_eff_e_emb", "lnN", processes::emb, channels::em_et, 1.02);
            builder.AddSyst("CMS_eff_m_emb", "lnN", processes::emb, channels::em_mt, 1.02);
        }
    }
    // Tau ID

    // MC

    const std::vector<std::string> detailed_tau_decaymodes = {
        "1prong0pizero", "1prong1pizero", "3prong0pizero", "3prong1pizero"
    };
    std::vector<std::string> tau_pt_suffixes;

    if (split_tau_id_and_es_by_pt) {
        tau_pt_suffixes = {"20to40", "40toInf"};
    } else {
        tau_pt_suffixes = {""}; // "CMS_eff_t_1prong0pizero_Run2018" vs "CMS_eff_t_1prong0pizero_20to40_Run2018"
    }

    // Tau ID + Tau Energy Scale
    for (const auto& dm : detailed_tau_decaymodes) {
        for (const auto& pt_suffix : tau_pt_suffixes) {
            
            std::string bin_name = dm;

            if (!pt_suffix.empty()){
                bin_name += "_" + pt_suffix;
            }

            // ---
            // Tau ID

            // Common Component (Correlated)
            // CMS_eff_t_1prong0pizero_20to40_Run2018
            std::string id_common = "CMS_eff_t_" + bin_name + "_Run$ERA";
            builder.AddSyst(id_common, "shape", processes::real_taus_mc, channels::lt);

            if (embedding) {
                if (correlate_emb) {
                    builder.AddSyst(id_common, "shape", processes::emb, channels::lt);
                }

                // Embedding Specific Component (Uncorrelated)
                // CMS_eff_t_emb_1prong0pizero20to40_Run2018
                std::string id_emb = "CMS_eff_t_emb_" + bin_name + "_Run$ERA";
                builder.AddSyst(id_emb, "shape", processes::emb, channels::lt);
            }

            // ---
            // Tau Energy Scale

            // Common Component (Correlated)
            // CMS_scale_t_1prong0pizero_20to40_Run2018
            std::string tes_common = "CMS_scale_t_" + bin_name + "_Run$ERA";
            builder.AddSyst(tes_common, "shape", processes::real_taus_mc, channels::lt);

            if (embedding) {
                if (correlate_emb) {
                    builder.AddSyst(tes_common, "shape", processes::emb, channels::lt);
                }

                // Embedding Specific Component (Uncorrelated)
                // CMS_scale_t_emb_1prong0pizero_20to40_Run2018
                std::string tes_emb = "CMS_scale_t_emb_" + bin_name + "_Run$ERA";
                builder.AddSyst(tes_emb, "shape", processes::emb, channels::lt);
            }

            if (jetfakes) {
                std::string tes_mc_jf  = "CMS_scale_t_" + bin_name + "_Run$ERA";
                std::string tes_emb_jf = "CMS_scale_t_emb_" + bin_name + "_Run$ERA";

                builder.AddSyst(tes_mc_jf, "shape", processes::jetFakes, channels::lt);

                if (embedding) {
                    builder.AddSyst(tes_emb_jf, "shape", processes::jetFakes, channels::lt);
                }
            } 
        }
    }

    for (const auto& _dm: tau_decaymodes){
        builder.AddSyst("CMS_eff_t_dm"s + _dm + "_Run$ERA", "shape", processes::real_taus_mc, channels::tt);
    }
    builder.AddSyst("CMS_eff_t_$CHANNEL_Run$ERA", "lnN", processes::real_taus_mc, channels::tt, 1.014);

    for (const auto& _dm : tau_decaymodes){
        // builder.AddSyst("CMS_eff_t_emb_dm"s + _dm + "_Run$ERA", "shape", processes::emb, channels::lt, 0.866); // add when available
        builder.AddSyst("CMS_eff_t_dm"s + _dm + "_Run$ERA", "shape", processes::emb, channels::tt, 0.5);
    }
    builder.AddSyst("CMS_eff_t_emb_$CHANNEL_Run$ERA", "lnN", processes::emb, channels::tt, 1.012);
    builder.AddSyst("CMS_eff_t_$CHANNEL_Run$ERA", "lnN", processes::emb, channels::tt, 1.007);

    builder.AddSyst("CMS_eff_t_Run$ERA", "lnN", {"W", "ZJ", "TTJ", "VVJ"}, channels::tt, 1.06);
    builder.AddSyst("CMS_eff_t_$CHANNEL_Run$ERA", "lnN", {"W", "ZJ", "TTJ", "VVJ"}, channels::tt, 1.02);

    // btag uncertainties

    for (const auto& src : {"btag_b_HF", "btag_c_CFerr1", "btag_c_CFerr2", "btag_j_LF"}) {
        builder.AddSyst("CMS_"s + src, "shape", processes::mc, channels::all);
    }

    for (const auto& src : {"btag_b_HFstats1", "btag_b_HFstats2", "btag_j_LFstats1", "btag_j_LFstats2"}) {
        builder.AddSyst("CMS_"s + src + "_Run$ERA", "shape", processes::mc, channels::all);
    }

    // electron and tau energy scales

    builder.AddSyst("CMS_scale_e", "shape", processes::mc, channels::em_et);
    builder.AddSyst("CMS_res_e", "shape", processes::mc, channels::em_et);
    builder.AddSyst("CMS_scale_e_emb", "shape", processes::emb, channels::em_et);

    // for (const auto& _dm : {"1prong", "1prong1pizero", "3prong", "3prong1pizero"}) {
    //     builder.AddSyst("CMS_scale_t_"s + _dm + "_Run$ERA", "shape", processes::tes_affected_processes, channels::lt_tt);
    //     builder.AddSyst("CMS_scale_t_"s + _dm + "_Run$ERA", "shape", processes::emb, channels::lt_tt, 0.5); // Correlated
    //     if (embedding){
    //         builder.AddSyst("CMS_scale_t_emb_"s + _dm + "_Run$ERA", "shape", JoinStr({processes::emb, {"jetFakes"}}), channels::lt_tt, 0.866);
    //     }
    // }

    // jes/jer Uncertainties

    if (!regional_jec) {
        builder.AddSyst("CMS_scale_j_Run$ERA", "shape", processes::mc_gte1j, channels::all, 0.71);
        builder.AddSyst("CMS_scale_j", "shape", processes::mc_gte1j, channels::all, 0.71);
    } else {
        for (const auto& src : {"Absolute", "BBEC1", "EC2", "HF"}){
            builder.AddSyst("CMS_scale_j_"s + src, "shape", processes::mc_gte1j, channels::all);
            builder.AddSyst("CMS_scale_j_"s + src + "_Run$ERA", "shape", processes::mc_gte1j, channels::all);
        }
        builder.AddSyst("CMS_scale_j_RelativeSample_Run$ERA", "shape", processes::mc_gte1j, channels::all);
        builder.AddSyst("CMS_scale_j_FlavorQCD", "shape", processes::mc_gte1j, channels::all);
        builder.AddSyst("CMS_scale_j_RelativeBal", "shape", processes::mc_gte1j, channels::all);
    }

    builder.AddSyst("CMS_res_j_Run$ERA", "shape", processes::mc_gte1j, channels::all);

    if (era == 2018){
        builder.AddSyst("CMS_scale_j_HEMIssue_Run$ERA", "shape", processes::mc_gte1j, channels::all);
    }

    // met energy scale and recoil
    //Z and W processes are only included due to the EWK fraction. Make sure that there is no contribution to the shift from the DY or Wjets samples.
    builder.AddSyst("CMS_scale_met_unclustered_Run$ERA", "shape", JoinStr({processes::htt, processes::hww, processes::z, processes::ttbar, processes::w, processes::vv}), channels::all);
    builder.AddSyst("CMS_htt_boson_scale_met_Run$ERA", "shape", JoinStr({processes::htt, processes::hww, processes::z, processes::w}), channels::all);
    builder.AddSyst("CMS_htt_boson_res_met_Run$ERA", "shape", JoinStr({processes::htt, processes::hww, processes::z, processes::w}), channels::all);
    

    // Uncertainties: Background normalizations
    // Notes:
    // - TODO: Remeasure QCD extrapolation factors for SS and ABCD methods?
    //          Current values are measured by KIT.
    // - TODO: Adapt for fake factor and embedding
    // - TODO: W uncertainties: Do we need lnN uncertainties based on the Ersatz
    //          study in Run1 (found in HIG-16043 uncertainty model)
    // - TODO: References?

    builder.AddSyst("CMS_htt_vvXsec", "lnN", processes::vv, channels::all, 1.05);
    builder.AddSyst("CMS_htt_tjXsec", "lnN", processes::ttbar, channels::all, 1.06);
    builder.AddSyst("CMS_htt_wjXsec", "lnN", processes::w, channels::all, 1.04);
    builder.AddSyst("CMS_htt_zjXsec", "lnN", processes::z, channels::all, 1.02);

    builder.AddSyst("CMS_ExtrapSSOS_$CHANNEL_Run$ERA", "lnN", processes::qcd, channels::et, 1.05);
    builder.AddSyst("CMS_ExtrapSSOS_$CHANNEL_Run$ERA", "lnN", processes::qcd, channels::mt, 1.03);
    builder.AddSyst("CMS_ExtrapABCD_$CHANNEL_Run$ERA", "lnN", processes::qcd, channels::tt, 1.03);
    
    for (const auto& src : {"0jet", "1jet", "2jet"}){
        for (const auto& src2 : {"rate", "shape", "shape2"}){
            builder.AddSyst("CMS_htt_qcd_"s + src + "_" + src2 + "_Run$ERA", "shape", processes::qcd, channels::em);
        }
    }
    builder.AddSyst("CMS_htt_qcd_iso", "shape", processes::qcd, channels::em);
    
    // Uncertainty: Drell-Yan LO->NLO reweighting

    if (era == 2016){
        builder.AddSyst("CMS_htt_dyShape_Run$ERA", "shape", processes::z, channels::all, 0.10);
    }
    else {
        builder.AddSyst("CMS_htt_dyShape", "shape", processes::z, channels::all, 0.10);
    }

    // Uncertainty: TT shape reweighting

    builder.AddSyst("CMS_htt_ttbarShape", "shape", processes::ttbar, channels::all);

    // Uncertainty: Electron/muon to tau fakes and ZL energy scale

    builder.AddSyst("CMS_ZLShape_mt_Run$ERA", "shape", {"ZL"}, channels::mt);

    for (int i = 1; i <= 5; ++i){
        builder.AddSyst("CMS_fake_m_WH"s + std::to_string(i) + "_Run$ERA", "shape", {"ZL"}, channels::mt);
    }

    builder.AddSyst("CMS_fake_e_BA_Run$ERA", "shape", {"ZL"}, channels::et);
    builder.AddSyst("CMS_fake_e_EC_Run$ERA", "shape", {"ZL"}, channels::et);

    // TODO: Add corresponding fake_m_{BA,EC} ?

    // Uncertainty: PileUp

    builder.AddSyst("CMS_PileUp", "shape", processes::mc, channels::all);

    // Uncertainty: Jet to tau fakes

    builder.AddSyst("CMS_htt_fake_j_Run$ERA", "shape", processes::jetFakes_bkg_mc, channels::lt_tt);

    // Uncertainty: Embedded events
    // Embedded Normalization: No Lumi, Zjxsec information used, instead derived from data using dimuon selection efficiency
    // TTbar contamination in embedded events: 10% shape uncertainty of assumed ttbar->tautau event shape

    builder.AddSyst("CMS_htt_doublemutrg_Run$ERA", "lnN", processes::emb, channels::all, 1.04);
    builder.AddSyst("CMS_htt_emb_ttbar_Run$ERA", "shape", processes::emb, channels::all);
  
    // jetFakes uncertainties

    if (jetfakes){
        if (use_ml_ff_scheme){
            const std::vector<std::string> ff_base_sources = {
                "ff_QCD",
                "ff_Wjets",
                // "ff_ttbar", done on MC so no MC subtraction uncertainty here
                // ---
                "ff_QCDStat",
                "ff_WjetsStat",
                "ff_ttbarStat",
                // ---
                "ff_QCDNormalization",
                "ff_WjetsNormalization",
                "ff_ttbarNormalization",
                // ---
                "fractions_QCD",
                "fractions_Wjets",
                "fractions_ttbar",
                // ---
                "fractions_QCDStat",
                "fractions_WjetsStat",
                "fractions_ttbarStat",
                // ---
                "ff_total_sub_syst",
                // ---
                "QCD_DR_SR_correction",
                "QCD_DR_SR_correctionStat",
                // ---
                // "Wjets_DR_SR_correction", done on MC so no MC subtraction uncertainty here
                "Wjets_DR_SR_correctionStat",
                // ---
                "QCD_non_closure_CorrStat1Sigma",
                "QCD_non_closure_CorrSystMCShift",
                "QCD_non_closure_CorrSystBandAsym",
                // ---
                "Wjets_non_closure_CorrStat1Sigma",
                "Wjets_non_closure_CorrSystMCShift",
                "Wjets_non_closure_CorrSystBandAsym",
                // ---
                // "ttbar_non_closure_CorrSystMCShift", done on MC so no MC subtraction uncertainty here
                "ttbar_non_closure_CorrStat1Sigma",
                "ttbar_non_closure_CorrSystBandAsym",

            };

            for (const auto& unc : ff_base_sources) {
                builder.AddSyst("CMS_"s + unc + "_$CHANNEL_Run$ERA", "shape", processes::jetFakes, channels::lt);
            }
        } else {

            const std::vector<std::string> ff_base_sources = {
                "QCDFFunc",
                "WjetsFFunc",
                "ttbarFFunc",
                // ---
                "QCDFFmcSubUnc",
                "WjetsFFmcSubUnc",
                // ---
                "process_fractionsfracQCDUnc",
                "process_fractionsfracWjetsUnc",
                "process_fractionsfracTTbarUnc",
                // ---
                "ff_total_sub_syst",
                // ---
                "QCD_DR_SR_CorrStat1Sigma",
                "QCD_DR_SR_CorrSystMCShift",
                "QCD_DR_SR_CorrSystBandAsym",
                // ---
                "Wjets_DR_SR_CorrStat1Sigma",
                "Wjets_DR_SR_CorrSystBandAsym",
                // ---
                "QCD_non_closure_CorrStat1Sigma",
                "QCD_non_closure_CorrSystMCShift",
                "QCD_non_closure_CorrSystBandAsym",
                // ---
                "Wjets_non_closure_CorrStat1Sigma",
                "Wjets_non_closure_CorrSystMCShift",
                "Wjets_non_closure_CorrSystBandAsym",
                // ---
                "ttbar_non_closure_CorrStat1Sigma",
                "ttbar_non_closure_CorrSystBandAsym",
                // ---
            };

            for (const auto& unc : ff_base_sources) {
                builder.AddSyst("CMS_"s + unc + "_$CHANNEL_Run$ERA", "shape", processes::jetFakes, channels::lt);
            }
        }
    }

    // Uncertainty: Theory uncertainties

    builder.AddSyst("BR_htt_THU", "lnN", processes::htt, channels::all, 1.0117);
    builder.AddSyst("BR_htt_PU_mq", "lnN", processes::htt, channels::all, 1.0099);
    builder.AddSyst("BR_htt_PU_alphas", "lnN", processes::htt, channels::all, 1.0061);
    builder.AddSyst("BR_hww_THU", "lnN", processes::hww, channels::all, 1.0098);
    builder.AddSyst("BR_hww_PU_mq", "lnN", processes::hww, channels::all, 1.0097);
    builder.AddSyst("BR_hww_PU_alphas", "lnN", processes::hww, channels::all, 1.0063);

    if (!ggh_wg1){
        builder.AddSyst("QCDScale_ggH", "lnN", JoinStr({processes::ggH, processes::ggHToWW}), channels::all, 1.039);
    }
    if (!qqh_wg1){
        builder.AddSyst("QCDScale_qqH", "lnN", JoinStr({processes::qqH, processes::qqHToWW}), channels::all, 1.005);
    }

    if (ggh_wg1) {
        const std::vector<std::string> thu_ggh = {"Mig01", "Mig12", "Mu", "PT60", "PT120", "Res", "VBF2j", "VBF3j", "qmtop"};
        for(const auto& src : thu_ggh){
            builder.AddSyst("THU_ggH_" + src, "shape", processes::ggH, channels::all);
        }
    }

    if (qqh_wg1) {
        const std::vector<std::string> thu_qqh = {"TOT", "PTH200", "Mjj60", "Mjj120", "Mjj350", "Mjj700", "Mjj1000", "Mjj1500", "25", "JET01"};
        for(const auto& src : thu_qqh){
            builder.AddSyst("THU_qqH_" + src, "shape", processes::qqH, channels::all);
        }
    }

    builder.AddSyst("LHE_alphaS", "shape", JoinStr({processes::ggH, processes::qqH}), channels::all);
    builder.AddSyst("LHE_pdf", "shape", JoinStr({processes::ggH, processes::qqH}), channels::all);
    builder.AddSyst("LHE_scale_muF_norm", "shape", JoinStr({processes::ggH, processes::qqH}), channels::all);
    builder.AddSyst("LHE_scale_muR_norm", "shape", JoinStr({processes::ggH, processes::qqH}), channels::all);
    builder.AddSyst("PS_scale_Fsr_norm", "shape", processes::ggH, channels::all);
    builder.AddSyst("PS_scale_Isr_norm", "shape", processes::ggH, channels::all);

    for (const auto& src : {"muF_", "muR_"}){
        builder.AddSyst("ggH_scale_0jet_"s + src, "shape", processes::ggH, channels::all);
        builder.AddSyst("ggH_scale_1jet_lowpt_"s + src, "shape", processes::ggH, channels::all);
        builder.AddSyst("ggH_scale_2jet_lowpt_"s + src, "shape", processes::ggH, channels::all);
        builder.AddSyst("ggH_scale_highpt_"s + src, "shape", processes::ggH, channels::all);
        builder.AddSyst("ggH_scale_vbf_"s + src, "shape", processes::ggH, channels::all);
        builder.AddSyst("ggH_scale_very_highpt_"s + src, "shape", processes::ggH, channels::all);

        builder.AddSyst("vbf_scale_0jet_"s + src, "shape", processes::qqH, channels::all);
        builder.AddSyst("vbf_scale_1jet_"s + src, "shape", processes::qqH, channels::all);
        builder.AddSyst("vbf_scale_highmjj_highpt_"s + src, "shape", processes::qqH, channels::all);
        builder.AddSyst("vbf_scale_highmjj_lowpt_"s + src, "shape", processes::qqH, channels::all);
        builder.AddSyst("vbf_scale_lowmjj_"s + src, "shape", processes::qqH, channels::all);
    }
}
} // namespace ch
