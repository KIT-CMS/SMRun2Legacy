import logging

import CombineHarvester.CombineTools.ch as ch
from CombineHarvester.SMRun2Legacy.custom_logging import setup_logging


logger = setup_logging(logger=logging.getLogger(__name__))


class Channels:
    def __init__(self) -> None:
        self.all = ["et", "mt", "tt", "em"]
        self.lt_tt = ["et", "mt", "tt"]
        self.mt_tt = ["mt", "tt"]
        self.et_tt = ["et", "tt"]
        self.em_et = ["em", "et"]
        self.em_mt = ["em", "mt"]
        self.em_lt = ["em", "mt", "et"]
        self.lt = ["et", "mt"]
        self.et = ["et"]
        self.mt = ["mt"]
        self.tt = ["tt"]
        self.em = ["em"]

    def __str__(self) -> str:
        _string = "Channels:\n"
        for attr in dir(self):
            if not attr.startswith('_') and attr != 'all':
                _string += f"  {attr}: {getattr(self, attr)}\n"
        return _string

    def __repr__(self) -> str:
        return self.__str__()


class Processes:
    def __init__(self) -> None:
        self.ggH = [
            # STXS stage 0
            "ggH_htt",
            # STXS stage 1.1
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
            # STXS stage 1.2 syst
            "ggH125-ggh_htautau_bin101to104_selection",
            "ggH125-ggh_htautau_bin105to106_selection",
            "ggH125-ggh_htautau_bin107to109_selection",
            "ggH125-ggh_htautau_bin110to116_selection",
            "ggH125-ggh_htautau_bin101to104_selection125",
            "ggH125-ggh_htautau_bin105to106_selection125",
            "ggH125-ggh_htautau_bin107to109_selection125",
            "ggH125-ggh_htautau_bin110to116_selection125",

        ]
        self.ggZH_had = [
            # STXS stage 0
            "ggZH_had_htt",
            # STXS stage 1.1
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
        ]
        self.qqH = [
            # STXS stage 0
            "qqH_htt",
            # STXS stage 1
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
            "qqH_GE2J_MJJ_GT700_PTH_0_200_PTHJJ_GT25_htt",
            # STXS stage 1.2 syst
            "qqH125-qqh_htautau_bin201to210_selection",  # inclusive
            "qqH125-qqh_htautau_bin201to202_selection",
            "qqH125-qqh_htautau_bin203to210_selection",
            "qqH125-qqh_htautau_bin201to202_selection",
            "qqH125-qqh_htautau_bin203to210_selection",
            "qqH125-qqh_htautau_bin201to210_selection125",  # inclusive
            "qqH125-qqh_htautau_bin201to202_selection125",
            "qqH125-qqh_htautau_bin203to210_selection125",
            "qqH125-qqh_htautau_bin201to202_selection125",
            "qqH125-qqh_htautau_bin203to210_selection125",
        ]
        self.VH_had = [
            # STXS stage 0
            "WH_had_htt",
            "ZH_had_htt",
            # STXS stage 1
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
        ]
        self.VH = [
            # STXS stage 0
            "WH_lep_htt", "ZH_lep_htt", "ggZH_lep_htt", "ttH_htt",
            "WH_htt", "ZH_htt",
            # STXS stage 1
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
        ]
        self.ggHToWW = ["ggH_hww"]
        self.qqHToWW = ["qqH_hww"]
        self.signals = self.ggH + self.ggZH_had + self.qqH + self.VH_had + self.VH

        self.htt = self.ggH + self.ggZH_had + self.qqH + self.VH_had + self.VH
        self.hww = self.ggHToWW + self.qqHToWW + ["WH_hww", "ZH_hww"]

        self.z = ["ZTT", "ZL", "ZJ"]
        self.ttbar = ["TTT", "TTL", "TTJ", "TT"]
        self.vv = ["VVT", "VVJ", "VVL", "VV", "ST"]
        self.w = ["W"]
        self.real_tau_bkg_mc = ["ZTT", "TTT", "TTL", "VVT", "VVL"]
        self.jetFakes_bkg_mc = ["W", "TTJ", "ZJ", "VVJ"]
        self.qcd = ["QCD"]
        self.emb = ["EMB"]
        self.jetFakes = ["jetFakes"]

        self.mc_bkg = self.z + self.ttbar + self.vv + self.w
        self.mc = self.htt + self.hww + self.mc_bkg
        self.real_taus_mc = self.htt + self.real_tau_bkg_mc
        self.tes_affected_processes = self.signals + self.real_tau_bkg_mc + self.jetFakes

        # mc_gte1j not migrated, since potential error

    def __str__(self) -> str:
        _string = "Processes:\n"
        for attr in dir(self):
            if not attr.startswith('_') and attr != 'all':
                _string += f"  {attr}: {getattr(self, attr)}\n"
        return _string

    def __repr__(self) -> str:
        return self.__str__()


tau_decaymodes = {"1prong0pizero", "1prong1pizero", "3prong0pizero", "3prong1pizero"}


def add_systematics(
    cb: ch.CombineHarvester,
    era: str,
    jetfakes: bool = True,
    embedding: bool = True,
    split_tau_id_and_es_by_pt: bool = True,
    use_ml_ff_scheme: bool = True,
    correlate_emb: bool = True,
    regional_jec: bool = False,
    ggh_wg1: bool = True,
    qqh_wg1: bool = True,
    shape_systematics: bool = True,
) -> None:
    logger.debug(f"calling add_systematics {locals()}")

    channels, processes = Channels(), Processes()

    logger.debug(f"Defined channels: {channels}")
    logger.debug(f"Defined processes: {processes}")

    def add_syst(name, syst_type, processes, channels, value=1.0):
        if syst_type == "shape" and not shape_systematics:
            logger.debug(f"Skipping shape systematic {name} (shape systematics disabled)")
            return

        if isinstance(value, (float, int)):
            syst_map = ch.SystMap()(float(value))
        else:
            syst_map = value

        logger.debug(f"Adding systematic {name} of type {syst_type} to processes {processes} in channels {channels} with value/map {value}")
        cb.cp().process(processes).channel(channels).AddSyst(cb, name, syst_type, syst_map)

    # Lumi
    # TODO: Run3 entries below reuse the Run2 2018 value as a placeholder
    # (https://twiki.cern.ch/twiki/bin/view/CMS/TWikiLUM requires CERN SSO and
    # couldn't be fetched) -- replace with the official Run3 numbers once available.
    lumi_unc = {
        "2016": 1.012, "2016preVFP": 1.012, "2016postVFP": 1.012,
        "2017": 1.0082,
        "2018": 1.0084,
        "2022preEE": 1.0084, "2022postEE": 1.0084,
        "2023preBPix": 1.0084, "2023postBPix": 1.0084,
        "2024": 1.0084, "2025": 1.0084,
    }
    lumi_unc_corr = {
        "2016": 1.006, "2016preVFP": 1.006, "2016postVFP": 1.006,
        "2017": 1.009,
        "2018": 1.020,
        "2022preEE": 1.020, "2022postEE": 1.020,
        "2023preBPix": 1.020, "2023postBPix": 1.020,
        "2024": 1.020, "2025": 1.020,
    }
    # lumi_13TeV_1718 models the partial correlation specific to the 2017/2018
    # Run2 measurement and does not apply to other eras -- intentionally left
    # Run2-only, so it stays a no-op (1.0) elsewhere.
    lumi_unc_1718 = {"2016": 1.0,    "2017": 1.006,  "2018": 1.002}

    # lumi
    if int(era[:4]) < 2022:
        add_syst("lumi_13TeV_Run$ERA", "lnN", processes.mc, channels.all, lumi_unc.get(era, 1.0))
        add_syst("lumi_13TeV_correlated", "lnN", processes.mc, channels.all, lumi_unc_corr.get(era, 1.0))
        add_syst("lumi_13TeV_1718", "lnN", processes.mc, channels.all, lumi_unc_1718.get(era, 1.0))
    else:
        add_syst("lumi_13.6TeV_Run$ERA", "lnN", processes.mc, channels.all, lumi_unc.get(era, 1.0))
        add_syst("lumi_13.6TeV_correlated", "lnN", processes.mc, channels.all, lumi_unc_corr.get(era, 1.0))

    # Prefiring "https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1ECALPrefiringWeightRecipe", Note: assumed uncorrelated accross years
    if era != "2018":
        add_syst("CMS_prefiring", "shape", processes.mc, channels.all)

    # Trigger
    add_syst("CMS_eff_trigger_et_Run$ERA", "shape", processes.mc, channels.et)
    add_syst("CMS_eff_trigger_mt_Run$ERA", "shape", processes.mc, channels.mt)
    add_syst("CMS_eff_trigger_em_Run$ERA", "lnN", processes.mc, channels.em, 1.02)

    if embedding:
        add_syst("CMS_eff_trigger_emb_et_Run$ERA", "shape", processes.emb, channels.et)
        add_syst("CMS_eff_trigger_emb_mt_Run$ERA", "shape", processes.emb, channels.mt)
        add_syst("CMS_eff_trigger_emb_em_Run$ERA", "lnN", processes.emb, channels.em, 1.02)

    # TODO: adjust those to the new naming scheme
    for df in tau_decaymodes:
        add_syst(f"CMS_eff_trigger_tt_dm_{df}_Run$ERA", "shape", processes.mc, channels.tt)
        add_syst(f"CMS_eff_trigger_tt_dm_{df}_Run$ERA", "shape", processes.emb, channels.tt, 0.5)
        add_syst(f"CMS_eff_trigger_emb_tt_dm_{df}_Run$ERA", "shape", processes.emb, channels.tt, 0.866)

    # 3% in Tau ID SF with different anti-l fake WP
    if correlate_emb:
        add_syst("CMS_eff_t_wp_Run$ERA", "lnN", processes.htt + processes.emb + processes.real_tau_bkg_mc, channels.mt_tt, 1.03)
    else:
        # Decorrelated: MC gets the standard name, EMB gets an emb-specific name
        add_syst("CMS_eff_t_wp_Run$ERA", "lnN", processes.htt + processes.real_tau_bkg_mc, channels.mt_tt, 1.03)
        add_syst("CMS_eff_t_emb_wp_Run$ERA", "lnN", processes.emb, channels.mt_tt, 1.03)

    # Lepton ID
    add_syst("CMS_eff_e", "lnN", processes.mc, channels.em_et, 1.02)
    add_syst("CMS_eff_m", "lnN", processes.mc, channels.em_mt, 1.02)
    if embedding:
        if correlate_emb:
            add_syst("CMS_eff_e", "lnN", processes.emb, channels.em_et, 1.01)
            add_syst("CMS_eff_m", "lnN", processes.emb, channels.em_mt, 1.01)
            add_syst("CMS_eff_e_emb", "lnN", processes.emb, channels.em_et, 1.017)
            add_syst("CMS_eff_m_emb", "lnN", processes.emb, channels.em_mt, 1.017)
        else:
            add_syst("CMS_eff_e_emb", "lnN", processes.emb, channels.em_et, 1.02)
            add_syst("CMS_eff_m_emb", "lnN", processes.emb, channels.em_mt, 1.02)

    # Tau ID + Tau ES
    for dm in tau_decaymodes:
        for tau_pt_suffix in ["20to40", "40toInf"] if split_tau_id_and_es_by_pt else [""]:
            bin_name = f"{dm}_{tau_pt_suffix}" if tau_pt_suffix else dm  # "CMS_eff_t_1prong0pizero_Run2018" vs "CMS_eff_t_1prong0pizero_20to40_Run2018"

            # ID
            id_common = f"CMS_eff_t_{bin_name}_Run$ERA"
            add_syst(id_common, "shape", processes.real_taus_mc, channels.lt)
            if embedding:
                if correlate_emb:
                    add_syst(id_common, "shape", processes.emb, channels.lt)
                add_syst(f"CMS_eff_t_emb_{bin_name}_Run$ERA", "shape", processes.emb, channels.lt)

            # ES
            tes_common = f"CMS_scale_t_{bin_name}_Run$ERA"
            add_syst(tes_common, "shape", processes.real_taus_mc, channels.lt)
            if embedding:
                if correlate_emb:
                    add_syst(tes_common, "shape", processes.emb, channels.lt)
                add_syst(f"CMS_scale_t_emb_{bin_name}_Run$ERA", "shape", processes.emb, channels.lt)

            if jetfakes:
                add_syst(f"CMS_scale_t_{bin_name}_Run$ERA", "shape", processes.jetFakes, channels.lt)
                if embedding:
                    add_syst(f"CMS_scale_t_emb_{bin_name}_Run$ERA", "shape", processes.jetFakes, channels.lt)

    for _dm in tau_decaymodes:
        add_syst(f"CMS_eff_t_dm_{_dm}_Run$ERA", "shape", processes.real_taus_mc, channels.tt)
    add_syst("CMS_eff_t_$CHANNEL_Run$ERA", "lnN", processes.real_taus_mc, channels.tt, 1.014)

    for _dm in tau_decaymodes:
        # add_syst("CMS_eff_t_emb_dm"s + _dm + "_Run$ERA", "shape", processes.emb, channels.lt, 0.866) // add when available
        add_syst(f"CMS_eff_t_dm_{_dm}_Run$ERA", "shape", processes.emb, channels.tt, 0.5)

    add_syst("CMS_eff_t_emb_$CHANNEL_Run$ERA", "lnN", processes.emb, channels.tt, 1.012)
    add_syst("CMS_eff_t_$CHANNEL_Run$ERA", "lnN", processes.emb, channels.tt, 1.007)

    add_syst("CMS_eff_t_Run$ERA", "lnN", {"W", "ZJ", "TTJ", "VVJ"}, channels.tt, 1.06)
    add_syst("CMS_eff_t_$CHANNEL_Run$ERA", "lnN", {"W", "ZJ", "TTJ", "VVJ"}, channels.tt, 1.02)

    # btag uncertainties
    for src in {"btag_b_HF", "btag_c_CFerr1", "btag_c_CFerr2", "btag_j_LF"}:
        add_syst(f"CMS_{src}", "shape", processes.mc, channels.all)

    for src in {"btag_b_HFstats1", "btag_b_HFstats2", "btag_j_LFstats1", "btag_j_LFstats2"}:
        add_syst(f"CMS_{src}_Run$ERA", "shape", processes.mc, channels.all)

    # electron and tau energy scales
    add_syst("CMS_scale_e", "shape", processes.mc, channels.em_et)
    add_syst("CMS_res_e", "shape", processes.mc, channels.em_et)
    add_syst("CMS_scale_e_emb", "shape", processes.emb, channels.em_et)

    # jes/jer Uncertainties
    if regional_jec:
        for src in {"Absolute", "BBEC1", "EC2", "HF"}:
            add_syst(f"CMS_scale_j_{src}", "shape", processes.mc, channels.all)
            add_syst(f"CMS_scale_j_{src}_Run$ERA", "shape", processes.mc, channels.all)
        add_syst("CMS_scale_j_RelativeSample_Run$ERA", "shape", processes.mc, channels.all)
        add_syst("CMS_scale_j_FlavorQCD", "shape", processes.mc, channels.all)
        add_syst("CMS_scale_j_RelativeBal", "shape", processes.mc, channels.all)
    else:
        add_syst("CMS_scale_j_Run$ERA", "shape", processes.mc, channels.all, 0.71)
        add_syst("CMS_scale_j", "shape", processes.mc, channels.all, 0.71)

    add_syst("CMS_res_j_Run$ERA", "shape", processes.mc, channels.all)

    if era == "2018":
        add_syst("CMS_scale_j_HEMIssue_Run$ERA", "shape", processes.mc, channels.all)

    # met energy scale and recoil
    # Z and W processes are only included due to the EWK fraction. Make sure that there is no contribution to the shift from the DY or Wjets samples.
    add_syst("CMS_scale_met_unclustered_Run$ERA", "shape", processes.htt + processes.hww + processes.z + processes.ttbar + processes.w + processes.vv, channels.all)
    add_syst("CMS_htt_boson_scale_met_Run$ERA", "shape", processes.htt + processes.hww + processes.z + processes.w, channels.all)
    add_syst("CMS_htt_boson_res_met_Run$ERA", "shape", processes.htt + processes.hww + processes.z + processes.w, channels.all)

    # Uncertainties: Background normalizations
    # Notes:
    # - TODO: Remeasure QCD extrapolation factors for SS and ABCD methods?
    #          Current values are measured by KIT.
    # - TODO: Adapt for fake factor and embedding
    # - TODO: W uncertainties: Do we need lnN uncertainties based on the Ersatz
    #          study in Run1 (found in HIG-16043 uncertainty model)
    # - TODO: References?

    add_syst("CMS_htt_vvXsec", "lnN", processes.vv, channels.all, 1.05)
    add_syst("CMS_htt_tjXsec", "lnN", processes.ttbar, channels.all, 1.06)
    add_syst("CMS_htt_wjXsec", "lnN", processes.w, channels.all, 1.04)
    add_syst("CMS_htt_zjXsec", "lnN", processes.z, channels.all, 1.02)

    add_syst("CMS_ExtrapSSOS_$CHANNEL_Run$ERA", "lnN", processes.qcd, channels.et, 1.05)
    add_syst("CMS_ExtrapSSOS_$CHANNEL_Run$ERA", "lnN", processes.qcd, channels.mt, 1.03)
    add_syst("CMS_ExtrapABCD_$CHANNEL_Run$ERA", "lnN", processes.qcd, channels.tt, 1.03)

    for src in {"0jet", "1jet", "2jet"}:
        for src2 in {"rate", "shape", "shape2"}:
            add_syst(f"CMS_htt_qcd_{src}_{src2}_Run$ERA", "shape", processes.qcd, channels.em)
    add_syst("CMS_htt_qcd_iso", "shape", processes.qcd, channels.em)

    # Uncertainty: Drell-Yan LO->NLO reweighting
    if era == "2016":
        add_syst("CMS_htt_dyShape_Run$ERA", "shape", processes.z, channels.all, 0.10)
    else:
        add_syst("CMS_htt_dyShape", "shape", processes.z, channels.all, 0.10)

    # Uncertainty: TT shape reweighting
    add_syst("CMS_htt_ttbarShape", "shape", processes.ttbar, channels.all)

    # Uncertainty: Electron/muon to tau fakes and ZL energy scale
    add_syst("CMS_ZLShape_mt_Run$ERA", "shape", {"ZL"}, channels.mt)

    for i in range(1, 6):
        add_syst(f"CMS_fake_m_WH{i}_Run$ERA", "shape", {"ZL"}, channels.mt)

    add_syst("CMS_fake_e_BA_Run$ERA", "shape", {"ZL"}, channels.et)
    add_syst("CMS_fake_e_EC_Run$ERA", "shape", {"ZL"}, channels.et)

    # TODO: Add corresponding fake_m_{BA,EC} ?

    # Uncertainty: PileUp
    add_syst("CMS_PileUp", "shape", processes.mc, channels.all)

    # Uncertainty: Jet to tau fakes
    add_syst("CMS_htt_fake_j_Run$ERA", "shape", processes.jetFakes_bkg_mc, channels.lt_tt)

    # Uncertainty: Embedded events
    # Embedded Normalization: No Lumi, Zjxsec information used, instead derived from data using dimuon selection efficiency
    # TTbar contamination in embedded events: 10% shape uncertainty of assumed ttbar->tautau event shape
    add_syst("CMS_htt_doublemutrg_Run$ERA", "lnN", processes.emb, channels.all, 1.04)
    add_syst("CMS_htt_emb_ttbar_Run$ERA", "shape", processes.emb, channels.all)

    # jetFakes uncertainties
    if jetfakes:
        if use_ml_ff_scheme:
            for unc in [
                *[f"ff_{p}" for p in ["QCD", "Wjets"]],  # ttbar, done on MC so no MC subtraction uncertainty here
                *[f"ff_{p}Stat" for p in ["QCD", "Wjets", "ttbar"]],
                *[f"fractions_{p}" for p in ["QCD", "Wjets", "ttbar"]],
                *[f"fractions_{p}Stat" for p in ["QCD", "Wjets", "ttbar"]],
                *[f"{p}_DR_SR_correction" for p in ["QCD"]],  # Wjets, done on MC so no MC subtraction uncertainty here
                *[f"{p}_DR_SR_correctionStat" for p in ["QCD", "Wjets"]],
                *[f"{p}_non_closure_CorrStat1Sigma" for p in ["QCD", "Wjets", "ttbar"]],
                *[f"{p}_non_closure_CorrSystBandAsym" for p in ["QCD", "Wjets", "ttbar"]],
                *[f"{p}_non_closure_CorrSystMCShift" for p in ["QCD", "Wjets"]],  # ttbar, done on MC so no MC subtraction uncertainty here
                # --- total or per process
                # "ff_total_sub_syst",
                *[f"{p}Normalization" for p in ["ff_QCD", "ff_Wjets", "ff_ttbar"]],
            ]:
                add_syst(f"CMS_{unc}_$CHANNEL_Run$ERA", "shape", processes.jetFakes, channels.lt)
        else:
            for unc in [
                *[f"{p}FFunc" for p in ["QCD", "Wjets", "ttbar"]],
                *[f"{p}FFuncSubUnc" for p in ["QCD", "Wjets"]],  # ttbar, done on MC so no MC subtraction uncertainty here
                *[f"process_fractionsfrac{p}Unc" for p in ["QCD", "Wjets", "TTbar"]],  # TODO: fix naming convention for ttbar
                *[f"{p}_DR_SR_CorrStat1Sigma" for p in ["QCD", "Wjets"]],
                *[f"{p}_DR_SR_CorrSystBandAsym" for p in ["QCD", "Wjets"]],
                *[f"{p}_DR_SR_CorrSystMCShift" for p in ["QCD"]],  # Wjets, done on MC so no MC subtraction uncertainty here
                *[f"{p}_non_closure_CorrStat1Sigma" for p in ["QCD", "Wjets", "ttbar"]],
                *[f"{p}_non_closure_CorrSystBandAsym" for p in ["QCD", "Wjets", "ttbar"]],
                *[f"{p}_non_closure_CorrSystMCShift" for p in ["QCD", "Wjets"]],  # ttbar, done on MC so no MC subtraction uncertainty here
                "ff_total_sub_syst",
            ]:
                add_syst(f"CMS_{unc}_$CHANNEL_Run$ERA", "shape", processes.jetFakes, channels.lt)

    # Theory uncertainties
    add_syst("BR_htt_THU", "lnN", processes.htt, channels.all, 1.0117)
    add_syst("BR_htt_PU_mq", "lnN", processes.htt, channels.all, 1.0099)
    add_syst("BR_htt_PU_alphas", "lnN", processes.htt, channels.all, 1.0061)
    add_syst("BR_hww_THU", "lnN", processes.hww, channels.all, 1.0098)
    add_syst("BR_hww_PU_mq", "lnN", processes.hww, channels.all, 1.0097)
    add_syst("BR_hww_PU_alphas", "lnN", processes.hww, channels.all, 1.0063)

    if ggh_wg1:
        for src in {"Mig01", "Mig12", "Mu", "PT60", "PT120", "Res", "qqh2j", "qqh3j", "qmtop"}:
            add_syst(f"THU_ggH_{src}", "shape", processes.ggH, channels.all)
    else:
        add_syst("QCDScale_ggH", "lnN", processes.ggH + processes.ggHToWW, channels.all, 1.039)

    if qqh_wg1:
        for src in {"TOT", "PTH200", "Mjj60", "Mjj120", "Mjj350", "Mjj700", "Mjj1000", "Mjj1500", "25", "JET01"}:
            add_syst(f"THU_qqH_{src}", "shape", processes.qqH, channels.all)
    else:
        add_syst("QCDScale_qqH", "lnN", processes.qqH + processes.qqHToWW, channels.all, 1.005)

    add_syst("LHE_alphaS", "shape", processes.ggH + processes.qqH, channels.all)
    add_syst("LHE_pdf", "shape", processes.ggH + processes.qqH, channels.all)
    add_syst("LHE_scale_muF_norm", "shape", processes.ggH + processes.qqH, channels.all)
    add_syst("LHE_scale_muR_norm", "shape", processes.ggH + processes.qqH, channels.all)
    add_syst("PS_scale_Fsr_norm", "shape", processes.ggH, channels.all)
    add_syst("PS_scale_Isr_norm", "shape", processes.ggH, channels.all)

    for src in {"muF_", "muR_"}:
        add_syst(f"ggH_scale_0jet_{src}", "shape", processes.ggH, channels.all)
        add_syst(f"ggH_scale_1jet_lowpt_{src}", "shape", processes.ggH, channels.all)
        add_syst(f"ggH_scale_2jet_lowpt_{src}", "shape", processes.ggH, channels.all)
        add_syst(f"ggH_scale_highpt_{src}", "shape", processes.ggH, channels.all)
        add_syst(f"ggH_scale_qqh_{src}", "shape", processes.ggH, channels.all)
        add_syst(f"ggH_scale_very_highpt_{src}", "shape", processes.ggH, channels.all)

        add_syst(f"qqh_scale_0jet_{src}", "shape", processes.qqH, channels.all)
        add_syst(f"qqh_scale_1jet_{src}", "shape", processes.qqH, channels.all)
        add_syst(f"qqh_scale_highmjj_highpt_{src}", "shape", processes.qqH, channels.all)
        add_syst(f"qqh_scale_highmjj_lowpt_{src}", "shape", processes.qqH, channels.all)
        add_syst(f"qqh_scale_lowmjj_{src}", "shape", processes.qqH, channels.all)
