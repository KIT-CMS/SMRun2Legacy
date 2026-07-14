import logging

from CombineHarvester.SMRun2Legacy.custom_logging import setup_logging


logger = setup_logging(logger=logging.getLogger(__name__))


def get_backgrounds(
    channel: str,
    embedding: bool = True,
    jetfakes: bool = True,
) -> list:
    logger.debug(f"Calling get_backgrounds with {locals()}")

    # Run3 non-embedding path: explicit per-channel lists matching the nngrouped
    # shape file processes.  W, ZJ, TTJ, VVJ are folded into jetFakesMC (->jetFakes).
    if not embedding:
        if channel in ("mt", "et"):
            # genuine-tau MC backgrounds + lepton fakes + MC-based jet fakes
            bkgs = ["ZTT", "ZL", "ZJ", "TTT", "TTL", "TTJ", "VVT", "VVL", "VVJ", "W"]
        elif channel == "tt":
            # tt: lepton-fake contribution (ZL, TTL, VVL) negligible but kept for completeness
            bkgs = ["ZTT", "ZL", "ZJ", "TTT", "TTL", "TTJ", "VVT", "VVL", "VVJ", "W"]
        else:
            bkgs = ["ZTT", "ZL", "ZJ", "TTT", "TTL", "TTJ", "VVT", "VVL", "VVJ", "W"]

        if jetfakes:
            # MC-based fake estimation: W/ZJ/TTJ/VVJ folded into jetFakes
            bkgs = [b for b in bkgs if b not in ["W", "ZJ", "TTJ", "VVJ"]]
            bkgs.append("jetFakes")
            logger.info("Non-embedding + MC jetFakes: folded W/ZJ/TTJ/VVJ into jetFakes")

        logger.info(f"Run3 non-embedding background processes for {channel}: {bkgs}")
        return bkgs

    # Run2 embedding path (original behaviour preserved)
    bkgs = ["W", "ZTT", "QCD", "ZL", "ZJ", "TTT", "TTL", "TTJ", "VVJ", "VVT", "VVL"]
    logger.debug(f"Initial background processes for channel {channel}: {bkgs}")

    if channel == "em":
        bkgs = ["W", "ZTT", "TTT", "VVT", "QCD", "ZL", "TTL", "VVL", "ttH_htt", "ggH_hww", "qqH_hww", "WH_hww", "ZH_hww"]
        logger.info(f"Updated background processes for em channel: {bkgs}")

    bkgs = [bkg for bkg in bkgs if bkg not in ["ZTT", "TTT", "VVT"]]
    bkgs.append("EMB")
    logger.info("Replaced [ZTT, TTT, VVT] with [EMB]")

    if jetfakes:
        bkgs = [bkg for bkg in bkgs if bkg not in ["QCD", "W", "VVJ", "TTJ", "ZJ"]]
        bkgs.append("jetFakes")
        logger.info("Replaced [QCD, W, VVJ, TTJ, ZJ] with [jetFakes]")

    logger.info(f"Used background processes for {channel}: {bkgs}")
    return bkgs

    logger.info(f"Used background processes for {channel}: {bkgs}")
    return bkgs


def get_signals(stxs_version: str = "stxs_stage0") -> list:
    logger.debug(f"Calling get_signals with {locals()}")
    if stxs_version == "stxs_stage0":
        signals = ["ggH_htt", "qqH_htt"]
        logger.info(f"Using STXS stage 0 signals: {signals}")
        return signals

    elif stxs_version == "stxs_stage1p1":
        raise NotImplementedError("STXS stage 1.1 signals not implemented in this migration yet.")

    elif stxs_version == "stxs_stage1p2_syst":
        signals = [
            "qqH125-qqh_htautau_bin201to210_selection",
            "ggH125-ggh_htautau_bin101to104_selection",
            "ggH125-ggh_htautau_bin105to106_selection",
            "ggH125-ggh_htautau_bin107to109_selection",
            "ggH125-ggh_htautau_bin110to116_selection",
        ]
        logger.info(f"Using STXS stage 1.2 signals: {signals}")
        return signals

    else:
        raise ValueError(f"Unknown STXS version: {stxs_version}")
