def get_backgrounds(
    channel: str,
    embedding: bool = True,
    jetfakes: bool = True,
) -> list:
    bkgs = ["W", "ZTT", "QCD", "ZL", "ZJ", "TTT", "TTL", "TTJ", "VVJ", "VVT", "VVL"]
    
    if channel == "em":
        bkgs = ["W", "ZTT", "TTT", "VVT", "QCD", "ZL", "TTL", "VVL", "ttH_htt", "ggH_hww", "qqH_hww", "WH_hww", "ZH_hww"]

    if embedding:
        bkgs = [bkg for bkg in bkgs if bkg not in ["ZTT", "TTT", "VVT"]]
        bkgs.append("EMB")
        
    if jetfakes:
        bkgs = [bkg for bkg in bkgs if bkg not in ["QCD", "W", "VVJ", "TTJ", "ZJ"]]
        bkgs.append("jetFakes")
        
    return bkgs


def get_signals(stxs_version: str = "stxs_stage0") -> list:
    if stxs_version == "stxs_stage0":
        return ["ggH_htt", "qqH_htt"]
    
    elif stxs_version == "stxs_stage1p1":
        raise NotImplementedError("STXS stage 1.1 signals not implemented in this migration yet.")
        
    elif stxs_version == "stxs_stage1p2_syst":
        return [
            "qqH125-vbf_htautau_bin201to202_selection",
            "qqH125-vbf_htautau_bin203to210_selection",
            "ggH125-ggh_htautau_bin101to104_selection",
            "ggH125-ggh_htautau_bin105to106_selection",
            "ggH125-ggh_htautau_bin107to109_selection",
            "ggH125-ggh_htautau_bin110to116_selection",
        ]

    else:
        raise ValueError(f"Unknown STXS version: {stxs_version}")
