import logging
import sys
from itertools import count

from CombineHarvester.SMRun2Legacy.custom_logging import setup_logging

logger = setup_logging(logger=logging.getLogger(__name__))


def _get_dnn_class_mapping() -> dict:
    """
    Lazily import DNN_CLASS_MAPPING from the smhtt_ul analysis config.
    The config path is read from the SMHTT_UL_PATH environment variable,
    which must be set before calling get_categories with stxs_stage0_run3.
    """
    import os
    smhtt_path = os.environ.get("SMHTT_UL_PATH")
    if smhtt_path is None:
        raise EnvironmentError(
            "SMHTT_UL_PATH environment variable is not set. "
            "Pass --smhtt-path to make_datacards.py."
        )
    config_dir = os.path.join(smhtt_path, "config", "shapes")
    if config_dir not in sys.path:
        sys.path.insert(0, config_dir)
    from category_selection import DNN_CLASS_MAPPING
    return DNN_CLASS_MAPPING


def get_categories(
    channel: str,
    categorization: str,
    gof_category_name: str = "gof",
    embedding: bool = True,
    jetfakes: bool = True,
) -> list:
    logger.debug(f"Calling get_categories with {locals()}")

    counter = count(10)  # Start category IDs from 11
    background_categories = []

    if categorization == "gof":
        categories = [(300, gof_category_name)]
        logger.info(f"Using GOF categorization with category: {categories[0]}")
        return categories

    # stxs_stage0_run3: pure DNN categorization for Run3 — DNN score bins serve
    # as both signal and control regions.  embedding / jetfakes flags are irrelevant.
    if categorization == "stxs_stage0_run3":
        # mt/et: ggh(0) qqh(1) | dyjets_tt(2) dyjets_ll(3) jetFakes(4) ttbar(5) diboson(6)
        # tt:    ggh(0) qqh(1) | dyjets_tt(2) jetFakes(3) bkg_rest(4)
        # Category names match convert_to_synced_shapes.py output: {channel}_{category}
        # (same naming convention as the Run2 categories below, e.g. {channel}_jetFakes).
        dnn_class_mapping = _get_dnn_class_mapping()
        if channel not in dnn_class_mapping:
            raise NotImplementedError(f"stxs_stage0_run3 not defined for channel '{channel}'")
        mapping = dnn_class_mapping[channel]
        signal_dnn = [
            (info["index"] + 1, f"{channel}_{name}")
            for name, info in sorted(mapping.items(), key=lambda x: x[1]["index"])
            if info["is_signal"]
        ]
        bkg_dnn = [
            (info["index"] + 1, f"{channel}_{name}")
            for name, info in sorted(mapping.items(), key=lambda x: x[1]["index"])
            if not info["is_signal"]
        ]
        logger.info(f"Using stxs_stage0_run3: signal={signal_dnn}, background={bkg_dnn}")
        return signal_dnn + bkg_dnn

    # stxs_stage1p2_run3: placeholder for Run3 STXS stage 1.2 categorization
    if categorization == "stxs_stage1p2_run3":
        raise NotImplementedError("stxs_stage1p2_run3 is not yet implemented.")

    # All Run2 categorizations below use embedding / jetFakes control regions
    if embedding:
        background_categories.append((next(counter), f"{channel}_embedding"))
        logger.info(f"Added {background_categories[-1]}")
    else:
        logger.info(f"Skipping embedding category for channel {channel} (non-embedding mode)")

    if jetfakes:
        background_categories.append((next(counter), f"{channel}_jetFakes"))
        logger.info(f"Added {background_categories[-1]}")
    else:
        logger.info(f"Skipping jetFakes category for channel {channel}")

    if channel == "mt":
        background_categories.extend(
            [
                (next(counter), f"{channel}_ttbar"),
                (next(counter), f"{channel}_dyjets"),
                (next(counter), f"{channel}_diboson"),
            ],
        )
        logger.info(f"Added mt-specific background categories: {background_categories[-3:]}")
    elif channel == "et":
        background_categories.extend(
            [
                (next(counter), f"{channel}_ttbar"),
                (next(counter), f"{channel}_dyjets"),
                (next(counter), f"{channel}_diboson"),
                (next(counter), f"{channel}_wjets"),
            ],
        )
        logger.info(f"Added et-specific background categories: {background_categories[-4:]}")
    elif channel == "tt":
        background_categories.extend(
            [
                (next(counter), f"{channel}_ttbar"),
                (next(counter), f"{channel}_dyjets"),
            ],
        )
        logger.info(f"Added tt-specific background categories: {background_categories[-2:]}")
    else:
        raise NotImplementedError(f"Channel-specific background categories for {channel} not implemented.")

    if categorization == "stxs_stage0":
        signal_categories = [(1, f"{channel}_qqh"), (2, f"{channel}_ggh")]
        logger.info(f"Using stxs_stage0 with signal categories: {signal_categories}")
        return signal_categories + background_categories

    elif categorization == "stxs_stage1p2_run2":
        signal_categories = [
            (100, f"{channel}_qqh_bin201to210"),
            (101, f"{channel}_ggh_bin101to104"),
            (102, f"{channel}_ggh_bin105to106"),
            (103, f"{channel}_ggh_bin107to109"),
            (104, f"{channel}_ggh_bin110to116"),
        ]
        logger.info(f"Using stxs_stage1p2_run2 with signal categories: {signal_categories}")
        return signal_categories + background_categories

    elif categorization == "stxs_stage1p2_run2_bkg_only":
        logger.info(f"Using stxs_stage1p2_run2_bkg_only: {background_categories}")
        return background_categories

    else:
        raise ValueError(f"Unknown categorization: {categorization}")
