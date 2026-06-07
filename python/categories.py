import logging
from itertools import count

from CombineHarvester.SMRun2Legacy.custom_logging import setup_logging


logger = setup_logging(logger=logging.getLogger(__name__))


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

    if embedding:
        background_categories.append((next(counter), f"{channel}_embedding"))
        logger.info(f"Added {background_categories[-1]}")
    else:
        raise NotImplementedError("Non-embedding background categories not implemented in this migration yet.")

    if jetfakes:
        background_categories.append((next(counter), f"{channel}_jetFakes"))
        logger.info(f"Added {background_categories[-1]}")
    else:
        raise NotImplementedError("Non-jet-fake background categories not implemented in this migration yet.")

    if channel == "mt":
        background_categories.extend(
            [
                (next(counter), f"{channel}_ttbar"),
                (next(counter), f"{channel}_dyjets"),
                (next(counter), f"{channel}_diboson"),
            ],
        )
        logger.info(f"Added mt-specific background categories: {background_categories[-3:]}")
    else:
        raise NotImplementedError(f"Channel-specific background categories for {channel} not implemented in this migration yet.")

    if categorization == "stxs_stage0":
        signal_categories = [(1, f"{channel}_vbf"), (2, f"{channel}_ggh")]
        logger.info(f"Using STXS stage 0 categorization with signal categories: {signal_categories}")
        return signal_categories + background_categories

    elif categorization == "stxs_stage1p2_syst":
        signal_categories = [
            (100, f"{channel}_vbf_bin201to210"),
            (101, f"{channel}_ggh_bin101to104"),
            (102, f"{channel}_ggh_bin105to106"),
            (103, f"{channel}_ggh_bin107to109"),
            (104, f"{channel}_ggh_bin110to116"),
        ]
        logger.info(f"Using STXS stage 1p2 syst categorization with signal categories: {signal_categories}")
        return signal_categories + background_categories

    elif categorization == "stxs_stage1p2_syst_bkg_only":
        logger.info(f"Using STXS stage 1p2 syst background-only categorization with background categories: {background_categories}")
        return background_categories

    else:
        raise ValueError(f"Unknown categorization: {categorization}")
