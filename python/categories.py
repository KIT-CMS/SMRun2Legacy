from itertools import count

def get_categories(
    channel: str,
    categorization: str = "gof",
    gof_category_name: str = "gof",
    embedding: bool = True,
    jetfakes: bool = True,
) -> list:
    counter = count(11) # Start category IDs from 11
    background_categories = []

    if embedding:
        background_categories.append((next(counter), f"{channel}_embedding"))

    if jetfakes:
        background_categories.append((next(counter), f"{channel}_jetFakes"))

    if channel == "mt":
        background_categories.extend(
            [
                (next(counter), f"{channel}_ttbar"),
                (next(counter), f"{channel}_dyjets"),
                (next(counter), f"{channel}_diboson"),
            ],
        )

    if categorization == "gof":
        return [(300, gof_category_name)]
        
    elif categorization == "stxs_stage0":
        return [(1, f"{channel}_vbf"), (2, f"{channel}_ggh")] + background_categories
        
    elif categorization == "stxs_stage1p2_syst":
        return [
            (100, f"{channel}_vbf_bin201to202"),
            (101, f"{channel}_vbf_bin203to210"),
            (102, f"{channel}_ggh_bin101to104"),
            (103, f"{channel}_ggh_bin105to106"),
            (104, f"{channel}_ggh_bin107to109"),
            (105, f"{channel}_ggh_bin110to116"),
        ] + background_categories
        
    else:
        raise ValueError(f"Unknown categorization: {categorization}")
