import ctypes
import math
from typing import Literal, Tuple, Optional

import ROOT

import CombineHarvester.CombineTools.ch as ch


def filter_zero_yield_processes(cb: ch.CombineHarvester) -> None:
    print("[INFO] Filtering processes with null yield...")
    def filter_proc(p):
        if (null_yield := not (p.rate() > 0.0)):
            print(f"[WARNING] Removing process {p.process()} in bin {p.bin()} with null yield")
            cb.FilterSysts(lambda s: s.bin() == p.bin() and s.process() == p.process() and s.era() == p.era())
        return null_yield

    cb.FilterProcs(filter_proc)


def fix_negative_bins(cb: ch.CombineHarvester) -> None:
    print("[INFO] Fixing negative bins to zero...")
    def zero_negative_bins_th1(hist):
        if not hist:
            return False
        has_negative = False
        for i in range(1, hist.GetNbinsX() + 1):
            if hist.GetBinContent(i) < 0.0:
                has_negative = True
                hist.SetBinContent(i, 0.0)
        return has_negative

    cb.ForEachProc(lambda p: zero_negative_bins_th1(p.shape()))

    def fix_syst(s):
        if (s.type() == "shape") and (hist_up := s.shape_u()):
            zero_negative_bins_th1(hist_up)
        if (s.type() == "shape") and (hist_down := s.shape_d()):
            zero_negative_bins_th1(hist_down)
                
    cb.ForEachSyst(fix_syst)


def filter_zero_yield_systs(cb: ch.CombineHarvester) -> None:
    print("[INFO] Filtering shape systematics with null yield...")
    def hist_integral(hist):
        if not hist:
            return None
        return hist.Integral(1, hist.GetNbinsX())

    def filter_syst(s):
        if s.type() == "shape":
            shape_u, shape_d = s.shape_u(), s.shape_d()
            if not shape_u or not shape_d:
                print(f"  [WARNING] Removing systematic {s.name()} on {s.process()} with missing shape")
                return True

            yield_u, yield_d = hist_integral(shape_u), hist_integral(shape_d)

            if yield_u is None or yield_d is None:
                return True

            if yield_u == 0.0 or yield_d == 0.0:
                print(f"  [WARNING] Removing systematic {s.name()} on {s.process()} with null yield in shift")
                return True
        return False
    cb.FilterSysts(filter_syst)


def load_systematic_shapes(cb: ch.CombineHarvester, root_file: str, channel: str) -> None:
    print(f"[INFO] Loading shape systematics for {channel} from {root_file}")
    tfile = ROOT.TFile.Open(root_file, "READ")

    def has_shapes(s):
        mass, process = s.mass(), s.process()
        process = f"{process}{mass}" if mass and mass != "*" else process
        base = f"{s.bin()}/{process}_{s.name()}"
        return bool(tfile.Get(base + "Up")) and bool(tfile.Get(base + "Down"))

    def filter_missing(s):
        if s.channel() != channel:
            return False
        if s.type() not in ("shape",):
            return False
        if not has_shapes(s):
            print(f"  [WARNING] Removing systematic {s.name()} on {s.process()} with missing shape")
            return True
        return False

    cb.FilterSysts(filter_missing)
    tfile.Close()

    cb.cp().channel([channel]).backgrounds().ExtractShapes(root_file, "$BIN/$PROCESS", "$BIN/$PROCESS_$SYSTEMATIC")
    cb.cp().channel([channel]).signals().ExtractShapes(root_file, "$BIN/$PROCESS$MASS", "$BIN/$PROCESS$MASS_$SYSTEMATIC")


def convert_shapes_to_lnN(cb: ch.CombineHarvester) -> None:
    # Checks shape systematics. If the normalization shift is smaller than the statistical
    # uncertainty of the template, it replaces the shape systematic with a symmetrized lnN systematic.

    print("[INFO] Checking shape systematics for lnN conversion...")    
    count = {"lnN": 0, "all": 0}

    def check_and_convert(s):
        name = s.name()
        if any(substring in name for substring in["scale", "CMS_htt_boson_reso_met", "res_j", "res_e"]):
            count["all"] += 1
            shape_u, shape_d = s.shape_u(), s.shape_d()
            
            if not shape_u or not shape_d:
                return
                
            nbins = shape_u.GetNbinsX()
            err_u_ref, err_d_ref = ctypes.c_double(0.0), ctypes.c_double(0.0)
            
            yield_u = shape_u.IntegralAndError(1, nbins, err_u_ref)
            yield_d = shape_d.IntegralAndError(1, nbins, err_d_ref)
            
            err_u, err_d = err_u_ref.value, err_d_ref.value
            
            if yield_u <= 0.0 or yield_d <= 0.0:
                return

            value_u, value_d = s.value_u(), s.value_d()

            # Is the shift smaller than the statistical uncertainty of the shift?
            if abs(value_u - 1.0) + abs(value_d - 1.0) < (err_u / yield_u) + (err_d / yield_d):
                count["lnN"] += 1
                print(f"  [WARNING] Replacing systematic by lnN: {name} (bin: {s.bin()}, proc: {s.process()})")

                s.set_type("lnN")
                up_is_larger = (value_u > value_d)
                
                if value_u < 1.0:
                    value_u = 1.0 / value_u
                if value_d < 1.0:
                    value_d = 1.0 / value_d

                if up_is_larger:  # Symmetrize
                    value_u = math.sqrt(value_u * value_d)
                    value_d = 1.0 / value_u
                else:
                    value_d = math.sqrt(value_u * value_d)
                    value_u = 1.0 / value_d
                    
                s.set_value_u(value_u)
                s.set_value_d(value_d)

    cb.cp().ForEachSyst(check_and_convert)
    print(f"[WARNING] Turned {count['lnN']} of {count['all']} checked systematics into lnN.")


def replace_with_asimov(cb: ch.CombineHarvester) -> None:
    print("[INFO] Replacing observation with Asimov dataset...")

    def is_empty_shape(hist) -> bool:
        return not hist or (hist.GetNbinsX() == 1 and hist.Integral() == 0.0)

    for category in cb.cp().bin_set():
        category_bin = cb.cp().bin([category])
        bkg_shape = category_bin.backgrounds().GetShape()
        sig_shape = category_bin.signals().GetShape()
        has_bkg = not is_empty_shape(bkg_shape)
        has_sig = not is_empty_shape(sig_shape)

        if not has_bkg and not has_sig:
            print(f"  [WARNING] No signal and no background available in bin {category}")
            continue

        template = bkg_shape if has_bkg else sig_shape
        asimov_shape = template.Clone()
        asimov_shape.SetDirectory(0)
        asimov_shape.Reset("ICES")

        if has_bkg:
            asimov_shape.Add(bkg_shape)
        else:
            print(f"  [WARNING] No background available in bin {category}")

        if has_sig:
            asimov_shape.Add(sig_shape)
        else:
            print(f"  [WARNING] No signal available in bin {category}")

        category_bin.ForEachObs(lambda obs: obs.set_shape(asimov_shape, True))


def apply_nn_rebinning(
    cb: ch.CombineHarvester,
    strategy: Literal["total_bkg", "per_process", "total_bkg__per_process"] = "total_bkg",
    threshold: float = 10.0,
    secondary_threshold: Optional[float] = 1.0,
) -> None:
    # Merges from left (0 -> peak) and right (1 -> peak) until conditions are met.

    for category in cb.cp().bin_set():
        category_bin, bkg_shapes = cb.cp().bin([category]), []

        if strategy == "total_bkg":
            bkg_shapes.append(category_bin.backgrounds().GetShape())

            def is_valid(yields):
                return yields[0] >= threshold

        elif strategy == "per_process":
            for p_name in category_bin.backgrounds().process_set():
                bkg_shapes.append(category_bin.cp().process([p_name]).GetShape())

            def is_valid(yields):
                has_relevant = False
                for y, shape in zip(yields, bkg_shapes):
                    if shape.Integral() >= threshold:
                        has_relevant = True
                        if y < threshold:
                            return False
                return has_relevant

        elif strategy == "total_bkg__per_process":
            if secondary_threshold is None:
                raise ValueError("secondary_threshold required for 'total_bkg__per_process'")
            for p_name in category_bin.backgrounds().process_set():
                shape = category_bin.cp().process([p_name]).GetShape()
                if shape.Integral() >= secondary_threshold:
                    bkg_shapes.append(shape)

            def is_valid(yields):
                if not bkg_shapes:  # no relevant process
                    return False
                if sum(yields) < threshold:
                    return False
                for y in yields:
                    if y < secondary_threshold:
                        return False
                return True

        else:
            raise ValueError(f"Unknown rebinning strategy: {strategy}")

        if not bkg_shapes or not bkg_shapes[0]:
            continue

        ref_shape = bkg_shapes[0]
        nbins = ref_shape.GetNbinsX()
        peak_bin = category_bin.backgrounds().GetShape().GetMaximumBin()
        peak_bin = max(1, min(peak_bin, nbins))  # Protect boundaries

        edges = set()
        edges.add(ref_shape.GetBinLowEdge(1))  # Absolute min (0.0)
        edges.add(ref_shape.GetBinLowEdge(nbins + 1))  # Absolute max (1.0)

        current_yields = [0.0] * len(bkg_shapes)  # 0.0 (1) to peak (peak_bin - 1)
        for i in range(1, peak_bin):
            for j, shape in enumerate(bkg_shapes):
                current_yields[j] += shape.GetBinContent(i)

            if is_valid(current_yields):
                edges.add(ref_shape.GetBinLowEdge(i + 1))
                current_yields = [0.0] * len(bkg_shapes) # Reset yields for next bin

        current_yields = [0.0] * len(bkg_shapes)  # 1.0 (nbins) to peak (peak_bin + 1)
        for i in range(nbins, peak_bin, -1):
            for j, shape in enumerate(bkg_shapes):
                current_yields[j] += shape.GetBinContent(i)

            if is_valid(current_yields):
                edges.add(ref_shape.GetBinLowEdge(i))
                current_yields = [0.0] * len(bkg_shapes) # Reset yields for next bin

        final_edges = sorted(list(edges))
        cb.cp().bin([category]).VariableRebin(final_edges)  # leftovers absorbed by the peak since no internal edge was placed


def scale_higgs_mass_processes(
    cb: ch.CombineHarvester, 
    apply_scaling: bool = False, 
    scale_map: Tuple[Tuple[str, float], ...] = (
        ("ggH.*htt", 0.984),
        ("qqH.*htt", 0.987),
        ("WH.*htt", 0.979),
        ("ZH.*htt", 0.982),
        ("ggH.*hww", 1.025),
        ("qqH.*hww", 1.028),
        ("WH.*hww", 1.020),
        ("ZH.*hww", 1.022),
    ),
) -> None:
    if not apply_scaling:
        return
    
    print("[INFO] Scaling Higgs mass processes...")
    for proc_rgx, scale_factor in scale_map:
        cb.cp().process_rgx([proc_rgx]).ForEachProc(lambda p: p.set_rate(p.rate() * scale_factor))


def scale_2016_lumi(
    cb: ch.CombineHarvester, 
    era: int, 
    apply_scaling: bool = False, 
    scale_factor: float = 1.0128
) -> None:
    if str(era) == "2016" and apply_scaling:
        print(f"[INFO] Updating nominal lumi for 2016 MC by a factor of {scale_factor}...")
        cb.cp().process(["EMB", "QCD", "jetFakes"], False).ForEachProc(lambda p: p.set_rate(p.rate() * scale_factor))  # Grab everything NOT data-driven
