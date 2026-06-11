import argparse
import glob
import json
import logging
import os
import sys

cmssw_base = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../"))  # CMSSW_BASE/src/CombineHarvester/SMRun2Legacy/scripts

if os.path.exists(cmssw_python := os.path.join(cmssw_base, "python")) and cmssw_python not in sys.path:
    sys.path.insert(0, cmssw_python)

for p in [os.path.dirname(os.path.abspath(__file__)), os.getcwd(), '']:
    while p in sys.path:
        sys.path.remove(p)

import CombineHarvester.CombineTools.ch as ch
from CombineHarvester.SMRun2Legacy.categories import get_categories
from CombineHarvester.SMRun2Legacy.custom_logging import setup_logging
from CombineHarvester.SMRun2Legacy.processes import (get_backgrounds,
                                                     get_signals)
from CombineHarvester.SMRun2Legacy.systematics import add_systematics
from CombineHarvester.SMRun2Legacy.tools import (apply_nn_rebinning,
                                                 convert_shapes_to_lnN,
                                                 filter_zero_yield_processes,
                                                 filter_zero_yield_systs,
                                                 fix_negative_bins,
                                                 load_systematic_shapes,
                                                 replace_with_asimov)


def str2bool(v):
    if isinstance(v, bool):
        return v
    if str(v).lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif str(v).lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


parser = argparse.ArgumentParser()
parser.add_argument("--era", type=int, default=2018)
parser.add_argument("--channels", type=str, default="mt", help="Comma-separated list of channels")
parser.add_argument("--categories", type=str, default="gof", help="Categorization scheme (gof, stxs_stage0, stxs_stage1p2_syst etc.)")
parser.add_argument("--gof-category-name", type=str, default="gof")
parser.add_argument("--stxs-signals", type=str, default="stxs_stage0")

parser.add_argument("--real-data", type=str2bool, default=False, help="Use data_obs instead of generating Asimov")
parser.add_argument("--embedding", type=str2bool, default=True, help="Use embedded samples")
parser.add_argument("--jetfakes", type=str2bool, default=True, help="Use jetFakes (FF method)")
parser.add_argument("--bbb", type=str2bool, default=True, help="Apply autoMCStats (Bin-by-Bin) uncertainties")

parser.add_argument("--ggh-wg1", type=str2bool, default=True)
parser.add_argument("--qqh-wg1", type=str2bool, default=True)
parser.add_argument("--split-tau-id-and-es-by-pt", type=str2bool, default=True)
parser.add_argument("--use-ml-ff-scheme", type=str2bool, default=True)
parser.add_argument("--correlate-emb", type=str2bool, default=True)
parser.add_argument("--regional-jec", type=str2bool, default=True)
parser.add_argument("--convert-shapes-to-lnN", type=str2bool, default=True, help="Symmetrize noise shapes to lnN")

parser.add_argument("--rebinning-strategy", type=str, default="none", choices=["none", "total_bkg", "per_process", "total_bkg__per_process", "combine"])

parser.add_argument("--rebinning-threshold", type=float, default=10.0)

parser.add_argument("--rebinning-using-combine", type=str2bool, default=False)
parser.add_argument("--rebinning-using-combine-threshold", type=float, default=10.0)
parser.add_argument("--rebinning-using-combine-uncert-fraction", type=float, default=0.1)
parser.add_argument("--rebinning-using-combine-mode", type=int, default=1)

parser.add_argument("--custom-binning-file", type=str, default=None, help="Path to JSON file containing custom binning definitions")

parser.add_argument("--base-path", type=str, default=os.path.join(cmssw_base, "src/CombineHarvester/SMRun2Legacy/shapes"))
parser.add_argument("--input-folder-mt", type=str, default="shapes")
parser.add_argument("--input-folder-et", type=str, default="shapes")
parser.add_argument("--input-folder-tt", type=str, default="shapes")
parser.add_argument("--input-folder-em", type=str, default="shapes")
parser.add_argument("--output-folder", type=str, default="output_cards")
parser.add_argument("--postfix", type=str, default="-ML")

parser.add_argument("--log-level", type=str, default="INFO")

args = parser.parse_args()
logger = setup_logging(logger=logging.getLogger(__name__), level=getattr(logging, args.log_level.upper()))
logger.info(f"Starting datacard generation with arguments: {args}")

if __name__ == "__main__":
    cb = ch.CombineHarvester()

    channels, era, masses = args.channels.split(','), str(args.era), ["125"]

    for channel in channels:
        categories = get_categories(
            channel=channel,
            categorization=args.categories,
            gof_category_name=args.gof_category_name,
            embedding=args.embedding,
            jetfakes=args.jetfakes,
        )
        bkgs = get_backgrounds(channel=channel, embedding=args.embedding, jetfakes=args.jetfakes)
        sigs = get_signals(stxs_version=args.stxs_signals)

        logger.info(f"Initializing channel {channel} with {len(categories)} categories: {categories}")
        cb.AddObservations(["*"], ["htt"], [era], [channel], categories)
        cb.AddProcesses(["*"], ["htt"], [era], [channel], bkgs, categories, False)
        cb.AddProcesses(masses, ["htt"], [era], [channel], sigs, categories, True)

    input_folders = {
        "mt": args.input_folder_mt,
        "et": args.input_folder_et,
        "tt": args.input_folder_tt,
        "em": args.input_folder_em,
    }

    base_path, root_files = args.base_path, {}
    for channel in channels:
        input_folder = input_folders.get(channel, "shapes")
        root_file = os.path.join(base_path, input_folder, f"htt_{channel}.inputs-sm-Run{era}{args.postfix}.root")
        root_files[channel] = root_file
        logger.info(f"Extracting shapes for {channel} from {root_file}")

        cb.cp().channel([channel]).backgrounds().ExtractShapes(root_file, "$BIN/$PROCESS", "$BIN/$PROCESS_$SYSTEMATIC")
        cb.cp().channel([channel]).signals().ExtractShapes(root_file, "$BIN/$PROCESS$MASS", "$BIN/$PROCESS$MASS_$SYSTEMATIC")

    fix_negative_bins(cb)
    filter_zero_yield_processes(cb)

    if not args.real_data:
        replace_with_asimov(cb)

    logger.info("Adding Systematics...")
    add_systematics(
        cb,
        era=int(era),
        jetfakes=args.jetfakes,
        embedding=args.embedding,
        split_tau_id_and_es_by_pt=args.split_tau_id_and_es_by_pt,
        use_ml_ff_scheme=args.use_ml_ff_scheme,
        correlate_emb=args.correlate_emb,
        regional_jec=args.regional_jec,
        ggh_wg1=args.ggh_wg1,
        qqh_wg1=args.qqh_wg1,
    )

    for channel in channels:
        load_systematic_shapes(cb, root_files[channel], channel)

    filter_zero_yield_systs(cb)

    ch.SetStandardBinNames(cb, "$ANALYSIS_$CHANNEL_$BINID_$ERA")

    custom_binnings = None
    if args.custom_binning_file:
        with open(args.custom_binning_file, "r") as f:
            custom_binnings = json.load(f)

    if args.categories != "gof":
        if args.rebinning_strategy == "none" and custom_binnings is not None:
            for category, edges in custom_binnings.items():
                if category in cb.cp().bin_set():
                    logger.info(f"Applying dedicated custom binning to {category}: {edges}")
                    cb.cp().bin([category]).VariableRebin(edges)

        elif args.rebinning_strategy in {"total_bkg", "per_process", "total_bkg__per_process", "combine"}:
            if args.rebinning_strategy == "combine":
                logger.info("Applying rebinning using Combine's Rebin method with custom overrules...")
                rebinner = (
                    ch.AutoRebin()
                    .SetBinThreshold(args.rebinning_threshold)
                    .SetBinUncertFraction(args.rebinning_using_combine_uncert_fraction)
                    .SetRebinMode(args.rebinning_using_combine_mode)
                    .SetPerformRebin(True)
                    .SetVerbosity(1)
                )
                
                # Apply per-category: overrule or run Combine's AutoRebin
                for category in cb.cp().bin_set():
                    if custom_binnings is not None and category in custom_binnings:
                        logger.info(f"Overruling category {category} with dedicated custom binning: {custom_binnings[category]}")
                        cb.cp().bin([category]).VariableRebin(custom_binnings[category])
                    else:
                        category_cb = cb.cp().bin([category])
                        rebinner.Rebin(category_cb, cb)
            else:
                apply_nn_rebinning(
                    cb,
                    strategy=args.rebinning_strategy,
                    threshold=args.rebinning_threshold,
                    custom_binnings=custom_binnings,
                )
        else:
            raise ValueError(f"Invalid rebinning strategy: {args.rebinning_strategy}")

    if args.convert_shapes_to_lnN:
        convert_shapes_to_lnN(cb)

    output_dir = args.output_folder
    logger.info(f"Writing datacards to {output_dir}")
    os.makedirs(output_dir, exist_ok=True)

    if args.bbb:
        logger.info("Adding AutoMCStats...")
        cb.SetAutoMCStats(cb, 0.0)

    writer = ch.CardWriter(
        os.path.join(output_dir, "$TAG/$MASS/$BIN.txt"),
        os.path.join(output_dir, f"$TAG/common/htt_input_{era}.root")
    )
    writer.WriteCards("cmb", cb)
    for channel in channels:
        writer.WriteCards(channel, cb.cp().channel([channel]))

    if args.bbb:
        for card in glob.glob(os.path.join(output_dir, "**/*.txt"), recursive=True):
            with open(card, "r", encoding="utf-8") as handle:
                lines = [line for line in handle if "autoMCStats" not in line]
            lines.append(" * autoMCStats 0.0\n")
            with open(card, "w", encoding="utf-8") as handle:
                handle.writelines(lines)

    logger.info("Done!")
