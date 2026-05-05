import os
import sys
import argparse
import glob
import logging

cmssw_base = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../../../"))  # CMSSW_BASE/src/CombineHarvester/SMRun2Legacy/scripts

if os.path.exists(cmssw_python := os.path.join(cmssw_base, "python")) and cmssw_python not in sys.path:
    sys.path.insert(0, cmssw_python)

for p in[os.path.dirname(os.path.abspath(__file__)), os.getcwd(), '']:
    while p in sys.path:
        sys.path.remove(p)

import CombineHarvester.CombineTools.ch as ch

from CombineHarvester.SMRun2Legacy.processes import get_backgrounds, get_signals
from CombineHarvester.SMRun2Legacy.categories import get_categories
from CombineHarvester.SMRun2Legacy.systematics import add_systematics
from CombineHarvester.SMRun2Legacy.tools import (
    filter_zero_yield_processes,
    filter_zero_yield_systs,
    fix_negative_bins,
    convert_shapes_to_lnN,
    apply_nn_rebinning,
    scale_higgs_mass_processes,
    scale_2016_lumi,
    replace_with_asimov,
    load_systematic_shapes,
)
from CombineHarvester.SMRun2Legacy.custom_logging import setup_logging


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

parser.add_argument("--rebinning-strategy", type=str, default="none", choices=["none", "total_bkg", "per_process", "total_bkg__per_process"])
parser.add_argument("--rebinning-threshold", type=float, default=10.0)
parser.add_argument("--apply-mass-scaling", type=str2bool, default=False)
parser.add_argument("--apply-2016-lumi-scaling", type=str2bool, default=False)

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
        cb.AddProcesses(["*"], ["htt"], [era],[channel], bkgs, categories, False)
        cb.AddProcesses(masses, ["htt"],[era], [channel], sigs, categories, True)

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

    mh_125_00_to_125_11_map = (("ggH.*htt", float("nan")), ("qqH.*htt", float("nan")))
    scale_higgs_mass_processes(cb, apply_scaling=args.apply_mass_scaling, scale_map=mh_125_00_to_125_11_map)
    scale_2016_lumi(cb, era=era, apply_scaling=args.apply_2016_lumi_scaling) 

    fix_negative_bins(cb)
    filter_zero_yield_processes(cb)

    if not args.real_data:
        replace_with_asimov(cb)

    logger.info(f"Adding Systematics...")
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

    if args.categories != "gof" and args.rebinning_strategy != "none":
        apply_nn_rebinning(cb, strategy=args.rebinning_strategy, threshold=args.rebinning_threshold)

    if args.convert_shapes_to_lnN:
        convert_shapes_to_lnN(cb)

    output_dir = args.output_folder
    logger.info(f"Writing datacards to {output_dir}")
    os.makedirs(output_dir, exist_ok=True)
    
    ch.SetStandardBinNames(cb, "$ANALYSIS_$CHANNEL_$BINID_$ERA")

    if args.bbb:
        logger.info(f"Adding AutoMCStats...")
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
        
    logger.info(f"Done!")
