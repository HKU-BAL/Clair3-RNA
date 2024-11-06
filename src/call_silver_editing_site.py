import sys
import shlex
import logging
from subprocess import PIPE
from os.path import isfile
from argparse import ArgumentParser, SUPPRESS
from collections import Counter, defaultdict

import shared.param_p as param
from shared.interval_tree import bed_tree_from, is_region_in
from shared.utils import subprocess_popen, file_path_from, IUPAC_base_to_num_dict as BASE2NUM, region_from, \
    reference_sequence_from, str2bool, vcf_candidates_from, log_error

logging.getLogger().setLevel(logging.INFO)
BASES = ('A', 'T')
flanking_base_num = param.flankingBaseNum
sliding_window_size = no_of_positions = 2 * flanking_base_num + 1
BASE2NUMBER = dict(zip(
    "ACGTURYSWKMBDHVN-",
    (0, 1, 2, 3, 3, 0, 1, 1, 0, 2, 0, 1, 0, 0, 0, 0, 4)
))
channel = param.channel
channel_size = len(channel)
BASE2INDEX = dict(zip(channel, tuple(range(channel_size))))


def phredscore2raw_score(qual):
    return ord(qual) - 33


def get_positive_candidates_from_file(fn):
    candidates = set()
    if isfile(fn):
        with open(fn, 'r') as fp:
            for line in fp:
                contig, pos,_ , _ = line.strip().split()
                candidates.add(int(pos))
    return candidates

def evc_base_from(base):
    if base == 'N':
        return 'A'
    elif base == 'n':
        return 'a'
    elif base in 'ACGTacgt':
        return base
    elif base.isupper():
        return 'A'
    else:
        return 'a'

def get_base_list(columns, min_bp = 0):
    if len(columns) < 5:
        print("less than 5")
        return Counter(), []
    min_bq_cut = min_bp
    pileup_bases = columns[4]
    base_idx = 0
    base_list = []
    bq_list = [ord(qual) - 33 for qual in columns[5]]
    while base_idx < len(pileup_bases):
        base = pileup_bases[base_idx]
        if base == '+' or base == '-':
            base_idx += 1
            advance = 0
            while True:
                num = pileup_bases[base_idx]
                if num.isdigit():
                    advance = advance * 10 + int(num)
                    base_idx += 1
                else:
                    break
            try:
                base_list[-1][1] = base + pileup_bases[base_idx: base_idx + advance]
            except:
                print(columns)
                print(pileup_bases)
                print(base_idx)
                print(advance)
                print(f"lengh of base_list: {len(base_list)}")
                print(f"base_list: {base_list}")
                print(f"length of bq_list: {len(bq_list)}")

            base_idx += advance - 1

        elif base in "ACGTNacgtn#*":
            base_list.append([base, ""])
        elif base == '^':
            base_idx += 1
        base_idx += 1
    upper_base_counter = Counter([''.join(item).upper() for item, bq in zip(base_list, bq_list) if bq >= min_bq_cut])
    return upper_base_counter, base_list


def CreateTensorPileup(args):
    """
    Create pileup tensor for pileup model training or calling.
    Use slide window to scan the whole candidate regions, keep all candidates over specific minimum allelic frequency
    and minimum depth, use samtools mpileup to store pileup info for pileup tensor generation. Only scan candidate
    regions once, we could directly get all variant candidates directly.
    """
    fasta_file_path = args.ref_fn
    ctg_name = args.ctgName
    samtools_execute_command = args.samtools
    rna_bam_fn = args.rna_bam_fn
    dna_bam_fn = args.dna_bam_fn
    minimum_dna_dp = args.min_dna_dp
    maximum_dna_af = args.max_dna_af
    minimum_rna_dp = args.min_rna_dp
    minimum_rna_af = args.min_rna_af
    out_put_path = args.output_path

    extend_bed = args.extend_bed
    is_extend_bed_file_given = extend_bed is not None
    min_mapping_quality = param.min_mq
    min_base_quality = param.min_bq
    dsrna_tree = None
    if args.dsrna_bed:
        dsrna_tree= bed_tree_from(bed_file_path=args.dsrna_bed, contig_name=ctg_name)
    threshold = args.dsrna_threshold
    ref_regions = []
    reads_regions = []

    tree, bed_start, bed_end = bed_tree_from(bed_file_path=extend_bed,
                                             contig_name=ctg_name,
                                             return_bed_region=True)
    ctg_start = bed_start if is_extend_bed_file_given else None
    ctg_end = bed_end if is_extend_bed_file_given else None
    is_ctg_name_given = ctg_name is not None
    is_ctg_range_given = is_ctg_name_given and ctg_start is not None and ctg_end is not None

    if is_ctg_name_given:
        reads_regions.append(region_from(ctg_name=ctg_name))
        ref_regions.append(region_from(ctg_name=ctg_name))
        reference_start = 1

    reference_sequence = reference_sequence_from(
        samtools_execute_command=samtools_execute_command,
        fasta_file_path=fasta_file_path,
        regions=ref_regions
    )

    if reference_sequence is None or len(reference_sequence) == 0:
        sys.exit(log_error("[ERROR] Failed to load reference sequence from file ({}).".format(fasta_file_path)))
    mq_option = ' --min-MQ {}'.format(min_mapping_quality)
    bq_option = ' --min-BQ {}'.format(min_base_quality)
    flags_option = ' --excl-flags {}'.format(param.SAMTOOLS_VIEW_FILTER_FLAG)
    bed_option = ' -l {}'.format(extend_bed) if is_extend_bed_file_given else ""
    rna_samtools_mpileup_process = subprocess_popen(
        shlex.split(
            "{} mpileup  {} -r {} --reverse-del".format(samtools_execute_command,
                                                        rna_bam_fn,
                                                        " ".join(reads_regions), )
            + mq_option + bq_option + bed_option + flags_option  ))

    output_fn = out_put_path + "/" + ctg_name + "_edit_site"
    with open (output_fn, 'w') as output_fp:
        for row in rna_samtools_mpileup_process.stdout:
            alt_base = None
            rna_alt_base = None
            columns = row.strip().split('\t')
            pos = int(columns[1])
            reference_base = reference_sequence[pos - reference_start].upper()
            if reference_base not in BASES:
                continue
            if reference_base == "A":
                alt_base = "G"
                rna_alt_base = "G"
            elif reference_base == "T":
                alt_base = "C"
                rna_alt_base = "C"
            rna_counter, rna_base_list = get_base_list(columns, param.min_bq)
            rna_depth = len(rna_base_list)
            if rna_depth < minimum_rna_dp:
                continue
            rna_alt_cnt = rna_counter[rna_alt_base]
            rna_af = rna_alt_cnt / rna_depth
            is_in_dsrna = dsrna_tree and is_region_in(dsrna_tree,ctg_name, pos-1, pos+1)
            if is_in_dsrna:
                if not (rna_af > 0.4):
                    continue
            else:
                if not (rna_af > minimum_rna_af):
                    continue
            dna_samtools_mpileup_process = subprocess_popen(
                shlex.split(
                    "{} mpileup  {} -r {}:{}-{} --reverse-del".format(samtools_execute_command,
                                                                dna_bam_fn,
                                                                ctg_name,  pos, pos )
                    + mq_option + bq_option + bed_option + flags_option ))
            dna_output = dna_samtools_mpileup_process.stdout.read().rstrip()
            dna_columns = dna_output.split('\t')
            dna_base_counter, dna_base_list = get_base_list(dna_columns, param.min_bq)
            dna_depth = len(dna_base_list)
            dna_alt_depth = dna_base_counter[alt_base]
            if dna_depth < minimum_dna_dp:
                continue
            dna_af = dna_alt_depth / dna_depth
            if dna_af > maximum_dna_af:
                continue
            output_fp.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(ctg_name, pos, reference_base, alt_base, rna_depth, rna_af, dna_depth, dna_af, str(is_in_dsrna)))
def main():
    parser = ArgumentParser(description="Generate variant candidate tensors using pileup")


    parser.add_argument('--ctgName', type=str, default=None,
                        help="Contig name to process")
    parser.add_argument('--output_path', type=str, default=".",
                        help="Output path to store the tensor")
    parser.add_argument('--rna_bam_fn', type=str, default="input.bam", required=True,
                        help="Sorted RNA BAM file input, required")
    
    parser.add_argument('--dna_bam_fn', type=str, default="input.bam", required=True,
                        help="Sorted DNA BAM file input, required")
    parser.add_argument('--min_dna_dp', type=int, default=8,
                        help="Minimum depth for DNA")
    
    parser.add_argument('--max_dna_af', type=float, default=0.25,
                        help="Minimum AF for DNA")
    
    parser.add_argument('--min_rna_dp', type=int, default=8,
                        help="Minimum depth for RNA")
    
    parser.add_argument('--min_rna_af', type=float, default=0.75,
                        help="Minimum AF for RNA")
    
    parser.add_argument('--ref_fn', type=str, default="ref.fa", required=True,
                        help="Reference fasta file input, required")

    parser.add_argument('--extend_bed', type=str, default=None,
                        help="Call variant only in the provided regions. Will take an intersection if --ctgName and/or (--ctgStart, --ctgEnd) are set")

    parser.add_argument('--samtools', type=str, default="samtools",
                        help="Path to the 'samtools', samtools version >= 1.10 is required. default: %(default)s")
    
    parser.add_argument('--dsrna_bed', type=str, default=None,
                        help="Call variant only in the provided regions. Will take an intersection if --ctgName and/or (--ctgStart, --ctgEnd) are set")
    parser.add_argument('--dsrna_threshold', type=float, default=0.1,
                        help="Threshold for DSRNA")


    args = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)
    CreateTensorPileup(args)


if __name__ == "__main__":
    main()
