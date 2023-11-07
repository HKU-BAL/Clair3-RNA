# BSD 3-Clause License
#
# Copyright 2023 The University of Hong Kong, Department of Computer Science
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import sys
import os
import subprocess
import concurrent.futures
import shlex

from collections import Counter
from argparse import ArgumentParser, SUPPRESS
from collections import defaultdict

from shared.vcf import VcfReader
import shared.param_p as param

from shared.interval_tree import bed_tree_from, is_region_in

from shared.utils import subprocess_popen, file_path_from, region_from, \
    reference_sequence_from, str2bool, str_none

def get_base_list(columns, args=None):
    if len(columns) < 5:
        return Counter(), []
    min_bq_cut = args.min_bq_cut if args is not None else 0
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
                base_list[-1][1] = base + pileup_bases[base_idx: base_idx + advance]  # add indel seq
            except:
                print(columns)
            base_idx += advance - 1

        elif base in "ACGTNacgtn#*":
            base_list.append([base, ""])
        elif base == '^':  # start of read, next base is mq, update mq info
            base_idx += 1
        # skip $, the end of read
        base_idx += 1
    upper_base_counter = Counter([''.join(item).upper() for item, bq in zip(base_list, bq_list) if bq >= min_bq_cut])
    return upper_base_counter, base_list

#
# def extract_base(POS):
#
#     pos = POS.pos
#     ref_base = POS.reference_bases
#     alt_base = POS.alternate_bases[0].upper()
#     ctg_name = POS.ctg_name
#
#     match_alt_base = alt_base
#     if len(ref_base) == 1 and len(alt_base) > 1:
#         match_alt_base = alt_base[0].upper() + '+' + alt_base[1:].upper()
#     if len(ref_base) > 1 and len(alt_base) == 1:
#         match_alt_base = ref_base[0].upper() + '-' + len(ref_base[1:]) * "N"
#
#     base_list = []
#     alt_count = 0
#
#
#     samtools_command_with_region = samtools_command + ' -r {}:{}-{}'.format(ctg_name, pos, pos)
#     output = subprocess.run(samtools_command_with_region, shell=True, stdout=subprocess.PIPE,
#                                   stderr=subprocess.PIPE, universal_newlines=True)
#     output = output.stdout.rstrip()
#     columns = output.split('\t')
#     base_counter, base_list = get_base_list(columns)
#     alt_count = base_counter[match_alt_base]
#
#     return [ctg_name, pos, len(base_list), alt_count]
#
#
# class INFO():
#     def __init__(self):
#         self.normal_base_counter = None
#         self.normal_depth = None
#         self.tumor_base_counter = None
#         self.tumor_depth = None
#         self.alt_base = None

#
# def parser_info(row):
#     columns = row.split('\t')
#     info_list = columns[8].split(':')
#     FORMAT = columns[9].split(':')
#     DP_index = info_list.index("DP")
#     NDP_index = info_list.index('NDP')
#     AF_index = info_list.index('AF')
#     NAF_index = info_list.index('NAF')
#     tumor_depth = int(FORMAT[DP_index])
#     normal_depth = int(FORMAT[NDP_index])
#     tumor_alt_depth = round(float(FORMAT[AF_index]) * tumor_depth)
#     normal_alt_depth = round(float(FORMAT[NAF_index]) * normal_depth)
#
#     return columns[0], columns[1], normal_depth, tumor_depth, normal_alt_depth, tumor_alt_depth


def cal_af(args, truth_variant_dict=None, input_variant_dict=None):
    ctg_name = args.ctg_name
    output_path = args.output_path
    chunk_id = args.chunk_id - 1 if args.chunk_id else None  # 1-base to 0-base
    chunk_num = args.chunk_num

    bed_fn = args.bed_fn

    bed_tree = bed_tree_from(bed_file_path=bed_fn, contig_name=ctg_name)

    if truth_variant_dict is None:
        truth_vcf_fn = args.truth_vcf_fn
        vcf_reader = VcfReader(vcf_fn=truth_vcf_fn,
                               ctg_name=ctg_name,
                               show_ref=False,
                               keep_row_str=True,
                               filter_tag=args.truth_filter_tag)
        vcf_reader.read_vcf()
        truth_variant_dict = vcf_reader.variant_dict

    results_dict = defaultdict()
    variant_dict = defaultdict()
    #make sure all position in the BED regions
    for k, v in truth_variant_dict.items():
        pos = k if ctg_name is not None else k[1]
        ctg = ctg_name if ctg_name is not None else k[0]
        # ref_base = v.alternate_bases
        if bed_fn is not None and not is_region_in(bed_tree, ctg, pos):
            continue
        variant_dict[k] = v

    if output_path is not None:
        output_dir = os.path.dirname(output_path)
        if not os.path.exists(output_dir):
            subprocess.run("mkdir -p {}".format(output_dir), shell=True)
        output_file = open(output_path, 'w')
    output_dir = os.path.dirname(output_path)

    min_mq = args.min_mq
    min_bq = args.min_bq

    pos_list = sorted(list(variant_dict.keys()))
    total_variants_size = len(pos_list)
    chunk_variants_size = total_variants_size // chunk_num if total_variants_size % chunk_num == 0 else total_variants_size // chunk_num + 1
    chunk_start_pos = chunk_id * chunk_variants_size
    chunk_variants_list = pos_list[chunk_start_pos: chunk_start_pos + chunk_variants_size]
    if len(chunk_variants_list) == 0:
        return
    ctg_start, ctg_end = min(chunk_variants_list) - 100, max(chunk_variants_list) + 100

    bed_path = os.path.join(output_dir, "bed", '{}_{}.bed'.format(ctg_name, chunk_id))
    if not os.path.exists(os.path.join(output_dir, 'bed')):
        output = subprocess.run("mkdir -p {}".format(os.path.join(output_dir, 'bed')), shell=True)
    output_bed = open(bed_path, 'w')
    for pos in sorted(list(chunk_variants_list)):
        output_bed.write('\t'.join([ctg_name, str(pos - 1), str(pos)]) + '\n')
    output_bed.close()

    read_region = "{}:{}-{}".format(ctg_name, ctg_start, ctg_end)

    phasing_option = "--output-extra HP " if args.phase_output else " "
    samtools_command = "{} mpileup --min-MQ {} --min-BQ {} --excl-flags 2316 {} -l {} -r {} {}".format(args.samtools,
                                                                                          min_mq,
                                                                                          min_bq,
                                                                                          phasing_option,
                                                                                          bed_path,
                                                                                          read_region,
                                                                                          args.bam_fn)


    # global samtools_command
    # samtools_command = (samtools_command + args.bam_fn) if args.bam_fn is not None else None

    total_num = 0
    print("[INFO] Total truth need to calculate AF: {}, ctg/chunk:{}/{}".format(len(chunk_variants_list), ctg_name, chunk_id))

    samtools_mpileup_process = subprocess_popen(
        shlex.split(samtools_command), stdin=sys.stdin, stderr=subprocess.PIPE)

    for row in samtools_mpileup_process.stdout:  # chr position N depth seq BQ read_name mapping_quality phasing_info
        columns = row.strip().split('\t')
        pos = int(columns[1])

        if pos not in variant_dict:
            continue

        POS = variant_dict[pos]
        pos = POS.pos
        ref_base = POS.reference_bases
        alt_base = POS.alternate_bases[0].upper()
        ctg_name = POS.ctg_name

        match_alt_base = alt_base
        if len(ref_base) == 1 and len(alt_base) > 1:
            match_alt_base = alt_base[0].upper() + '+' + alt_base[1:].upper()
        if len(ref_base) > 1 and len(alt_base) == 1:
            match_alt_base = ref_base[0].upper() + '-' + len(ref_base[1:]) * "N"

        base_list = []
        HAP_LIST = [0, 0, 0]
        ALL_HAP_LIST = [0, 0, 0]

        tumor_base_counter, base_list = get_base_list(columns)
        alt_count = tumor_base_counter[match_alt_base]

        # if len(columns) >= 7:
        #     phasing_info = columns[6].split(',')
        #     for hap_idx, (b, hap) in enumerate(zip(base_list, phasing_info)):
        #         if hap not in '12':
        #             hap = 0
        #         ALL_HAP_LIST[int(hap)] += 1
        #         if ''.join(b).upper() == match_alt_base:
        #             HAP_LIST[int(hap)] += 1
        #
        # HAP_LIST = " ".join([str(i) for i in HAP_LIST])
        # ALL_HAP_LIST = " ".join([str(i) for i in ALL_HAP_LIST])

        result = [ctg_name, pos, len(base_list), alt_count]

        results_dict[pos] = result
        total_num += 1
        # if total_num % 2000 == 0 and total_num > 0:
        #     print("[INFO] Total processed positions: {}".format(total_num))

    for pos in chunk_variants_list:
        base_list = []
        normal_alt_count = 0
        base_list = []
        tumor_alt_count = 0
        HAP_LIST = [0, 0, 0]
        ALL_HAP_LIST = [0, 0, 0]
        HAP_LIST = " ".join([str(i) for i in HAP_LIST])
        ALL_HAP_LIST = " ".join([str(i) for i in ALL_HAP_LIST])
        result = [ctg_name, pos, len(base_list), alt_count]

        if pos in results_dict:
            result = results_dict[pos]

            if output_path is not None:
                output_file.write(' '.join(str(item) for item in result) + '\n')

    if output_path is not None:
        output_file.close()
        return

def main():
    parser = ArgumentParser(description="Calculate AF distribution in tumor and normal BAMs")


    parser.add_argument('--bam_fn', type=str, default=None,
                        help="Sorted tumor BAM file input, required")

    parser.add_argument('--input_vcf_fn', type=str, default=None,
                        help="Input VCF filename")

    parser.add_argument('--bed_fn', type=str, default=None,
                        help="Input BED filename")

    parser.add_argument('--truth_vcf_fn', type=str, default=None,
                        help="Input truth VCF filename")

    parser.add_argument('--input_filter_tag', type=str, default=None,
                        help="Filter variants with tag from the input VCF")

    parser.add_argument('--truth_filter_tag', type=str, default=None,
                        help="Filter variants with tag from the truth VCF")

    parser.add_argument('--ctg_name', type=str, default=None,
                        help="Output VCF filename, required")

    parser.add_argument('--output_path', type=str, default=None,
                        help="Output VCF filename, required")

    parser.add_argument('--threads', type=int, default=4,
                        help="Max #threads to be used")

    parser.add_argument('--phase_output', type=str2bool, default=False,
                        help="Output phasing INFO")

    parser.add_argument('--samtools', type=str, default="samtools",
                        help="Absolute path to the 'samtools', samtools version >= 1.10 is required. Default: %(default)s")

    # options for advanced users
    parser.add_argument('--min_mq', type=int, default=param.min_mq,
                        help="EXPERIMENTAL: If set, reads with mapping quality with <$min_mq are filtered. Default: %(default)d")

    parser.add_argument('--min_bq', type=int, default=param.min_bq,
                        help="EXPERIMENTAL: If set, bases with base quality with <$min_bq are filtered. Default: %(default)d")

    parser.add_argument('--min_bq_cut', type=int, default=0,
                        help="EXPERIMENTAL: Minimal base quality cut-off")

    ## The number of chucks to be divided into for parallel processing
    parser.add_argument('--chunk_num', type=int, default=None,
                        help=SUPPRESS)

    ## The chuck ID to work on
    parser.add_argument('--chunk_id', type=int, default=None,
                        help=SUPPRESS)

    ## The chuck ID to work on
    parser.add_argument('--debug', type=str2bool, default=False,
                        help=SUPPRESS)

    global args
    args = parser.parse_args()

    # if args.debug:
    #     check_equal_output()
    # else:
    cal_af(args)


if __name__ == "__main__":
    main()

# /autofs/bal36/zxzheng/env/miniconda3/envs/clair3/bin/pypy3 /autofs/bal36/zxzheng/somatic/Clair-somatic/clairs.py cal_af_distribution_parallel --debug 1