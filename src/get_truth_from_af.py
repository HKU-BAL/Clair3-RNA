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

from collections import Counter
from argparse import ArgumentParser
from collections import defaultdict

from shared.vcf import VcfReader
import shared.param as param
from shared.utils import str2bool


def find_low_confident_variants(args):
    ctg_name = args.ctg_name
    output_path = args.output_path
    truth_log_fn = args.truth_log_fn.split(',')
    debug = args.debug

    low_af_truth = set()
    low_tumor_count = 0
    no_tumor_alt_count = 0

    truth_vcf_fn = args.truth_vcf_fn
    vcf_reader = VcfReader(vcf_fn=truth_vcf_fn,
                           ctg_name=ctg_name,
                           show_ref=False,
                           keep_row_str=True,
                           filter_tag=args.truth_filter_tag)
    vcf_reader.read_vcf()
    truth_variant_dict = vcf_reader.variant_dict

    select_truth = set()
    #in debug mode, cal all AF with depth > 10
    debug_af_dict = defaultdict(float)
    debug_count = 0
    result_dict = defaultdict()
    homo_mismatch = set()
    het_mismatch = set()

    for idx, low_af_path in enumerate(truth_log_fn):
        fp = open(low_af_path).readlines()
        fp = [item.rstrip().split(' ') for item in fp]
        for row in fp:
            ctg, pos, normal_cov, tumor_cov, normal_alt, tumor_alt, *hap_info = row
            if ctg_name is not None and ctg != ctg_name:
                continue
            key = (ctg, int(pos)) if ctg_name is None else int(pos)
            result_dict[key] = ctg, pos, normal_cov, tumor_cov, normal_alt, tumor_alt, hap_info

    for k, v in result_dict.items():
        ctg_name, pos, normal_cov, tumor_cov, normal_alt, tumor_alt = v[:6]
        key = int(pos) if args.ctg_name is not None else (ctg_name, int(pos))
        if int(tumor_alt) == 0 or int(tumor_cov) == 0:
            low_af_truth.add(key)
            no_tumor_alt_count += 1
            continue
        elif int(tumor_alt) / float(tumor_cov) <= args.min_af or int(tumor_alt) <= args.min_alt_coverage:
            low_tumor_count += 1
            low_af_truth.add(key)
            continue
        gt = truth_variant_dict[k].genotype
        is_homo = sum(gt) == 2
        is_het = sum(gt) == 1
        af = int(tumor_alt) / max(float(tumor_cov), 1.0)
        if af < 0.2 and is_homo:
            homo_mismatch.add(key)
        if af > 0.8 and is_het:
            het_mismatch.add(key)

        if debug and float(tumor_cov) >= 20:
            if k not in truth_variant_dict:
                continue
            gt = truth_variant_dict[k].genotype
            is_homo = sum(gt) == 2
            is_het = sum(gt) == 1


            af = int(tumor_alt) / float(tumor_cov)

            if af > 0.9 and is_het:
                debug_count += 1
                select_truth.add(k)
            debug_af_dict[key] = af




    variant_dict = defaultdict()

    if output_path is not None:
        output_dir = os.path.dirname(output_path)
        if not os.path.exists(output_dir):
            subprocess.run("mkdir -p {}".format(output_dir), shell=True)
        output_file = open(output_path, 'w')

    for k, v in truth_variant_dict.items():
        if k in low_af_truth:
            continue
        if k in homo_mismatch or k in het_mismatch:
            continue

        output_file.write(v.row_str)

    output_file.close()
    ctg_info = args.ctg_name + '/' if args.ctg_name is not None else ""
    pro = max(0.0, round(len(low_af_truth) / float(len(truth_variant_dict)+1),2))
    print("[INFO] Total truth/Low AF truth/PRO: {}{}/{}/{}".format(ctg_info, len(truth_variant_dict), len(truth_variant_dict)-len(low_af_truth), pro))
    print("[INFO] Total HOM/HET mismatch: {}{}/{}".format(ctg_info, len(homo_mismatch), len(het_mismatch)))
    if debug:
        print("[INFO] Debug count:{}{}/{}/{}".format(ctg_info, debug_count, len(debug_af_dict), round(debug_count/float(len(debug_af_dict)),2)))

    if output_path is not None:
        output_file.close()


def main():
    parser = ArgumentParser(description="Calculate AF distribution in tumor and normal BAMs")

    parser.add_argument('--truth_log_fn', type=str, default=None,
                        help="Sorted tumor BAM file input, required")

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

    parser.add_argument('--samtools', type=str, default="samtools",
                        help="Absolute path to the 'samtools', samtools version >= 1.10 is required. Default: %(default)s")

    # options for advanced users
    parser.add_argument('--min_mq', type=int, default=param.min_mq,
                        help="EXPERIMENTAL: If set, reads with mapping quality with <$min_mq are filtered. Default: %(default)d")

    parser.add_argument('--min_bq', type=int, default=param.min_bq,
                        help="EXPERIMENTAL: If set, bases with base quality with <$min_bq are filtered. Default: %(default)d")

    parser.add_argument('--min_bq_cut', type=int, default=0,
                        help="EXPERIMENTAL: Minimal base quality cut-off")

    parser.add_argument('--min_af', type=float, default=0.05,
                        help="EXPERIMENTAL: Minimum VAF for a variant to be included in bechmarking")

    parser.add_argument('--min_alt_coverage', type=int, default=2,
                        help="EXPERIMENTAL: Minimum alt base count for a variant to be included in bechmarking")

    parser.add_argument('--min_coverage', type=int, default=4,
                        help="EXPERIMENTAL: Minimum coverage for a variant to be included in bechmarking")

    parser.add_argument('--debug', type=str2bool, default=0,
                        help="Add hetero candidates into training")

    global args
    args = parser.parse_args()

    find_low_confident_variants(args)


if __name__ == "__main__":
    main()

