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

import shlex
import os
import sys
import math
import argparse
import random
import shlex
import math
import subprocess

from os.path import isfile, abspath
from sys import exit, stderr
from subprocess import check_output, PIPE, Popen
from subprocess import PIPE
# from shared.vcf import VcfWriter
# from shared.bed import BedWriter
# from collections import defaultdict
from argparse import ArgumentParser, SUPPRESS
# from collections import Counter, defaultdict
#
# from shared.interval_tree import bed_tree_from, is_region_in
# import shared.param as param
# from shared.utils import AltInfos, str2bool


major_contigs = {"chr" + str(a) for a in list(range(1, 23))}.union(
    {str(a) for a in list(range(1, 23))})
major_contigs_order = ["chr" + str(a) for a in list(range(1, 23))] + [str(a) for a in
                                                                                   list(range(1, 23))]

def subprocess_popen(args, stdin=None, stdout=PIPE, stderr=stderr, bufsize=8388608):
    return Popen(args, stdin=stdin, stdout=stdout, stderr=stderr, bufsize=bufsize, universal_newlines=True)


def get_candidates(args):

    ctg_name = args.ctg_name
    # bed_fn = args.bed_fn
    bam_fn = args.bam_fn
    # ref_fn = args.ref_fn
    high_confident_bed_fn = args.high_confident_bed_fn
    mosdepth = args.mosdepth
    output_dir = args.output_dir
    # input_vcf_fn = args.input_vcf_fn
    excl_flags = args.excl_flags
    min_mq = args.min_mq
    min_coverage = args.min_coverage
    threads = args.threads
    bedtools = args.bedtools

    truth_vcf_fn = args.truth_vcf_fn
    depth_output = os.path.join(output_dir, "coverage",)
    if not os.path.exists(output_dir):
        output = subprocess.run("mkdir -p {}".format(output_dir), shell=True)
    # cal mosdepth for the bam
    mos_depth_command = "{}".format(mosdepth)
    mos_depth_command += ' --flag ' + str(excl_flags) if excl_flags is not None else ""
    mos_depth_command += ' --mapq ' + str(min_mq) + ' ' if min_mq is not None else ""
    mos_depth_command += ' --threads ' + str(threads) + ' ' if threads is not None else ""
    mos_depth_command += ' --chrom ' + str(ctg_name) + ' ' if ctg_name is not None else ""
    mos_depth_command += ' ' + depth_output
    mos_depth_command += ' ' + bam_fn

    print('[INFO] Run mosdepth depth calculation command: {}'.format(mos_depth_command))
    subprocess.run(mos_depth_command, shell=True)

    depth_output_fn = depth_output + '.per-base.bed.gz'
    #mosdepth output
    unzip_command = "gzip -fdc {} | awk '$4 >= {}'".format(depth_output_fn, args.min_coverage)

    # depth_file_path_process = subprocess_popen(shlex.split("gzip -fdc %s" % (depth_output_fn)))
    # depth_file_path_output = depth_file_path_process.stdout
    #
    # for row in depth_file_path_output:
    #     columns = row.rstrip().split('\t')
    #     if len(columns) < 4:
    #         continue
    #     if ctg_name is not None and columns[0] != ctg_name:
    #         continue
    #     ctg = columns[0]
    #     #exclude HLA|decoy|random|alt|chrUn|chrEBV
    #
    #     start = int(columns[1])
    #     end = int(columns[2])
    #     coverage = int(columns[3])
    #     if coverage < args.min_coverage:
    #         continue
    # depth_file_path_output.close()
    # depth_file_path_process.wait()

    # stdin = depth_file_path_process.stdout

    bed_output_path = os.path.join(output_dir, 'coverage.bed')
    bedtools_merge_command = "{} | {} merge -d 1 -c 4 -o mean -i - > {}".format(unzip_command, bedtools, bed_output_path)
    print('[INFO] Run Extract regions exceeding minimum coverage command: {}'.format(bedtools_merge_command))
    subprocess.run(bedtools_merge_command, shell=True)

    # bedtools_merge_process.close()

    #intersect with high-confident BED
    output_bed_path = os.path.join(output_dir, 'final.bed')
    bedtools_intersect_command = "{} intersect".format(bedtools)
    bedtools_intersect_command += " -a " + bed_output_path
    bedtools_intersect_command += " -b " + high_confident_bed_fn
    bedtools_intersect_command += " > " + output_bed_path
    print('[INFO] Run bedtools BED intersection command: {}'.format(bedtools_intersect_command))
    subprocess.run(bedtools_intersect_command, shell=True)


    file_directory = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    main_entry = os.path.join(file_directory, "clair3_rna.py")

    chunk_num = 15
    major_contigs_order = ['chr20']
    chunk_list = []
    chunk_list_path = os.path.join(args.output_dir, 'CHUNK_LIST')

    with open(chunk_list_path, 'w') as output_file:
        for contig_name in major_contigs_order:
            chunk_num = args.chunk_num
            for chunk_id in range(1, chunk_num + 1):
                output_file.write(contig_name + ' ' + str(chunk_id) + ' ' + str(chunk_num) + '\n')
                chunk_list.append((contig_name, chunk_id, chunk_num))

    parallel_command = "{} -C ' ' -j{} {} {} cal_truth_vcf_af_distribution".format(args.parallel, threads, args.pypy, main_entry)
    parallel_command += " --output_path {}".format(os.path.join(output_dir, 'vcf_output')) + "/truths_{1}_{2} "
    parallel_command += " --bam_fn " + str(args.bam_fn)
    parallel_command += " --bed_fn " + str(output_bed_path)
    parallel_command += " --truth_vcf_fn " + str(truth_vcf_fn)
    parallel_command += " --threads " + str(args.threads)
    parallel_command += " --min_mq " + str(args.min_mq)
    parallel_command += " --chunk_id {2}"
    parallel_command += " --chunk_num " + str(chunk_num)
    parallel_command += " --samtools " + str(args.samtools)
    parallel_command += " --ctg_name {1}"
    parallel_command += " :::: " + str(chunk_list_path)
    parallel_command += ' && ' + args.pypy + ' ' + main_entry + ' concat_files'
    parallel_command += ' --input_dir ' + os.path.join(output_dir, 'vcf_output')
    parallel_command += ' --input_prefix ' + "truths_"
    parallel_command += ' --output_fn truths '

    print(parallel_command)

    subprocess.run(parallel_command, shell=True)

    # # bed_tree = bed_tree_from(bed_file_path=output_bed_path, contig_name=ctg_name)
    #
    # if truth_vcf_fn is not None:
    #     if not os.path.exists(truth_vcf_fn):
    #         sys.exit("[ERROR] Cannot find truth vcf fn:{}, exit!".format(truth_vcf_fn))
    #         return




def main():
    parser = ArgumentParser(description="Get BED with sufficient RNA read support")

    parser.add_argument('--output_dir', type=str, default=None,
                        help="Output directory")

    parser.add_argument('--platform', type=str, default='ont',
                        help="Select the sequencing platform of the input. Default: %(default)s")

    parser.add_argument('--bam_fn', type=str, default="input.bam",
                        help="Sorted BAM file input, required")

    parser.add_argument('--ref_fn', type=str, default=None,
                        help="Reference fasta file input, required")

    parser.add_argument('--truth_vcf_fn', type=str, default=None,
                        help="Reference fasta file input, required")

    parser.add_argument('--ctg_name', type=str, default=None,
                        help="The name of sequence to be processed")

    parser.add_argument('--high_confident_bed_fn', type=str, default=None,
                        help="Call variant only in the provided regions")

    parser.add_argument('--pypy', type=str, default="pypy3",
                        help="Path to the 'parallel', parallel version >= 1.10 is required. default: %(default)s")

    parser.add_argument('--parallel', type=str, default="parallel",
                        help="Path to the 'parallel', parallel version >= 1.10 is required. default: %(default)s")

    parser.add_argument('--samtools', type=str, default="samtools",
                        help="Path to the 'samtools', samtools version >= 1.10 is required. default: %(default)s")

    parser.add_argument('--mosdepth', type=str, default="mosdepth",
                        help="Path to the 'mosdepth', samtools version >= 1.10 is required. default: %(default)s")

    parser.add_argument('--bedtools', type=str, default="bedtools",
                        help="Path to the 'mosdepth', samtools version >= 1.10 is required. default: %(default)s")

    parser.add_argument('--excl_flags', type=int, default=2316,
                        help="Path to the 'mosdepth', samtools version >= 1.10 is required. default: %(default)s")

    parser.add_argument('--min_mq', type=int, default=5,
                        help="EXPERIMENTAL: If set, reads with mapping quality with <$min_mq are filtered. Default: %(default)d")

    parser.add_argument('--min_bq', type=int, default=0,
                        help="EXPERIMENTAL: If set, bases with base quality with <$min_bq are filtered. Default: %(default)d")

    parser.add_argument('--threads', type=int, default=None,
                        help="Max #threads to be used.")

    parser.add_argument('--chunk_num', type=int, default=15,
                        help="Max #threads to be used.")

    parser.add_argument('--min_alt_coverage', type=int, default=2,
                        help="EXPERIMENTAL: Minimum alt base count for a variant to be included in bechmarking")

    parser.add_argument('--min_coverage', type=int, default=4,
                        help="EXPERIMENTAL: Minimum coverage for a variant to be included in bechmarking")




    args = parser.parse_args()

    get_candidates(args)


if __name__ == "__main__":
    main()
