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

import os
import subprocess
import shlex
from sys import stdin, exit
from argparse import ArgumentParser
from collections import defaultdict

from shared.utils import log_error, log_warning, file_path_from, subprocess_popen, str2bool, str_none
from shared.vcf import vcf_header

major_contigs_order = ["chr" + str(a) for a in list(range(1, 23)) + ["X", "Y"]] + [str(a) for a in
                                                                                   list(range(1, 23)) + ["X", "Y"]]



def MarkLowQual(row, quality_score_for_pass, qual):
    if row == '':
        return row

    if quality_score_for_pass and qual <= quality_score_for_pass:
        row = row.split("\t")
        row[6] = "LowQual"
        return '\t'.join(row)
    return row

def mark_rediportal(row, item):
    tag_by_rediportal = False
    if row == '' or "Germline" in row or "RefCall" in row:
        return row, tag_by_rediportal
    red_ref_base, red_alt_base, db_filter = item[:3]

    columns = row.split('\t', maxsplit=8)
    ref_base, alt_base = columns[3], columns[4]
    if red_ref_base == ref_base and red_alt_base == alt_base:
        columns[6] = "RNAEditing"
        tag_by_rediportal = True

    return '\t'.join(columns), tag_by_rediportal

def compress_index_vcf(input_vcf):
    # use bgzip to compress vcf -> vcf.gz
    # use tabix to index vcf.gz
    proc = subprocess.run('bgzip -f {}'.format(input_vcf), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    proc = subprocess.run('tabix -f -p vcf {}.gz'.format(input_vcf), shell=True, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)

def output_header(output_fn, reference_file_path, cmd_fn=None, sample_name='SAMPLE'):
    output_file = open(output_fn, "w")
    # header_str = get_header(reference_file_path=reference_file_path, cmd_fn=cmd_fn, sample_name=sample_name)
    # output_file.write(header_str)
    output_file.close()

def print_calling_step(output_fn=""):
    merge_output = os.path.join(os.path.dirname(output_fn), 'merge_output.vcf.gz')
    pileup_output = os.path.join(os.path.dirname(output_fn), 'pileup.vcf.gz')


def sort_vcf_from_stdin(args):
    """
    Sort vcf file according to variants start position and contig name.
    """

    row_count = 0
    header = []
    contig_dict = defaultdict(defaultdict)
    no_vcf_output = True
    for row in stdin:
        row_count += 1
        if row[0] == '#':
            if row not in header:
                header.append(row)
            continue
        # use the first vcf header
        columns = row.strip().split(maxsplit=3)
        ctg_name, pos = columns[0], columns[1]
        contig_dict[ctg_name][int(pos)] = row
        no_vcf_output = False
    if row_count == 0:
        print(log_warning("[WARNING] No vcf file found, please check the setting"))
    if no_vcf_output:
        print(log_warning("[WARNING] No variant found, please check the setting"))

    contigs_order = major_contigs_order + list(contig_dict.keys())
    contigs_order_list = sorted(contig_dict.keys(), key=lambda x: contigs_order.index(x))
    with open(args.output_fn, 'w') as output:
        output.write(''.join(header))
        for contig in contigs_order_list:
            all_pos = sorted(contig_dict[contig].keys())
            for pos in all_pos:
                output.write(contig_dict[contig][pos])


def sort_vcf_from(args):
    """
    Sort vcf file from providing vcf filename prefix.
    """
    output_fn = args.output_fn
    input_dir = args.input_dir
    vcf_fn_prefix = args.vcf_fn_prefix
    vcf_fn_suffix = args.vcf_fn_suffix
    sample_name = args.sample_name
    ref_fn = args.ref_fn
    contigs_fn = args.contigs_fn
    compress_vcf = args.compress_vcf
    QUAL = args.qual
    output_no_tagging_fn = args.output_no_tagging_fn
    show_ref = args.show_ref

    print("[INFO] Sorting VCFs...")

    if not os.path.exists(input_dir):
        exit(log_error("[ERROR] Input directory: {} not exists!").format(input_dir))
    all_files = os.listdir(input_dir)

    if vcf_fn_prefix is not None:
        all_files = [item for item in all_files if item.startswith(vcf_fn_prefix)]
        if len(all_files) == 0:
            output_header(output_fn=output_fn, reference_file_path=ref_fn, cmd_fn=args.cmd_fn, sample_name=sample_name)
            print (log_warning(
                "[WARNING] No vcf file found with prefix:{}/{}, output empty vcf file".format(input_dir,vcf_fn_prefix)))
            compress_index_vcf(output_fn)
            print_calling_step(output_fn=output_fn)
            return

    if vcf_fn_suffix is not None:
        all_files = [item for item in all_files if item.endswith(vcf_fn_suffix)]
        if len(all_files) == 0:
            output_header(output_fn=output_fn, reference_file_path=ref_fn, cmd_fn=args.cmd_fn, sample_name=sample_name)
            print (log_warning(
                "[WARNING] No vcf file found with suffix:{}/{}, output empty vcf file".format(input_dir,vcf_fn_prefix)))
            compress_index_vcf(output_fn)
            print_calling_step(output_fn=output_fn)
            return

    all_contigs_list = []
    if contigs_fn and os.path.exists(contigs_fn):
        with open(contigs_fn) as f:
            all_contigs_list = [item.rstrip() for item in f.readlines()]
    else:
        exit(log_error("[ERROR] Cannot find contig file {}. Exit!").format(contigs_fn))

    contigs_order = major_contigs_order + all_contigs_list
    contigs_order_list = sorted(all_contigs_list, key=lambda x: contigs_order.index(x))

    rediportal_variant_dict = defaultdict()
    if args.tag_variant_using_readiportal:
        print("[INFO] Reading readiportal source file ...")
        readiportal_database_filter_tag = set(args.readiportal_database_filter_tag.split(':')) if args.readiportal_database_filter_tag is not None else None

        #allow gzip format
        if args.readiportal_source_fn is None or not os.path.exists(args.readiportal_source_fn):
            print('[WARNING] Enabled tagging variant using readiportal, but --readiportal_source_fn {} file not found, skip tagging!'.format(args.readiportal_source_fn))
        else:
            unzip_process = subprocess_popen(shlex.split("gzip -fdc %s" % (args.readiportal_source_fn)))
            for row_idx, row in enumerate(unzip_process.stdout):
                if row_idx == 0:
                    continue

                columns = row.rstrip().split('\t', maxsplit=6)
                if len(contigs_order_list) and columns[0] not in contigs_order_list:
                    continue
                try:
                    key = (columns[0], int(columns[1]))
                except:
                    print(columns[1])
                    continue
                db_filter = columns[5]

                if readiportal_database_filter_tag is not None and db_filter not in readiportal_database_filter_tag:
                    continue

                ref_base, alt_base = columns[2:4]
                rediportal_variant_dict[key] = (ref_base, alt_base, db_filter)

            unzip_process.stdout.close()
            unzip_process.wait()
    row_count = 0
    header = []
    no_vcf_output = True
    need_write_header = True

    tag_by_rediportal_num = 0

    output = open(output_fn, 'w')
    if args.tag_variant_using_readiportal:
        output_no_tagging = open(output_no_tagging_fn, 'w')

    for contig in contigs_order_list:
        contig_dict = defaultdict(str)
        contig_vcf_fns = [fn for fn in all_files if contig in fn]
        for vcf_fn in contig_vcf_fns:
            file = os.path.join(input_dir, vcf_fn)

            fn = open(file, 'r')
            for row in fn:
                row_count += 1
                if row[0] == '#':
                    if row not in header:
                        header.append(row)
                    continue
                # use the first vcf header
                columns = row.strip().split(maxsplit=6)
                ctg_name, pos = columns[0], columns[1]
                # skip vcf file sharing same contig prefix, ie, chr1 and chr11
                if ctg_name != contig:
                    break

                pos = int(columns[1])
                qual = float(columns[5])
                ref_base, alt_base = columns[3], columns[4]
                is_reference = (alt_base == "." or ref_base == alt_base)

                if not show_ref and is_reference:
                    continue
                if not is_reference:
                    row = MarkLowQual(row, QUAL, qual)
                key = (ctg_name, int(pos))
                if key in rediportal_variant_dict:
                    row, tag_by_rediportal = mark_rediportal(row, rediportal_variant_dict[key])
                    tag_by_rediportal_num += int(tag_by_rediportal)
                contig_dict[int(pos)] = row
                no_vcf_output = False
            fn.close()

        if need_write_header and len(header):
            output.write(''.join(header))
            need_write_header = False
            if args.tag_variant_using_readiportal:
                output_no_tagging.write(''.join(header))
        all_pos = sorted(contig_dict.keys())
        for pos in all_pos:
            output.write(contig_dict[pos])
            if args.tag_variant_using_readiportal:
                output_no_tagging.write(contig_dict[pos].replace('RNAEditing', 'PASS'))

    output.close()
    if args.tag_variant_using_readiportal:
        output_no_tagging.close()

    if row_count == 0:
        print(log_warning("[WARNING] No vcf file found, output empty vcf file"))
        output_header(output_fn=output_fn, reference_file_path=ref_fn, sample_name=sample_name)
        if compress_vcf:
            compress_index_vcf(output_fn)
        print_calling_step(output_fn=output_fn)
        return
    if no_vcf_output:
        output_header(output_fn=output_fn, reference_file_path=ref_fn, sample_name=sample_name)
        print(log_warning("[WARNING] No variant found, output empty vcf file"))
        if compress_vcf:
            compress_index_vcf(output_fn)
        return

    if compress_vcf:
        compress_index_vcf(output_fn)
        if args.tag_variant_using_readiportal:
            compress_index_vcf(output_no_tagging_fn)

    if args.tag_variant_using_readiportal:
        print('[INFO] Dataset size:{}, total variants tagged by REDIportal dataset: {}'.format(len(rediportal_variant_dict), tag_by_rediportal_num))

    print("[INFO] Finished VCF sorting!")


def main():
    parser = ArgumentParser(description="Sort a VCF file according to contig name and starting position")

    parser.add_argument('--output_fn', type=str, default=None, required=True,
                        help="Output VCF filename, required")

    parser.add_argument('--input_dir', type=str, default=None,
                        help="Input directory")

    parser.add_argument('--vcf_fn_prefix', type=str, default=None,
                        help="Input vcf filename prefix")

    parser.add_argument('--vcf_fn_suffix', type=str, default='.vcf',
                        help="Input vcf filename suffix")

    parser.add_argument('--ref_fn', type=str, default=None,
                        help="Reference fasta file input")

    parser.add_argument('--sample_name', type=str, default="SAMPLE",
                        help="Define the sample name to be shown in the VCF file, optional")

    parser.add_argument('--contigs_fn', type=str, default=None,
                        help="Contigs file with all processing contigs")

    parser.add_argument('--compress_vcf', type=str2bool, default=False,
                        help="Compress and index the output VCF")

    parser.add_argument('--show_ref', type=str2bool, default=False,
                        help="Compress and index the output VCF")

    parser.add_argument('--cmd_fn', type=str_none, default=None,
                        help="If defined, added command line into VCF header")

    parser.add_argument('--qual', type=int, default=2,
                        help="If set, variants with >$qual will be marked 'PASS', or 'LowQual' otherwise, optional")

    parser.add_argument('--output_no_tagging_fn', type=str, default=None,
                        help="Output VCF filename without tagging")

    parser.add_argument('--tag_variant_using_readiportal', type=str2bool, default=None,
                        help="If defined, added command line into VCF header")

    parser.add_argument('--readiportal_source_fn', type=str_none, default=None,
                        help="If defined, added command line into VCF header")

    parser.add_argument('--readiportal_database_filter_tag', type=str, default=None,
                        help="If defined, added command line into VCF header")

    args = parser.parse_args()
    if args.input_dir is None:
            sort_vcf_from_stdin(args)
    else:
        # default entry
        sort_vcf_from(args)


if __name__ == "__main__":
    main()
