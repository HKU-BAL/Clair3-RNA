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
import sys
import shlex
import logging
import argparse

from sys import stderr
from subprocess import Popen
from argparse import ArgumentParser
from subprocess import PIPE
from collections import defaultdict
logging.basicConfig(format='%(message)s', level=logging.INFO)

from shared.utils import str2bool, str_none
from shared.vcf import VcfReader, VcfWriter
from shared.interval_tree import bed_tree_from, is_region_in

def subprocess_popen(args, stdin=None, stdout=PIPE, stderr=stderr, bufsize=8388608):
    return Popen(args, stdin=stdin, stdout=stdout, stderr=stderr, bufsize=bufsize, universal_newlines=True)

def covert_gt(gt):
    if gt == None:
        return None
    if gt == '.':
        return -1
    gt = gt.replace('|', '/').split('/')
    gt = sum([int(item) for item in gt])
    return gt

def str2bool(v):
    if v is None:
        return v
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'ture', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'flase', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def metrics(query_fp, query_tp, truth_fn, truth_tp):
    # https://github.com/Illumina/hap.py/blob/master/doc/happy.md
    precision = query_tp / (query_tp + query_fp)
    recall = truth_tp / (truth_tp + truth_fn)
    f1_score = 2 * precision * recall / (precision + recall)
    return round(precision, 6), round(recall, 6), round(f1_score, 6)


def Cal(args):
    happy_vcf_fn = args.happy_vcf_fn
    ctg_name = args.ctg_name
    output_fn = args.output_fn
    skip_genotyping = args.skip_genotyping
    truths_info_fn = args.truths_info_fn
    input_vcf_fn = args.input_vcf_fn
    bed_fn = args.bed_fn


    min_coverage = args.min_coverage
    min_alt_coverage = args.min_alt_coverage
    min_af = args.min_af

    truths_dict = defaultdict()
    input_variant_dict = defaultdict()

    if args.bed_fn is not None:
        bed_tree = bed_tree_from(bed_file_path=bed_fn, contig_name=ctg_name)

    if truths_info_fn is not None:
        if not os.path.exists(truths_info_fn):
            sys.exit('[ERROR] Truths info file not found:{}, exit!'.format(truths_info_fn))
        with open(truths_info_fn) as f:
            for row in f:
                ctg, pos, cov, alt_cov = row.split()[:4]
                cov = float(cov)
                alt_cov = float(alt_cov)
                af = alt_cov / cov if cov > 0 else 0.0
                key = (ctg, int(pos)) if ctg is None else int(pos)
                truths_dict[key] = (cov, alt_cov, af)

    low_af_set = set()
    if input_vcf_fn is not None:
        input_vcf_reader = VcfReader(vcf_fn=input_vcf_fn,
                                     ctg_name=ctg_name,
                                     # ctg_start=args.ctg_start,
                                     # ctg_end=args.ctg_end,
                                     show_ref=False,
                                     keep_row_str=True,
                                     skip_genotype=skip_genotyping,
                                     filter_tag=args.input_filter_tag,
                                     keep_af=True,
                                     # min_qual=args.min_qual,
                                     # max_qual=args.max_qual,
                                     # discard_multi=args.discard_multi,
                                     discard_indel=False)
        input_vcf_reader.read_vcf()
        input_variant_dict = input_vcf_reader.variant_dict
        for k in list(input_variant_dict.keys()):
            v = input_variant_dict[k]
            DP = v.DP
            AD = v.AD.split(',')
            if min_coverage is not None:
                if DP < min_coverage:
                    del input_variant_dict[k]
                if len(AD) == 2 and min_alt_coverage is not None:
                    if AD[1] < min_alt_coverage:
                        del input_variant_dict[k]
                        low_af_set.add(k)

    if output_fn:
        output_file = open(output_fn, 'w')
    else:
        output_file = None
    happy_vcf_unzip_process = subprocess_popen(shlex.split("gzip -fdc %s" % (happy_vcf_fn)))

    truth_all_tp, query_all_tp, query_all_fp, truth_all_fn = 0, 0, 0, 0
    truth_snp_tp, query_snp_tp, query_snp_fp, truth_snp_fn = 0, 0, 0, 0
    truth_indel_tp, query_indel_tp, query_indel_fp, truth_indel_fn = 0, 0, 0, 0
    truth_ins_tp, query_ins_tp, query_ins_fp, truth_ins_fn = 0, 0, 0, 0
    truth_del_tp, query_del_tp, query_del_fp, truth_del_fn = 0, 0, 0, 0

    snp_query_fp_set = set()
    snp_query_tp_set = set()
    snp_truth_fn_set = set()
    snp_truth_tp_set = set()

    indel_query_fp_set = set()
    indel_query_tp_set = set()
    indel_truth_fn_set = set()
    indel_truth_tp_set = set()

    happy_variant_dict = defaultdict()
    for row in happy_vcf_unzip_process.stdout:
        if row[0] == '#':
            continue
        columns = row.strip().split()

        ctg, pos = columns[0], int(columns[1])
        key = (ctg, int(pos)) if ctg_name is None else int(pos)

        if ctg_name is not None and ctg != ctg_name:
            continue

        if key in low_af_set:
            continue
        if key in truths_dict:
            dp, ad, af = truths_dict[key]
            if min_coverage is not None and dp < min_coverage:
                continue
            if min_alt_coverage is not None and ad < min_alt_coverage:
                continue

        # happy_variant_dict[key] = 1

    # for k in input_variant_dict:
    #     ctg = ctg_name if ctg_name is not None else k[0]
    #     pos = k if ctg_name is not None else k[1]
    #     if not is_region_in(bed_tree, ctg_name, pos):
    #         continue
    #     if k not in happy_variant_dict:
    #         print(k)



        FORMAT, TRUTH, QUERY = columns[8], columns[9], columns[10]
        FORMAT = FORMAT.split(':')
        TRUTH = TRUTH.split(':')
        QUERY = QUERY.split(':')

        ft_dict = dict(zip(FORMAT, TRUTH))
        fq_dict = dict(zip(FORMAT, QUERY))

        # hap.py vcf header
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        ##FORMAT=<ID=BD,Number=1,Type=String,Description="Decision for call (TP/FP/FN/N)">
        ##FORMAT=<ID=BK,Number=1,Type=String,Description="Sub-type for decision (match/mismatch type). (Loose match distance is 30)">
        ##FORMAT=<ID=BI,Number=1,Type=String,Description="Additional comparison information">
        ##FORMAT=<ID=QQ,Number=1,Type=Float,Description="Variant quality for ROC creation">
        ##INFO=<ID=Regions,Number=.,Type=String,Description="Tags for regions.">
        ##FORMAT=<ID=BVT,Number=1,Type=String,Description="High-level variant type (SNP|INDEL).">
        ##FORMAT=<ID=BLT,Number=1,Type=String,Description="High-level location type (het|homref|hetalt|homalt|nocall).">

        t_BD = ft_dict['BD'] if 'BD' in ft_dict else None
        t_BI = ft_dict['BI'] if 'BI' in ft_dict else None
        t_BVT = ft_dict['BVT'] if 'BVT' in ft_dict else None
        t_gt = ft_dict['GT'] if 'GT' in ft_dict else None
        t_bk = ft_dict['BK'] if 'BK' in ft_dict else None

        q_BD = fq_dict['BD'] if 'BD' in fq_dict else None
        q_BI = fq_dict['BI'] if 'BI' in fq_dict else None
        q_BVT = fq_dict['BVT'] if 'BVT' in fq_dict else None
        q_gt = fq_dict['GT'] if 'GT' in ft_dict else None
        q_bk = fq_dict['BK'] if 'BK' in ft_dict else None

        if not t_BD or not t_BI or not t_BVT or not q_BD or not q_BI or not q_BVT \
                or not t_gt or not q_gt or not t_bk or not q_bk:
            sys.exit("[ERROR] Happy format not match, exit!")

        query_fp = q_BD == 'FP'
        query_tp = q_BD == 'TP'
        truth_fn = t_BD == 'FN'
        truth_tp = t_BD == 'TP'

        if skip_genotyping:
            #alternative match, only genotye mismatch
            if q_bk == 'am' and t_bk == 'am':
                query_fp = False
                query_tp = True
                truth_fn = False
                truth_tp = True

        is_query_snp_fp = (q_BVT == 'SNP') and query_fp
        is_query_snp_tp = (q_BVT == 'SNP') and query_tp
        is_truth_snp_fn = (t_BVT == 'SNP') and truth_fn
        is_truth_snp_tp = (t_BVT == 'SNP') and truth_tp

        is_query_indel_fp = (q_BVT == 'INDEL') and query_fp
        is_query_indel_tp = (q_BVT == 'INDEL') and query_tp
        is_truth_indel_fn = (t_BVT == 'INDEL') and truth_fn
        is_truth_indel_tp = (t_BVT == 'INDEL') and truth_tp

        query_snp_fp = query_snp_fp + 1 if is_query_snp_fp else query_snp_fp
        query_snp_tp = query_snp_tp + 1 if is_query_snp_tp else query_snp_tp
        truth_snp_fn = truth_snp_fn + 1 if is_truth_snp_fn else truth_snp_fn
        truth_snp_tp = truth_snp_tp + 1 if is_truth_snp_tp else truth_snp_tp


        query_indel_fp = query_indel_fp + 1 if is_query_indel_fp else query_indel_fp
        query_indel_tp = query_indel_tp + 1 if is_query_indel_tp else query_indel_tp
        truth_indel_fn = truth_indel_fn + 1 if is_truth_indel_fn else truth_indel_fn
        truth_indel_tp = truth_indel_tp + 1 if is_truth_indel_tp else truth_indel_tp

        if is_query_snp_fp:
            snp_query_fp_set.add(key)
        if is_query_snp_tp:
            snp_query_tp_set.add(key)
        if is_truth_snp_fn:
            snp_truth_fn_set.add(key)
        if is_truth_snp_tp:
            snp_truth_tp_set.add(key)

        if is_query_indel_fp:
            indel_query_fp_set.add(key)
        if is_query_indel_tp:
            indel_query_tp_set.add(key)
        if is_truth_indel_fn:
            indel_truth_fn_set.add(key)
        if is_truth_indel_tp:
            indel_truth_tp_set.add(key)


        is_query_ins_fp = q_BI[0] == 'i' and is_query_indel_fp
        is_query_ins_tp = q_BI[0] == 'i' and is_query_indel_tp
        is_truth_ins_fn = t_BI[0] == 'i' and is_truth_indel_fn
        is_truth_ins_tp = t_BI[0] == 'i' and is_truth_indel_tp

        is_query_del_fp = q_BI[0] == 'd' and is_query_indel_fp
        is_query_del_tp = q_BI[0] == 'd' and is_query_indel_tp
        is_truth_del_fn = t_BI[0] == 'd' and is_truth_indel_fn
        is_truth_del_tp = t_BI[0] == 'd' and is_truth_indel_tp

        query_ins_fp = query_ins_fp + 1 if is_query_ins_fp else query_ins_fp
        query_ins_tp = query_ins_tp + 1 if is_query_ins_tp else query_ins_tp
        truth_ins_fn = truth_ins_fn + 1 if is_truth_ins_fn else truth_ins_fn
        truth_ins_tp = truth_ins_tp + 1 if is_truth_ins_tp else truth_ins_tp

        query_del_fp = query_del_fp + 1 if is_query_del_fp else query_del_fp
        query_del_tp = query_del_tp + 1 if is_query_del_tp else query_del_tp
        truth_del_fn = truth_del_fn + 1 if is_truth_del_fn else truth_del_fn
        truth_del_tp = truth_del_tp + 1 if is_truth_del_tp else truth_del_tp

    # truth_all_tp = truth_snp_tp + truth_indel_tp
    # truth_all_fn = truth_snp_fn + truth_indel_fn
    # query_all_fp = query_snp_fp + query_indel_fp
    # query_all_tp = query_snp_tp + query_indel_tp
    #
    # # p->precision, r->recall, f1->f1_score
    # # a->overall, s->snp, id->indel, i->insertion, d->deletion
    # ap, ar, af1 = metrics(query_fp=query_all_fp, query_tp=query_all_tp, truth_fn=truth_all_fn, truth_tp=truth_all_tp)
    # sp, sr, sf1 = metrics(query_fp=query_snp_fp, query_tp=query_snp_tp, truth_fn=truth_snp_fn, truth_tp=truth_snp_tp)
    # idp, idr, idf1 = metrics(query_fp=query_indel_fp, query_tp=query_indel_tp, truth_fn=truth_indel_fn, truth_tp=truth_indel_tp)
    # ip, ir, if1 = metrics(query_fp=query_ins_fp, query_tp=query_ins_tp, truth_fn=truth_ins_fn, truth_tp=truth_ins_tp)
    # dp, dr, df1 = metrics(query_fp=query_del_fp, query_tp=query_del_tp, truth_fn=truth_del_fn, truth_tp=truth_del_tp)
    #
    # print (''.join([item.ljust(20) for item in ["VariantType", 'TRUTH.FP', 'TRUTH.FN', 'TRUTH.TP','QUERY.TP', 'METRIC.Precision', 'METRIC.Recall', 'METRIC.F1_Score']]), file=output_file)
    # print (''.join([str(item).ljust(20) for item in ["Overall", query_all_fp, truth_all_fn, truth_all_tp, query_all_tp, ap, ar, af1]]), file=output_file)
    # print (''.join([str(item).ljust(20) for item in ["SNP", query_snp_fp, truth_snp_fn, truth_snp_tp, query_snp_tp, sp, sr, sf1]]),file=output_file)
    # print (''.join([str(item).ljust(20) for item in ["INDEL", query_indel_fp, truth_indel_fn, truth_indel_tp, query_indel_tp, idp, idr, idf1]]), file=output_file)
    # print (''.join([str(item).ljust(20) for item in ["INS", query_ins_fp, truth_ins_fn, truth_ins_tp, query_ins_tp, ip, ir, if1]]), file=output_file)
    # print (''.join([str(item).ljust(20) for item in ["DEL", query_del_fp, truth_del_fn, truth_del_tp, query_del_tp, dp, dr, df1]]), file=output_file)
    # print('\n', file=output_file)



    #
    # # print log_happy output
    # pass_row = []
    # snp_row = []
    # indel_row = []
    # if args.log_happy and os.path.exists(args.log_happy):
    #     log_happy = open(args.log_happy)
    #     for row in log_happy.readlines():
    #         if 'PASS' not in row:
    #             continue
    #         pass_row.append(row)
    #
    #     for row in pass_row:
    #
    #         if 'INDEL' in row:
    #             row = row.split()
    #             tp, fn, fp = row[3], row[4], row[6]
    #             precision, recall, f1 = row[11], row[10], row[13]
    #             indel_row = [fp, fn, tp, precision, recall, f1]
    #         if 'SNP' in row:
    #             row = row.split()
    #             tp, fn, fp = row[3], row[4], row[6]
    #             precision, recall, f1 = row[11], row[10], row[13]
    #             snp_row = [fp, fn, tp, precision, recall, f1]
    #     print('Double check with happy log:', file=output_file)
    #     print(' '.join(['%.6f' % item for item in [sp, sr, sf1, idp, idr, idf1]] + [str(item) for item in
    #                                                                                 [query_snp_fp, truth_snp_fn,
    #                                                                                  truth_snp_tp, query_indel_fp,
    #                                                                                  truth_indel_fn, truth_indel_tp]]),
    #           file=output_file)
    #     print(' '.join(snp_row[3:]) + ' ' + ' '.join(indel_row[3:]) + ' ' + ' '.join(snp_row[:3]) + ' ' + ' '.join(
    #         indel_row[:3]), file=output_file)
    #     print('\n', file=output_file)
    #
    # print(' '.join([str(item) for item in [ap, ar, af1, sp, sr, sf1, idp, idr, idf1, ip, ir, if1, dp, dr, df1]]),
    #       file=output_file)
    # print(' '.join([str(item) for item in
    #                 [query_all_fp, truth_all_fn, truth_all_tp, query_snp_fp, truth_snp_fn, truth_snp_tp, query_indel_fp,
    #                  truth_indel_fn, truth_indel_tp, query_ins_tp, truth_ins_tp, truth_ins_tp, query_del_fp,
    #                  truth_del_fn, truth_del_tp]]), file=output_file)

    if output_fn:
        output_file.close()


def main():
    parser = ArgumentParser(description="Overall Metrics of hap.py output")

    parser.add_argument('--happy_vcf_fn', type=str, default=None,
                        help="Path to the happy vcf output file")

    parser.add_argument('--log_happy', type=str, default=None,
                        help="Path to the happy vcf output file")

    parser.add_argument('--bed_fn', type=str, default=None,
                        help="Path to the happy vcf output file")

    parser.add_argument('--ctg_name', type=str, default=None,
                        help="The name of sequence to be processed")

    parser.add_argument('--output_fn', type=str, default=None,
                        help="Filename of the metrics output")

    parser.add_argument('--min_af', type=float, default=None,
                        help="EXPERIMENTAL: Minimum VAF for a variant to be included in bechmarking")

    parser.add_argument('--min_alt_coverage', type=int, default=2,
                        help="EXPERIMENTAL: Minimum alt base count for a variant to be included in bechmarking")

    parser.add_argument('--min_coverage', type=int, default=4,
                        help="EXPERIMENTAL: Minimum coverage for a variant to be included in bechmarking")

    parser.add_argument('--truths_info_fn', type=str, default=None,
                        help="")

    parser.add_argument('--input_vcf_fn', type=str, default=None,
                        help="")

    ## Only benchmark 'HighConf' tag in seqc VCF
    parser.add_argument('--skip_genotyping', type=str2bool, default=True,
                        help="Skip calculating VCF genotype")

    parser.add_argument('--input_filter_tag', type=str_none, default=None,
                        help="Filter variants with tag from the input VCF")

    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    Cal(args)


if __name__ == "__main__":
    main()
