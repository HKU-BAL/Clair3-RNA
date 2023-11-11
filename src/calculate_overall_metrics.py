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
from argparse import ArgumentParser, SUPPRESS
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


def output_best_cut_off(fp_qual_dict, tp_qual_dict, fn_count, truth_tp_count, use_int_cut_off=True,add_tp_fn=False):
    results = []
    if use_int_cut_off:
        qual_list = set([int(q) for q in list(fp_qual_dict.values()) + list(tp_qual_dict.values())])
    else:
        qual_list = [item / 100.0 for item in range(0, 101)]

    for qual in qual_list:
        fp = sum([1 for k, v in fp_qual_dict.items() if v >= qual])
        tp = sum([1 for k, v in tp_qual_dict.items() if v >= qual])
        fn = fn_count + len(tp_qual_dict) - tp
        pass_qual_tp_count= truth_tp_count - (len(tp_qual_dict) - tp)
        # snv_pre, snv_rec, snv_f1 = metrics(tp=tp_snv, fp=fp_snv, fn=fn_snv)
        pre, rec, f1 = metrics(query_fp=fp, query_tp=tp, truth_fn=fn, truth_tp=pass_qual_tp_count)
        tp_fn = tp + fn
        results.append([qual, pre, rec, f1, tp, fp, fn, tp_fn])

    results = sorted(results, key=lambda x: x[3], reverse=True)
    return results


def Cal(args):
    happy_vcf_fn = args.happy_vcf_fn
    ctg_name = args.ctg_name
    output_fn = args.output_fn
    skip_genotyping = args.skip_genotyping
    truths_info_fn = args.truths_info_fn
    input_vcf_fn = args.input_vcf_fn
    bed_fn = args.bed_fn

    min_qual = args.min_qual
    min_coverage = args.min_coverage
    min_alt_coverage = args.min_alt_coverage
    min_af = args.min_af
    debug = args.debug

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

                alt_cov_list = [float(item) for item in alt_cov.split(',')]

                af = [float(item / cov) if cov > 0 else 0.0 for item in alt_cov_list]
                key = (ctg, int(pos)) if ctg is None else int(pos)
                truths_dict[key] = (cov, alt_cov_list, af)
                if args.bed_fn is not None and not is_region_in(bed_tree, ctg, pos):
                    print('Not in BED')

    low_confident_set = set()
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
                                     parser_dp=True,
                                     parser_ad=True,
                                     min_qual=min_qual,
                                     # max_qual=args.max_qual,
                                     # discard_multi=args.discard_multi,
                                     discard_indel=False)
        input_vcf_reader.read_vcf()
        input_variant_dict = input_vcf_reader.variant_dict
        for k in list(input_variant_dict.keys()):
            v = input_variant_dict[k]
            DP = v.DP
            AD = None
            if v.AD is not None:
                AD = v.AD.split(',')[1:] # ref,alt1,alt2
            if min_coverage is not None:
                if DP < min_coverage:
                    low_confident_set.add(k)
                    del input_variant_dict[k]
                    continue

                if min_alt_coverage is not None and AD is not None:
                    for ad in AD:
                        if int(ad) < min_alt_coverage:
                            del input_variant_dict[k]
                            low_confident_set.add(k)
                            break

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

    snp_query_fp_qual_dict = defaultdict(float)
    snp_query_tp_qual_dict = defaultdict(float)
    snp_truth_fn_set = set()
    snp_truth_tp_set = set()

    indel_query_fp_set = set()
    indel_query_tp_set = set()
    indel_truth_fn_set = set()
    indel_truth_tp_set = set()

    fail_to_pass_filtering_count = 0
    for row in happy_vcf_unzip_process.stdout:
        if row[0] == '#':
            continue
        columns = row.strip().split()

        ctg, pos = columns[0], int(columns[1])
        key = (ctg, int(pos)) if ctg_name is None else int(pos)

        if ctg_name is not None and ctg != ctg_name:
            continue

        if key in low_confident_set:
            continue
        fail_to_pass_filtering = False
        if key in truths_dict:
            dp, ad_list, af_list = truths_dict[key]

            if min_coverage is not None and dp < min_coverage:
                fail_to_pass_filtering = True
                fail_to_pass_filtering_count += 1
                continue
            for ad, af in zip(ad_list, af_list):

                if min_alt_coverage is not None and ad < min_alt_coverage:
                    fail_to_pass_filtering = True
                    break
                if af < min_af:
                    fail_to_pass_filtering = True
                    break
        if fail_to_pass_filtering:
            fail_to_pass_filtering_count += 1
            continue

        FORMAT, TRUTH, QUERY = columns[8], columns[9], columns[10]
        FORMAT = FORMAT.split(':')
        TRUTH = TRUTH.split(':')
        QUERY = QUERY.split(':')

        ft_dict = dict(zip(FORMAT, TRUTH))
        fq_dict = dict(zip(FORMAT, QUERY))

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
            #switch fp and fn to tp
            if q_bk == 'am' and t_bk == 'am' and t_BD != "UNK" and q_BD != "UNK":
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
            qual = input_variant_dict[key].qual if key in input_variant_dict else None
            snp_query_fp_qual_dict[key] = float(qual)
        if is_query_snp_tp:
            qual = input_variant_dict[key].qual if key in input_variant_dict else None
            snp_query_tp_qual_dict[key] = float(qual)
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

    truth_all_tp = truth_snp_tp + truth_indel_tp
    truth_all_fn = truth_snp_fn + truth_indel_fn
    query_all_fp = query_snp_fp + query_indel_fp
    query_all_tp = query_snp_tp + query_indel_tp

    # p->precision, r->recall, f1->f1_score
    # a->overall, s->snp, id->indel, i->insertion, d->deletion
    ap, ar, af1 = metrics(query_fp=query_all_fp, query_tp=query_all_tp, truth_fn=truth_all_fn, truth_tp=truth_all_tp)
    sp, sr, sf1 = metrics(query_fp=query_snp_fp, query_tp=query_snp_tp, truth_fn=truth_snp_fn, truth_tp=truth_snp_tp)
    idp, idr, idf1 = metrics(query_fp=query_indel_fp, query_tp=query_indel_tp, truth_fn=truth_indel_fn, truth_tp=truth_indel_tp)
    ip, ir, if1 = metrics(query_fp=query_ins_fp, query_tp=query_ins_tp, truth_fn=truth_ins_fn, truth_tp=truth_ins_tp)
    dp, dr, df1 = metrics(query_fp=query_del_fp, query_tp=query_del_tp, truth_fn=truth_del_fn, truth_tp=truth_del_tp)

    print("[INFO] Total low confident variants in input/truths:{}/{}".format(len(low_confident_set),fail_to_pass_filtering_count))
    print (''.join([item.ljust(20) for item in ["VariantType", 'TRUTH.FP', 'TRUTH.FN', 'TRUTH.TP','QUERY.TP', 'METRIC.Precision', 'METRIC.Recall', 'METRIC.F1_Score']]), file=output_file)
    print (''.join([str(item).ljust(20) for item in ["Overall", query_all_fp, truth_all_fn, truth_all_tp, query_all_tp, ap, ar, af1]]), file=output_file)
    print (''.join([str(item).ljust(20) for item in ["SNP", query_snp_fp, truth_snp_fn, truth_snp_tp, query_snp_tp, sp, sr, sf1]]),file=output_file)
    print (''.join([str(item).ljust(20) for item in ["INDEL", query_indel_fp, truth_indel_fn, truth_indel_tp, query_indel_tp, idp, idr, idf1]]), file=output_file)
    print (''.join([str(item).ljust(20) for item in ["INS", query_ins_fp, truth_ins_fn, truth_ins_tp, query_ins_tp, ip, ir, if1]]), file=output_file)
    print (''.join([str(item).ljust(20) for item in ["DEL", query_del_fp, truth_del_fn, truth_del_tp, query_del_tp, dp, dr, df1]]), file=output_file)
    print('\n', file=output_file)

    if args.output_best_f1_score:

        results = output_best_cut_off(snp_query_fp_qual_dict, snp_query_tp_qual_dict, len(snp_truth_fn_set), truth_tp_count=len(snp_truth_tp_set), use_int_cut_off=args.use_int_cut_off)
        best_match = results[0].copy()
        best_match[0] = 'SNV(Best F1)'
        print(
            ''.join(
                [str(item).ljust(13) if idx >= 4 or idx == 0 else ('%.4f' % item).ljust(13) for idx, item in
                 enumerate(best_match)]),)
            # file=output_file)

        if args.debug:
            print("")
            for result in results:
                print(''.join(
                    [str(item).ljust(13) if idx >= 4 or idx == 0 else ('%.4f' % item).ljust(13) for idx, item in
                     enumerate(result)]),)
                    # file=output_file)

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

    parser.add_argument('--min_alt_coverage', type=int, default=None,
                        help="EXPERIMENTAL: Minimum alt base count for a variant to be included in bechmarking")

    parser.add_argument('--min_coverage', type=int, default=None,
                        help="EXPERIMENTAL: Minimum coverage for a variant to be included in bechmarking")

    parser.add_argument('--truths_info_fn', type=str, default=None,
                        help="Truth infomation file acquried form `cal_truth_vcf_af_distribution` submodule")

    parser.add_argument('--input_vcf_fn', type=str, default=None,
                        help="Input VCF file name")

    ## Only benchmark 'HighConf' tag in seqc VCF
    parser.add_argument('--skip_genotyping', type=str2bool, default=True,
                        help="Skip calculating VCF genotype")

    parser.add_argument('--min_qual', type=float, default=None,
                        help="Minimum quality score")

    parser.add_argument('--input_filter_tag', type=str_none, default=None,
                        help="Filter variants with tag from the input VCF")

    parser.add_argument('--output_best_f1_score', action='store_true',
                        help=SUPPRESS)

    parser.add_argument('--use_int_cut_off', type=str2bool, default=True,
                        help=SUPPRESS)

    parser.add_argument('--debug', action='store_true',
                        help=SUPPRESS)

    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    Cal(args)


if __name__ == "__main__":
    main()
