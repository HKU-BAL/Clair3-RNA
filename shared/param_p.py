#Clair3-RNA pileup parameters
caller_name=REPO_NAME="Clair3-RNA"
version = "0.0.2"
import re
from itertools import accumulate

zstd='zstd'
default_optimizer = "Radam"
default_loss_function = "FocalLoss"
support_platform = {'ont', 'hifi','ilmn'}
min_af = 0.08
min_af_dict = {'ont':0.15, 'hifi':min_af, 'ilmn':min_af }
#as three platform training data vary in depth distribution, we recommend below max_depth base on max training data depth for calling
max_depth = 144
max_depth_dict = {'ont':max_depth, 'hifi':max_depth, 'ilmn':max_depth}
maximum_variant_length_that_need_infer = 50
maximum_variant_length_that_need_infer_include_long_indel = 100000
cal_precise_long_indel_af = False
long_indel_distance_proportion = 0.1
min_mq = 5
min_bq = 0
min_coverage = 2
tensorflow_threads = 4

#GVCF parameters
base_err = 0.001
gq_bin_size = 5

#Pileup input feature list
#           0    1    2    3    4    5    6    7     8    9    10   11  12   13    14  15   16    17
channel = ('A', 'C', 'G', 'T', 'I', 'I1', 'D', 'D1', '*', 'a', 'c', 'g','t', 'i', 'i1','d', 'd1','#')
channel_size = len(channel)
flankingBaseNum = 16
no_of_positions = 2 * flankingBaseNum + 1
ont_input_shape = input_shape = [no_of_positions, channel_size]
label_shape = [21, 3, no_of_positions, no_of_positions]
label_size = sum(label_shape)
label_shape_cum = list(accumulate(label_shape))
expandReferenceRegion = 1000
SAMTOOLS_VIEW_FILTER_FLAG = 2316
partition_size = 500000
region_size =1000
phasing_window_size = 30000
extend_bp=10

#Training hyperparameters
chunk_size = 200
trainBatchSize = 2000
predictBatchSize = 200
initialLearningRate = 1e-3
l2RegularizationLambda = 1e-7
trainingDatasetPercentage = 0.9
maxEpoch = 30
OPERATION_SEED = None
RANDOM_SEED = None

model_name_platform_dict = {
    'ont_r9_guppy_cdna': 'ont_guppy_cdna',
    'ont_r9_guppy_drna': 'ont_guppy_drna002',
    'hifi_sequel2': 'hifi_sequel2_pbmm2',
    'hifi_mas': 'hifi_mas_pbmm2'
}

min_thred_qual = {'ont': 8,
                  'hifi': 2}

snp_min_af = 0.08
indel_min_af = 0.15
min_coverage = 4
CHUNK_SIZE = 5000000
qual_cut_off = 2
readiportal_database_filter_tag = "A,D:A,R:A,R,D"
