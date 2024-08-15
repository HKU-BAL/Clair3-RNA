# Parameters
PLATFORM="ont_guppy_drna002"
THREADS=4
INPUT_DIR="${HOME}/ont_quick_demo"
INPUT_DIR="/autofs/bal33/zxzheng/rna/ont_quick_demo"
OUTPUT_DIR="${INPUT_DIR}/output"

mkdir -p ${INPUT_DIR}
mkdir -p ${OUTPUT_DIR}

# Download quick demo data
# GRCh38_no_alt reference
wget -P ${INPUT_DIR} -nc http://www.bio8.cs.hku.hk/clair3_rna/quick_demo/ont/GRCh38_no_alt_chr1.fa
wget -P ${INPUT_DIR} -nc http://www.bio8.cs.hku.hk/clair3_rna/quick_demo/ont/GRCh38_no_alt_chr1.fa.fai
# RNA BAM
wget -P ${INPUT_DIR} -nc http://www.bio8.cs.hku.hk/clair3_rna/quick_demo/ont/HG004_chr1_demo.bam
wget -P ${INPUT_DIR} -nc http://www.bio8.cs.hku.hk/clair3_rna/quick_demo/ont/HG004_chr1_demo.bam.bai

# Truth VCF and BED
wget -P ${INPUT_DIR} -nc http://www.bio8.cs.hku.hk/clair3_rna/quick_demo/ont/HG004_GRCh38_1_22_v4.2.1_benchmark_chr1.vcf.gz
wget -P ${INPUT_DIR} -nc http://www.bio8.cs.hku.hk/clair3_rna/quick_demo/ont/HG004_GRCh38_1_22_v4.2.1_benchmark_chr1.vcf.gz.tbi
wget -P ${INPUT_DIR} -nc http://www.bio8.cs.hku.hk/clair3_rna/quick_demo/ont/HG004_GRCh38_1_22_v4.2.1_benchmark_chr1.bed

REF_FILE_PATH="${INPUT_DIR}/GRCh38_no_alt_chr1.fa"
BAM_FILE_PATH="${INPUT_DIR}/HG004_chr1_demo.bam"
BASELINE_VCF_FILE_PATH="${INPUT_DIR}/HG004_GRCh38_1_22_v4.2.1_benchmark_chr1.vcf.gz"
BASELINE_BED_FILE_PATH="${INPUT_DIR}/HG004_GRCh38_1_22_v4.2.1_benchmark_chr1.bed"
RNA_BED_FILE_PATH="${OUTPUT_DIR}/final.bed"
OUTPUT_VCF_FILE_PATH="${OUTPUT_DIR}/output.vcf.gz"

cd ${OUTPUT_DIR}

# Run Clair3-RNA in one command
docker run -it \
-v ${INPUT_DIR}:${INPUT_DIR} \
-v ${OUTPUT_DIR}:${OUTPUT_DIR} \
hkubal/clair3-rna:latest \
/opt/bin/run_clair3_rna \
    --bam_fn ${BAM_FILE_PATH} \
    --ref_fn ${REF_FILE_PATH} \
    --output_dir ${OUTPUT_DIR} \
    --platform ${PLATFORM} \
    --threads ${THREADS} \
    --region chr1:816000-828000 \
    --tag_variant_using_readiportal

# Acquire intersected GIAB high-confidence BED with RNA coverage >=4
docker run -it \
-v ${INPUT_DIR}:${INPUT_DIR} \
-v ${OUTPUT_DIR}:${OUTPUT_DIR} \
hkubal/clair3-rna:latest \
pypy3 /opt/bin/clair3_rna.py get_rna_bed \
    --output_dir ${OUTPUT_DIR} \
    --bam_fn ${BAM_FILE_PATH} \
    --threads 40 \
    --min_mq 5 \
    --min_bq 0 \
    --ctg_name chr1 \
    --high_confident_bed_fn ${BASELINE_BED_FILE_PATH} \
    --truth_vcf_fn ${BASELINE_VCF_FILE_PATH}

# Run hap.py for evaluation
docker run \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
    ${BASELINE_VCF_FILE_PATH} \
    ${OUTPUT_VCF_FILE_PATH} \
    -o ${OUTPUT_DIR}/happy \
    -V \
    --no-json \
    --no-write-counts \
    -r ${REF_FILE_PATH} \
    -f ${RNA_BED_FILE_PATH} \
    --threads ${THREADS} \
    --no-roc \
    --pass-only \
    --engine=vcfeval \
    -l chr1:816000-828000

# Get overall metrics with coverage cut-off
docker run -it \
-v ${INPUT_DIR}:${INPUT_DIR} \
-v ${OUTPUT_DIR}:${OUTPUT_DIR} \
hkubal/clair3-rna:latest \
pypy3 /opt/bin/clair3_rna.py calculate_overall_metrics \
    --happy_vcf_fn ${OUTPUT_DIR}/happy.vcf.gz \
    --input_vcf_fn ${OUTPUT_VCF_FILE_PATH} \
    --output_fn ${OUTPUT_DIR}/METRICS \
    --min_coverage 4 \
    --min_alt_coverage 2 \
    --min_af 0.05 \
    --skip_genotyping False \
    --ctg_name chr1 \
    --input_filter_tag 'PASS' \
    --truths_info_fn ${OUTPUT_DIR}/vcf_output/truths

