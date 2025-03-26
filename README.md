<div align="center">
    <img src="images/icon.png" width = "200" title="Clair3-RNA">
    <img src="images/fjsClaiRR.jpg" width = "200" title="aka. ClaiRR. Image credits to Fritz Sedlezeck.">
</div>

# Clair3-RNA - A deep learning-based small variant caller for long-read RNA sequencing data

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

Contact: Ruibang Luo, Zhenxian Zheng  
Email: {rbluo,zxzheng}@cs.hku.hk

----

## Introduction

Clair3-RNA is a small variant caller for **long-read RNA sequencing (lrRNA-seq)** data. Clair3-RNA supports ONT R10.4.1 and R9.4.1 complementary DNA sequencing (cDNA) and direct RNA sequencing (dRNA). dRNA sequencing support the ONT latest [SQK-RNA004 kit](https://community.nanoporetech.com/docs/prepare/library_prep_protocols/direct-rna-sequencing-sqk-rna004/v/drs_9195_v4_revd_20sep2023) data for variant calling. Clair3-RNA also supports PacBio Sequel and PacBio MAS-Seq RNA sequencing data. Clair3-RNA reached a ~95% F1-score for ONT dRNA using SQK-RNA004 kit and ~96% F1-score using PacBio Iso-Seq and MAS-Seq, respectively, with at least ten supporting reads and disregarding the zygosity. With read phased, the performance reached ~97% for ONT and ~98% for PacBio.

A preprint describing Clair3-RNA's algorithms and results is at [bioRxiv](https://doi.org/10.1101/2024.11.17.624050).

For germline small variant calling, please use [Clair3](https://github.com/HKU-BAL/Clair3). 

For somatic small variant calling using a tumor-normal pair, please try [ClairS](https://github.com/HKU-BAL/ClairS).

For somatic small variant calling using tumor sample only, please try [ClairS-TO](https://github.com/HKU-BAL/ClairS-TO).

----

## Contents

* [Pre-trained Models](#pre-trained-models)
- [Installation](#installation)
  - [Option 1. Docker pre-built image](#option-1--docker-pre-built-image)
  - [Option 2. Singularity](#option-2-singularity)
  - [Option 3. Build an anaconda virtual environment](#option-3-build-an-anaconda-virtual-environment)
  - [Option 4. Docker Dockerfile](#option-4-docker-dockerfile)
- [Quick Demo](#quick-demo)
- [Usage](#usage)


----

## Latest Updates
*v0.2.2 (Mar 13, 2025)*: 1. Added a new ONT R10 cDNA model (`ont_r10_dorado_cdna`) for R10.4.1 chemistry. 2. Added the  `--enable_padding_in_splice_junction_regions` option to enable padding in pileup input tensor in splice junction regions and exon boundaries. 3. Added chromosomes X and Y to the default contigs when the `--include_all_ctgs` option is not specified.

*v0.2.1 (Dec 2, 2024)*: 1. Fixed a bug that misses some variants when `--print_ref_calls` is enabled ([#6](https://github.com/HKU-BAL/Clair3-RNA/issues/6)). 2. Added the  `--enable_variant_calling_at_sequence_head_and_tail` option to enable variant calling at the head and tail 16bp of each sequence. Use with caution because alignments are less reliable in the regions, and there would be insufficient context to be fed to the neural network for reliable calling.

*v0.2.0 (Nov 18, 2024)* : 1. Added a new pileup phasing model (enable by using `--enable_phasing_model` opiton) for ONT dRNA004(`ont_dorado_drna004`), PacBio Iso-Seq(`hifi_sequel2_minimap2`), and PacBio MAS-Seq(`hifi_mas_minimap2`), the SNP performance improved by ~2% and Indel performance improved by ~6%. 2. Fixed some formatting issues in the calling workflow.

*v0.1.0 (Aug 15, 2024)* : 1. Added a new ONT dRNA004 direct RNA sequencing model (`ont_dorado_drna004`) for SQK-RNA004 kit. 2. Added new PacBio Sequel (`hifi_sequel2_minimap2`) and Revio (`hifi_mas_minimap2`) model to support minimap2 alignment. 3. Enhance model training techniques to boost performance by incorporating strategies such as managing low-coverage sites, verifying variant zygosity, filtering RNA editing sites, etc.  4. Renamed all ONT and PacBio model names, check [here](https://github.com/HKU-BAL/Clair3-RNA?tab=readme-ov-file#pre-trained-models) for more details.

*v0.0.1 (Nov 27, 2023)*: Initial release for early access.

---

## Quick Demo

- Oxford Nanopore (ONT) data as input, see [ONT Quick Demo](docs/ont_quick_demo.md).
- PacBio HiFi data as input, see [PacBio HiFi Quick Demo](docs/pacbio_hifi_quick_demo.md).

### Quick start

After following [installation](#installation), you can run Clair3-RNA with one command:

```bash
./run_clair3_rna -B input.bam -R ref.fa -o output -t 8 -p ont_dorado_drna004
## Final output file: output/output.vcf.gz
```

Check [Usage](#Usage) for more options.

----

------

## Pre-trained Models

Clair-RNA was trained using GIAB RNA sequencing data. All models were trained with chr20 excluded (including only chr1-19, 21, 22). 

|         Platform         |       Chemistry/Kit/Instruments       | Basecaller | Aligner |           Support  phasing model (`--enable_phasing_model`)?           |           Option (`-p/--platform`)            |   Reference   | Training samples |
| :----------------------: |:-------------------------------------:| :--------: | :----------------------: |:---------------------------------------------:| ---------------- | :--------------: | ---------------- |
| ONT | SQK-RNA004 kit, direct RNA sequencing | Dorado | minimap2 |             Yes             |             `ont_dorado_drna004`              | GRCh38 | HG002 |
| ONT | R10.4.1, complementary DNA sequencing | Dorado | minimap2 |             Yes           |             `ont_r10_dorado_cdna`           | GRCh38 | HG002 |
| ONT | SQK-RNA002 kit, direct RNA sequencing | Guppy | minimap2 |                            |              `ont_guppy_drna002`              | GRCh38 | HG002 |
|           ONT            | R9.4.1, complementary DNA sequencing  |   Guppy    | minimap2 |                               |               `ont_guppy_cdna`                | GRCh38 | HG002       |
|       PacBio HiFi       |        Sequel with Iso-Seq kit        |     -      | pbmm2/minimap2 | Yes for `hifi_sequel2_minimap2` | `hifi_sequel2_pbmm2`, `hifi_sequel2_minimap2` | GRCh38 | HG002            |
|       PacBio HiFi    |        Revio with MAS-Seq kit         |     -      | pbmm2/minimap2 |     Yes for `hifi_mas_minimap2`     |     `hifi_mas_pbmm2`, `hifi_mas_minimap2`     | GRCh38 | HG002            |


------


## Installation

### Option 1.  Docker pre-built image

A pre-built docker image is available at [DockerHub](https://hub.docker.com/r/hkubal/clairs). 

**Caution**: Absolute path is needed for both `INPUT_DIR` and `OUTPUT_DIR` in docker. 

```bash
docker run -it \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hkubal/clair3-rna:latest \
  /opt/bin/run_clair3_rna \
  --bam_fn ${INPUT_DIR}/input.bam \      ## use your input bam file name here
  --ref_fn ${INPUT_DIR}/ref.fa \         ## use your reference file name here
  --threads ${THREADS} \                 ## maximum threads to be used
  --platform ${PLATFORM} \               ## options: {ont_dorado_drna004, ont_r10_dorado_cdna, ont_guppy_drna002, ont_guppy_cdna, hifi_sequel2_pbmm2, hifi_sequel2_minimap2, hifi_mas_pbmm2, hifi_sequel2_minimap2}
  --tag_variant_using_readiportal \      ## optional, tag variants using REDIportal dataset
  --enable_phasing_model \               ## optional, enable calling using phasing model
  --output_dir ${OUTPUT_DIR}             ## output path prefix 
```

Check [Usage](#Usage) for more options.

### Option 2. Singularity

**Caution**: 1. Absolute path is needed for both `INPUT_DIR` and `OUTPUT_DIR` in singularity. 2. Please add `--conda_prefix /opt/conda/envs/clair3_rna` to specify the conda environment path.

```bash
INPUT_DIR="[YOUR_INPUT_FOLDER]"        # e.g. /home/user1/input (absolute path needed)
OUTPUT_DIR="[YOUR_OUTPUT_FOLDER]"      # e.g. /home/user1/output (absolute path needed)
mkdir -p ${OUTPUT_DIR}

conda config --add channels defaults
conda create -n singularity-env -c conda-forge singularity -y
conda activate singularity-env

# singularity pull docker pre-built image
singularity pull docker://hkubal/clair3-rna:latest

# run the sandbox like this afterward
singularity exec \
  -B ${INPUT_DIR},${OUTPUT_DIR} \
  clair3-rna_latest.sif \
  /opt/bin/run_clair3_rna \
  --bam_fn ${INPUT_DIR}/input.bam \            ## use your input bam file name here
  --ref_fn ${INPUT_DIR}/ref.fa \               ## use your reference file name here
  --threads ${THREADS} \                       ## maximum threads to be used
  --platform ${PLATFORM} \                     ## options: {ont_dorado_drna004, ont_r10_dorado_cdna, ont_guppy_drna002, ont_guppy_cdna, hifi_sequel2_pbmm2, hifi_sequel2_minimap2, hifi_mas_pbmm2, hifi_sequel2_minimap2}
  --tag_variant_using_readiportal \            ## optional, tag variants using REDIportal dataset
  --enable_phasing_model \                     ## optional, enable calling using phasing model 
  --output_dir ${OUTPUT_DIR} \                 ## output path prefix
  --conda_prefix /opt/conda/envs/clair3_rna
```

### Option 3. Build an anaconda virtual environment

**Anaconda install**:

Please install anaconda using the official [guide](https://docs.anaconda.com/anaconda/install) or using the commands below:

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x ./Miniconda3-latest-Linux-x86_64.sh 
./Miniconda3-latest-Linux-x86_64.sh
```

**Install Clair3-RNA using anaconda step by step:**

```bash
# create and activate an environment named clair3_rna
# install pypy and packages in the environment
conda create -n clair3_rna -c conda-forge -c bioconda clair3 mosdepth bedtools -y
source activate clair3_rna

git clone https://github.com/HKU-BAL/Clair3-RNA.git
cd Clair3-RNA

# make sure in conda environment
# download pre-trained models
echo ${CONDA_PREFIX}
mkdir -p ${CONDA_PREFIX}/bin/clair3_rna_models
wget http://www.bio8.cs.hku.hk/clair3_rna/models/clair3_rna_models.tar.gz
tar -zxvf clair3_rna_models.tar.gz -C ${CONDA_PREFIX}/bin/clair3_rna_models/

./run_clair3_rna --help
```

### Option 4. Docker Dockerfile

This is the same as option 1 except that you are building a docker image yourself. Please refer to option 1 for usage. 

```bash
git clone https://github.com/HKU-BAL/Clair3-RNA.git
cd Clair-RNA

# build a docker image named hkubal/clairs:latest
# might require docker authentication to build docker image
docker build -f ./Dockerfile -t hkubal/clair3-rna:latest .

# run the docker image like option 1
docker run -it hkubal/clair3-rna:latest /opt/bin/clair3_rna --help
```

------

## Usage

### General Usage

```bash
./run_clair3_rna \
  --bam_fn ${INPUT_DIR}/input.bam \          ## use your input bam file name here
  --ref_fn ${INPUT_DIR}/ref.fa \             ## use your reference file name here
  --threads ${THREADS} \                     ## maximum threads to be used
  --platform ${PLATFORM} \                   ## options: {ont_dorado_drna004, ont_r10_dorado_cdna, ont_guppy_drna002, ont_guppy_cdna, hifi_sequel2_pbmm2, hifi_sequel2_minimap2, hifi_mas_pbmm2, hifi_sequel2_minimap2}
  --tag_variant_using_readiportal \          ## optional, tag variants using REDIportal dataset
  --enable_phasing_model \                   ## optional, enable calling using phasing model
  --output_dir ${OUTPUT_DIR}                 ## output path prefix
 
## Final output file: ${OUTPUT_DIR}/output.vcf.gz
```

### Options

**Required parameters:**

```bash
  -B BAM_FN, --bam_fn BAM_FN
                        RNA BAM file input. The input file must be samtools indexed.
  -R REF_FN, --ref_fn REF_FN
                        FASTA reference file input. The input file must be samtools indexed.
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        VCF output directory.
  -t THREADS, --threads THREADS
                        Max #threads to be used.
  -p PLATFORM, --platform PLATFORM
                        Select the sequencing platform of the input. Possible options: {ont_dorado_drna004, ont_r10_dorado_cdna, ont_guppy_drna002, ont_guppy_cdna, hifi_sequel2_pbmm2, hifi_sequel2_minimap2, hifi_mas_pbmm2, hifi_sequel2_minimap2}.
```

**Miscellaneous parameters:**

```bash
  -P PILEUP_MODEL_PATH, --pileup_model_path PILEUP_MODEL_PATH
                        Specify the path prefix to your own pileup model. Including ${pileup_model_path}.data-00000-of-00001, ${pileup_model_path}.index.
  -c CTG_NAME, --ctg_name CTG_NAME
                        The name of the contigs to be processed. Split by ',' for multiple contigs. Default: call in chr{1..22} and {1..22}.
  -r REGION, --region REGION
                        A region to be processed. Format: `ctg_name:start-end` (start is 1-based).
  -b BED_FN, --bed_fn BED_FN
                        Path to a BED file. Call variants only in the provided BED regions.
  -G GENOTYPING_MODE_VCF_FN, --genotyping_mode_vcf_fn GENOTYPING_MODE_VCF_FN
                        VCF file input containing candidate sites to be genotyped. Variants will only be called at the sites in the VCF file if provided.
  -q QUAL, --qual QUAL  If set, variants with >QUAL will be marked as PASS, or LowQual otherwise.
  --enable_phasing_model
                        Enable phasing with whatshap or longphase. Usually leads to performance improvement when coverage is sufficient. Default: False.  
  --snp_min_af SNP_MIN_AF
                        Minimal SNP AF required for a variant to be called. Decrease SNP_MIN_AF might increase a bit of sensitivity, but in trade of precision, speed and accuracy. Default: 0.08.
  --indel_min_af INDEL_MIN_AF
                        Minimum Indel AF required for a candidate variant to be called. Default: ont:0.15,hifi:0.08.
  --min_coverage MIN_COVERAGE
                        Minimal coverage required for a variant to be called. Default: 4.
  --tag_variant_using_readiportal
                        Tag variants uisng REDIportal dataset, If set, called variants that are also in the readiportal dataset will be marked as "RNAEditing". Default: disable.
  --readiportal_database_filter_tag READIPORTAL_DATABASE_FILTER_TAG
                        Use only editing sites with these tags in the readiportal dataset, split by ":" for multiple tags. Default: using sites supported by two or more sources - "A,D:A,R:A,R,D".
  --readiportal_reference_genome_version READIPORTAL_REFERENCE_GENOME_VERSION
                        Select the reference genome version in the readiportal dataset. Possible options: {grch38, grch37}. Default: "grch38".
  --chunk_size CHUNK_SIZE
                        The size of each chuck for parallel processing. Default: 5000000.
  -s SAMPLE_NAME, --sample_name SAMPLE_NAME
                        Define the sample name to be shown in the VCF file. Default: SAMPLE.
  --output_prefix OUTPUT_PREFIX
                        Prefix for output VCF filename. Default: output.
  --remove_intermediate_dir
                        Remove intermediate directory before finishing to save disk space.
  --include_all_ctgs    Call variants on all contigs, otherwise call in chr{1..22} and {1..22}.
  --print_ref_calls     Show reference calls (0/0) in VCF file.
  -d, --dry_run         Print the commands that will be ran.
  --min_mq MIN_MQ       Minimal mapping quality required for an alignment to be considered. Default: 5.
  --phased_pileup_model_path PHASED_PILEUP_MODEL_PATH
                        Specify the path prefix to your own pileup phasing model. Including ${phased_pileup_model_path}.data-00000-of-00001, ${phased_pileup_model_path}.index.
  --python PYTHON       Absolute path of python, python3 >= 3.9 is required.
  --pypy PYPY           Absolute path of pypy3, pypy3 >= 3.6 is required.
  --samtools SAMTOOLS   Absolute path of samtools, samtools version >= 1.10 is required.
  --parallel PARALLEL   Absolute path of parallel, parallel >= 20191122 is required.
  --longphase LONGPHASE
                        Absolute path of longphase, longphase >= 1.7 is required.
  --whatshap WHATSHAP   Absolute path of whatshap, whatshap >= 1.0 is required.
  --use_longphase_for_intermediate_phasing USE_LONGPHASE_FOR_INTERMEDIATE_PHASING
                        Use longphase for intermediate phasing. Default:False.
  --use_longphase_for_intermediate_haplotagging USE_LONGPHASE_FOR_INTERMEDIATE_HAPLOTAGGING
                        Use longphase for intermediate haplotagging. Default:False.
  --enable_variant_calling_at_sequence_head_and_tail
                        EXPERIMENTAL: Enable variant calling at the head and tail 16bp of each sequence. Default: disable.
```

#### Call variants in one or multiple chromosomes using the `-C/--ctg_name` parameter

```bash
./run_clair3_rna -B input.bam -R ref.fa -o output -t 8 -p ont_dorado_drna004 -C chr21,chr22
```

#### Call variants in one specific region using the `-r/--region` parameter

```bash
./run_clair3_rna -B input.bam -R ref.fa -o output -t 8 -p ont_dorado_drna004 -r chr20:1000000-2000000
```

#### Call variants at interested variant sites (genotyping) using the `-G/--genotyping_mode_vcf_fn` parameter

```bash
./run_clair3_rna -B input.bam -R ref.fa -o output -t 8 -p ont_dorado_drna004 -G input.vcf
```

#### Call variants in the BED regions using the `-b/--bed_fn` parameter

We highly recommended using BED file to define multiple regions of interest like:

```bash
./run_clair3_rna -B input.bam -R ref.fa -o output -t 8 -p ont_dorado_drna004 -b input.bed
```

------

#### Tag RNA editing site using REDIportal dataset

RNA undergoes editing by ADAR (adenosine deaminases acting on RNA), resulting in Adenosine-to-inosine (A-to-I) changes. These A-to-I changes can be observed in RNA-seq datasets as A-to-G and T-to-C changes, which do not represent genuine RNA variants. To address this, we provide users with the option to utilize external datasets such as [REDIportal](http://srv00.recas.ba.infn.it/atlas/) to annotate RNA editing sites. In Clair3-RNA's VCF output, variants that are also RNA editing sites reported in REDIportal can be tagged. These sites will be marked as `RNAEditing` instead of `PASS` in the `FILTER` column when the `--tag_variant_using_readiportal` option is enabled.

**Caution**: `--tag_variant_using_readiportal` option currently works for GRCh38 and GRCh37 reference genome only, use can specify the reference genome version by using option `--readiportal_reference_genome_version={grch38, grch37}`.

