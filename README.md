
# Clair3-RNA - long-read short variant caller for RNA sequencing data

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

Contact: Ruibang Luo, Zhenxian Zheng  
Email: {rbluo,zxzheng}@cs.hku.hk

----

## Introduction

Clair3-RNA is a small variant caller for RNA long-read data. Clair3-RNA supports ONT R9 chemistry with complementary DNA sequencing (cDNA) and direct RNA sequencing (dRNA). It also support PacBio Sequel and PacBio MAS-Seq RNA sequencing data.

For germline small variant calling with DNA long-read, please use [Clair3](https://github.com/HKU-BAL/Clair3). 

For somatic small variant calling with DNA long-read, please try [ClairS](https://github.com/HKU-BAL/ClairS).

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

## Quick Demo

- Oxford Nanopore (ONT) data as input, see [ONT Quick Demo](docs/ont_quick_demo.md).
- PacBio HiFi data as input, see [PacBio HiFi Quick Demo](docs/pacbio_hifi_quick_demo.md).

### Quick start

After following [installation](#installation), you can run Clair3-RNA with one command:

```bash
./run_clair3_rna -B input.bam -R ref.fa -o output -t 8 -p ont_r9_cdna
## Final output file: output/output.vcf.gz
```

Check [Usage](#Usage) for more options.

----

------

## Pre-trained Models

Clair-RNA was trained using GIAB RNA sequencing data. All models were trained with chr20 excluded (including only chr1-19, 21, 22). 

|         Platform         |        Chemistry/Instruments        | Basecaller | Option (`-p/--platform`) |   Reference   | Training samples |
| :----------------------: | :----------------------------------: | :--------: | :----------------------: | :-----------: | ---------------- |
|           ONT            |    R9.4.1, complementary DNA sequencing     |   Guppy    |   `ont_r9_guppy_cdna`    | GRCh38_no_alt | HG002            |
|           ONT            | R9.4.1, direct RNA sequencing |   Guppy    |   `ont_r9_guppy_drna`    | GRCh38_no_alt | HG002            |
|       PacBio HIFI        |       Sequel with Iso-Seq kit        |     -      |      `hifi_sequel2`      | GRCh38_no_alt | HG002            |
|       PacBio HIFI        |        Revio with MAS-Seq kit        |     -      |        `hifi_mas`        | GRCh38_no_alt | HG002            |


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
  --platform ${PLATFORM} \               ## options: {ont_r9_guppy_cdna, ont_r9_guppy_drna, hifi_sequel2, hifi_mas}
  --output_dir ${OUTPUT_DIR}             ## output path prefix 
```

Check [Usage](#Usage) for more options.

### Option 2. Singularity

**Caution**: Absolute path is needed for both `INPUT_DIR` and `OUTPUT_DIR` in singularity. 

```bash
INPUT_DIR="[YOUR_INPUT_FOLDER]"        # e.g. /home/user1/input (absolute path needed)
OUTPUT_DIR="[YOUR_OUTPUT_FOLDER]"      # e.g. /home/user1/output (absolute path needed)
mkdir -p ${OUTPUT_DIR}

conda config --add channels defaults
conda create -n singularity-env -c conda-forge singularity -y
conda activate singularity-env

# singularity pull docker pre-built image
singularity pull docker://hkubal/clair3_rna:latest

# run the sandbox like this afterward
singularity exec \
  -B ${INPUT_DIR},${OUTPUT_DIR} \
  clair3_rna_latest.sif \
  hkubal/clair3_rna:latest \
  /opt/bin/run_clair3_rna \
  --bam_fn ${INPUT_DIR}/input.bam \            ## use your input bam file name here
  --ref_fn ${INPUT_DIR}/ref.fa \               ## use your reference file name here
  --threads ${THREADS} \                       ## maximum threads to be used
  --platform ${PLATFORM} \                     ## options: {ont_r9_guppy_cdna, ont_r9_guppy_drna, hifi_sequel2, hifi_mas}
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
# install pypy and packages in the environemnt
conda create -n clair3_rna -c conda-forge -c bioconda clair3 mosdepth bedtools -y
source activate clair3_rna

git clone https://github.com/HKU-BAL/Clair3-RNA.git
cd Clair-RNA

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
  --platform ${PLATFORM} \                   ## options: {ont_r9_guppy_cdna, ont_r9_guppy_drna, hifi_sequel2, hifi_mas}
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
                        Select the sequencing platform of the input. Possible options: {ont_r9_guppy_cdna, ont_r9_guppy_drna, hifi_sequel2, hifi_mas}.
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
  --python PYTHON       Absolute path of python, python3 >= 3.9 is required.
  --pypy PYPY           Absolute path of pypy3, pypy3 >= 3.6 is required.
  --samtools SAMTOOLS   Absolute path of samtools, samtools version >= 1.10 is required.
  --parallel PARALLEL   Absolute path of parallel, parallel >= 20191122 is required.
  --min_mq MIN_MQ       Minimal mapping quality required for an alignment to be considered. Default: 5.
```

#### Call variants in one or mutiple chromosomes using the `-C/--ctg_name` parameter

```bash
./run_clair3_rna -B input.bam -R ref.fa -o output -t 8 -p ont_r9_cdna -C chr21,chr22
```

#### Call variants in one specific region using the `-r/--region` parameter

```bash
./run_clair3_rna -B input.bam -R ref.fa -o output -t 8 -p ont_r9_cdna -r chr20:1000000-2000000
```

#### Call variants at interested variant sites (genotyping) using the `-G/--genotyping_mode_vcf_fn` parameter

```bash
./run_clair3_rna -B input.bam -R ref.fa -o output -t 8 -p ont_r9_cdna -G input.vcf
```

#### Call variants in the BED regions using the `-b/--bed_fn` parameter

We highly recommended using BED file to define multiple regions of interest like:

```bash
./run_clair3_rna -B input.bam -R ref.fa -o output -t 8 -p ont_r9_cdna -b input.bed
```

------

#### Tag RNA editing site using REDIportal dataset

RNA undergoes editing by ADAR (adenosine deaminases acting on RNA), resulting in Adenosine-to-inosine (A-to-I) changes. These A-to-I changes can be observed in RNA-seq datasets as A-to-G and T-to-C changes, which do not represent genuine RNA variants. To address this, we provide users with the option to utilize external datasets such as [REDIportal](http://srv00.recas.ba.infn.it/atlas/) to annotate RNA editing sites. In Clair3-RNA's VCF output, variants that are also RNA editing sites reported in REDIportal can be tagged. These sites will be marked as `RNAEditing` instead of `PASS` in the `FILTER` column when the `--tag_variant_using_readiportal` option is enabled.

