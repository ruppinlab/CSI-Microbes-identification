# CSI-Microbes Identification

This repository contains part of the workflows for reproducing the results from the Science Advances paper [Identification of intracellular bacteria from multiple single-cell RNA-seq platforms using CSI-Microbes]([https://www.biorxiv.org/content/10.1101/2020.05.14.096230v3](https://www.science.org/doi/10.1126/sciadv.adj7402)) by Welles Robinson, Josh Stone, Fiorella Schischlik, Billel Gasmi, Michael Kelly, Charlie Seibert, Kimia Dadkhah, E. Michael Gertz, Joo Sang Lee, Kaiyuan Zhu, Lichun Ma, Xin Wang, S. Cenk Sahinalp, Rob Patro, Mark D.M. Leiserson, Curtis Harris, Alejandro A. Schäffer, and Eytan Ruppin. This repository contains the workflows to identify microbial reads from 10x and Smart-seq2 scRNA-seq datasets. These microbial reads can then be analyzed using the [CSI-Microbes-analysis repository](https://github.com/ruppinlab/CSI-Microbes-analysis). The code in this repository was written by Welles Robinson and alpha-tested by Alejandro Schaffer.

## Requirements

The workflows in this repository are set-up to be run specifically on Biowulf, which is a Linux cluster provided by the NIH for intramural use. This workflow expects that conda has been installed. For instructions on how to install conda, see [conda install documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/install/).

## Software Installation

It should take 5 minutes to install the software, which involves downloading the codebase from GitHub using the below instructions. The below instructions assume that you have an ssh key associated with your GitHub account. If you do not, you can generate a new ssh key and associate it with your GitHub username by following [these instructions](https://docs.github.com/en/github/authenticating-to-github/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent).

```
git clone git@github.com:ruppinlab/CSI-Microbes-identification.git
```

The above command may return an error such as

```
Permission denied (publickey).
fatal: Could not read from remote repository.
Please make sure you have the correct access rights and the repository exists.
```

In this case, your ssh key may not be associated with your GitHub account. Please follow the instructions above and re-try. The below git submodule commands assume that your ssh has been associated with your GitHub account.

This repository uses two [git submodules](https://git-scm.com/book/en/v2/Git-Tools-Submodules) ([RNA-snakemake-rules](https://github.com/ruppinlab/RNA-snakemake-rules) and [pathogen-discovery-rules](https://github.com/ruppinlab/pathogen-discovery-rules)), which need to be initialized and fetched using the following commands.

```
cd CSI-Microbes-identification
git submodule init
git submodule update
```

## Software Dependencies

CSI-Microbes-identification depends on the following software packages that are loaded via the Biowulf module system: snakemake (6.0.5)<sup>[REF](#Snakemake)</sup>, sratoolkit (2.10.9)<sup>[REF](#SRAToolkit)</sup>, cellranger (5.0.1)<sup>[REF](#CellRanger)</sup>, samtools (1.11)<sup>[REF](#SAMtools)</sup>, bedtools (2.29.2)<sup>[REF](#BedTools)</sup>, picard (latest=2.25.0)<sup>[REF](#Picard)</sup> and  and the following conda software packages from the conda-forge, bioconda and defaults channels: fastp (0.20.1)<sup>[REF](#Fastp)</sup>, STAR (2.7.8a)<sup>[REF](#STAR), and pysam (1.16.0)<sup>[REF](#pysam)</sup>.

Currently, we support three distinct approaches for quantifying the number of reads assigned to one or more microbial taxa: PathSeq, CAMMiQ and SRPRISM. Running CSI-Microbes-identification with PathSeq requires the above packages as well as GATK (4.1.8.1) <sup>[REF](#PathSeq)</sup>. Running CSI-Microbes-identification with SRPRISM requires the above packages as well as SRPRISM (3.1.2)<sup>[REF](#SRPRISM)</sup>, which must be installed via [https://github.com/ncbi/SRPRISM](https://github.com/ncbi/SRPRISM). Running CSI-Microbes-identification with CAMMiQ requires the above packages as well as CAMMiQ (0.1)<sup>[REF](#CAMMiQ)</sup>, which must be installed via [https://github.com/algo-cancer/CAMMiQ](https://github.com/algo-cancer/CAMMiQ), and Gurobi (tested using both 9.0.0 and 9.1.0 due to license limitations) <sup>[REF](#Gurobi)</sup> and gcc (7.4.0). GATK (4.1.8.1) <sup>[REF](#PathSeq)</sup> is also required to run CSI-Microbes-identification with CAMMiQ for now although this is due to how we structured the workflow and not an inherent dependence of CAMMiQ.

## Database Dependencies

The workflows assume that the files needed by PathSeq and CAMMiQ are pre-built and their location is specified in `config/PathSeq-config.yaml`. The index used by PathSeq in this project is ~41 GB while the indices used by CAMMiQ are > 200 GB so we do not distribute them with the rest of the package although they are available upon request.

The existing projects assume that these database files will be in the top-level directory so the first step is to create the directory `CSI-Microbes-identification/data` (if it doesn't already exist) using

```
mkdir data
```

The PathSeq files can be divided into two groups: the microbial reference files and the host (human) reference files. To download the host reference files (which can take a while), use the below commands

```
cd data
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/pathseq/pathseq_host.tar.gz
tar -xf pathseq_host.tar.gz
cd ..
```

 The microbial reference fasta file (`microbev1.fa`), VecScreen hits (`microbev1-vecscreen-combined-matches.bed`) and taxonomy hierarchy file (`microbev1_release201_taxonomy.db`) used in our paper are available from [zenodo](https://zenodo.org/record/5604433). The microbial reference fasta file (`microbev1.fa`) can be used to build the additional microbial files required by PathSeq (`microbev1.fa.fai`, `microbev1.dict` and `microbev1.fa.img`) using the below command.

 ```
 ./scripts/run-build-PathSeq-microbe-files.sh
 ```

 Please note that building the `microbev1.fa.img` file will take ~1 day.

### Construction of PathSeq files

An example of how to build the PathSeq index (and other required files) is available in `build-PathSeq-microbes-files/Snakefile`.

### Construction of CAMMiQ files

The fastest way to obtain the index (of unique and doubly unique substrings) for CAMMiQ is to run

```
./cammiq --build --both -f <MAP_FILE> -D <FASTA_DIR> -k <int> -L <int> -Lmax <int> -h <int> -i <INDEX_FILES> -t <int>
```

Required arguments:

```
-f <MAP_FILE> gives a list of reference genomes in fasta format e.g., all/selected complete genomes in RefSeq for the bacterial, archaeal, and viral domains (can be downloaded with `CAMMiQ-download`), which constitute CAMMiQ's database, possibly along with NCBI's taxonomic information.
-D <FASTA_DIR> should contain the list of (fasta) file names given in <MAP_FILE>.
-k <int> specifies the minimum length of a unique and doubly-unique substring to be considered in CAMMiQ index. Default value is 26.
-L <int> specifies the potential read length in a query supported by CAMMiQ index. Default value is 100.
-Lmax <int> specifies the maximum length of a unique or doubly-unique substring to be considered in CAMMiQ index. Default value is 50.
-h <int> specifies the length of the common prefixes of the unique or doubly-unique substrings to be hashed. Default value is 26, and the value of h is required to be less than or equal to k.
-i <INDEX_FILES> include two file names (with directory): .bin1 specifies the index consisting of unique substrings and *.bin2 specifies the index specifies the index consisting of doubly unique substrings. The default file names (if ```-i``` is not specified) are index_u.bin1 and index_d.bin2.
-t <int> specifies the number of threads used during CAMMiQ's index construction. Note that CAMMiQ uses OpenMP during its index construction, which, by default, is 'auto-threaded' (i.e. attempts to use all available CPUs on a computer).
```

The input lines in `<MAP_FILE>` should contain at least the 4 tab-delimited fields in the below example

| File names | Genome IDs  | NCBI taxonomic IDs | Organism names |
| ---------------------- | ------------ | ------- | ----------------------------------- |
| GCF_000010525.1_ASM1052v1_genomic.fna | 1 | 7 | Azorhizobium caulinodans ORS 571 |
| GCF_000007365.1_ASM736v1_genomic.fna | 2 | 9 | Buchnera aphidicola str. Sg (Schizaphis graminum) |
| GCF_000218545.1_ASM21854v1_genomic.fna | 3 | 11 | Cellulomonas gilvus ATCC 13127 |
| GCF_000020965.1_ASM2096v1_genomic.fna | 4 | 14 | Dictyoglomus thermophilum H-6-12 |
| GCF_000012885.1_ASM1288v1_genomic.fna | 5 | 19 | Pelobacter carbinolicus DSM 2380 |

Other detailed format requirement for ```<MAP_FILE>``` can be found at [https://github.com/algo-cancer/CAMMiQ](https://github.com/algo-cancer/CAMMiQ);



## Running CSI-Microbes-identification

Currently, CSI-Microbes-identification can be run with three distinct approaches for quantifying the number of reads assigned to one or more microbial genomes: PathSeq, CAMMiQ and SRPRISM. PathSeq and CAMMiQ are used to quantify the number of reads assigned to a large number of microbial genomes while SRPRISM is an error-tolerant and ambiguity-character-tolerant aligner used to align reads against a single microbial genome and identify the genome location of reads. PathSeq uses a (more computationally expensive) alignment-based approach and reports matches across all levels of the taxonomy while CAMMiQ uses a much faster _k_-mer-based approach and reports matches at the specificied taxonomic level of interest.

The code in this pipeline is designed to analyze scRNA-seq datasets generated using either Smart-Seq2 (or similar plate-based approaches) or 10x although it should be straightforward to extend this pipeline to analyze scRNA-seq datasets from additional approaches.

### Test for CSI-Microbes-identification using PathSeq

We have created two small test directories: `test-10x` and `test-SS2` for users to quickly test that they can properly run CSI-Microbes-identification using PathSeq.

#### 10x test run

The 10x test dataset consists of all reads from 12 cells sequenced by Ben-Moshe2019<sup>[REF](#BenMoshe2019)</sup>. The input fastq files were generated by subsetting the CellRanger BAM file and are located at `test-10x/FASTQ/raw/Pt0/` (Yes, they were added to git. No, I do not feel good about it.). To run CSI-Microbes-identification on this dataset, use the below command, which submits the job running the Snakemake instance to the server.

```
./scripts/run-test-10x.sh
```

The job that runs the Snakemake instance will also submit additional jobs to the server. In total, these jobs should last ~20 minutes (depending on waiting times for job allocation). The key output directories are placed in the directory `output/PathSeq`. This directory should contain 13 directories: one directory for the microbial reads found in the sample (`output/PathSeq/{patient}-{sample}/`) and one directory for each of the 12 cells processed (`output/PathSeq/{patient}-{sample}-{cell}/`). These directories contain the `pathseq.txt` file that contains the microbial reads identified in this cell. To check the results of your run, you can compare against expected results in `expected_output/PathSeq` using the following command (from within the `test-10x` directory)

```
diff -r output/PathSeq/ expected_output/PathSeq/
```

which should yield a list of files that exist only in `output/PathSeq/` but no differences between the files that exist in both directories.

#### Smart-seq2 test run

The Smart-seq2 test dataset consists of all reads from 15 cells sequenced by Aulicino2018<sup>[REF](#Aulicino2018)</sup>. To run CSI-Microbes-identification on this dataset, use the below command to run the script, which submits the job running the Snakemake instance to the server.

```
./scripts/run-test-SS2.sh
```

The job that runs the Snakemake instance will also submit additional jobs to the server. In total, these jobs should last 1.5 hours (50 minutes of this is spent building the STAR index). The key output directories are placed in the directory `output/PathSeq`. This directory should contain 16 directories: one directory for the microbial reads found in the plate (`output/PathSeq/{patient}-{sample}-{plate}/`) and one directory for each of the 15 cells processed (`output/PathSeq/{patient}-{sample}-{plate}-{cell}/`). These directories contain the `pathseq.txt` file that contains the microbial reads identified in this cell. To check the results of your run, you can compare against expected results in `expected_output/PathSeq` using the following command (from within the `test-SS2` directory)

```
diff -r output/PathSeq/ expected_output/PathSeq/
```

which should yield a list of files that exist only in `output/PathSeq/` but no differences between the files that exist in both directories.


### Comparing CSI-Microbes, INVADE-seq and SAHMI on Robinson2023-10x

As part of the manuscript, we compare CSI-Microbes with two alternative approaches for identifying microbial reads from scRNA-seq (SAHMI<sup>[REF](#Ghaddar2022)</sup> and INVADE-seq<sup>[REF](#GaleanoNino2022)</sup>). SAHMI, CSI-Microbes and INVADE-seq all download the same Robinson2023-10x FASTQ files from SRA so you should run one tool first (which will download the files) and then after the files have been downloaded, run the remaining tools (which will use the already downloaded files).

#### Running SAHMI on Robinson2023-10x

[SAHMI](https://github.com/sjdlabgroup/SAHMI) involves running Kraken2Uniq first. To build the required input files for Kraken2Uniq, you need to navigate to the directory `build-Kraken2-microbe-files` and following the below steps (which require singularity to be installed). 

```
cd build-Kraken2-microbe-files
mkdir data
cp ../data/microbev1.fa data/microbev1.fa
./scripts/run-snakemake.sh
```

This should create a directory named `microbev1`, which contains all the Kraken2Uniq files. You'll need to copy the directory containing the Kraken2Uniq files to the Robinson2023-10x directory.

```
cd Robinson2023-10x
cp -r ../build-Kraken2-microbe-files/microbev1 microbev1
```

Next, you need to download our fork of the SAHMI codebase, which includes small changes to the code to get it running with our database (I suggest to do this in the same directory where the CSI-Microbes-identification codebase was cloned; otherwise you'll need to change any references to SAHMI within Robinson2023-10x/run-SAHMI.smk). 

```
git clone git@github.com:ruppinlab/SAHMI.git
```

Finally, you need to navigate to Robinson2023-10x and run SAHMI using the below commands

```
cd Robinson2023-10x
./scripts/run-SAHMI.sh
```

#### Running INVADE-seq on Robinson2023-10x

Similar to CSI-Microbes, INVADE-seq uses PathSeq and expects the PathSeq dependencies to exist (see above for how to download/build them). Next, you need to download the INVADE-seq GitHub repo (I suggest to do this in the same directory where the CSI-Microbes-identification codebase was cloned; otherwise you'll need to change any references to INVADE-seq within Robinson2023-10x/run-INVADE-seq.smk).

```
git clone git@github.com:FredHutch/Galeano-Nino-Bullman-Intratumoral-Microbiota_2022.git
```

Next, you can navigate to Robinson2023-10x and run INVADE-seq using the below commands

```
cd Robinson2023-10x
./scripts/run-INVADE-seq.sh
```

#### Running CSI-Microbes on Robinson2023-10x

CSI-Microbes uses PathSeq and expects the PathSeq dependencies to exist (see above for how to download/build them). Next, you can navigate to Robinson2023-10x and run CSI-Microbes using the below commands

```
cd Robinson2023-10x
./scripts/run-CSI-Microbes.sh
```

#### Running CSI-Microbes on Robinson2023-plexWell

CSI-Microbes uses PathSeq and expects the PathSeq dependencies to exist (see above for how to download/build them). Next, you can navigate to Robinson2023-SS2 and run CSI-Microbes using the below commands

```
cd Robinson2023-SS2
./scripts/run-snakemake.sh
```

### Reproducing other results in the manuscript

#### Reproducing results from Ben-Moshe2019 10x dataset

In our paper, we analyze a 10x dataset generated by `Ben-Moshe2019`<sup>[REF](#BenMoshe2019)</sup>, which performed 10x 3' v2 scRNA-seq (read length=58bp) on ~3,500 immune cells exposed to _Salmonella SL1344_ strain and ~3,500 control cells.

To reproduce our analysis using CSI-Microbes-identification with PathSeq on Ben-Moshe2019 (which can take ~18 hours), run the below command.

```
./scripts/run-Ben-Moshe2019-PathSeq.sh
```

The expected output of this pipeline is one pathseq.txt file for each cell. A cell is specified by the unique combination of patient, sample and cell barcode. The cell metadata is specified in the `data/units.tsv` file. For Ben-Moshe2019, you can see that there is one patient (Pt0), two samples (GSM3454528, the control cells, and GSM3454529, the exposed cells) and 7,000 cells. Moreover, the reads from each sample were sequenced in four lanes, which are represented in `data/samples.tsv`.

CSI-Microbes-identification groups cells together by their sample, which should come from a single GEM (Gelbeads in EMulsion) well, which can contain multiple sequencing runs or multiple lanes from the same sequencing run. One sample is processed through CellRanger count, the aligned reads are filtered and the unaligned reads are trimmed and cleaned. Next, the cleaned unaligned are mapped to microbial genomes through PathSeq, which produces a score file (`output/PathSeq/{patient}-{sample}/pathseq.txt`) and BAM file (`output/PathSeq/{patient}-{sample}/pathseq.bam`) for all the reads in the sample. Next, the reads in the PathSeq BAM are merged with the CellRanger BAM file to add the cell barcode and UMI tags to the microbe annotations. Finally, we split the sample-level PathSeq BAM files into a cell barcode-level PathSeq BAM file and score it to create a cell-specific pathseq.txt file (`output/PathSeq/{patient}-{sample}-{cell}/pathseq.txt`).

CSI-Microbes-identification using PathSeq, SRPRISM and CAMMiQ share many rules and will re-use data generated by each other but do not run them in parallel.

To run CSI-Microbes-identification using SRPRISM on Ben-Moshe2019, run the below command.

```
./scripts/run-Ben-Moshe2019-SRPRISM.sh
```

The key output file from this analysis is `output/SRPRISM/Pt0/GSM3454529/CB-UMI-count-SL1344.tsv`, which is a table where the rows are cell barcodes and the column is the number of UMIs found.

To run CSI-Microbes-identification using CAMMiQ on Ben-Moshe, run the below command.

```
./scripts/run-Ben-Moshe2019-CAMMiQ.sh
```

The key output file from this analysis is `output/CAMMiQ/read_cnts_genus.txt`, which is a table where the rows are NCBI tax ids and the rows are cells.


#### Reproducing results from Aulicino2018

In our paper, we analyze a Smart-seq2 dataset generated by `Aulicino2018`<sup>[REF](#Aulicino2018)</sup>, which performed Smart-seq2 scRNA-seq (read length=58bp) on 342 monocyte-derived dendritic cells exposed to either the _Salmonella LT2_ or _Salmonella D23580_ strain.

To reproduce our analysis using CSI-Microbes-identification with PathSeq on Aulicino2018 (which can take ~8 hours), run the below command.

```
./scripts/run-Aulicino2018-PathSeq.sh
```

The expected output of this pipeline is one pathseq.txt file for each cell. A cell is specified by the unique combination of patient, sample, plate and cell. The cell metadata is specified in the `data/units.tsv` file. For Aulicino2018, you can see that there is one patient (Pt0), one samples (S0), 4 plates (P1, P2, P3, P4) and 342 cells.

CSI-Microbes-identification groups cells together by their plate. Each plate of cells is processed through STARSolo, the aligned reads are filtered and the unaligned reads are trimmed and cleaned. Next, the cleaned unaligned are mapped to microbial genomes through PathSeq, which produces a score file (`output/PathSeq/{patient}-{sample}-{plate}/pathseq.txt`) and BAM file (`output/PathSeq/{patient}-{sample}-{plate}/pathseq.bam`) for all the reads in the sample. Finally, we split the plate-level PathSeq BAM files into a cell-level PathSeq BAM file and score it to create a cell-specific pathseq.txt file (`output/PathSeq/{patient}-{sample}-{plate}-{cell}/pathseq.txt`).

CSI-Microbes-identification using PathSeq, SRPRISM and CAMMiQ share many rules and will re-use data generated by each other but do not run them in parallel.

To run CSI-Microbes-identification using SRPRISM on Aulicino2018, run the below command.

```
./scripts/run-Aulicino2018-SRPRISM.sh
```

The key output file from this analysis is a .gff file (`output/SRPRISM/Pt0/S0-{plate}-{cell}/{genome}-paired-count.gff` (where genome is either LT2 or D23580)) for each cell, which can be parsed by CSI-Microbes-analysis to determine the number of reads per cell.

To run CSI-Microbes-identification using CAMMiQ on Aulicino2018, run the below command.

```
./scripts/run-Aulicino2018-CAMMiQ.sh
```

The key output file from this analysis is `output/CAMMiQ/read_cnts_genus.txt`, which is a table where the rows are NCBI tax ids and the rows are cells.


#### Reproducing results from Pelka2021

In our paper, we analyze a 10x dataset generated by `Pelka2021`<sup>[REF](#Pelka2021)</sup>, which performed 10x scRNA-seq on a large number of colorectal carcinoma tumors. The raw FASTQ files required to reproduce these results are available via dbGaP, which requires an application.

To reproduce our analysis using CSI-Microbes-identification with PathSeq on Pelka2021, run the below command (once the raw FASTQ files have been obtained)

```
./scripts/run-Pelka2021-PathSeq.sh
```


#### Reproducing results from Zhang2021

In our paper, we analyze a 10x dataset generated by `Zhang2021`<sup>[REF](#Zhang2021)</sup>, which performed 10x scRNA-seq on a large number of colorectal carcinoma tumors.

To reproduce our analysis using CSI-Microbes-identification with PathSeq on Zhang2021, run the below command (once the raw FASTQ files have been obtained)

```
./scripts/run-Zhang2021-PathSeq.sh
```

### Running CSI-Microbes-identification for additional datasets

To run CSI-Microbes-identification using PathSeq on additional datasets, you will need to set up some data files and some instructions for running the pipeline in the form of scripts. Below are the directories and files that are needed to add an example dataset `Example2023` for either 10x or Smart-seq2. Many of these can be copied from `Ben-Moshe2019` for 10x datasets or `Aulicino2018` for Smart-seq2 datasets.

```
Example2023/
-   config/
    -   cluster.json - specifies the cluster requirements for particular rules
    -   PathSeq-config.yaml - specifies parameters and files (such as host genome files)
-   scripts/
  -   run-snakemake.sh - code for running the Snakemake instance that runs local rules and submits jobs to the cluster
-   Snakefile - contains rules for downloading the data and includes .smk files that contain rules that are reused
scripts/
-   run-Example2023-PathSeq.sh - submits Example2023/scripts/run-snakemake.sh to the Biowulf cluster
```

For a 10x dataset, you will also need to include

```
Example2023/
-   data/
  -   samples.tsv - must contain columns named patient, sample and lane (column order is irrelevant and additional columns may be included)
  -   units.tsv - must contain columns named patient, sample and barcode (column order is irrelevant and additional columns may be included); the barcode must match the "CB" tag from the BAM outputted by CellRanger (check to see whether you need to include the -1 at the end); usually the cell-barcode and cell-type annotations are published by the original authors (when they are not, I have successfully requested them via email)
```

For a Smart-seq2 dataset, you will need to include the same files as above but with slightly different specifications

```
Example2023/
- data/
  - samples.tsv - must contain columns named patient and sample (column order is irrelevant and additional columns may be included)
  - units.tsv - must contain columns named patient, sample, plate and cell (column order is irrelevant and additional columns may be included);
```


## FAQs

### How to monitor progress?

This workflow creates one slurm job named `run-snakemake.sh`, which run an instance of snakemake. The jobid of `run-snakemake.sh` is the output of the `./scripts/run-Ben-Moshe2019-PathSeq.sh` (or other submission scripts in the `scripts` directory). For example, if `./scripts/run-Ben-Moshe2019-PathSeq.sh` returns 6586175 then `Ben-Moshe2019/slurm-6586175.out` contains the snakemake output.

Snakemake will first build a DAG of the jobs that need to be run and will output this information in the top of `Ben-Moshe2019/slurm-6586175.out` as shown below.

```
Job counts:
        count   jobs
        2       FastqToBam
        2       PathSeqPipelineSpark
        7000    PathSeqScoreSpark
        2       add_CB_UB_tags_to_PathSeq_BAM
        1       all
        2       cellranger_count
        8       compress_rename_SRA_FASTQ_files
        2       convert_to_fastq
        8       download_FASTQ_from_SRA
        2       filter_aligned_reads
        2       filter_vector_contaminant_reads
        2       get_query_names_for_vector_contaminants
        2       identify_reads_with_vector_contamination
        2       run_fastp
        2       sort_by_query_name
        7000    split_PathSeq_BAM_by_CB_UB
        2       trim_reads
        14041
```

The job counts do not correspond to the number of slurm jobs that will be submitted as many of these jobs will be batched together (group jobs) or run as `localrules` on the slurm job running the snakemake instance. For example, all 7,000 `PathSeqScoreSpark` and `split_PathSeq_BAM_by_CB_UB` jobs, which are very short jobs run for each cell, will be run on the slurm job running the snakemake instance. Snakemake will submit computationally intensive jobs to SLURM as appropriate. To monitor progress, you can check the slurm output file. For example, if `./scripts/run-Ben-Moshe2019-PathSeq.sh` returns 6586175 then you can check `Ben-Moshe2019/slurm-6586175.out` for progress. To obtain the percentage of the jobs completed (which may not correspond to the percentage of time), you can use the following command: `grep % Ben-Moshe2019/slurm-6586175.out`. The workflow is done when the snakemake instance is no longer running. To check whether the workflow completed successfully, you can either check the status of the slurm job (6586175 in this case) using the HPC dashboard or tail the slurm output file. If you get an output like below with 100% of the jobs done, then the run was successful.

```
Finished job 0.
14041 of 14041 steps (100%) done
```

To calculate the number of jobs to be done without submitting them, you can perform a dry-run using Snakemake (which may take a long time) using the below commands. The below commands are also potentially useful for debugging.

```
sinteractive
cd Ben-Moshe2019
module load snakemake
snakemake -n --quiet
```


### What if I don't have access to the ccr partition?

CSI-Microbes-identification is intended to be run on distributed compute farms that use slurm, such as NIH’s Biowulf system. On such systems, there may be multiple queues for jobs. On Biowulf, the relevant queues are called 'ccr' and 'norm'. All biowulf users have access to the 'norm' partition. This workflow assumes additional access to the 'ccr' partition. This assumption is hardcoded in files such as `scripts/run-Ben-Moshe2019.sh` script and `Ben-Moshe2019/config/cluster.json`.

If you do not have access to the ccr partition, you should change `scripts/run-Ben-Moshe2019-PathSeq.sh` from

```
cd Ben-Moshe2019
sbatch \
--time=1-00:00:00 \
--cpus-per-task=4 \
--mem=10g \
--partition=norm,ccr \
scripts/run-snakemake.sh
```

to

```
cd Ben-Moshe2019
sbatch \
--time=1-00:00:00 \
--cpus-per-task=4 \
--mem=10g \
--partition=norm \
scripts/run-snakemake.sh
```

and change the first entry in `Ben-Moshe2019/config/cluster.json` from

```
"__default__" :
{
    "time" : "08:00:00",
    "mem": "16g",
    "nthreads" : 4,
    "gres": 10,
    "partition": "norm,ccr"
}
```
to
```
"__default__" :
{
    "time" : "08:00:00",
    "mem": "16g",
    "nthreads" : 4,
    "gres": 10,
    "partition": "norm"
}
```

### What if I don't have access to Biowulf?

Biowulf is the NIH's linux cluster that uses the slurm workload manager. It should be relatively straightforward to run CSI-Microbes-identification to another linux cluster that uses the slurm workload manager. Users should ensure that there is a module system that includes snakemake (6.0.5)<sup>[REF](#Snakemake)</sup>, sratoolkit (2.10.9)<sup>[REF](#SRAToolkit)</sup>, cellranger (5.0.1)<sup>[REF](#CellRanger)</sup>, samtools (1.11)<sup>[REF](#SAMtools)</sup>, bedtools (2.29.2)<sup>[REF](#BedTools)</sup>, and picard (latest=2.25.0)<sup>[REF](#Picard)</sup>.

Users may also need to change the names of the partitions to the partitions used by their server (see above question for an example). To speed up and effectively parallelize the PathSeq step, it is important to use nodes with at least 200 GB of local storage because we copy the PathSeq files to the local node before running PathSeq. In our experience, running PathSeq on multiple nodes using the same reference files will be much slower because of network latency as well competition for access to the one PathSeq reference files between the nodes (which will also sometimes cause errors).

Currently, the example 10x analyses use an HG38 genome that exists on biowulf. Therefore, the `CellRanger`: `genome_dir` value in `config/PathSeq-config.yaml` will need to updated as well.

### What are the expected output files?

The expected output files from CSI-Microbes-identification are pathseq.txt files, which are output in `output/PathSeq`. For example, the pathseq file for cell barcode TTTCCTCTCCACTGGG-1 from sample GSM3454529 (exposed to _Salmonella_) is located at `output/PathSeq/Pt0-GSM3454529-TTTCCTCTCCACTGGG-1/pathseq.txt`. These output files are used as input to [CSI-Microbes-analysis](https://github.com/ruppinlab/CSI-Microbes-analysis), which computes the differential abundance of microbes across cell-types.


### How to file a problem report?

Send a detailed e-mail to Welles Robinson (wir963@gmail.com). You may wish to record the commands you issued and the responses you received within an instance of the unix command `script` and then send the contents of the script as part of your report. A run of a program or pipeline such as CSI-Microbes-Identification is an experiment. Please report what you observe and why the output you received  differs from what you expected without speculating as to the cause of the discrepancy.

### How to resume a run that did not complete?

If your run does not complete due to a failed job, you can kick off the pipeline again (after making the necessary changes) and Snakemake will resume from the last completed job. If there is an issue with the snakemake instance (ex. the node running the snakemake instance is killed), then it may be necessary to [unlock the directory](https://snakemake.readthedocs.io/en/stable/project_info/faq.html#how-does-snakemake-lock-the-working-directory) before kicking off the pipeline again.


## References

### Software Tools

This pipeline leverages the many open-source tools that are listed below.

<a id="CellRanger"></a> CellRanger [https://github.com/10XGenomics/cellranger](https://github.com/10XGenomics/cellranger)

<a id="Fastp"></a> Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). Fastp: An ultra-fast all-in-one FASTQ preprocessor. Bioinformatics, 34(17), i884–i890. [https://doi.org/10.1093/bioinformatics/bty560](https://doi.org/10.1093/bioinformatics/bty560)

<a id="STAR"></a> Dobin, A., Davis, C. A., Schlesinger, F., Drenkow, J., Zaleski, C., Jha, S., Batut, P., Chaisson, M., & Gingeras, T. R. (2013). STAR: Ultrafast universal RNA-seq aligner. Bioinformatics, 29(1), 15–21. [https://doi.org/10.1093/bioinformatics/bts635](https://doi.org/10.1093/bioinformatics/bts635)

<a id="Gurobi"></a> Gurobi Optimization, LLC (2021). Gurobi Optimizer Reference Manual. [http://www.gurobi.com](http://www.gurobi.com)

<a id="Snakemake"></a> Köster, J., & Rahmann, S. (2012). Snakemake-a scalable bioinformatics workflow engine. Bioinformatics, 28(19), 2520–2522. [https://doi.org/10.1093/bioinformatics/bts480](https://doi.org/10.1093/bioinformatics/bts480)

<a id="SAMtools"></a> Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., & Durbin, R. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25(16), 2078–2079. [https://doi.org/10.1093/bioinformatics/btp352](https://doi.org/10.1093/bioinformatics/btp352)

<a id="SRPRISM"></a> Morgulis, A., & Agarwala, R. (2020). SRPRISM (Single Read Paired Read Indel Substitution Minimizer): an efficient aligner for assemblies with explicit guarantees. GigaScience, 9(4), giaa023. [https://doi.org/10.1093/gigascience/giaa023](https://doi.org/10.1093/gigascience/giaa023)

<a id="Picard"></a> Picard [http://broadinstitute.github.io/picard/](http://broadinstitute.github.io/picard/)

<a id="Pysam"></a> pysam [https://github.com/pysam-developers/pysam](https://github.com/pysam-developers/pysam)

<a id="BedTools"></a> Quinlan, A. R., & Hall, I. M. (2010). BEDTools: A flexible suite of utilities for comparing genomic features. Bioinformatics, 26(6), 841–842. [https://doi.org/10.1093/bioinformatics/btq033](https://doi.org/10.1093/bioinformatics/btq033)

<a id="SRAToolkit"></a> SRAToolkit [https://github.com/ncbi/sra-tools](https://github.com/ncbi/sra-tools)

<a id="PathSeq"></a> Walker, M. A., Pedamallu, C. S., Ojesina, A. I., Bullman, S., Sharpe, T., Whelan, C. W., & Meyerson, M. (2018). GATK PathSeq: a customizable computational tool for the discovery and identification of microbial sequences in libraries from eukaryotic hosts. Bioinformatics (Oxford, England), 34(24), 4287–4289. [https://doi.org/10.1093/bioinformatics/bty501](https://doi.org/10.1093/bioinformatics/bty501)

<a id="CAMMiQ"></a> Zhu, K., Robinson, W., Schaffer, A. A., Xu, J., Ruppin, E., Ergun, F., Ye, Y., & Sahinalp, S. C. (2020). Strain Level Microbial Detection and Quantification with Applications to Single Cell Metagenomics. Abstract to appear in Proceedings of the 25th International Conference on Research in Computational Molecular Biology. [https://www.biorxiv.org/content/10.1101/2020.06.12.149245v2](https://www.biorxiv.org/content/10.1101/2020.06.12.149245v2)

### Publications

This repository processes data from the following publications.

<a id="Aulicino2018"></a> Aulicino, A. et al. Invasive Salmonella exploits divergent immune evasion strategies in infected and bystander dendritic cell subsets. Nat. Commun. 9, 4883 (2018).

<a id="BenMoshe2019"></a> Bossel Ben-Moshe, N. et al. Predicting bacterial infection outcomes using single cell RNA-sequencing analysis of human immune cells. Nat. Commun. 10, 3266 (2019).

<a id="GaleanoNino2022"></a> Galeano Niño, J.L., Wu, H., LaCourse, K.D. et al. Effect of the intratumoral microbiota on spatial and cellular heterogeneity in cancer. Nature 611, 810–817 (2022).

<a id="Ghaddar2022"></a> Ghaddar, B., et al. Tumor microbiome links cellular programs and immunity in pancreatic cancer. Cancer Cell Volume 40, Issue 10 (2022).

<a id="Pelka2021"></a> Pelka, K. et al. Spatially organized multicellular immune hubs in human colorectal cancer. Cell, (2021).

<a id="Zhang2021"></a> Zhang, X. et al. Dissecting esophageal squamous-cell carcinoma ecosystem by single-cell transcriptomic analysis. Nat.Commun. 12, 5291 (2021).
