# CSI-Microbes Identification

This repository contains part of the workflows for reproducing the results from the bioarxiv paper [Identifying the Landscape of Intratumoral Microbes via a Single Cell Transcriptomic Analysis](https://www.biorxiv.org/content/10.1101/2020.05.14.096230v1) by Welles Robinson, Fiorella Schischlik, Michael Gertz, Alejandro Schaffer and Eytan Ruppin.  This repository contains the workflows to identify microbial reads from 10x and Smart-seq2 scRNA-seq datasets. These microbial reads can then be analyzed using the [CSI-Microbes-analysis repository](https://github.com/ruppinlab/CSI-Microbes-analysis). The code in this repository was written by Welles Robinson and alpha-tested by Alejandro Schaffer.

## Set-up

The workflows in this repository are set-up to be run specifically on the NIH Biowulf cluster. This workflow expects that conda has been installed. For instructions on how to install conda, see [conda install documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/install/).

Next, the GitHub repository needs to be cloned as below. The below instructions assume that you have an ssh key associated with your GitHub account. If you do not, you can generate a new ssh key and associate it with your GitHub username by following [these instructions](https://docs.github.com/en/github/authenticating-to-github/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent).

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

CSI-Microbes-identification depends on the following software packages that are loaded via the Biowulf module system: snakemake (6.0.5)<sup>[REF](#Snakemake)</sup>, sratoolkit (2.10.9)<sup>[REF](#SRAToolkit)</sup>, cellranger (5.0.1)<sup>[REF](#CellRanger)</sup>, samtools (1.11)<sup>[REF](#SAMtools)</sup>, bedtools (2.29.2)<sup>[REF](#BedTools)</sup>, picard (latest=2.25.0)<sup>[REF](#Picard)</sup> and GATK (4.1.8.1) <sup>[REF](#PathSeq)</sup>. CSI-Microbes-identification also depends on the following conda software packages from the conda-forge channel: fastp (0.20.1)<sup>[REF](#Fastp)</sup> and pysam (1.16.0)<sup>[REF](#pysam)</sup>. CSI-Microbes-identification for Smart-seq2 also depends on STAR (2.7.6a_patch_2020-11-16)<sup>[REF](#STAR)</sup>, which must be downloaded from GitHub (there is currently no rule to download the executable and the current path to the STAR executable is hardecoded and requires access to the Biowulf directory `/data/Robinson-SB`).

## Database Dependencies

Currently, CSI-Microbes-identification requires the path to PathSeq database files to be hardcoded in the `config/PathSeq-config.yaml`. The code to generate these files is not ready and so currently execution of the pipeline depends on access to the Biowulf directory `/data/Robinson-SB`.

## Running CSI-Microbes-identification

The output of the CSI-Microbes-identification pipeline is microbial abundance files. Currently, we support three distinct approaches for quantifying the read abundance of one or more microbes: PathSeq, CAMMiQ and SRPRISM. PathSeq and CAMMiQ both map reads against a large number of microbial genomes. PathSeq uses a (more computationally expensive) alignment-based approach and reports matches across all levels of the taxonomy while CAMMiQ uses a much faster _k_-mer-based approach and reports matches at the specificied taxonomic level of interest. SRPRISM is a read alignment tool that we use to align reads against a single microbial genome and identify the genome location of reads. For simplicity, we will focus on using PathSeq.

The code in this pipeline is designed to analyze scRNA-seq datasets generated using either Smart-Seq2 (or similar plate-based approaches) or 10x although it should be straightforward to extend this pipeline to analyze scRNA-seq datasets from additional approaches.

### Test Run CSI-Microbes-identification

We have created two small test directories: `test-10x` and `test-SS2` for users to quickly test that they can properly run the code base.

#### Running 10x test dataset

The input fastq files are located at `test-10x/FASTQ/raw/Pt0/` (Yes, they were added to git. No, I do not feel good about it.). The input fastq files contain all reads from 12 cells from Ben-Moshe2019. To run this directory, use the below command, which submits the job running the Snakemake instance to the server.

```
./scripts/run-test-10x.sh
```

The job that runs the Snakemake instance will also submit additional jobs to the server. In total, these jobs should last ~20 minutes (depending on waiting times for job allocation). The key output directories are placed in the directory `output/PathSeq`. This directory should contain 13 directories: one directory for the microbial reads found in the sample (`output/PathSeq/{patient}-{sample}/`) and one directory for each of the 12 cells processed (`output/PathSeq/{patient}-{sample}-{cell}/`). These directories contain the `pathseq.txt` file that contains the microbial reads identified in this cell. To check the results of your run, you can compare against expected results in `expected_output/PathSeq` using the following command (from within the `test-10x` directory)

```
diff -r output/PathSeq/ expected_output/PathSeq/
```

which should yield a list of files that exist only in `output/PathSeq/` but no differences between the files that exist in both directories.

#### Running Smart-seq2 test dataset

### Running CSI-Microbes-identification for Ben-Moshe2019 10x dataset

The gold-standard 10x dataset analyzed in our paper is `Ben-Moshe2019`<sup>[REF](#BenMoshe2019)</sup>, which performed 10x 3' v2 scRNA-seq (read length=58bp) on ~3,500 immune cells exposed to _Salmonella SL1344_ strain and ~3,500 control cells. To run Ben-Moshe2019 (which can take ~18 hours), run the below command.

```
./scripts/run-Ben-Moshe2019.sh
```

The expected output of this pipeline is one pathseq.txt file for each cell. A cell is specified by the unique combination of patient, sample and cell barcode. The cell metadata is specified in the `data/units.tsv` file. For Ben-Moshe2019, you can see that there is one patient (Pt0), two samples (GSM3454528, the control cells, and GSM3454529, the exposed cells) and 7,000 cells. Moreover, the reads from each sample were sequenced in four lanes, which are represented in `data/samples.tsv`.

CSI-Microbes-identification groups cells together by their sample, which should come from a single GEM (Gelbeads in EMulsion) well, which can contain multiple sequencing runs or multiple lanes from the same sequencing run. One sample is processed through CellRanger count, the aligned reads are filtered and the unaligned reads are trimmed and cleaned. Next, the cleaned unaligned are mapped to microbial genomes through PathSeq, which produces a score file (`output/PathSeq/{patient}-{sample}/pathseq.txt`) and BAM file (`output/PathSeq/{patient}-{sample}/pathseq.bam`) for all the reads in the sample. Next, the reads in the PathSeq BAM are merged with the CellRanger BAM file to add the cell barcode and UMI tags to the microbe annotations. Finally, we split the sample-level PathSeq BAM files into a cell barcode-level PathSeq BAM file and score it to create a cell-specific pathseq.txt file (`output/PathSeq/{patient}-{sample}-{cell}/pathseq.txt`).

### Running CSI-Microbes-identification for additional 10x datasets

To run CSI-Microbes-identification on additional 10x datasets, you will need to set up some data files and some instructions for running the pipeline in the form of scripts. Below are the directories and files that are needed to add an example dataset `Robinson2021`.

```
Robinson2021/
- config/
  - cluster.json - specifies the cluster requirements for particular rules
  - PathSeq-config.yaml - specifies parameters and files (such as host genome files)
- data/
  - samples.tsv - must contain columns named patient, sample and lane (column order is irrelevant and additional columns may be included)
  - units.tsv - must contain columns named patient, sample and barcode (column order is irrelevant and additional columns may be included); the barcode must match the "CB" tag from the BAM outputted by CellRanger; usually the cell-barcode and cell-type annotations are published by the original authors (when they are not, I have successfully requested them via email)
- scripts/
  - run-snakemake.sh - code for running the Snakemake instance that runs local rules and submits jobs to the cluster
- Snakefile - contains rules for downloading the data and includes .smk files that contain rules that are reused
scripts/
- run-Robinson2021-microbes.sh - submits Robinson2021/scripts/run-snakemake.sh to the Biowulf cluster
```


## FAQs

### How to monitor progress?

This workflow creates one slurm job named `run-snakemake.sh`, which run an instance of snakemake. The jobid of `run-snakemake.sh` is the output of the `./scripts/run-Ben-Moshe2019.sh` (or other submission scripts in the `scripts` directory). For example, if `./scripts/run-Ben-Moshe2019.sh` returns 6586175 then `Ben-Moshe2019/slurm-6586175.out` contains the snakemake output.

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

The job counts do not correspond to the number of slurm jobs that will be submitted as many of these jobs will be batched together (group jobs) or run as `localrules` on the slurm job running the snakemake instance. For example, all 7,000 `PathSeqScoreSpark` and `split_PathSeq_BAM_by_CB_UB` jobs, which are very short jobs run for each cell, will be run on the slurm job running the snakemake instance. Snakemake will submit computationally intensive jobs to SLURM as appropriate. To monitor progress, you can check the slurm output file. For example, if `./scripts/run-Ben-Moshe2019.sh` returns 6586175 then you can check `Ben-Moshe2019/slurm-6586175.out` for progress. To obtain the percentage of the jobs completed (which may not correspond to the percentage of time), you can use the following command: `grep % Ben-Moshe2019/slurm-6586175.out`. The workflow is done when the snakemake instance is no longer running. To check whether the workflow completed successfully, you can either check the status of the slurm job (6586175 in this case) using the HPC dashboard or tail the slurm output file. If you get an output like below with 100% of the jobs done, then the run was successful.

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

This workflow assumes access to the ccr and norm partition. This assumption is hardcoded in files such as `scripts/run-Ben-Moshe2019.sh` script and `Ben-Moshe2019/config/cluster.json`.

If you do not have access to the ccr partition, you should change `scripts/run-Ben-Moshe2019.sh` from

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

### What are the expected output files?

The expected output files from CSI-Microbes-identification are pathseq.txt files, which are output in `output/PathSeq`. For example, the pathseq file for cell barcode TTTCCTCTCCACTGGG-1 from sample GSM3454529 (exposed to _Salmonella_) is located at `output/PathSeq/Pt0-GSM3454529-TTTCCTCTCCACTGGG-1/pathseq.txt`. These output files are used as input to [CSI-Microbes-analysis](https://github.com/ruppinlab/CSI-Microbes-analysis), which computes the differential abundance of microbes across cell-types.

<!-- ## Aulicino2018

There are two key sets of output files generated by this pipeline for Aulicino2018. The first set of files are the read count files resulting from directly mapping the unaligned reads (the reads that do not map to the human genome) against the _Salmonella_ strain genome used for infection (_LT2_ or _D23580_) using SRPRISM.

The second set of files are the read counts files resulting from running PathSeq with a large database of microbial genomes on the unaligned reads (the reads that do not map to the human genome).


## Ben-Moshe2019

[Ben-Moshe2019](https://www.nature.com/articles/s41467-019-11257-y) performed 10x 3' v2 scRNA-seq (read length=58bp) on ~3,500 immune cells exposed to _Salmonella SL1344_ strain and ~3,500 control cells.

There are two key sets of output files generated by this pipeline for Ben-Moshe2019. The first set of files are the read count files resulting from directly mapping the unaligned reads (the reads that do not map to the human genome) against the _Salmonella SL1344_ strain genome used for infection.

The second set of files are the read counts files resulting from running PathSeq with a large database of microbial genomes on the unaligned reads (the reads that do not map to the human genome).

Currently, the pipeline is set-up to generate the second set of files. To generate these files, from the top-level directory (CSI-Microbes-analysis), you should run the below command, which uses snakemake to submit the necessary jobs to the Biowulf cluster

```
./scripts/run-Ben-Moshe2019.sh
```

## Lee2020

[Lee2020](https://www.nature.com/articles/s41588-020-0636-z) performed 10x 3' v2 scRNA-seq (read length=98 bp) on ~90,000 cells from patients with colorectal carcinomas (although we only analyze the 6 patients where raw reads are available).

## Paulson2018

[Paulson2018](https://www.nature.com/articles/s41467-018-06300-3) performed scRNA-seq (10x 3' scRNA-seq for patient 2586-4 (read length=98bp) and 5' scRNA-seq for patient 9245-3 (read length=91bp)) on two patients with Merkel cell carcinoma. -->

### How to file a problem report?

Send a detailed e-mail to Welles Robinson (Wells.Robinson@nih.gov). You may wish to record the commands you issued and the responses you received within an instance of the unix script and then send the contents of the script as part of your report. A run of a program or pipeline such as CSI-Microbes-Identification is an experiment. Please report what you observe and why the output you received  differs from what you expected without speculating as to the cause of the discrepancy.

### How to resume a run that did not complete?

If your run does not complete due to a failed job, you can kick off the pipeline again (after making the necessary changes) and Snakemake will resume from the last completed job. If there is an issue with the snakemake instance (ex. the node running the snakemake instance is killed), then it may be necessary to [unlock the directory](https://snakemake.readthedocs.io/en/stable/project_info/faq.html#how-does-snakemake-lock-the-working-directory) before kicking off the pipeline again.


## References

### Software Tools

This pipeline leverages the many open-source tools that are listed below.

<a id="CellRanger"></a> CellRanger [https://github.com/10XGenomics/cellranger](https://github.com/10XGenomics/cellranger)

<a id="Fastp"></a> Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). Fastp: An ultra-fast all-in-one FASTQ preprocessor. Bioinformatics, 34(17), i884–i890. [https://doi.org/10.1093/bioinformatics/bty560](https://doi.org/10.1093/bioinformatics/bty560)

<a id="STAR"></a> Dobin, A., Davis, C. A., Schlesinger, F., Drenkow, J., Zaleski, C., Jha, S., Batut, P., Chaisson, M., & Gingeras, T. R. (2013). STAR: Ultrafast universal RNA-seq aligner. Bioinformatics, 29(1), 15–21. [https://doi.org/10.1093/bioinformatics/bts635](https://doi.org/10.1093/bioinformatics/bts635)

<a id="Snakemake"></a> Köster, J., & Rahmann, S. (2012). Snakemake-a scalable bioinformatics workflow engine. Bioinformatics, 28(19), 2520–2522. [https://doi.org/10.1093/bioinformatics/bts480](https://doi.org/10.1093/bioinformatics/bts480)

<a id="SAMtools"></a> Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., & Durbin, R. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25(16), 2078–2079. [https://doi.org/10.1093/bioinformatics/btp352](https://doi.org/10.1093/bioinformatics/btp352)

<a id="SRPRISM"></a> Morgulis, A., & Agarwala, R. (2020). SRPRISM (Single Read Paired Read Indel Substitution Minimizer): an efficient aligner for assemblies with explicit guarantees. GigaScience, 9(4), 1–12. [https://doi.org/10.1093/gigascience/giaa023](https://doi.org/10.1093/gigascience/giaa023)

<a id="Picard"></a> Picard [http://broadinstitute.github.io/picard/](http://broadinstitute.github.io/picard/)

<a id="Pysam"></a> pysam [https://github.com/pysam-developers/pysam](https://github.com/pysam-developers/pysam)

<a id="BedTools"></a> Quinlan, A. R., & Hall, I. M. (2010). BEDTools: A flexible suite of utilities for comparing genomic features. Bioinformatics, 26(6), 841–842. [https://doi.org/10.1093/bioinformatics/btq033](https://doi.org/10.1093/bioinformatics/btq033)

<a id="SRAToolkit"></a> SRAToolkit [https://github.com/ncbi/sra-tools](https://github.com/ncbi/sra-tools)

<a id="PathSeq"></a> Walker, M. A., Pedamallu, C. S., Ojesina, A. I., Bullman, S., Sharpe, T., Whelan, C. W., & Meyerson, M. (2018). GATK PathSeq: a customizable computational tool for the discovery and identification of microbial sequences in libraries from eukaryotic hosts. Bioinformatics (Oxford, England), 34(24), 4287–4289. [https://doi.org/10.1093/bioinformatics/bty501](https://doi.org/10.1093/bioinformatics/bty501)

<a id="CAMMiQ"></a> Zhu, K., Robinson, W., Schaffer, A. A., Xu, J., Ruppin, E., Ergun, F., Ye, Y., & Sahinalp, C. (2020). Strain Level Microbial Detection and Quantification with Applications to Single Cell Metagenomics. BioRxiv, 2020.06.12.149245. [https://www.biorxiv.org/content/10.1101/2020.06.12.149245v2](https://www.biorxiv.org/content/10.1101/2020.06.12.149245v2)

### Publications

This repository processes data from the following publications.

<a id="Aulicino2018"></a>Aulicino, A., Rue-Albrecht, K. C., Preciado-Llanes, L., Napolitani, G., Ashley, N., Cribbs, A., Koth, J., Lagerholm, B. C., Ambrose, T., Gordon, M. A., Sims, D., & Simmons, A. (2018). Invasive Salmonella exploits divergent immune evasion strategies in infected and bystander dendritic cell subsets. Nature Communications, 9(1). [https://doi.org/10.1038/s41467-018-07329-0](https://doi.org/10.1038/s41467-018-07329-0)

<a id="BenMoshe2019"></a> Bossel Ben-Moshe, N., Hen-Avivi, S., Levitin, N., Yehezkel, D., Oosting, M., Joosten, L. A. B., Netea, M. G., & Avraham, R. (2019). Predicting bacterial infection outcomes using single cell RNA-sequencing analysis of human immune cells. Nature Communications, 10(1), 1–16. [https://doi.org/10.1038/s41467-019-11257-y](https://doi.org/10.1038/s41467-019-11257-y)
