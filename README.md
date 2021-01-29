# CSI-Microbes Identification

This repository contains the workflows to identify microbial reads from scRNA-seq data.

## Set-up

This workflow is only set-up to be run on the NIH Biowulf cluster. This workflow expects that conda has been installed. For instructions on how to install conda, see [conda install documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/install/).

Next, the GitHub repository needs to be cloned as below. The below instructions assume that you have an ssh key associated with your GitHub account. If you do not, you can generate a new ssh key and associate it with your GitHub username by following [these instructions](https://docs.github.com/en/github/authenticating-to-github/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent).

```
git clone git@github.com:ruppinlab/CSI-Microbes-identification.git
```

The above command may return an error like

```
Permission denied (publickey).
fatal: Could not read from remote repository.
Please make sure you have the correct access rights and the repository exists.
```

In this case, your ssh key may not be associated with your GitHub account. Please follow the instructions above and re-try. The below git submodule commands assume that your ssh has been associated with your GitHub account.

This repository uses [git submodules](https://git-scm.com/book/en/v2/Git-Tools-Submodules), which need to be initialized and fetched.

```
cd CSI-Microbes-identification
git submodule init
git submodule update
```

## Running CSI-Microbes-identification

The output of the CSI-Microbes-identification pipeline is microbial abundance files. Currently, we support three distinct approaches for quantifying the read abundance of one or more microbes: PathSeq, CAMMiQ and SRPRISM. PathSeq and CAMMiQ both map reads against a large number of microbial genomes. PathSeq uses a (more computationally expensive) alignment-based approach and reports matches across all levels of the taxonomy while CAMMiQ uses a much faster _k_-mer-based approach and reports matches at the specificied taxonomic level of interest. SRPRISM is a read alignment tool that we use to align reads against a single microbial genome and identify the genome location of reads. For simplicity, we will focus on using PathSeq.

The code in this pipeline is designed to analyze scRNA-seq datasets generated using either Smart-Seq2 (or similar plate-based approaches) or 10x although it should be straightforward to extend this pipeline to analyze scRNA-seq datasets from additional approaches.

### Running CSI-Microbes-identification for Ben-Moshe2019 10x dataset

The example 10x datasets analyzed in our paper is `Ben-Moshe2019`, which performed 10x 3' v2 scRNA-seq (read length=58bp) on ~3,500 immune cells exposed to _Salmonella SL1344_ strain and ~3,500 control cells. To run Ben-Moshe2019 (which can take ~18 hours), run the below command.

```
./scripts/run-Ben-Moshe2019.sh
```

The expected output of this pipeline is one pathseq.txt file for each cell. A cell is specified by the unique combination of patient, sample and cell barcode. The cell metadata is specified in the `data/units.tsv` file. For Ben-Moshe2019, you can see that there is one patient (Pt0), two samples (GSM3454528, the control cells, and GSM3454529, the exposed cells) and 7,000 cells. However, the reads from each sample were sequenced in four lanes, which are represented in `data/samples.tsv`.

CSI-Microbes-identification groups cells together by their sample, which should come from a single GEM well (which can contain multiple sequencing runs or multiple lanes from the same sequencing run). One sample is processed through CellRanger count, the aligned reads are filtered and the unaligned reads are trimmed and cleaned. Next, the cleaned unaligned are mapped to microbial genomes through PathSeq, which produces a score file (`output/PathSeq/{patient}-{sample}/pathseq.txt`) and BAM file (`output/PathSeq/{patient}-{sample}/pathseq.bam`) for all the reads in the sample. Next, the reads in the PathSeq BAM are merged with the CellRanger BAM file to add the cell barcode and UMI tags to the microbe annotations. Finally, we split the sample-level PathSeq BAM files into a cell barcode-level PathSeq BAM file and score it to create a cell-specific pathseq.txt file (`output/PathSeq/{patient}-{sample}-{cell}/pathseq.txt`).

### Running CSI-Microbes-identification for additional 10x datasets

To run CSI-Microbes-identification on additional 10x datasets, you should create a separate top-level directory (first-author last name publication year). Many files including `config/cluster.json`, which specifies the requirements for each job submitted to the cluster, `config/PathSeq-config.yaml`, which specifies rule parameters and host genome files, `scripts/run-snakemake.sh`, which contains the code for the Snakemake instance that runs localrules and submits jobs to the cluster) and `Snakemake` can likely be copied verbatim from `Ben-Moshe2019`. It is important to generate appropriate `data/units.tsv` and `data/samples.tsv` files. The `data/samples.tsv` file must contain the following columns: patient, sample and lane. The `data/units.tsv` must contain the following columns: patient, sample and barcode. The barcode must be identical to the barcode reported by CellRanger in the CB tag in the output BAM or else the pipeline will report all cells have zero microbial reads. It is likely that additional datasets will require their own rules for downloading the FASTQ files (in my experience, there are even differences between datasets from the same host ex. SRA).


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

The job counts do not correspond to the number of slurm jobs that will be submitted as many of these jobs will be batched together (group jobs) or run as `localrules` on the slurm job running the snakemake instance. For example, all 7,000 `PathSeqScoreSpark` and `split_PathSeq_BAM_by_CB_UB` jobs, which are very short jobs run for each cell, will be run on the slurm job running the snakemake instance. Snakemake will submit computationally intensive jobs to SLURM as appropriate. To monitor progress, you can check the slurm output file. For example, if `./scripts/run-Ben-Moshe2019.sh` returns 6586175 then you can check `Ben-Moshe2019/slurm-6586175.out` for progress. To obtain the percentage of the jobs completed (which may not correspond to the percentage of time), you can use the following command: `grep % Ben-Moshe2019/slurm-6586175.out`.

To calculate the number of jobs to be done without submitting them, you can perform a dry-run using Snakemake (which may take a long time) using the below commands. The below commands are also potentially useful for debugging.

```
sinteractive
cd Ben-Moshe2019
module load snakemake
snakemake -n --quiet
```


### What if I don't have access to the ccr partition?

This workflow assumes access to the ccr and norm partition. This assumption is hardcoded in the `scripts/run-Ben-Moshe2019.sh` script and `Ben-Moshe2019/config/cluster.json`.

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

## Aulicino2018

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

[Paulson2018](https://www.nature.com/articles/s41467-018-06300-3) performed scRNA-seq (10x 3' scRNA-seq for patient 2586-4 (read length=98bp) and 5' scRNA-seq for patient 9245-3 (read length=91bp)) on two patients with Merkel cell carcinoma.
