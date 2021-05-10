# include main Snakefile
include: "Snakefile"

# include rules
include: "../pathogen-discovery-rules/rules/CAMMiQ-10x.smk"

rule:
    input:
        expand(CAMMIQ_COUNT_FILE, tax_level="genus")
