# include main Snakefile
include: "Snakefile"

# include rules
include: "../pathogen-discovery-rules/rules/CAMMiQ-SS2.smk"

rule CAMMiQ:
    input:
        expand(CAMMIQ_COUNT_FILE, tax_level="genus")
