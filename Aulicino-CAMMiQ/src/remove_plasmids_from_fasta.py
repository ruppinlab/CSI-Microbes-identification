from Bio import SeqIO

seq_records = list(SeqIO.parse(snakemake.input[0], "fasta"))

non_plasmid_seq_records = []

for seq_record in seq_records:
    if "plasmid" not in seq_record.description:
        non_plasmid_seq_records.append(seq_record)

SeqIO.write(non_plasmid_seq_records, snakemake.output[0], "fasta")
