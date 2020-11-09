import pandas as pd
import pysam

output = []

for input in snakemake.input:
    output.append(pysam.view("-c", input))

df = pd.DataFrame(output, index=snakemake.params["cells"])
df.to_csv(snakemake.output[0], header=False, sep="\t")
