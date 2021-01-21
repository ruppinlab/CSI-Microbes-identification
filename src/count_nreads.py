import pandas as pd
import pysam

output = []
print(snakemake.params)
for input in snakemake.input:
    output.append(int(pysam.view("-c", input)))
df = pd.DataFrame(output, index=snakemake.params)
df.to_csv(snakemake.output[0], header=False, sep="\t")
