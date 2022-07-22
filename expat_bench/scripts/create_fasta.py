import seqcurate as sc
import pandas as pd

sc.write_to_fasta(pd.read_csv(snakemake.input[0]), snakemake.output[0])



