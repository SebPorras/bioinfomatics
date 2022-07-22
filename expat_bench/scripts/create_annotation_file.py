import pandas as pd
print (snakemake.wildcards)

align_df = pd.read_csv(snakemake.input.csv)

annotation_cols =[x.strip() for x in open(snakemake.input.annotation_cols).read().splitlines()]

if annotation_cols[0] != 'Entry':

    annotation_cols = ['Entry'] + annotation_cols

print (annotation_cols)

subset_df = align_df[[x for x in annotation_cols if x in align_df.columns]]

subset_df = subset_df.fillna("None")

subset_df.to_csv(snakemake.output.tsv, sep='\t', index=False)
