import pandas as pd

df = pd.read_csv("expression.tsv", sep="\t", index_col=0)
print(df.shape)  # Should be (genes, samples)
print(df.head())  # Should match the image pattern
