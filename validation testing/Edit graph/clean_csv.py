import pandas as pd
# Remove first edit record with edit ID 0
df = pd.read_csv("genomic_edits.csv")
df = df.drop(0)
df.to_csv("genomic_edits.csv", index=False)

