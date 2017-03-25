import pandas as pd
import numpy as np

Graham_Objects_path = "C:/Users/Christopher/Documents/UNI/Year 4/Project/AGN-code/Graham_Period_Data.txt"
Graham_Objects_df = pd.read_table(Graham_Objects_path, sep = "\t", header=0, index_col = 0)
Graham_Objects_df.columns = ['Object', 'Period']

Graham_Objects_df.sort_values("Period", inplace = True)
Graham_Objects_df.drop_duplicates("Period", inplace = True)

print(len(Graham_Objects_df))
print(Graham_Objects_df.Period)
