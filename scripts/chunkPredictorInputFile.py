import sys
import pandas as pd
import numpy as np
import os

def chunk_predictor_input_file(xls, chunksize):
    df = pd.read_excel(xls)
    rowcount = len(df)
    chunks = int(np.ceil(rowcount/chunksize))
    dfs = np.array_split(df, chunks)
    
    chunk = 0
    for df in dfs:
        chunk += 1
        tmpfn = f"predict_input{chunk}.xlsx"
        df.to_excel(tmpfn, index=False)
        os.rename("predict_input" + str(chunk) + ".xlsx", "predict_input" + str(chunk) + ".xls")
        
if __name__ == "__main__":
    xls = sys.argv[1]
    chunksize = int(sys.argv[2])
    chunk_predictor_input_file(xls, chunksize)