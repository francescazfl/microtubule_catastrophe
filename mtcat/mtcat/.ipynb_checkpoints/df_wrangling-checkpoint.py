import os
import pandas as pd



def load_mtcat(
    csv_name="gardner_mt_catastrophe_only_tubulin.csv", numb_head=9, if_display=True
):
    """Load the gardner mt catastrophe only tubulin dataset"""
    data_path = "../data/"
    df = pd.read_csv(os.path.join(data_path, csv_name), header=numb_head)
    if if_display:
        display(df)
    return df



def df_cleanup(df, if_display=True):
    """Clean up the data to tidy form"""
    df.columns = df.columns.str.rstrip(" uM")
    df = pd.melt(df, value_name="time to catastrophe (s)").dropna()
    df["variable"] = df["variable"].astype(int)
    df = (
        df.rename(columns={"variable": "concentration (μM)"})
        .sort_values(by=["concentration (μM)"])
        .reset_index(drop=True)
    )
    if if_display:
        display(df)
    return df



def get1conc(df, conc):
    """Get one chosen concentration"""
    return df.loc[df["concentration (μM)"] == conc]



