import pandas as pd
import os 

def combine_site_files(site_files, out_file=None):
    df_list = [pd.read_csv(f) for f in site_files]
    combined_df = pd.concat(df_list)
    if out_file:
        combined_df.to_csv(out_file)
    return combined_df

def main():
    lake_ids = ["120020150", "107072210"]
    doy_files = [f"2_process/out/doy_{i}.csv" for i in lake_ids]
    combine_site_files(doy_files, "2_process/out/combined_doy.csv")

if __name__ == '__main__':
    main()