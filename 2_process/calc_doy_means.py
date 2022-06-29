import pandas as pd
import os

def read_pred_csv(csv_path):
    df = pd.read_csv(csv_path,
                     parse_dates = ['date'],
                     infer_datetime_format = True)
    return df

def calc_doy_means(df, site_id, out_file=None):
    doy_means = df.groupby([df.date.dt.month, df.date.dt.day]).mean()
    doy_means = doy_means.reset_index(drop=True)
    doy_means.index.name = "doy"
    doy_means["site_id"] = site_id
    if out_file:
        doy_means.to_csv(out_file)
    return doy_means

def main(out_file, in_file, lake_id):
    out_dir = os.path.dirname(out_file)
    if not os.path.exists(out_dir):
	       os.makedirs(out_dir)
    df = read_pred_csv(in_file)
    calc_doy_means(df, lake_id, f"2_process/out/doy_{lake_id}.csv")

if __name__ == '__main__':
    out_file = snakemake.output['out_file']
    in_file = snakemake.input['in_file']
    lake_id = snakemake.wildcards['lake_id']
    main(out_file, in_file, lake_id)
