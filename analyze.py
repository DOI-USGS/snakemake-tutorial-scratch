import os
import zipfile
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sciencebasepy import SbSession

def sb_get(item_id, sb_data_file=None, destination_dir='.'):
    """
    Download an item from ScienceBase by ID. 
    Written by @AndyMcAliley

    :param item_id: ScienceBase ID of item to download
    :param sb_data_file: Name of file to download. If None, download all files. (Default value = None)
    :param destination_dir: directory to save to (Default value = '.')
    :returns: ScienceBase JSON response

    """
    sb_session = SbSession()
    if not (sb_session.ping()['result'] == 'OK'):
        raise ConnectionError('ScienceBase ping unsuccessful')
    # If data are public, no need to log in
    # if not sb_session.is_logged_in():
        # sb_session.login()

    item_json = sb_session.get_item(item_id)
    if sb_data_file is None:
        response = sb_session.get_item_files(item_json, destination=destination_dir)
    else:
        all_file_info = sb_session.get_item_file_info(item_json)
        file_found = False
        for file_info in all_file_info:
            if file_info['name'] == sb_data_file:
                file_found = True
                response = sb_session.download_file(file_info['url'], file_info['name'], destination_dir)
        if not file_found:
            raise FileNotFoundError(f'{sb_data_file} not found on ScienceBase')

    return response

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


def combine_site_files(site_files, out_file=None):
    df_list = [pd.read_csv(f) for f in site_files]
    combined_df = pd.concat(df_list)
    if out_file:
        combined_df.to_csv(out_file)
    return combined_df


def plot_doy_means(combined_doy_means, out_file, depths = (0, 1, 2, 5)):
    df_combined = pd.read_csv(combined_doy_means).set_index(["doy", "site_id"])
    depth_cols = [f"temp_{d}" for d in depths]
    df_sel_depths = df_combined[depth_cols]
    df_long = df_sel_depths.reset_index().melt(var_name="depth",
                                               value_name = "temperature (C)",
                                               id_vars=["doy", "site_id"])
    sns.relplot(data=df_long, x="doy", y="temperature (C)", col="site_id",
                hue="depth", col_wrap=2, kind='line')
    plt.tight_layout()
    plt.savefig(out_file)
    return out_file


if __name__ == "__main__":
    sb_item = "5e5d0bb9e4b01d50924f2b36"
    sb_file = "pgdl_predictions_04_N45.50-48.00_W92.00-93.00.zip"
    out_dir = "out_data/"
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    sb_get(sb_item, sb_file, out_dir)
    file_path = os.path.join(out_dir, sb_file)

    with zipfile.ZipFile(file_path, 'r') as zip_ref:
        zip_ref.extractall(out_dir)

    lake_ids = ["120020150", "107072210"]
    for i in lake_ids:
        pred_csv_file = os.path.join("out_data/tmp",
                                     f"pgdl_nhdhr_{i}_temperatures.csv")
        df = read_pred_csv(pred_csv_file)
        calc_doy_means(df, i, f"doy_{i}.csv")
    doy_files = [f"doy_{i}.csv" for i in lake_ids]
    combine_site_files(doy_files, "combined_doy.csv")
    plot_doy_means("combined_doy.csv", "doy_plot.png")


