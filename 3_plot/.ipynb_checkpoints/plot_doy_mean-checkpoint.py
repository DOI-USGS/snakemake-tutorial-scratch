import pandas as pd 
import os 
import seaborn as sns
import matplotlib.pyplot as plt

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

def main():
    plot_doy_means("2_process/out/combined_doy.csv", "3_plot/out/doy_plot.png")  
    
if __name__ == '__main__':
    main()