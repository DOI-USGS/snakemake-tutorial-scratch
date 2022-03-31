rule all:
    """
    Collect main outputs of the workflow
    """
    input:
        'pgdl_predictions_04_N45.50-48.00_W92.00-93.00.zip'
        'pgdl_nhdhr_{lake_id1}_temperature.csv'
        'pgdl_nhdhr_{lake_id2}_temperature.csv'
        'pgdl_nhdhr_{lake_id3}_temperature.csv'
        'doy_{lake_id1}.csv'
        'doy_{lake_id2}.csv'
        'doy_{lake_id3}.csv'
        'combined_doy.csv'
        'doy_plot.png'

rule download_data_sb:
    """
    Retreive data from ScienceBase
    """
    output:
        '/out_data/pgdl_nhdhr_{lake_id1}_temperature.csv'
    shell:
        python analyze.py sb_get
