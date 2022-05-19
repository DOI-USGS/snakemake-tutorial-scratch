import zipfile
import os

def unzip_file(file_path, out_dir):
    with zipfile.ZipFile(file_path, 'r') as zip_ref:
        zip_ref.extractall(out_dir)

def main(sb_file, out_dir):
    if not os.path.exists(out_dir):
	    os.makedirs(out_dir)
    unzip_file(sb_file, out_dir)

if __name__ == '__main__':
    sb_file = snakemake.input['sb_file']
    out_dir = snakemake.params['out_dir']
    main(sb_file, out_dir)
