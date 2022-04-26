import zipfile
import os

def unzip_file(file_path, out_dir, sb_file):
    file_path = os.path.join(out_dir, sb_file)
    with zipfile.ZipFile(file_path, 'r') as zip_ref:
        zip_ref.extractall(out_dir)
        
def main():
    sb_file = "pgdl_predictions_04_N45.50-48.00_W92.00-93.00.zip"
    out_dir = "1_fetch/out/"
    file_path = os.path.join(out_dir, sb_file)
    unzip_file(file_path, out_dir, sb_file)
        
if __name__ == '__main__':
    main()