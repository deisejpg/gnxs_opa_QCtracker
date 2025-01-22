from pathlib import Path
import subprocess
import sys
import csv

'''
  A few functions to download specific files from an Azure blob  
  keeping the folder structure from the origin to handle uninformative  
  file names.  

  Input files:
  list_of_folders_to_download.csv>  - list of paths to folders of interest  
              in the format "GNXS-3/reports/AssayDev_sample1_OncominePrecision_GX5_DNA_FUSION_w2.6.0_J01SP_J01SP"  

  <path_destination_root> - path where the folders will be copied to following  
              the format: "C:\User\Bioinformatics\7_Data_orders\om_po1\reports"
  
  Optional:
  [id_csv] - optional, provide a csv file with Case ID, sample ID and  
              Screened sample ID to update the folder name in the format:  
              CaseID	SampleID	ScreenedSampleID
              100312	495019PB	495019PB


  *> download_rename.log - append screen output in a log file for good  
              record keeping
'''


def download_files(file_type, azcopy_path, base_url, destination_root, folders):
    """ 
        Downloads specified file type from Azure blob storage into corresponding folders.
    
        Args:
            file_type (str): The type of file to download (e.g., 'Snvindel.tsv', 'cnv.tsv').
            azcopy_path (str): Path to the AzCopy executable.
            base_url (str): Base URL of the Azure blob storage container.
            destination_root (str): Local root directory to store downloaded files.
            folders (list): List of folder names to process (each corresponds to a subfolder in Azure).
        Returns:
            None
    """

    destination_root = Path(destination_root)

    for folder in folders:
        # Construct the full Azure URL for the specific file in the folder
        source_url = f"{base_url}/{folder}/{file_type}"

        # Construct destination folder path
        destination_path = destination_root / Path(folder).name

        # Ensure the destination folder exists
        destination_path.mkdir(parents = True, exist_ok=True)

        # Construct the AzCopy command
        command = [
            azcopy_path, "copy", source_url, str(destination_path)
        ]

        try:
            # Run the AzCopy command as a subprocess
            print(f"Downloading {file_type} from {source_url} to {destination_path}.")
            subprocess.run(command, check=True)
            print(f"Successfully downloaded {file_type} to {destination_path}'")
        except subprocess.CalledProcessError as e:
            print(f"Error downloading {file_type} from {source_url}: {e}")

def read_folders_from_file(file_path):
    """
        Reads folder names from CSV file.
        Args:
            file_path (str): Path to the CSV file containing folder names.
        Returns:
            list: list of folder names.
    """
    file_path = Path(file_path)
    folders = []
    with file_path.open(mode='r') as csv_file:
        csv_reader = csv.reader(csv_file)
        for row in csv_reader:
            if row: # Avoid emtpy lines
                folders.append(row[0])
    return folders


def rename_folders(folder_root, id_csv):
    """
        Renames folders based on IDs provided in a CSV file.
        Args:
            folder_root (str): Path to the root directory containing folders to rename.
            id_csv (str) : Path to the CSV file with caseID, SampleID, and screened Sample ID.
        Returns:
        None
    """
    # convert folder_root to a Path object
    folder_root = Path(folder_root)
    #Read ID mappings from the csv file
    id_mapping = {}
    with open(id_csv, mode='r') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        for row in csv_reader:
            sample_id = row['ScreenedSampleID']
            if sample_id not in id_mapping:
                id_mapping[sample_id] = []
            id_mapping[sample_id].append({ 'CaseID': row['CaseID'], 
                                        'SampleID': row['SampleID']
                                    })
    # Rename folders
    for folder in folder_root.iterdir():
        if folder.is_dir():
            folder_name = folder.name
            parts = folder_name.split('_')
            if len(parts) >1:
                sample_id = parts[1]
                if sample_id in id_mapping:
                    # handles multiple mpping for the same sampleID
                    for mapping in id_mapping[sample_id]:
                        case_id = mapping['CaseID']
                        screened_sample_id = mapping['SampleID']
                    # create the new folder name
                    new_folder_name = f"Company_{case_id}_{sample_id}_{screened_sample_id}_" + '_'.join(parts[2:])
                    new_folder_path = folder_root / new_folder_name
                    # update the name of the folder
                    folder.rename(new_folder_path)
                    print(f"Renamed {folder_name} to {new_folder_name}")

if __name__ == "__main__":
    if len(sys.argv) <3:
        print("Usage: python azure_tsv_download.py <list_of_folders_to_download.csv> <path_destination_root> [id_csv] *> download_rename.log")
        sys.exit(1)
    # Parse command-line arguments
    folders_file = sys.argv[1]
    destination_root = sys.argv[2]
    id_csv = sys.argv[3] if len(sys.argv) > 3 else None
    # argv[1] example of the first line of the file "GNXS-3/reports/AssayDev_sampleX_OncominePrecision_GX5_DNA_FUSION_w2.6.0_J01SP_J01SP"
    # argv[2] example 'L:\MUTATION Screened Data\Bioinformatics\7_Product_Data_orders\12345_om'
    # argv[3] example of the header of the file "CaseID,SampleID,ScreenedSampleID"

    # Define other parameters needed
    azcopy_path_exe = r'C:/Users/dgoncalves/Downloads/azcopy_windows_amd64_10.27.1/azcopy.exe'
    base_url = "https://name.blob.core.net/tfs-data"

    # Read the list of folders from the csv file
    folders = read_folders_from_file(folders_file)

    # Download Snvindel.tsv files
    download_files("Snvindel.tsv", azcopy_path_exe, base_url, destination_root, folders)

    # Download cnv.tsv files 
    download_files("Cnv.tsv", azcopy_path_exe, base_url, destination_root, folders)

    # Download Fusion.tsv files
    download_files("Fusion.tsv", azcopy_path_exe, base_url, destination_root, folders)
    
    # Rename folders if id_csv is provided
    if id_csv:
        rename_folders(destination_root, id_csv)
