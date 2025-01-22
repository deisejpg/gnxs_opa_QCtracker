import pandas as pd
from pathlib import Path

# Setting up base directory and the input directory
base_dir = Path('./').resolve()
input_dir = Path('input_R2024')

# List of all Excel files in the input directory
qctrackers_path = list(input_dir.glob('*.xlsx'))

#Extracting IDs from file names
file_ids = []
for file_path in qctrackers_path:
	file_name = file_path.name
	file_id = file_name.split('_')[0]
	file_ids.append(file_id)

all_data = pd.DataFrame()

for file_path, file_id in zip(qctrackers_path, file_ids):
	df = pd.read_excel(file_path, sheet_name = 'QC Tracker')

	# Add a new column for the file ID
	df['Requisition_ID'] = file_id

	# Concatenate this DataFrame with the main one
	all_data = pd.concat([all_data, df], ignore_index=True)


# Reorganize column names
all_data_columns_names = all_data.columns.tolist()
reordered_columns = ['Requisition_ID']
reordered_columns.append(all_data_columns_names)

reordered_columns = ['Requisition_ID','Specimen ID',
		'QCTracker best PPG', 'DNA Quality Score',       
       'RNA Quality Score', 'SNV Summary', 'CNV Summary', 'Fusion Summary',
       'AssayID', 'DNA Mapped Reads (>=500,000)',
       'DNA Mean Read Length (>=50)', 'Uniformity of base coverage (>=90%)',
       'Percent end-to-end Reads (>=50%)', 'MAPD Value (<=0.5)',
       'MAPD Outcome', 'RNA Mapped Reads (>=20,000)',
       'RNA Mean Read Length (>=35)', 'RNA Expression Ctrls Value (>=5)',
       'RNA Expression Ctrls Outcome', 'SNV Indel ID', 'SNV Locus',
       'SNV Allele Ref', 'SNV Allele Alt', 'SNV Type', 'SNV Cosmic',
       'SNV Variant Class', 'SNV Gene Class', 'SNV Allele Frequency',
       'SNV AA Change', 'CNV ID', 'CNV Locus', 'CNV Gene', 'CNV Variant Class',
       'CNV Copy Number', 'CNV Ratio', 'CNV Gene Class', 'CNV Call',
       'Fusion Variant ID', 'Fusion Locus', 'Fusion Variant Class',
       'Fusion Gene Class', 'Fusion Read Counts', 'Fusion Type',
       'Fusion Driver Gene', 'Fusion Mol Cov', 'Fusion Imbalance Score',
       'Fusion Predicted Break-point Range', 'Fusion Read Counts Per Million',
       'Percent Loading', 'Key Signal', 'Raw Read Accuracy',
       'DNA Mean AQ20 Read Length', 'CF-1 Avg ReadsperLane',
       'CF-1 BaseCallAccuracy', 'CF-1 AQ20 MeanReadLength',
       'Total Assigned Amplicon Reads', 'Uniformity of Amplicon Coverage (%)',
       'Amplicons With 1+ Read (%)', 'Amplicons With 100+ Reads (%)',
       'Amplicons With 500+ Reads (%)', 'Amplicons Without Strand Bias (%)',
       'Amplicons reading end-to-end', 'Avg Base Cov Depth',
       'Target Base Cov at 1x (%)', 'Target Base Cov at 100x (%)',
       'Target Base Cov at 500x (%)', 'Target bases Without Strand Bias (%)',
       'DNA_bamID', 'RNA_bamID', 'DNA_barcode', 'RNA_barcode', 'GNXS',
       'Analysis Date', 'Completion Date', 'Genexus Software Version',
       'Assay Version', 'PPG1', 'PPG2', 'PPG3', 'PPG4', 'PPG5', 'PPG6']

all_data = all_data[reordered_columns]

# After the loop, all_data contains all concatenated data with file IDS
all_data.to_excel('2024_concatenated_QCTrackers.xlsx', index=False)