#!/usr/bin/env python3

'''
	This script reads some of the OPA output files and creates a .csv file.

    Important:
    1. With the current set up you have to be in the directory where the directory
    'reports' is, which is the directory that has all the 'AssayDev_*' sub-directories.
    I may implement the path to reports to be an argument to run the script in the
    future!

    Usage:
    python3 test_fetch_qc_genexus_02032023.py

    To do next:
    fix paths to run scripts in a different directory

'''

# Import needed libraries to deal with paths and to read and fetch data from .tsv files
from pathlib import Path
import pandas as pd 
import numpy as np

# Create an object to indicate where is your base directory
base_dir = Path("reports").resolve()

# Use list() from pathlib to list the paths of all files in 'reports/AssayDev_*/'
# These will be helpful to loop over all the files using the function implemented below
info_csvs = list(base_dir.rglob("Info.csv"))
cov_txts = list(base_dir.rglob("*.stats.cov.txt"))
snv_tsvs = list(base_dir.rglob("Snvindel.tsv"))
cnv_tsvs = list(base_dir.rglob("Cnv.tsv"))
fus_tsvs = list(base_dir.rglob("Fusion.tsv"))

# Function to fetch desired data from Info.csv
def get_csv_info(csv_path):
    csv_info_data = {}
    with open(csv_path, "rt") as csv_in:
        for content in "".join(csv_in.readlines()).split("\n\n"):
            if "Library Name" in content:
                lib_blocks = [content.split("\n")]
                for block in lib_blocks:
                    for line in block:
                        if "Library Type" in line:
                            lib_type = f"{line.split(',')[1]}_barcode"
                        elif "Barcode Id" in line:
                            barcode_id = line.split(",")[1]
                    csv_info_data[lib_type] = barcode_id
            if "MAPD" in content:
                dna_blocks = [content.split("\n")]
                for block in dna_blocks:
                    for line in block:
                        if "MAPD" in line:
                            csv_info_data["MAPD Value (<=0.5)"] = line.split(",")[1]
                            csv_info_data["MAPD Outcome"] = line.split(",")[3]
                        if "Mapped Reads" in line:
                            csv_info_data["DNA Mapped Reads (>=500,000)"] = line.split(",")[1]
                        if "Mean AQ20 Read Length (bp)" in line:
                            csv_info_data["DNA Mean AQ20 Read Length"] = line.split(",")[1]
                        if "Mean Read Length (bp)" in line:
                            csv_info_data["DNA Mean Read Length (>=50)"] = line.split(",")[1]
            if "RNA Expression Ctrls Detected" in content:
                rna_blocks = [content.split("\n")]
                for block in rna_blocks:
                    for line in block:
                        if "Mapped Fusion Reads" in line:
                            csv_info_data["RNA Mapped Reads (>=20,000)"] = line.split(",")[1]
                        if "Mean Read Length (bp)" in line:
                            csv_info_data["RNA Mean Read Length (>=35)"] = line.split(",")[1]
                        if "RNA Expression Ctrls Detected" in line:
                            csv_info_data["RNA Expression Ctrls Value (>=5)"] = line.split(",")[1]
                            csv_info_data["RNA Expression Ctrls Outcome"] = line.split(",")[3]
            if "Average Reads Per Lane" in content:
                cf1_blocks = [content.split("\n")]
                for block in cf1_blocks:
                    for line in block:
                        if "Average Reads Per Lane" in line:
                            csv_info_data["CF-1 Avg ReadsperLane"] = line.split(",")[1]
                        if "Base Call Accuracy" in line:
                            csv_info_data["CF-1 BaseCallAccuracy"] = line.split(",")[1]
                        if "Mean AQ20 Read Length (bp)" in line:
                            csv_info_data["CF-1 AQ20 MeanReadLength"] = line.split(",")[1]
            if "Percent Loading" in content:
                loading_blocks = [content.split("\n")]
                for block in loading_blocks:
                    for line in block:
                        if "Percent Loading" in line:
                            csv_info_data["Percent Loading"] = line.split(",")[1]
                        if "Key Signal" in line:
                            csv_info_data["Key Signal"] = line.split(",")[1]
                        if "Raw Read Accuracy" in line:
                            csv_info_data["Raw Read Accuracy"] = line.split(",")[1]
            for line in content.split("\n"):
                if line.startswith("Name,GNXS-"):
                    csv_info_data["GNXS"] = line.split(',')[1]
                if line.startswith("Completion Date,"):
                    csv_info_data["Completion Date"] = line.split(",")[1].split(" ")[0]
                if line.startswith("Date,"):
                    csv_info_data["Analysis Date"] = line.split(",")[1].split(" ")[0]
                if line.startswith("Sample Name") and not line.endswith("Control Kit"):
                    csv_info_data["Specimen ID"] = line.split(",")[1]
                if line.startswith("Name,Oncomine"):
                    csv_info_data["Assay Version"] = line.split(',')[1]
                if line.startswith("ts_dx_suite_version"):
                    csv_info_data["Genexus Software Version"]  = line.split(',')[1]
    bams = list(base_dir.rglob("*.bam"))
    for file in bams:
        bamfile = file.parts[-1]
        if bamfile.startswith("IonHDdual"):
            if bamfile.split("_rawlib")[0] == csv_info_data['DNA_barcode']:
                csv_info_data["DNA_bamID"] = bamfile
            if bamfile.split("_rawlib")[0] == csv_info_data['RNA_barcode']:
                csv_info_data["RNA_bamID"] = bamfile
    return csv_info_data

# Function to fetch data from *.stats.cov.txt
def get_txt_cov(cov_txt):
    txt_info_data = {}
    with open(cov_txt, "rt") as txt_in:
        for line in "".join(txt_in.readlines()).split("\n"):
            if "Percent end-to-end reads" in line:
                txt_info_data["Percent end-to-end Reads (>=50%)"] = line.strip(" ").split(":")[1].strip().split("%")[0]
            if "Total assigned amplicon reads" in line:
                txt_info_data["Total Assigned Amplicon Reads"] = line.strip(" ").split(":")[1].strip()
            if "Average base coverage depth" in line:
                txt_info_data["Avg Base Cov Depth"] = line.strip(" ").split(":")[1].strip()
            if "Target base coverage at 500x" in line:
                txt_info_data["Target Base Cov at 500x (%)"] = line.strip(" ").split(":")[1].strip().split("%")[0]
            if "Target base coverage at 100x" in line:
                txt_info_data["Target Base Cov at 100x (%)"] = line.strip(" ").split(":")[1].strip().split("%")[0]
            if "Target base coverage at 1x" in line:
                txt_info_data["Target Base Cov at 1x (%)"] = line.strip(" ").split(":")[1].strip().split("%")[0]
            if "Amplicons with at least 500 reads" in line:
                txt_info_data["Amplicons With 500+ Reads (%)"] = line.strip(" ").split(":")[1].strip().split("%")[0]
            if "Amplicons with at least 100 reads" in line:
                txt_info_data["Amplicons With 100+ Reads (%)"] = line.strip(" ").split(":")[1].strip().split("%")[0]
            if "Amplicons with at least 1 read" in line:
                txt_info_data["Amplicons With 1+ Read (%)"] = line.strip(" ").split(":")[1].strip().split("%")[0]
            if "Uniformity of amplicon coverage" in line:
                txt_info_data["Uniformity of Amplicon Coverage (%)"] = line.strip(" ").split(":")[1].strip().split("%")[0]
            if "Uniformity of base coverage" in line:
                txt_info_data["Uniformity of base coverage (>=90%)"] = line.strip(" ").split(":")[1].strip().split("%")[0]
            if "Amplicons with no strand bias" in line:
                txt_info_data["Amplicons Without Strand Bias (%)"] = line.strip(" ").split(":")[1].strip().split("%")[0]
            if "Target bases with no strand bias" in line:
                txt_info_data["Target bases Without Strand Bias (%)"] = line.strip(" ").split(":")[1].strip().split("%")[0]
            if "Amplicons reading end-to-end" in line:
                txt_info_data["Amplicons reading end-to-end"] = line.strip(" ").split(":")[1].strip()
    return txt_info_data

# Function to fetch data from Snvindel.tsv
def get_info_snv(snv_path):
    snv_data = {}
    snv_df = pd.read_csv(snv_path, sep='\t', skip_blank_lines=True)
    # strip blank spaces at the end of column names     
    snv_df.columns = snv_df.columns.to_series().apply(lambda x: x.strip())
    # filtering to get only the rows with calls
    snv_df_call = snv_df[(snv_df.Call.isin(["PRESENT (HETEROZYGOUS)", "PRESENT (HOMOZYGOUS)"]))]
    # removing empty spaces from the names of the columns; i.e., "Amino Acid" > "AminoAcid"
    snv_df_call.columns = snv_df_call.columns.str.replace(' ', '')
    # checking for alternative names of headers from different versions of the analysis software
    if 'AlleleFraction' and 'AAChange' in snv_df_call:
        snv_df_call.rename(columns={'AlleleFraction' : 'AlleleFrequency', 'AAChange' : 'AminoAcidChange'}, inplace = True)
    # copying the data frame into a new data frame to manipulate according to our needs
    snv_df_call_summary_copy = snv_df_call.copy()

    #added the following line for cases when there are no calls
    if not snv_df_call_summary_copy.empty:
        snv_df_call_summary_copy['AminoAcidChange'] = snv_df_call_summary_copy['AminoAcidChange'].map(lambda x: x.lstrip('p.'))
        # creating a new column on the new data frame to summarize the calls
        snv_df_call_summary_copy['SNV_summary_info'] = snv_df_call_summary_copy['Gene']+':'+snv_df_call_summary_copy['AminoAcidChange']+snv_df_call_summary_copy['OncomineVariantClass']+' (' + snv_df_call_summary_copy['AlleleFrequency'].astype(str) + ')'
        snv_df_call_summary_copy['SNV_summary_info'] = snv_df_call_summary_copy['SNV_summary_info'].astype(str).map(lambda x: x.replace('Hotspot', ''))
        snv_df_call_summary_copy['SNV_summary_info'] = snv_df_call_summary_copy['SNV_summary_info'].map(lambda x: x.replace('?', ''))

    # if the snv file has no calls, the data frame we just created will be empty
    if snv_df_call_summary_copy.empty:
        snv_data['SNV Summary'] = 'No Variant Detected'
        snv_data['SNV Indel ID'] = 'na'
        snv_data['SNV Locus'] = 'na'
        snv_data['SNV Allele Ref'] = 'na'
        snv_data['SNV Allele Alt'] = 'na'
        snv_data['SNV Type'] = 'na'
        snv_data['SNV Cosmic'] = 'na'
        snv_data['SNV Variant Class'] = 'na'
        snv_data['SNV Gene Class'] = 'na'
        snv_data['SNV Allele Frequency'] = 'na'
        snv_data['SNV AA Change'] = 'na'

    # if the snv file has calls, transform the data frame into a dictionary to get the information we want to report
    else:
        snv_dict = snv_df_call_summary_copy.to_dict()
        # for each gene call get the following information
        for i in snv_dict['Gene']:
            snv_data['SNV Summary'] = "; ".join(repr(x) for x in snv_dict['SNV_summary_info'].values()).replace("'", "")
            snv_data['SNV Indel ID'] = "; ".join(repr(x) for x in snv_dict['Gene'].values()).replace("'", "")
            snv_data['SNV Locus'] =  "; ".join(repr(x) for x in snv_dict['Locus'].values()).replace(":", "|").replace("'", "")
            snv_data['SNV Allele Ref'] =  "; ".join(repr(x) for x in snv_dict['Ref'].values()).replace("'", "")
            snv_data['SNV Allele Alt'] =  "; ".join(repr(x) for x in snv_dict['Alt'].values()).replace("'", "")
            snv_data['SNV Type'] = "; ".join(repr(x) for x in snv_dict['Type'].values()).replace("'", "")

            snv_data['SNV Cosmic'] = "; ".join(repr(x) for x in snv_dict['VariantID'].values()).replace("'", "")
            snv_data['SNV Variant Class'] = "; ".join(repr(x) for x in snv_dict['OncomineVariantClass'].values()).replace("'", "")
            snv_data['SNV Gene Class'] = ";".join(repr(x) for x in snv_dict['OncomineGeneClass'].values()).replace("'", "")

            snv_data['SNV Allele Frequency'] = "; ".join(repr(x) for x in snv_dict['AlleleFrequency'].values()).replace("'", "")
            snv_data['SNV AA Change'] = "; ".join(repr(x) for x in snv_dict['AminoAcidChange'].values()).replace("'", "")
    return(snv_data)

# Function to fetch data from Cnv.tsv
def get_info_cnv(cnv_path):
    cnv_data = {}
    cnv_df = pd.read_csv(cnv_path, sep='\t', skip_blank_lines=True)
    cnv_df_call = cnv_df[(cnv_df.Call.isin(["PRESENT", "PRESENT (GAIN)", "PRESENT (LOSS)"]))]
    cnv_df_call_summary = cnv_df_call.copy()
    cnv_df_call_summary['CNV_summary_info'] = cnv_df_call_summary['Variant ID'] + ' CNV:' + cnv_df_call_summary['CNV Ratio'].astype(str)
    if cnv_df_call_summary.empty:
        cnv_data['CNV Summary'] = 'No Variant Detected'
        cnv_data['CNV ID'] = 'na'
        cnv_data['CNV Locus'] = 'na'
        cnv_data['CNV Gene'] = 'na'
        cnv_data['CNV Variant Class'] = 'na'
        cnv_data['CNV Copy Number'] = 'na'
        cnv_data['CNV Ratio'] = 'na'
        cnv_data['CNV Gene Class'] = 'na'
        cnv_data['CNV Call'] = 'na'
    else:
        cnv_dict = cnv_df_call_summary.to_dict()
        for i in cnv_dict['Variant ID']:
            cnv_data['CNV Summary'] =  "; ".join(repr(x) for x in cnv_dict['CNV_summary_info'].values()).replace("'", "")
            cnv_data['CNV ID'] = "; ".join(repr(x) for x in cnv_dict['Variant ID'].values()).replace("'", "")
            cnv_data['CNV Locus'] = "; ".join(repr(x) for x in cnv_dict['Locus'].values()).replace("'", "")
            cnv_data['CNV Gene'] = "; ".join(repr(x) for x in cnv_dict['Gene'].values()).replace("'", "")
            cnv_data['CNV Variant Class'] = "; ".join(repr(x) for x in cnv_dict['Oncomine Variant Class'].values()).replace("'", "")
            cnv_data['CNV Copy Number'] = "; ".join(repr(x) for x in cnv_dict['Copy Number'].values()).replace("'", "")
            cnv_data['CNV Ratio'] = "; ".join(repr(x) for x in cnv_dict['CNV Ratio'].values()).replace("'", "")
            cnv_data['CNV Gene Class'] = "; ".join(repr(x) for x in cnv_dict['Oncomine Gene Class'].values()).replace("'", "")
            cnv_data['CNV Call'] = "; ".join(repr(x) for x in cnv_dict['Call'].values()).replace("'", "")
    return(cnv_data)

# Function to fetch data from Fusion.tsv
def get_info_fusion(fus_path):
    fus_data = {}
    fus_df = pd.read_csv(fus_path, sep='\t')
    fus_df.columns = [column.replace(" ", "_") for column in fus_df.columns]
    fus_df_call = fus_df[(fus_df.Call.isin(["PRESENT"]))]
    # filtering only target calls
    fus_df_var = fus_df[fus_df['Oncomine_Variant_Class'].notnull()]
    # creating a new data frame to manipulate the data according to our needs
    fus = fus_df_var.copy()
    fus['FUS_summary_info'] = fus['Variant_ID'] +' '+ fus['Oncomine_Variant_Class']+'; '+ 'Read Counts: ' + fus['Read_Counts'].astype(str) +'; Imbalance Score: ' + fus['Imbalance_Score'].astype(str)
    #fus['FUS_summary_info'] = fus['Variant_ID'] +' '+ fus['Oncomine_Variant_Class']+'; '+ fus['Predicted_Break-point_Range'] +'; Imbalance Score: ' + fus['Imbalance_Score'].astype(str) +'; '+ fus['Oncomine_Driver_Gene'] + '; Read Counts: ' + fus['Read_Counts'].astype(str)
    # adding a column to summarize fusion call
    #for item in fus['Oncomine_Variant_Class']:
        #if item == 'ExpressionImbalance' or item == 'Fusion' or item == 'RNAExonVariant':
        #fus['FUS_summary_info'] = fus['Genes_(Exons)'] +'; '+ fus['Predicted_Break-point_Range'] +'; Imbalance Score: ' + fus['Imbalance_Score'].astype(str) +'; '+ fus['Oncomine_Driver_Gene'] +'; '+ fus['Genes_(Exons)'] + '; Read Counts: ' + fus['Read_Counts'].astype(str)
        #if item == 'Fusion' or item == 'RNAExonVariant':
        #    fus['FUS_summary_info']  = fus['Oncomine_Driver_Gene'] +'; '+ fus['Genes_(Exons)'] + '; Read Counts: ' + fus['Read_Counts'].astype(str)
    if fus.empty:
        fus_data['Fusion Summary'] = 'No Fusion Detected'
        fus_data['Fusion Variant ID'] = 'na'
        fus_data['Fusion Locus'] = 'na'
        fus_data['Fusion Variant Class'] ='na'
        fus_data['Fusion Gene Class'] = 'na'
        #fus_data['Genes_(Exons)'] = 'na'
        fus_data['Fusion Read Counts'] = 'na'
        fus_data['Fusion Type'] = 'na'
        #fus_data['Call'] = 'na'
        #fus_data['Call_Details'] = 'na'
        fus_data['Fusion Driver Gene'] = 'na'
        fus_data['Fusion Mol Cov'] = 'na'
        fus_data['Fusion Imbalance Score'] = 'na'
        #fus_data['Imbalance_P-Value'] = 'na'
        fus_data['Fusion Predicted Break-point Range'] = 'na'
        #fus_data['Ratio_To_Wild_Type'] = 'na'
        #fus_data['Norm_Count_Within_Gene'] = 'na'
        fus_data['Fusion Read Counts Per Million'] = 'na'
    else:
        fus_dict = fus.to_dict()
        for i in fus_dict['Variant_ID']:
            fus_data['Fusion Summary'] = "; ".join(repr(x) for x in fus_dict['FUS_summary_info'].values()).replace("'", "")
            fus_data['Fusion Variant ID'] = "; ".join(repr(x) for x in fus_dict['Variant_ID'].values()).replace("'", "")
            fus_data['Fusion Locus'] = "; ".join(repr(x) for x in fus_dict['Locus'].values()).replace("'", "")
            fus_data['Fusion Variant Class'] = "; ".join(repr(x) for x in fus_dict['Oncomine_Variant_Class'].values()).replace("'", "")
            fus_data['Fusion Gene Class'] = "; ".join(repr(x) for x in fus_dict['Oncomine_Gene_Class'].values()).replace("'", "")
            #fus_data['Genes_(Exons)'] = "; ".join(repr(x) for x in fus_dict['Genes_(Exons)'].values()).replace("'", "")
            fus_data['Fusion Read Counts'] = "; ".join(repr(x) for x in fus_dict['Read_Counts'].values()).replace("'", "")
            fus_data['Fusion Type'] = "; ".join(repr(x) for x in fus_dict['Type'].values()).replace("'", "")
            #fus_data['Call'] = "; ".join(repr(x) for x in fus_dict['Call'].values()).replace("'", "")
            #fus_data['Call_Details'] = "; ".join(repr(x) for x in fus_dict['Call_Details'].values()).replace("'", "")
            fus_data['Fusion Driver Gene'] = "; ".join(repr(x) for x in fus_dict['Oncomine_Driver_Gene'].values()).replace("'", "")
            fus_data['Fusion Mol Cov'] = "; ".join(repr(x) for x in fus_dict['Mol_Cov._Mutant'].values()).replace("'", "")
            fus_data['Fusion Imbalance Score'] = "; ".join(repr(x) for x in fus_dict['Imbalance_Score'].values()).replace("'", "")
            #fus_data['Imbalance_P-Value'] = "; ".join(repr(x) for x in fus_dict['Imbalance_P-Value'].values()).replace("'", "")
            fus_data['Fusion Predicted Break-point Range'] = "; ".join(repr(x) for x in fus_dict['Predicted_Break-point_Range'].values()).replace("'", "")
            #fus_data['Ratio_To_Wild_Type'] = "; ".join(repr(x) for x in fus_dict['Ratio_To_Wild_Type'].values()).replace("'", "")
            #fus_data['Norm_Count_Within_Gene'] = "; ".join(repr(x) for x in fus_dict['Norm_Count_Within_Gene'].values()).replace("'", "")
            fus_data['Fusion Read Counts Per Million'] = "; ".join(repr(x) for x in fus_dict['Read_Counts_Per_Million'].values()).replace("'", "")
    return(fus_data)

# create a dictionary to store data from all Info.csv files from all AssayDev_* folders present in 'reports/'
assay_csv_data = {}
# loop over all Info.csv files in the paths present in the object 'info_csvs'
for csv in info_csvs:
    # use the get_csv_info() to fetch the desired information
    csv_data = get_csv_info(csv)
    # get name of folder ('AssayDev_*') to serve as key of our dictionary
    assay_id = csv.parts[-2]
    # add data form Info.csv file to the dictionary
    assay_csv_data[assay_id] = csv_data

# follows the same logic presented for Info.csv to fetch data of coverage, snv, cnv, and fusion
assay_txt_data = {}
for txt in cov_txts:
    txt_data = get_txt_cov(txt)
    assay_id = txt.parts[-2]
    assay_txt_data[assay_id] = txt_data

assay_snv_data = {}
for snv in snv_tsvs:
    snv_data = get_info_snv(snv)
    assay_id = snv.parts[-2]
    assay_snv_data[assay_id] = snv_data

assay_cnv_data = {}
for cnv in cnv_tsvs:
    cnv_data = get_info_cnv(cnv)
    assay_id = cnv.parts[-2]
    assay_cnv_data[assay_id] = cnv_data

assay_fus_data = {}
for fus in fus_tsvs:
    fus_data = get_info_fusion(fus)
    assay_id = fus.parts[-2]
    assay_fus_data[assay_id] = fus_data

# transform dictionaries into pandas data frames to save the output file as .xlsx
csv_df = pd.DataFrame.from_dict(assay_csv_data, orient='index')
txt_df = pd.DataFrame.from_dict(assay_txt_data, orient='index')
snv_df = pd.DataFrame.from_dict(assay_snv_data, orient='index')
cnv_df = pd.DataFrame.from_dict(assay_cnv_data, orient='index')
fus_df = pd.DataFrame.from_dict(assay_fus_data, orient='index')

csv_txt_snv_cnv_fus_df = pd.concat([csv_df, txt_df, snv_df, cnv_df, fus_df], axis=1)
full_df = csv_txt_snv_cnv_fus_df.rename_axis('AssayID')
full_df = full_df.reset_index()
full_df = full_df.replace('N/A','0')


# convert columns type to assess DNA and RNA quality
df_assess_quality = full_df.astype({'DNA Mapped Reads (>=500,000)':'float', 'DNA Mean Read Length (>=50)':'int','Uniformity of base coverage (>=90%)':'float',
                                'Percent end-to-end Reads (>=50%)':'float','MAPD Value (<=0.5)':'float',
                                'RNA Mapped Reads (>=20,000)':'float','RNA Mean Read Length (>=35)':'int','RNA Expression Ctrls Value (>=5)':'int'
                                })

# Information to score DNA:
# DNA Mapped Reads (>=500,000)
# DNA Mean Read Length (>=50)
# Uniformity of base coverage (>=90%)
# Percent end-to-end Reads (>=50%)
# MAPD Value (<=0.5)
# DNA Quality Score (0-1 poor, 2-3 questionable, 4-5 good)

df_assess_quality['DNA Score Count'] = 0
df_assess_quality.loc[df_assess_quality['DNA Mapped Reads (>=500,000)'] >= 500000, 'DNA Score Count'] += 1
df_assess_quality.loc[df_assess_quality['DNA Mean Read Length (>=50)'] >= 50, 'DNA Score Count'] += 1
df_assess_quality.loc[df_assess_quality['Uniformity of base coverage (>=90%)'] >= 90, 'DNA Score Count'] += 1
df_assess_quality.loc[df_assess_quality['Percent end-to-end Reads (>=50%)'] >= 50, 'DNA Score Count'] += 1
df_assess_quality.loc[df_assess_quality['MAPD Value (<=0.5)'] <= 0.5, 'DNA Score Count'] += 1

dna_conditions = [
            # questionable
            (df_assess_quality['DNA Score Count'] >= 2) & (df_assess_quality['DNA Score Count'] < 4),
            # good
            (df_assess_quality['DNA Score Count'] >= 4),
            # poor
            (df_assess_quality['DNA Score Count'] <= 1)
            ]
dna_choices = ['Questionable', 'Good', 'Poor']
df_assess_quality['DNA Quality Score'] = np.select(dna_conditions, dna_choices)

# Information to score RNA:
# RNA Mapped Reads >=20000
# RNA Mean Read Length >=35
# RNA Expression Ctrls Value >=5
# RNA Quality Score (0 poor, 1 questionable, 2-3 good)
df_assess_quality['RNA Score Count'] = 0
df_assess_quality.loc[df_assess_quality['RNA Mapped Reads (>=20,000)'] >= 20000, 'RNA Score Count'] += 1
df_assess_quality.loc[df_assess_quality['RNA Mean Read Length (>=35)'] >= 35, 'RNA Score Count'] += 1
df_assess_quality.loc[df_assess_quality['RNA Expression Ctrls Value (>=5)'] >= 5, 'RNA Score Count'] += 1

rna_conditions = [
            # questionable
            (df_assess_quality['RNA Score Count'] == 1), 
            # good
            (df_assess_quality['RNA Score Count'] >= 2),
            # poor
            (df_assess_quality['RNA Score Count'] == 0)
            ]
rna_choices = ['Questionable', 'Good', 'Poor']
df_assess_quality['RNA Quality Score'] = np.select(rna_conditions, rna_choices)

df_dna_rna_quality = df_assess_quality.sort_values(['DNA Quality Score', 'RNA Quality Score'], ascending=[True,True])
#df_assess_quality.sort_values(by='RNA Quality Score')

# Adding new PPG here
# PPG dictionary for SNV, CNV, and Fusion
ppgs_dict = {
        'PPG1': ['AKT3','ALK','CDK4','RAF1','ROS1','CD274'],
		'PPG2': ['AKT1','CDKN2A','CTNNB1','FGFR1','FGFR2','FGFR3','FGFR4',
            'HRAS','IDH1','IDH2','MET','RET','TACC3'],
		'PPG3': ['AR:','AR-','-AR' ,'AR CNV', 'AR Fusion','BRAF', 'BRAF Fusion','CHEK2:','EGFR','FLT3:','GNA11:','GNAQ:','GNAS:',
            'KIT:','KRAS','KRAS CNV','MAP2K1:','MAP2K2:','MTOR:','NRAS:','NTRK1','NTRK2','NTRK3'],
		'PPG4': ['ERBB2','ERBB3','ERBB4','ESR1','PIK3CA','PTEN',
            'NRG1','NUTM1'],
		'PPG5': ['AKT2','ARAF','PDGFRA','SMO','RSPO2','RSPO3'],
		'PPG6': ['TP53', 'No Variant Detected', 'No Fusion Detected'],
		}

for category, words_to_check in ppgs_dict.items():
    snv_matches = df_dna_rna_quality['SNV Summary'].str.contains('|'.join(words_to_check), case=False, regex=True)
    cnv_matches = df_dna_rna_quality['CNV Summary'].str.contains('|'.join(words_to_check), case=False, regex=True)
    fusion_matches = df_dna_rna_quality['Fusion Summary'].str.contains('|'.join(words_to_check), case=False, regex=True)
    
    combined_matches = snv_matches | cnv_matches | fusion_matches
    df_dna_rna_quality[category] = combined_matches

# Function to rank PPG columns
def rank_ppg(column_name):
    return int(column_name[-1])

# Function to get the highest ranked PPG column with a value of True for a given row
def get_highest_ranked_ppg(row):
    true_columns = [col for col in ppgs_dict.keys() if row[col]]
    if not true_columns:
        return None
    return min(true_columns, key=rank_ppg)

# Apply the function to each row and assign the result to a new column
df_dna_rna_quality['QCTracker best PPG'] = df_dna_rna_quality.apply(get_highest_ranked_ppg, axis=1)


df_reordered = df_dna_rna_quality[[
                        'Specimen ID','QCTracker best PPG','DNA Quality Score','RNA Quality Score','SNV Summary','CNV Summary','Fusion Summary',
                        
                        'AssayID','DNA Mapped Reads (>=500,000)','DNA Mean Read Length (>=50)','Uniformity of base coverage (>=90%)',
                        'Percent end-to-end Reads (>=50%)','MAPD Value (<=0.5)','MAPD Outcome',

                        'RNA Mapped Reads (>=20,000)','RNA Mean Read Length (>=35)','RNA Expression Ctrls Value (>=5)','RNA Expression Ctrls Outcome',

                        'SNV Indel ID','SNV Locus','SNV Allele Ref','SNV Allele Alt','SNV Type','SNV Cosmic',
                        'SNV Variant Class','SNV Gene Class','SNV Allele Frequency','SNV AA Change',

                        'CNV ID','CNV Locus','CNV Gene','CNV Variant Class','CNV Copy Number','CNV Ratio','CNV Gene Class','CNV Call',

                        'Fusion Variant ID','Fusion Locus','Fusion Variant Class','Fusion Gene Class','Fusion Read Counts',
                        'Fusion Type','Fusion Driver Gene','Fusion Mol Cov','Fusion Imbalance Score','Fusion Predicted Break-point Range',
                        'Fusion Read Counts Per Million',

                        'Percent Loading','Key Signal','Raw Read Accuracy','DNA Mean AQ20 Read Length','CF-1 Avg ReadsperLane',
                        'CF-1 BaseCallAccuracy','CF-1 AQ20 MeanReadLength','Total Assigned Amplicon Reads',
                        'Uniformity of Amplicon Coverage (%)','Amplicons With 1+ Read (%)','Amplicons With 100+ Reads (%)',
                        'Amplicons With 500+ Reads (%)','Amplicons Without Strand Bias (%)','Amplicons reading end-to-end',
                        'Avg Base Cov Depth','Target Base Cov at 1x (%)','Target Base Cov at 100x (%)',
                        'Target Base Cov at 500x (%)','Target bases Without Strand Bias (%)',

                        'DNA_bamID','RNA_bamID','DNA_barcode','RNA_barcode','GNXS','Analysis Date','Completion Date',
                        'Genexus Software Version','Assay Version','PPG1','PPG2','PPG3','PPG4','PPG5','PPG6',

                        ]]

# Getting reqid from the template file to give a more informative output file name
#mutUpload_template = 'MutationUploadTemplate.xlsx'
#batchinfo = pd.read_excel(mutUpload_template, sheet_name='BatchInfo')
#update the requisition name so we can use it to name the file
#batchinfo['ReqID'] = batchinfo['ReqID'].map(lambda x: x.replace(':', '-'))
#get the requisition name
#reqid =  batchinfo['ReqID'][0]
#Giving a name to the output file with requisition id
#outfile = f"{reqid}_OPA_QC_Tracker.xlsx"
outfile = f"tricom_2025_OPA_QC_Tracker.xlsx"

# Writing the pandas dataframe into an excel file
# start a writer because we want to add some formatting to the file
with pd.ExcelWriter(outfile, engine = 'xlsxwriter') as writer:
    # start an excel object with the pandas function 'to_excel()'
    df_reordered.to_excel(writer, sheet_name='QC Tracker', index = False)
    # start a workbook and the worksheet
    workbook = writer.book
    worksheet = writer.sheets['QC Tracker']

    # set up the formatting - cell color and font color
    # red
    format_poor = workbook.add_format({'bg_color': '#FFC7CE','font_color': '#9C0006'})
    format_poor_font = workbook.add_format({'font_color': '#FF0000'})
    # green
    format_good = workbook.add_format({'bg_color': '#C6EFCE','font_color': '#006100'})
    # yellow
    format_questionable = workbook.add_format({'bg_color': '#FFEB9C','font_color': '#9C6500'})
    
    # set up the rows and columns to have the formatting
    start_row = 1
    start_col = 2
    end_row = len(df_reordered)
    end_cold = 3

    worksheet.conditional_format(start_row, start_col, end_row, end_cold,{
        'type': 'text','criteria': 'containing','value': 'Poor','format': format_poor})

    worksheet.conditional_format(start_row, start_col, end_row, end_cold,{
        'type': 'text','criteria': 'containing','value': 'Good','format': format_good})

    worksheet.conditional_format(start_row, start_col, end_row, end_cold,{
        'type': 'text','criteria': 'containing','value': 'Questionable','format': format_questionable})

    #DNA metrics - highlight poor metrics
    dna_mapped_reads_start_row = 1
    dna_mapped_reads_start_col = 8
    dna_mapped_reads_end_row = len(df_reordered)
    dna_mapped_reads_end_cold = 8
    worksheet.conditional_format(dna_mapped_reads_start_row, dna_mapped_reads_start_col, dna_mapped_reads_end_row, dna_mapped_reads_end_cold,{
        'type': 'cell','criteria': '<','value': '500000','format': format_poor_font})

    dna_mean_read_length_start_row = 1
    dna_mean_read_length_start_col = 9
    dna_mean_read_length_end_row = len(df_reordered)
    dna_mean_read_length_end_cold = 9
    worksheet.conditional_format(dna_mean_read_length_start_row, dna_mean_read_length_start_col, dna_mean_read_length_end_row, dna_mean_read_length_end_cold,{
        'type': 'cell','criteria': '<','value': '50','format': format_poor_font})

    dna_uniformity_base_cov_start_row = 1
    dna_uniformity_base_cov_start_col = 10
    dna_uniformity_base_cov_end_row = len(df_reordered)
    dna_uniformity_base_cov_end_cold = 10
    worksheet.conditional_format(dna_uniformity_base_cov_start_row, dna_uniformity_base_cov_start_col, dna_uniformity_base_cov_end_row, dna_uniformity_base_cov_end_cold,{
        'type': 'cell','criteria': '<','value': '90','format': format_poor_font})

    dna_percent_etoe_start_row = 1
    dna_percent_etoe_start_col = 11
    dna_percent_etoe_end_row = len(df_reordered)
    dna_percent_etoe_end_cold = 11
    worksheet.conditional_format(dna_percent_etoe_start_row, dna_percent_etoe_start_col, dna_percent_etoe_end_row, dna_percent_etoe_end_cold,{
        'type': 'cell','criteria': '<','value': '50','format': format_poor_font})
    
    dna_mapd_start_row = 1
    dna_mapd_start_col = 12
    dna_mapd_end_row = len(df_reordered)
    dna_mapd_end_cold = 12
    worksheet.conditional_format(dna_mapd_start_row, dna_mapd_start_col, dna_mapd_end_row, dna_mapd_end_cold,{
        'type': 'cell','criteria': '>','value': '0.5','format': format_poor_font})

    dna_mapd_outcome_start_row = 1
    dna_mapd_outcome_start_col = 13
    dna_mapd_outcome_end_row = len(df_reordered)
    dna_mapd_outcome_end_cold = 13
    worksheet.conditional_format(dna_mapd_outcome_start_row, dna_mapd_outcome_start_col, dna_mapd_outcome_end_row, dna_mapd_outcome_end_cold,{
        'type': 'text','criteria': 'containing','value': 'Failed','format': format_poor_font})

    #RNA metrics - highlight poor metrics
    rna_mapped_reads_start_row = 1
    rna_mapped_reads_start_col = 14
    rna_mapped_reads_end_row = len(df_reordered)
    rna_mapped_reads_end_cold = 14
    worksheet.conditional_format(rna_mapped_reads_start_row, rna_mapped_reads_start_col, rna_mapped_reads_end_row, rna_mapped_reads_end_cold,{
        'type': 'cell','criteria': '<','value': '20000','format': format_poor_font})

    rna_mean_read_length_start_row = 1
    rna_mean_read_length_start_col = 15
    rna_mean_read_length_end_row = len(df_reordered)
    rna_mean_read_length_end_cold = 15
    worksheet.conditional_format(rna_mean_read_length_start_row, rna_mean_read_length_start_col, rna_mean_read_length_end_row, rna_mean_read_length_end_cold,{
        'type': 'cell','criteria': '<','value': '35','format': format_poor_font})

    rna_expr_ctrl_values_start_row = 1
    rna_expr_ctrl_values_start_col = 16
    rna_expr_ctrl_values_end_row = len(df_reordered)
    rna_expr_ctrl_values_end_cold = 16
    worksheet.conditional_format(rna_expr_ctrl_values_start_row, rna_expr_ctrl_values_start_col, rna_expr_ctrl_values_end_row, rna_expr_ctrl_values_end_cold,{
        'type': 'cell','criteria': '<','value': '5','format': format_poor_font})

    rna_expr_ctrl_outcome_start_row = 1
    rna_expr_ctrl_outcome_start_col = 17
    rna_expr_ctrl_outcome_end_row = len(df_reordered)
    rna_expr_ctrl_outcome_end_cold = 17
    worksheet.conditional_format(rna_expr_ctrl_outcome_start_row, rna_expr_ctrl_outcome_start_col, rna_expr_ctrl_outcome_end_row, rna_expr_ctrl_outcome_end_cold,{
        'type': 'text','criteria': 'containing','value': 'Failed','format': format_poor_font})

    worksheet.autofit()

print('QC tracker job completed')
