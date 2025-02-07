"""
This files contains helper functions to process quant files
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr, pearsonr
import pandas as pd
import datetime
import time
import os
import re
from collections import defaultdict
from scipy.interpolate import UnivariateSpline
# import seaborn as sns
from scipy.integrate import trapz
import sys
import csv


def create_image_name(name, format=".png"):

    crnt_tm = datetime.datetime.now()
    image_name = (name+"_" + str(crnt_tm.year) + "_" + str(crnt_tm.month) + "_" + str(crnt_tm.day) + "_"
                  + time.strftime("%H_%M_%S") + format)
    return image_name

def plot_EM_results(alpha_history, convergence_history, theta_history):
    # Plot convergence over iterations
    plt.figure(figsize=(12, 4))
    plt.subplot(1, 3, 1)
    plt.plot(convergence_history)
    plt.title('Convergence Over Iterations')
    plt.xlabel('Iteration')
    plt.ylabel('Convergence')

    # Plot alpha values over iterations
    plt.subplot(1, 3, 2)
    plt.plot(alpha_history)
    plt.title('Alpha Values Over Iterations')
    plt.xlabel('Iteration')
    plt.ylabel('Alpha Value')
    plt.legend()
    plt.savefig("figures/"+create_image_name("convergence_alpha"))


    for sample_key in theta_history:
        # Determine the consistent order of isoforms across all iterations
        isoforms = list(theta_history[sample_key][0].keys())  # Assuming the first Counter has all isoforms

        # Convert the Counters to a NumPy array
        theta_values = np.array([
            [counter[isoform] for isoform in isoforms]
            for counter in theta_history['sample1']
        ])

        # Now you can plot each column of theta_values across iterations
        plt.figure()
        for i, isoform in enumerate(isoforms):
            plt.plot(theta_values[:, i], label=isoform)
        plt.legend()
        plt.title(f'Isoform Values Across Iterations for {sample_key}')
        plt.xlabel('Iteration')
        plt.ylabel('Value')
        plt.savefig("figures/" + create_image_name("theta_"+sample_key))
        plt.show()
        plt.close()


def fraction_to_float(fraction_str):
    # Strip the input to remove any leading/trailing whitespace
    fraction_str = fraction_str.strip()

    # Try to convert the string directly to a float (works if there is no fraction)
    try:
        return float(fraction_str)
    except ValueError:
        # Split the string on the '/'
        fraction_parts = fraction_str.split('/')
        if len(fraction_parts) == 2:  # It's a fraction
            numerator = float(fraction_parts[0])
            denominator = float(fraction_parts[1])
            return numerator / denominator
        else:
            raise ValueError(f"Invalid fraction string: {fraction_str}")


def spearman_corr_SIRV(file_path1, sample):

    ground_truth = pd.read_csv('../../data/SIRV/SIRV_Shree/E2_molarity.csv')
    # Reset the index of our_quant so that the isoform names become a column
    our_quant = pd.read_csv(file_path1, sep="\t")
    our_quant_reset = our_quant.reset_index()

    # Rename the columns accordingly
    #our_quant_reset.columns = ['transcript_name', 'other_column1', 'est_count', 'tpm']

    # Clean the isoform names in our_quant to match the naming convention in df2
    our_quant_reset['cleaned_name'] = our_quant_reset['transcript_name'].str.replace(r'\(\+\)|\(\-\)', '', regex=True)

    # Initialize lists to store matched TPM and molarity values
    matched_tpm = []
    matched_molarity = []

    # Iterate over our_quant to find matching molarity values in ground_truth
    for index, row in our_quant_reset.iterrows():
        # Extract the cleaned isoform name and tpm value
        cleaned_name = row['cleaned_name']
        tpm_value = row['tpm']

        # Find the corresponding molarity value in ground_truth (assuming 'E2' column has the molarity)
        molarity_value = ground_truth.loc[ground_truth['Name'] == cleaned_name, 'E2'].values

        # If a matching isoform is found in ground_truth
        if len(molarity_value) > 0:
            # Append the tpm and molarity to the respective lists
            matched_tpm.append(tpm_value)
            matched_molarity.append(fraction_to_float(molarity_value[0]))  # Take the first match in case of multiple

    # Calculate Spearman's correlation using the matched lists
    correlation, p_value = spearmanr(matched_tpm, matched_molarity)

    # Output the results
    print(f'Spearman correlation coefficient {sample}: {correlation}')
    print(f'P-value {sample}: {p_value}')

def fraction_to_float_gen(value):
    try:
        return float(value)
    except ValueError:
        return None

def split_transcript_name(transcript):
    return transcript.split('|')[0]


def csv_tpm_processing(file_path1, file_path2, suffixes=('_quant', '_truth'), sep_by="\t"):
    # Load the datasets
    our_quant = pd.read_csv(file_path1, sep="\t")
    ground_truth = pd.read_csv(file_path2, sep=sep_by)

    # Apply the function to the 'transcript_name' column of your DataFrame
    our_quant['transcript_name'] = our_quant['transcript_name'].apply(split_transcript_name)
    ground_truth['transcript_name'] = ground_truth['transcript_name'].apply(split_transcript_name)

    # Find common isoforms
    common_isoforms = pd.merge(our_quant, ground_truth, on='transcript_name', suffixes=suffixes)

    return common_isoforms

def format_file_name(file_path1, file_path2):
    file_name1 = file_path1.split('/')[-1]
    file_name2 = file_path2.split('/')[-1]
    part1 = "_".join(file_name1.split('_')[3:13])
    part2 = "_".join(file_name2.split('_')[3:13])

    print()
    return part1, part2
def spearman_corr_generic(file_path1, file_path2, sep_by="\t"):
    common_isoforms = csv_tpm_processing(file_path1, file_path2, sep_by=sep_by)

    if 'tpm_quant' in common_isoforms and 'tpm_truth' in common_isoforms:
        common_isoforms['tpm_truth'] = common_isoforms['tpm_truth'].apply(fraction_to_float_gen)
        correlation, p_value = spearmanr(common_isoforms['tpm_quant'], common_isoforms['tpm_truth'])

        part1, part2 = format_file_name(file_path1, file_path2)

        formatted_output = (

            f"{part1} and {part2}.\nSpearman correlation: {correlation:.6f}"
        )
        # with open(log_file, 'a') as f:
        #     f.write(formatted_output + '\n')

        print(formatted_output)
        print(f"P-value for correlation: {p_value}")
    else:
        print("TPM columns missing or incorrectly named in one of the datasets.")
    
    return correlation


def spearman_pearson_corr_generic(file_path1, file_path2, sep_by="\t"):
    common_isoforms = csv_tpm_processing(file_path1, file_path2, sep_by=sep_by)

    if 'tpm_quant' in common_isoforms and 'tpm_truth' in common_isoforms:
        common_isoforms['tpm_truth'] = common_isoforms['tpm_truth'].apply(fraction_to_float_gen)
        spearman_corr, p_value = spearmanr(common_isoforms['tpm_quant'], common_isoforms['tpm_truth'])
        pearson_corr, p_value = pearsonr(common_isoforms['tpm_quant'], common_isoforms['tpm_truth'])

        part1, part2 = format_file_name(file_path1, file_path2)

        formatted_output = (

            f"{part1} and {part2}.\nSpearman correlation: {spearman_corr:.3f}\nPearson correlation: {pearson_corr:.3f}"
        )

        # log_file = '/'.join(file_path1.split('/')[0:-1])+'/corr.txt'
        # with open(log_file, 'a') as f:
        #     f.write(formatted_output + '\n')

        print(formatted_output)
        print(f"P-value for correlation: {p_value}")
    else:
        print("TPM columns missing or incorrectly named in one of the datasets.")
    
    return round(spearman_corr, 4), round(pearson_corr, 4)



def pair_files_exp4(directory, simulation, pair_type):
    # Dictionary to store the files based on downsampling percentage, length type, and additional metrics
    long_file_pairs = defaultdict(list)
    short_file_pairs = defaultdict(list)
    file_info = defaultdict(list)


    if simulation:
        file_pattern = re.compile(
    r'output_Simulation_VIGD_token_(\d+)_sample(\d+)_file(\d+)_ds(\d+)num(\d+)aln(\d+)(long|short)_file(\d+)_ds(\d+)num(\d+)aln(\d+)(long|short)_GDlr_(0\.\d+)_AlphaInitial_(\d+(?:\.\d+)?)_EMround_(\d+)_'
)
    else:
        file_pattern = re.compile(
        r'output_PacIllu_VIGD_token_(\d+)_sample(\d+)_file(\d+)_ds(\d+)num(\d+)aln(\d+)(long|short)_file(\d+)_ds(\d+)num(\d+)aln(\d+)(long|short)_GDlr_(0\.\d+)_AlphaInitial_(\d+(?:\.\d+)?)_EMround_(\d+)_'
    )
    # List all files in the directory
    for file in os.listdir(directory):
        match = file_pattern.search(file)
        if match:
            # Parse the matched groups
            (token, sample, file_num1, ds_percentage1, num1, aln_replica1, length1,
                file_num2, ds_percentage2, num2, aln_replica2, length2,
                GDlr, AlphaInitial, EMround
            ) = match.groups()
            # Determine which file to parse based on the sample number
            if sample == '1':
                ds_percentage, file_num, aln_replica, day, length = ds_percentage1, num1, aln_replica1[1], aln_replica1[0], length1
            elif sample == '2':
                ds_percentage, file_num, aln_replica, day, length = ds_percentage2, num2, aln_replica2[1], aln_replica1[0], length2
            
            if length == 'long':
                key = (ds_percentage, day, file_num, GDlr, AlphaInitial, EMround)
                long_file_pairs[key].append((file, token))
            else:
                key = (token)
                short_file_pairs[key].append((file, token))
    
    paired_files = []
    ### DO NOT ERASE
    # for key, files in long_file_pairs.items():
    #     paired_files.append((files[0][0], files[1][0]))
    
    ## (AT)
    # Create long read pairs
    if pair_type == 'replica':
        for key, files in long_file_pairs.items():
            if len(files) > 1:  # Ensure there are multiple replicas
                # Pair long read files for each replica
                for i in range(len(files)):
                    for j in range(i + 1, len(files)):
                        paired_files.append((files[i][0], files[j][0]))
                        
                        # Corresponding short read files for each replica
                        sr_file1 = short_file_pairs[files[i][1]][0][0]
                        sr_file2 = short_file_pairs[files[j][1]][0][0]
                        paired_files.append((sr_file1, sr_file2))

    # This will give you LR and SR file names those were trained together, eg: lr_01_replica1+sr_01_replica2,
    elif pair_type == 'within_trainee':
        for key, files in long_file_pairs.items():
            for i in range(len(files)):
                lr_file = files[i][0]
                sr_file = short_file_pairs[files[i][1]][0][0]
                paired_files.append((lr_file, sr_file))

    return paired_files

def pair_files_exp1(directory, simulation, pair_type):
    # Dictionary to store the files based on downsampling percentage, length type, and additional metrics
    long_file_pairs = defaultdict(list)
    short_file_pairs = defaultdict(list)
    file_info = defaultdict(list)


    if simulation:
        file_pattern = re.compile(
        r'output_Simulation_VIGD_token_(\d+)_sample(\d+)_file(\d+)_ds(\d+)num(\d+)aln(\d+)(long|short)_GDlr_(0\.\d+)_AlphaInitial_(\d+(?:\.\d+)?)_EMround_(\d+)_'
    )
    else:
        file_pattern = re.compile(
        r'output_PacIllu_VIGD_token_(\d+)_sample(\d+)_file(\d+)_ds(\d+)num(\d+)aln(\d+)(long|short)_GDlr_(0\.\d+)_AlphaInitial_(\d+(?:\.\d+)?)_EMround_(\d+)_'
    )
    # List all files in the directory
    for file in os.listdir(directory):
        match = file_pattern.search(file)
        if match:
            # Parse the matched groups
            (token, sample, file_num1, ds_percentage1, num1, aln_replica1, length1,
                GDlr, AlphaInitial, EMround
            ) = match.groups()
            # Determine which file to parse based on the sample number
            ds_percentage, file_num, aln_replica, day, length = ds_percentage1, num1, aln_replica1[1], aln_replica1[0], length1
            
            # if sample == '1':
            #     ds_percentage, file_num, aln_replica, day, length = ds_percentage1, num1, aln_replica1[1], aln_replica1[0], length1
            # elif sample == '2':
            #     ds_percentage, file_num, aln_replica, day, length = ds_percentage2, num2, aln_replica2[1], aln_replica1[0], length2
            
            if length == 'long':
                key = (day)
                long_file_pairs[key].append((file, token))
            else:
                key = (day)
                short_file_pairs[key].append((file, token))
    
    paired_files = []
    ### DO NOT ERASE
    # for key, files in long_file_pairs.items():
    #     paired_files.append((files[0][0], files[1][0]))
    
    ## (AT)
    # Create long read pairs
    if pair_type == 'replica':
        for key, files in long_file_pairs.items():
            if len(files) > 1:  # Ensure there are multiple replicas
                # Pair long read files for each replica
                for i in range(len(files)):
                    for j in range(i + 1, len(files)):
                        paired_files.append((files[i][0], files[j][0]))
    
    if pair_type == 'replica':
        for key1, files1 in short_file_pairs.items():
            if len(files1) > 1:  # Ensure there are multiple replicas
                # Pair long read files for each replica
                for i in range(len(files1)):
                    for j in range(i + 1, len(files1)):
                        paired_files.append((files1[i][0], files1[j][0]))
                        
                        # # Corresponding short read files for each replica
                        # sr_file1 = short_file_pairs[files[i][1]][0][0]
                        # sr_file2 = short_file_pairs[files[j][1]][0][0]
                        # paired_files.append((sr_file1, sr_file2))

    # This will give you LR and SR file names those were trained together, eg: lr_01_replica1+sr_01_replica2,
    # elif pair_type == 'within_trainee':
    #     for key, files in long_file_pairs.items():
    #         for i in range(len(files)):
    #             lr_file = files[i][0]
    #             sr_file = short_file_pairs[files[i][1]][0][0]
    #             paired_files.append((lr_file, sr_file))

    return paired_files


def pair_files_exp2(directory, simulation, pair_type):
    # Dictionary to store the files based on downsampling percentage, length type, and additional metrics
    long_file_pairs = defaultdict(list)
    short_file_pairs = defaultdict(list)
    file_info = defaultdict(list)


    if simulation:
        file_pattern = re.compile(
    r'output_Simulation_VIGD_token_(\d+)_sample(\d+)_file(\d+)_ds(\d+)num(\d+)aln(\d+)(long|short)_file(\d+)_ds(\d+)num(\d+)aln(\d+)(long|short)_GDlr_(0\.\d+)_AlphaInitial_(\d+(?:\.\d+)?)_EMround_(\d+)_'
)
    else:
        file_pattern = re.compile(
        r'output_PacIllu_VIGD_token_(\d+)_sample(\d+)_file(\d+)_ds(\d+)num(\d+)aln(\d+)(long|short)_file(\d+)_ds(\d+)num(\d+)aln(\d+)(long|short)_GDlr_(0\.\d+)_AlphaInitial_(\d+(?:\.\d+)?)_EMround_(\d+)_'
    )
    # List all files in the directory
    for file in os.listdir(directory):
        match = file_pattern.search(file)
        if match:
            # Parse the matched groups
            (token, sample, file_num1, ds_percentage1, num1, aln_replica1, length1,
                file_num2, ds_percentage2, num2, aln_replica2, length2,
                GDlr, AlphaInitial, EMround
            ) = match.groups()
            # Determine which file to parse based on the sample number
            ds_percentage, file_num, aln_replica, day, length = ds_percentage1, num1, aln_replica1[1], aln_replica1[0], length1
            
            # if sample == '1':
            #     ds_percentage, file_num, aln_replica, day, length = ds_percentage1, num1, aln_replica1[1], aln_replica1[0], length1
            # elif sample == '2':
            #     ds_percentage, file_num, aln_replica, day, length = ds_percentage2, num2, aln_replica2[1], aln_replica1[0], length2
            
            if length == 'long':
                key = (day)
                long_file_pairs[key].append((file, token))
            else:
                key = (day)
                short_file_pairs[key].append((file, token))
    
    paired_files = []
    ### DO NOT ERASE
    # for key, files in long_file_pairs.items():
    #     paired_files.append((files[0][0], files[1][0]))
    
    ## (AT)
    # Create long read pairs
    if pair_type == 'replica':
        for key, files in long_file_pairs.items():
            if len(files) > 1:  # Ensure there are multiple replicas
                # Pair long read files for each replica
                for i in range(len(files)):
                    for j in range(i + 1, len(files)):
                        paired_files.append((files[i][0], files[j][0]))
    
    # if pair_type == 'replica':
    #     for key1, files1 in short_file_pairs.items():
    #         if len(files1) > 1:  # Ensure there are multiple replicas
    #             # Pair long read files for each replica
    #             for i in range(len(files1)):
    #                 for j in range(i + 1, len(files1)):
    #                     paired_files.append((files1[i][0], files1[j][0]))
                        
                        # # Corresponding short read files for each replica
                        # sr_file1 = short_file_pairs[files[i][1]][0][0]
                        # sr_file2 = short_file_pairs[files[j][1]][0][0]
                        # paired_files.append((sr_file1, sr_file2))

    # This will give you LR and SR file names those were trained together, eg: lr_01_replica1+sr_01_replica2,
    # elif pair_type == 'within_trainee':
    #     for key, files in long_file_pairs.items():
    #         for i in range(len(files)):
    #             lr_file = files[i][0]
    #             sr_file = short_file_pairs[files[i][1]][0][0]
    #             paired_files.append((lr_file, sr_file))

    return paired_files

def pair_files_exp5(directory, simulation, pair_type):
    # Dictionary to store the files based on downsampling percentage, length type, and additional metrics
    long_file_pairs = defaultdict(list)
    short_file_pairs = defaultdict(list)
    file_info = defaultdict(list)

    ## (AT)
    # Updated regular expression to match the files with additional metrics
    # file_pattern = re.compile(
    #     r'output_PacIllu_VIGD_token_(\d+)_sample(\d+)_file(\d+)_ds(\d+)num(\d+)aln(\d+)(long|short)_file(\d+)_ds(\d+)num(\d+)aln(\d+)(long|short)_GDlr_(0\.\d+)_AlphaInitial_(\d+(?:\.\d+)?)_EMround_(\d+)_'
    # )
    if simulation:

        file_pattern = re.compile(
        r'output_Simulation_VIGD_token_(\d+)_sample(\d+)_file(\d+)_.*_file(\d+)_.*_GDlr_(0\.\d+)_AlphaInitial_(\d+(?:\.\d+)?)_EMround_(\d+)_'
    )
    else:
        file_pattern = re.compile(
        r'output_PacIllu_VIGD_token_(\d+)_sample(\d+)_file(\d+)_.*_file(\d+)_.*_GDlr_(0\.\d+)_AlphaInitial_(\d+(?:\.\d+)?)_EMround_(\d+)_'
    )
    # List all files in the directory
    for file in os.listdir(directory):
        match = file_pattern.search(file)
        if match:
            # Parse the matched groups
            (token, sample, file1, file2, GDlr, AlphaInitial, EMround) = match.groups()
            key = sample
            long_file_pairs[key].append((file, token))
            
    
    paired_files = []

    ### DO NOT ERASE
    # for key, files in long_file_pairs.items():
    #     paired_files.append((files[0][0], files[1][0]))
    
    ## (AT)
    # Create long read pairs
    if pair_type == 'replica':
        for key, files in long_file_pairs.items():
            if len(files) > 1:  # Ensure there are multiple replicas
                # Pair long read files for each replica
                for i in range(len(files)):
                    for j in range(i + 1, len(files)):
                        paired_files.append((files[i][0], files[j][0]))
                        
                        # # Corresponding short read files for each replica
                        # sr_file1 = short_file_pairs[files[i][1]][0][0]
                        # sr_file2 = short_file_pairs[files[j][1]][0][0]
                        # paired_files.append((sr_file1, sr_file2))

    # This will give you LR and SR file names those were trained together, eg: lr_01_replica1+sr_01_replica2,
    # elif pair_type == 'within_trainee':
    #     for key, files in long_file_pairs.items():
    #         for i in range(len(files)):
    #             lr_file = files[i][0]
    #             sr_file = short_file_pairs[files[i][1]][0][0]
    #             paired_files.append((lr_file, sr_file))

    return paired_files
    
# Function to calculate CV
def calculate_cv(data):
    # colm_name = 'long' if type == 'long' else 'short'
    colm_name = ''
    data_log = np.log(data[[f'tpm_{colm_name}Rep1', f'tpm_{colm_name}Rep2']] + 1)
    data['mean_abundance'] = data_log[[f'tpm_{colm_name}Rep1', f'tpm_{colm_name}Rep2']].mean(axis=1)
    data['std_abundance'] = data_log[[f'tpm_{colm_name}Rep1', f'tpm_{colm_name}Rep2']].std(axis=1)
    CV_ig = data['std_abundance'] / data['mean_abundance']
    data['CV'] = CV_ig

    CV_ig_squared = CV_ig ** 2
    IM = np.sqrt(CV_ig_squared.mean())

    sorted_u_ig = data['mean_abundance'].sort_values()
    sorted_CV_ig = data['CV'].loc[sorted_u_ig.index]

    ACVC = np.trapz(sorted_CV_ig, x=sorted_u_ig)
    result = f"IM {IM}, ACVC {ACVC}\n"
    # with open(log_file, 'a') as f:
    #     f.write(result)
    print(result)

    return data, ACVC, IM

def calculate_acvc(cv_values, abundance):
    return trapz(cv_values, abundance)

def calculate_im_acvc(rep1, rep2, directory):

    rep1_data = csv_tpm_processing(directory+rep1, directory+rep2, suffixes=('_Rep1', '_Rep2'))

    rep1_data, ACVC, IM = calculate_cv(rep1_data)
    rep1_data_sorted = rep1_data.sort_values('mean_abundance')

    part1, part2 = format_file_name(rep1, rep2)

    # Split the path into parts
    path_parts = directory.split('/')
    # Remove the last two parts ('files' and 'output_files')
    new_path_parts = path_parts[:-3]
    # Append 'figures' to the path
    new_path_parts.append('figures')
    # Join the parts to form the new path
    figure_dir = '/'.join(new_path_parts)
    
    # (AT)
    # fig, axes = plt.subplots(1, 2, figsize=(16, 8))

    # sns.boxplot(y='CV', data=rep1_data_sorted, palette="Set2", ax=axes[0], boxprops=dict(alpha=0.5))
    # axes[0].set_title('Box Plot for CV')
    # axes[0].set_ylabel('CV')
    # axes[0].grid(True, which='both', linestyle='--', linewidth=0.5)
    # axes[0].set_ylim(0, 1.5)

    # axes[1].plot(rep1_data_sorted['mean_abundance'], rep1_data_sorted['CV'], 'o')
    # axes[1].set_xlabel('Transcript abundance (log2(TPM+1))')
    # axes[1].set_ylabel('CV')
    # axes[1].set_title('CV vs Transcript Abundance')
    # axes[1].set_ylim(0, 1.5)
    # axes[1].set_xlim(0, 10)
    # plt.tight_layout()
    # plt.title(f"stats for {part1} and {part2}")
    # timestamp = time.strftime("_%Y_%m_%d__%H_%M_%S")
    # # (AT)
    # # plt.savefig(os.path.join(figure_dir, create_image_name(f"IM_CV_{type}" + timestamp + '.png')))
    # plt.show()
    # plt.close()
    return round(ACVC, 4), round(IM, 4)


def merge_csv_files(file1, file2, output_dir):
    # Load the first CSV file
    df1 = pd.read_csv(file1, sep="\t")

    # Load the second CSV file
    df2 = pd.read_csv(file2, sep="\t")

    # Identify common transcripts
    # Find common isoforms
    # common_isoforms = pd.merge(df1, df2, on='transcript_name')

    common_transcripts = pd.merge(df1, df2, left_on='transcript_name', right_on='transcript_name',
                                  suffixes=('_sample1', '_sample2'))

    # Select relevant columns
    merged_df = common_transcripts[['transcript_name', 'tpm_sample1', 'tpm_sample2']]

    # Rename columns to match the desired format
    merged_df.columns = ['ID', 'sample1', 'sample2']

    # Save the merged dataframe to a new CSV file
    output_file = output_dir+'merged'+time.strftime("_%Y_%m_%d__%H_%M_%S")+'.tsv'
    merged_df.to_csv(output_file, index=False)

    print(f'Merged file saved as {output_file}')

def csv_row_parsing_exp4(simulation, pair):

    ## (AT)
    # file_pattern = re.compile(
    #         r'output_PacIllu_VIGD_token_(\d+)_sample(\d+)_file(\d+)_ds(\d+)num(\d+)aln(\d+)(long|short)_file(\d+)_ds(\d+)num(\d+)aln(\d+)(long|short)_GDlr_(0\.\d+)_AlphaInitial_(\d+(?:\.\d+)?)_EMround_(\d+)_'
    #     )
    if simulation:
        file_pattern = re.compile(
        r'output_Simulation_VIGD_token_(\d+)_sample(\d+)_file(\d+)_ds(\d+)num(\d+)aln(\d+)(long|short)_file(\d+)_ds(\d+)num(\d+)aln(\d+)(long|short)_GDlr_(0\.\d+)_AlphaInitial_(\d+(?:\.\d+)?)_EMround_(\d+)_'
    )
    else:
        file_pattern = re.compile(
    r'output_PacIllu_VIGD_token_(\d+)_sample(\d+)_file(\d+)_ds(\d+)num(\d+)aln(\d+)(long|short)_file(\d+)_ds(\d+)num(\d+)aln(\d+)(long|short)_GDlr_(0\.\d+)_AlphaInitial_(\d+(?:\.\d+)?)_EMround_(\d+)_'
)
    i = 0
    for file in pair:
        match = file_pattern.search(file)
        (token, sample, file_num1, ds_percentage1, num1, aln_replica1, length1,
                file_num2, ds_percentage2, num2, aln_replica2, length2,
                GDlr, AlphaInitial, EMround
            ) = match.groups()
        if i == 0:
            if sample == '1':
                file_num1, replica1, dayF1, length, ds_prctF1 = num1, aln_replica1[1], aln_replica1[0], length1, ds_percentage1
            elif sample == '2':
                file_num1, replica1, dayF1, length, ds_prctF1 = num2, aln_replica2[1], aln_replica2[0], length2, ds_percentage2
            tokenF1, sampleF1 = token, sample
            file_name1 = f'ds{ds_prctF1}num{file_num1}aln{dayF1}{replica1}{length}'
        else:
            if sample == '1':
                file_num2, replica2, dayF2, length, ds_prctF2 = num1, aln_replica1[1], aln_replica1[0], length1, ds_percentage1
            elif sample == '2':
                file_num2, replica2, dayF2, length, ds_prctF2 = num2, aln_replica2[1], aln_replica2[0], length2, ds_percentage2
            tokenF2, sampleF2 = token, sample
            file_name2 = f'ds{ds_prctF2}num{file_num2}aln{dayF2}{replica2}{length}'
        i+=1
    return file_name1, file_name2, GDlr, AlphaInitial, EMround, length, file_num1, file_num2, replica1, replica2, dayF1, dayF2, tokenF1, tokenF2, sampleF1, sampleF2, ds_prctF1, ds_prctF2


def csv_row_parsing_exp1(simulation, pair):

    ## (AT)
    # file_pattern = re.compile(
    #         r'output_PacIllu_VIGD_token_(\d+)_sample(\d+)_file(\d+)_ds(\d+)num(\d+)aln(\d+)(long|short)_file(\d+)_ds(\d+)num(\d+)aln(\d+)(long|short)_GDlr_(0\.\d+)_AlphaInitial_(\d+(?:\.\d+)?)_EMround_(\d+)_'
    #     )
    if simulation:
        file_pattern = re.compile(
        r'output_Simulation_VIGD_token_(\d+)_sample(\d+)_file(\d+)_ds(\d+)num(\d+)aln(\d+)(long|short)_GDlr_(0\.\d+)_AlphaInitial_(\d+(?:\.\d+)?)_EMround_(\d+)_'
    )
    else:
        file_pattern = re.compile(
    r'output_PacIllu_VIGD_token_(\d+)_sample(\d+)_file(\d+)_ds(\d+)num(\d+)aln(\d+)(long|short)_GDlr_(0\.\d+)_AlphaInitial_(\d+(?:\.\d+)?)_EMround_(\d+)_'
)
    i = 0
    for file in pair:
        match = file_pattern.search(file)
        (token, sample, file_num1, ds_percentage1, num1, aln_replica1, length1,
                GDlr, AlphaInitial, EMround
            ) = match.groups()
        
        if i == 0:
            file_num1, replica1, dayF1, ds_prctF1 = num1, aln_replica1[1], aln_replica1[0], ds_percentage1
            # if sample == '1':
            #     file_num1, replica1, dayF1, ds_prctF1 = num1, aln_replica1[1], aln_replica1[0], ds_percentage1
            # elif sample == '2':
            #     file_num1, replica1, dayF1, ds_prctF1 = num2, aln_replica2[1], aln_replica2[0], ds_percentage2
            tokenF1, sampleF1 = token, sample
            file_name1 = f'ds{ds_prctF1}num{file_num1}aln{dayF1}{replica1}{length1}'
        else:
            file_num2, replica2, dayF2, length, ds_prctF2 = num1, aln_replica1[1], aln_replica1[0], length1, ds_percentage1
            # if sample == '1':
            #     file_num2, replica2, dayF2, length, ds_prctF2 = num1, aln_replica1[1], aln_replica1[0], length1, ds_percentage1
            # elif sample == '2':
            #     file_num2, replica2, dayF2, length, ds_prctF2 = num2, aln_replica2[1], aln_replica2[0], length2, ds_percentage2
            tokenF2, sampleF2 = token, sample
            file_name2 = f'ds{ds_prctF2}num{file_num2}aln{dayF2}{replica2}{length1}'
        i+=1
    return file_name1, file_name2, GDlr, AlphaInitial, EMround, length, file_num1, file_num2, replica1, replica2, dayF1, dayF2, tokenF1, tokenF2, sampleF1, sampleF2, ds_prctF1, ds_prctF2

def csv_row_parsing_exp5(simulation, pair):

    ## (AT)
    # file_pattern = re.compile(
    #         r'output_PacIllu_VIGD_token_(\d+)_sample(\d+)_file(\d+)_ds(\d+)num(\d+)aln(\d+)(long|short)_file(\d+)_ds(\d+)num(\d+)aln(\d+)(long|short)_GDlr_(0\.\d+)_AlphaInitial_(\d+(?:\.\d+)?)_EMround_(\d+)_'
    #     )
    if simulation:
        file_pattern = re.compile(
        r'output_Simulation_VIGD_token_(\d+)_sample(\d+)_file(\d+)_.*_file(\d+)_.*_GDlr_(0\.\d+)_AlphaInitial_(\d+(?:\.\d+)?)_EMround_(\d+)_'
    )
    else:
        file_pattern = re.compile(
        r'output_PacIllu_VIGD_token_(\d+)_sample(\d+)_file(\d+)_.*_file(\d+)_.*_GDlr_(0\.\d+)_AlphaInitial_(\d+(?:\.\d+)?)_EMround_(\d+)_'
    )
    i = 0
    for file in pair:
        match = file_pattern.search(file)
        if match:
            # Parse the matched groups
            (token, sample, file1, file2, GDlr, AlphaInitial, EMround) = match.groups()
          
        if i == 0:
            tokenF1, sampleF1 = token, sample
        else:
            tokenF2, sampleF2 = token, sample
        i+=1
    return GDlr, AlphaInitial, EMround, tokenF1, tokenF2, sampleF1, sampleF2
               


def main():
    experiment_file = 'exprmnt_2024_12_08__23_58_50'
    simulation = 1
    # "Different experiment setup, 1: for 1 sample, 2 for merged, 4 for multisample, 5 for merged multisample"
    experiment = 2




    main_dir = '/gpfs/commons/home/atalukder/RNA_Splicing/files/results'
    directory = os.path.join(main_dir, experiment_file, 'files/output_files/')
    # paired_files, 2 types
    # type == 'replica': This will give you LR and SR file names with same replica, eg: lr_01_replica1+lr_01_replica2,
    # type == 'within_trainee': This will give you LR and SR file names those were trained together, eg: lr_01_replica1+sr_01_replica2,
    if experiment == 4:
        paired_files = pair_files_exp4(directory, simulation, pair_type='replica')
    elif experiment == 1:
        paired_files = pair_files_exp1(directory, simulation, pair_type='replica')
    elif experiment == 2:
        paired_files = pair_files_exp2(directory, simulation, pair_type='replica')
    elif experiment == 5:
        paired_files = pair_files_exp5(directory, simulation, pair_type='replica')



    # paired_files = pair_files(directory, pair_type='within_trainee')
    timestamp = time.strftime("_%Y_%m_%d__%H_%M_%S")
    spearman_corr, ACVC, IM = 0, 0, 0

     # Prepare the CSV file for writing
    csv_file = directory + f'correlations_results{timestamp}.csv'

    # Open the CSV file for writing
    with open(csv_file, mode='w', newline='') as file:
        writer = csv.writer(file)
        
        # if experiment == 1:
        #     # Write the header row
        #     key_file_info = (['Spearman_Corr', 'Pearson_Corr', 'ACVC', 'IM', 'file_name', 'GDlr', 'AlphaInitial', 'EMround', 'length',
        #     'file_numF1', 'file_numF2', 'replicaF1', 'replicaF1', 'dayF1', 'dayF2', 'tokenF1', 'tokenF2', 'sampleF1', 'sampleF2','ds_prctF1','ds_prctF2'])
        # elif experiment == 5:
        #     key_file_info = (['Spearman_Corr', 'Pearson_Corr', 'ACVC', 'IM', 'GDlr', 'AlphaInitial', 'EMround', 'length',
        #     'file_numF1', 'file_numF2', 'replicaF1', 'replicaF1', 'dayF1', 'dayF2', 'tokenF1', 'tokenF2', 'sampleF1', 'sampleF2','ds_prctF1','ds_prctF2'])
        # else:
        #     # Write the header row
        #     key_file_info = (['Spearman_Corr', 'Pearson_Corr', 'ACVC', 'IM', 'GDlr', 'AlphaInitial', 'EMround', 
        #                       'tokenF1', 'tokenF2', 'sampleF1', 'sampleF2'])
        
        if experiment == 5:
            key_file_info = (['Spearman_Corr', 'Pearson_Corr', 'ACVC', 'IM', 'GDlr', 'AlphaInitial', 'EMround', 'length',
            'file_numF1', 'file_numF2', 'replicaF1', 'replicaF1', 'dayF1', 'dayF2', 'tokenF1', 'tokenF2', 'sampleF1', 'sampleF2','ds_prctF1','ds_prctF2'])
        else:
            # Write the header row
            key_file_info = (['Spearman_Corr', 'Pearson_Corr', 'ACVC', 'IM', 'file_name1','file_name2', 'GDlr', 'AlphaInitial', 'EMround', 'length',
            'file_numF1', 'file_numF2', 'replicaF1', 'replicaF1', 'dayF1', 'dayF2', 'tokenF1', 'tokenF2', 'sampleF1', 'sampleF2','ds_prctF1','ds_prctF2'])
       
            
          
        writer.writerow(key_file_info)
        # Iterate through the paired files
        for idx, pair in enumerate(paired_files):
            if experiment == 4 or experiment == 2:
                file_name1, file_name2, GDlr, AlphaInitial, EMround, length, file_num1, file_num2, replica1, replica2, dayF1, dayF2, tokenF1, tokenF2, sampleF1, sampleF2, ds_prctF1, ds_prctF2\
                    =csv_row_parsing_exp4(simulation, pair)
            elif experiment == 1:
                file_name1, file_name2, GDlr, AlphaInitial, EMround, length, file_num1, file_num2, replica1, replica2, dayF1, dayF2, tokenF1, tokenF2, sampleF1, sampleF2, ds_prctF1, ds_prctF2\
                    =csv_row_parsing_exp1(simulation, pair)
            elif experiment == 5:
                GDlr, AlphaInitial, EMround, tokenF1, tokenF2, sampleF1, sampleF2\
                    =csv_row_parsing_exp5(simulation, pair)

            #GDlr, AlphaInitial, EMround, tokenF1, tokenF2, sampleF1, sampleF2=csv_row_parsing_exp5(pair)
            
            # Perform calculations
            spearman_corr, pearson_corr = spearman_pearson_corr_generic(directory + pair[0], directory + pair[1])
            ACVC, IM = calculate_im_acvc(pair[0], pair[1], directory)

            if experiment == 5:
                row_writing = [spearman_corr, pearson_corr, ACVC, IM, GDlr, AlphaInitial, EMround, tokenF1, tokenF2, sampleF1, sampleF2]
            else:
                row_writing = [spearman_corr, pearson_corr, ACVC, IM, file_name1, file_name2, GDlr, AlphaInitial, EMround, length, file_num1, file_num2, replica1, replica2, dayF1, dayF2, tokenF1, tokenF2, sampleF1, sampleF2, ds_prctF1, ds_prctF2]
            
            # Write the row to the CSV file ##(AT)
            writer.writerow(row_writing)



# Example usage
if __name__ == "__main__":
    main()