# Code to divide JOLI prediction files into low, medium, and high isoform-abundance files

import pandas as pd
import os
import generate_result_stat as result_process
import generate_stat as stats
from util import ExperimentFileProcessor
import matplotlib.pyplot as plt
from gen_stats import get_stats
import numpy as np


def extract(data, predicted_path, output_path, gt_path, expt, gt):

    predicted = pd.read_csv(predicted_path, sep='\t')
    ground_truth = pd.read_csv(gt_path, sep='\t')

    # Calculate percentile thresholds
    percentile_33 = np.percentile(np.log(ground_truth['tpm']+1), 33.33)
    percentile_66 = np.percentile(np.log(ground_truth['tpm']+1), 66.66)
    
    print('Expt: ', expt)
    print('GT: ', gt_path.split('/')[-1])
    print(f'PERCENTILES: 33rd={percentile_33} | 66th={percentile_66}')

    # Create divisions in ground truth
    least_expressed = ground_truth[np.log(ground_truth['tpm']+1) <= percentile_33]
    moderately_expressed = ground_truth[(np.log(ground_truth['tpm']+1) > percentile_33) & 
                                        (np.log(ground_truth['tpm']+1) <= percentile_66)]
    most_expressed = ground_truth[np.log(ground_truth['tpm']+1) > percentile_66]

    # Create corresponding divisions in output file
    least_output = predicted[predicted['transcript_name'].isin(least_expressed['transcript_name'])]
    moderate_output = predicted[predicted['transcript_name'].isin(moderately_expressed['transcript_name'])]
    most_output = predicted[predicted['transcript_name'].isin(most_expressed['transcript_name'])]


    print(f'least: avg log tpm = {np.average(np.log(least_expressed["tpm"]+1))} | stdv log tpm = {np.std(np.log(least_expressed["tpm"]+1))}')
    print(f'least: num isoforms = {least_output.shape[0]}')
    print(f'moderate: avg log tpm = {np.average(np.log(moderately_expressed["tpm"]+1))} | stdv log tpm = {np.std(np.log(moderately_expressed["tpm"]+1))}')
    print(f'moderate: num isoforms = {moderate_output.shape[0]}')
    print(f'most: avg log tpm = {np.average(np.log(most_expressed["tpm"]+1))} | stdv log tpm = {np.std(np.log(most_expressed["tpm"]+1))}')
    print(f'most: num isoforms = {most_output.shape[0]}')
    print('\n')

    # Save to TSV files

    os.makedirs(output_path, exist_ok=True)

    # least_output.to_csv(os.path.join(output_path, f'least_expressed_exp{expt}_gt{gt}.tsv'), sep='\t', index=False)
    # moderate_output.to_csv(os.path.join(output_path, f'moderately_expressed_exp{expt}_gt{gt}.tsv'), sep='\t', index=False)
    # most_output.to_csv(os.path.join(output_path, f'most_expressed_exp{expt}_gt{gt}.tsv'), sep='\t', index=False)
    # print('3 TSVs saved')

def process_experiment_data(directory, experiment, data_type, file):
    # Initialize the file processor with the directory
    processor = ExperimentFileProcessor(directory)
    
    # Process files for the specified experiment and data type
    extracted_data = processor.process_files(experiment, data_type, file)
    
    # Print the extracted information
    # print(f"Extracted data for {experiment} ({data_type}):")
    return extracted_data

def process_file_arr(file_arr, directory, experiment, exp_data, groundTruth_main_dir, groundTruth_file, out_dir):
    
     for file in file_arr:
        ground_truth_arr = []
        extracted_data = process_experiment_data(directory, f"exp{experiment}", exp_data, file)
        if extracted_data:
            if extracted_data['result_aln_replica'][0] == '0':
                    # ground_truth_arr.append(groundTruth_file["PB_sample1"])
                    # ground_truth_arr.append(groundTruth_file["ill_sample1"])
                    ground_truth_arr.append(groundTruth_file["Ground_truth1"])
            if extracted_data['result_aln_replica'][0] == '2':
                    # ground_truth_arr.append(groundTruth_file["PB_sample2"])
                    # ground_truth_arr.append(groundTruth_file["ill_sample2"])
                    ground_truth_arr.append(groundTruth_file["Ground_truth2"])
            
            for ground_truth in ground_truth_arr:
                file_name = ground_truth    
                ground_truth = os.path.join(groundTruth_main_dir, ground_truth)

                
                gt = int(extracted_data['sample'])

                predicted_result = os.path.join(directory, file)
                if file_name.rsplit('.', 1)[0].rsplit('_')[0] == 'PB':
                    common_isoforms = result_process.csv_tpm_processing(file_path1=predicted_result, file_path2=ground_truth, sep_by=",")
                else:
                    common_isoforms =result_process.csv_tpm_processing(file_path1=predicted_result, file_path2=ground_truth)

                if 'tpm_quant' in common_isoforms and 'tpm_truth' in common_isoforms:
                    common_isoforms['tpm_truth'] = common_isoforms['tpm_truth'].apply(result_process.fraction_to_float_gen)

                #NEW
                if 'est_count' in common_isoforms and 'counts' in common_isoforms:
                    common_isoforms['counts'] = common_isoforms['counts'].apply(result_process.fraction_to_float_gen)

                output_path = '/'.join(directory.rstrip('/').split('/')[:-2])+out_dir

                extract(common_isoforms, predicted_result, output_path, ground_truth, experiment, gt)
                # print(f"ground_truth {file_name}")
        else:
            continue


def main():

    ######## parameters ##########
    main_result_dir = '/gpfs/commons/home/sraghavendra/Results' #'/gpfs/commons/home/atalukder/RNA_Splicing/files/results'

    experiments = ['exprmnt_2025_01_21__14_38_46', 'exprmnt_2025_01_22__14_57_21', # exp 1 LR, SR
                   'exprmnt_2025_01_21__14_41_09', # exp 4
                   'exprmnt_2025_01_22__14_58_51', # exp 2
                   'exprmnt_2025_01_22__15_00_09'] # exp 5 

    #Modify as Needed
    # experiment_file = 'exprmnt_2025_01_14__18_37_47'

    for experiment_file in experiments:
        print('EXPERIMENT: ', experiment_file)

        df_dir = '/out_dfs/'

        simulation = 1

        readme = os.path.join(main_result_dir, experiment_file, 'files/readme')
        with open(readme) as f:
            exp = int(f.readline().split(',')[1][4])

        experiment = exp  # "Different experiment setup, 1: for 1 sample, 2 for merged, 4 for multisample, 5 for merged multisample"
        groundTruth_main_dir = '/gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/ground_truths'
        groundTruth_file = {
        # "ill_sample1": "ill_sample1_gt.tsv",
        # "ill_sample2": "ill_sample2_gt.tsv",
        # "PB_sample1": "PB_sample1_gt.tsv",
        # "PB_sample2": "PB_sample1_gt.tsv",  
        "Ground_truth1": "sample1_gt.tsv",
        "Ground_truth2": "sample2_gt.tsv",
         }
        ######## parameters ##########


        predicted_theta_path_full = os.path.join(main_result_dir, experiment_file, 'files/output_files')
        file_arr  = os.listdir(predicted_theta_path_full)
        # ground_truth_path_full = os.path.join(main_dir, ground_truth_path)

        # LOCAL TESTING
        # simulation = 1
        # predicted_theta_path_full = r"C:\Users\shree\Shree\Isoform Quant\stats\exp2_MAP_same"
        # file_arr  = os.listdir(predicted_theta_path_full)
        # experiment = 2
        # groundTruth_main_dir = r'C:\Users\shree\Shree\Isoform Quant\ground_truths'
        # groundTruth_file = {
        #     # "ill_sample1": "ill_sample1_gt.tsv",
        #     # "ill_sample2": "ill_sample2_gt.tsv",
        #     # "PB_sample1": "PB_sample1_gt.tsv",
        #     # "PB_sample2": "PB_sample1_gt.tsv",  
        #     "Ground_truth1": "sample1_gt.tsv",
        #     "Ground_truth2": "sample2_gt.tsv",
        #     }
        # df_dir = "\\output_dfs\\"
        
        if simulation:
            exp_data = 'simulation'
        else:
            exp_data = 'real'
        
        process_file_arr(file_arr, predicted_theta_path_full, experiment, exp_data, groundTruth_main_dir, groundTruth_file, df_dir)
        # process_file_arr(file_arr, directory, experiment, exp_data, groundTruth_main_dir, groundTruth_file)



if __name__ == "__main__":
    main()


