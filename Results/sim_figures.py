# File to create Simulation scatter plots

import pandas as pd
import os
import generate_result_stat as result_process
import generate_stat as stats
from util import ExperimentFileProcessor
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, pearsonr
import numpy as np

# experiment_names = {1: 'Single Sample', 2: 'Merged Single Sample', 4: 'Multi-Sample', 5: 'Merged Multi-Sample'}
experiment_names = {1: 'JOLI SS', 2: 'JOLI SS (Simulation Merged)', 4: 'JOLI MS', 5: 'JOLI MS (Simulation Merged)'}
exp_len = {'long': 'LR', 'short': 'SR'}

def get_IM(data):

    data_log = np.log(data[[f'tpm_quant', f'tpm_truth']] + 1)
    data['mean_abundance'] = data_log[[f'tpm_quant', f'tpm_truth']].mean(axis=1)
    data['std_abundance'] = data_log[[f'tpm_quant', f'tpm_truth']].std(axis=1)
    CV_ig = data['std_abundance'] / data['mean_abundance']
    data['CV'] = CV_ig

    CV_ig_squared = CV_ig ** 2
    IM = np.sqrt(CV_ig_squared.mean())

    sorted_u_ig = data['mean_abundance'].sort_values()
    sorted_CV_ig = data['CV'].loc[sorted_u_ig.index]

    ACVC = np.trapz(sorted_CV_ig, x=sorted_u_ig)
    # result = f"IM {IM}, ACVC {ACVC}\n"
    # # with open(log_file, 'a') as f:
    # #     f.write(result)
    # print(result)

    return ACVC, IM

def plot_predictedVSgroundtruth_theta(common_isoforms, figure_path, experiment, sub_text, gt, l1, l2, d1, d2, color=True):

    predicted_theta = common_isoforms['tpm_quant']
    ground_truth_theta = common_isoforms['tpm_truth']

    x = common_isoforms['tpm_quant']+1
    y = common_isoforms['tpm_truth']+1

    log_sc, p_value = spearmanr(np.log(x), np.log(y))
    log_pc, p_value = pearsonr(np.log(x), np.log(y))
    acvc, im = get_IM(common_isoforms)
    
    # info_text+=f'\nlog-scale SC: {log_sc:.4f}'
    info_text = f'SC: {log_sc:.3f}\nPC: {log_pc:.3f}'#\nIM: {im:.4f}'

    # Calculate the range of data
    min_val = min(predicted_theta.min(), ground_truth_theta.min())
    max_val = max(predicted_theta.max(), ground_truth_theta.max())

    # Add some margin to the limits for better visualization
    margin = (max_val - min_val) * 0.05  # 5% margin
    lower_limit = min_val - margin
    upper_limit = max_val + margin

    plt.rcParams.update({'font.size': 17})  # Sets a base font size
    plt.rcParams['axes.titlesize'] = 26  # Sets the title font size
    plt.rcParams['axes.labelsize'] = 17  # Sets the axes labels font size

    # Plot
    plt.figure(figsize=(8, 6))
    plt.scatter(x, y, alpha=0.7)
    plt.xscale('log')
    plt.yscale('log')
    plt.plot([lower_limit, upper_limit], [lower_limit, upper_limit], linestyle='-', label=info_text) # Diagonal Line
    plt.xlim(lower_limit, upper_limit)
    plt.ylim(lower_limit, upper_limit)

    #If Colors are needed
    if color:
        plt.hexbin(x, y, 
                gridsize=100,
                cmap='jet',
                bins='log',
                xscale='log',
                yscale='log')
        
        plt.colorbar(label='counts')
    
    lengths = {1: l1, 2: l2}
    exp_length = lengths[gt]

    if experiment==1:
        exp_length = l1

    xlabel = False
    ylabel = False
    
    if experiment in [4,5]:
        xlabel=True
    if exp_length=='long' and experiment in [1,4]:
        ylabel=True


    if experiment in [1,4]:
        title = f'{experiment_names[experiment]} (Simulation {exp_len[exp_length]})'
    else:
        title = f'{experiment_names[experiment]}'


    plt.title(title)
    if xlabel:
        plt.xlabel('Estimated Abundances log(tpm+1)')
    if ylabel:
        plt.ylabel('Ground Truth Abundances log(tpm+1)')
    plt.legend(loc='upper left')
    plt.grid(True)

    # if sub_text:
    #     plt.annotate(sub_text, xy=(0.5, -0.15), xycoords='axes fraction', ha='center', va='center')
    
    #Check
    plt.tight_layout()

    # Create filename using samples, gt, and experiment number
    samples_info = []
    for sample in info_text.split('\n'):
        if 'Sample' in sample or 'Day' in sample:
            day = sample.split('Day ')[1].split(',')[0]
            length = sample.split(', ')[1].split(' ')[0].lower()
            # if length == 'short':
            #     length = 'short'
            # elif length == 'long':
            #     length = 'long'
            # else:
            #     length = 'mixed'
            samples_info.append(f"day{day}{length}")
    
    samples_info = "_".join(samples_info)
    
    
    #gt_info = sub_text.split('Sample ')[1].split(' ')[0] if sub_text else ''
    filename = f"exp{experiment}_{d1}{l1}_{d2}{l2}_gt{gt}_isoform_estimate_vs_truth.png"

    # Save or show plot
    os.makedirs(figure_path, exist_ok=True)
    plt.savefig(os.path.join(figure_path, filename))  # Save as a file
    print('saved')
    # plt.show()

# def plot_predictedVSgroundtruth_theta(common_isoforms, figure_path, experiment, info_text, sub_text, gt):
#     predicted_theta = common_isoforms['est_count']
#     ground_truth_theta = common_isoforms['counts']

#     x = common_isoforms['est_count']+1
#     y = common_isoforms['counts']+1

#     min_val = min(x.min(), y.min())
#     max_val = max(x.max(), y.max())

#     margin_factor = 0.1  # 10% margin
#     min_lim = min_val * (1 - margin_factor)
#     max_lim = max_val * (1 + margin_factor)

#     plt.figure(figsize=(8, 8))
    
#     # Plot the diagonal line with metrics in legend
#     #plt.plot([1e0, 1e5], [1e0, 1e5], label=info_text)
#     plt.plot([min_val, max_val], [min_val, max_val], label=info_text)
    
#     # Create hexbin plot with specific parameters
#     plt.hexbin(x, y, 
#                gridsize=100,
#                cmap='jet',
#                bins='log',
#                xscale='log',
#                yscale='log')
    
#     # Customize plot appearance
#     plt.colorbar(label='counts')
#     plt.xlim(min_val, max_val)
#     plt.ylim(min_val, max_val)
#     # plt.xlim(1e0, 1e5)
#     # plt.ylim(1e0, 1e5)
#     plt.grid(True) #Added grid

#     plt.xlabel("Estimated Abundances")
#     plt.ylabel("Ground Truth Abundances")
#     plt.title(f"{experiment_names[experiment]} Experiment")
#     plt.legend(loc='upper left')

#     # Create filename using samples, gt, and experiment number
#     samples_info = []
#     for sample in info_text.split('\n'):
#         if 'Sample' in sample or 'Day' in sample:
#             day = sample.split('Day ')[1].split(',')[0]
#             length = sample.split(', ')[1].split(' ')[0].lower()
#             samples_info.append(f"day{day}{length}")
    
#     samples_info = "_".join(samples_info)
    
    
#     #gt_info = sub_text.split('Sample ')[1].split(' ')[0] if sub_text else ''
#     filename = f"exp{experiment}_{samples_info}_gt{gt}_isoform_estimate_vs_truth.png"

#     # Save or show plot
#     os.makedirs(figure_path, exist_ok=True)
#     plt.savefig(os.path.join(figure_path, filename), dpi=300)  # Save as a file
#     print('saved')
#     # plt.show()



def process_experiment_data(directory, experiment, data_type, file):
    # Initialize the file processor with the directory
    processor = ExperimentFileProcessor(directory)
    
    # Process files for the specified experiment and data type
    extracted_data = processor.process_files(experiment, data_type, file)
    
    # Print the extracted information
    print(f"Extracted data for {experiment} ({data_type}):")
    return extracted_data

def process_file_arr(file_arr, directory, experiment, exp_data, groundTruth_main_dir, groundTruth_file, fig_dir, color):
    
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

                sc, pc, im = stats.get_stats(os.path.join(directory, file), ground_truth)

                # day1 = extracted_data['aln_replica1'][0]
                # day2 = extracted_data['aln_replica2'][0]
                # length1 = extracted_data['length1']
                # length2 = extracted_data['length2']
                
                gt = int(extracted_data['sample'])

                if experiment==1:
                    
                    day1 = extracted_data['aln_replica1'][0]
                    day2 = None
                    length1 = extracted_data['length1']
                    length2 = None

                    # info_text = f'Sample: Day {day1}, {length1} reads\nSC: {sc:.4f}\nPC: {pc:.4f}\nIM: {im:.4f}'
                    sub_text = None

                if experiment==2:
                    day1 = extracted_data['aln_replica1'][0]
                    day2 = extracted_data['aln_replica2'][0]
                    length1 = extracted_data['length1']
                    length2 = extracted_data['length2']

                    # info_text = f'Sample 1: Day {day1}, {length1} reads\nSample 2: Day {day2}, {length2} reads\nSC: {sc:.4f}\nPC: {pc:.4f}\nIM: {im:.4f}'
                    sub_text = f'Compared with Sample {gt} Ground Truth'

                if experiment==4:
                    day1 = extracted_data['aln_replica1'][0]
                    day2 = extracted_data['aln_replica2'][0]
                    length1 = extracted_data['length1']
                    length2 = extracted_data['length2']

                    # info_text = f'Sample 1: Day {day1}, {length1} reads\nSample 2: Day {day2}, {length2} reads\nSC: {sc:.4f}\nPC: {pc:.4f}\nIM: {im:.4f}'
                    sub_text = f'Compared with Sample {gt} Ground Truth'

                if experiment==5:
                    day1 = 0
                    day2 = 2
                    length1 = 'shortlong'
                    length2 = 'shortlong'    

                    # info_text = f'Sample 1: Day {day1}, {length1} reads\nSample 2: Day {day2}, {length2} reads\nSC: {sc:.4f}\nPC: {pc:.4f}\nIM: {im:.4f}'
                    sub_text = f'Compared with Sample {gt} Ground Truth'

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

                figure_path = '/'.join(directory.rstrip('/').split('/')[:-2])+fig_dir
                plot_predictedVSgroundtruth_theta(common_isoforms, figure_path, experiment, sub_text, gt, length1, length2, day1, day2, color)
                print(f"ground_truth {file_name}")
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

        fig_dir = '/figures/'
        color = 1

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

        if simulation:
            exp_data = 'simulation'
        else:
            exp_data = 'real'
        
        process_file_arr(file_arr, predicted_theta_path_full, experiment, exp_data, groundTruth_main_dir, groundTruth_file, fig_dir, color)
        # process_file_arr(file_arr, directory, experiment, exp_data, groundTruth_main_dir, groundTruth_file)


if __name__ == "__main__":
    main()