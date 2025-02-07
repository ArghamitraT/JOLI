# Code to create figures comparing correlation and error between low, medium, and high expression isoforms

import pandas as pd
import os
import generate_result_stat as result_process
import generate_stat as stats
from util import ExperimentFileProcessor
import matplotlib.pyplot as plt
from gen_stats import get_stats
import numpy as np

def plot_bar_graph(s1, s2, models, fig_name, e1, e2, measure, xlabel, ylabel, title, legend):
    # Data preparation
    
    metrics = ['sc', 'pc', 'im']

    # Set width of bars and positions of the bars
    bar_width = 0.15
    r1 = np.arange(0, len(models)/2, 0.5)
    r2 = [x + bar_width for x in r1]

    # Create the plot
    plt.figure(figsize=(4, 3))

    plt.rcParams['axes.titlesize'] = 16  # Sets the title font size

    # plt.rcParams.update({'font.size': 14})  # Sets a base font size
    # plt.rcParams['axes.titlesize'] = 16  # Sets the title font size
    # plt.rcParams['axes.labelsize'] = 14  # Sets the axes labels font size


    # Create bars
    plt.bar(r1, s1, width=bar_width, label=e1, color='salmon')
    plt.bar(r2, s2, width=bar_width, label=e2, color='lightseagreen')

    # Customize the plot
    if xlabel:
        plt.xlabel('Isoform Expression Level')
    if ylabel:
        plt.ylabel(measure)
    if title:
        plt.title(title)
    #plt.xticks([r + bar_width for r in range(len(models))], models)
    plt.xticks([r + bar_width for r in r1], models)
    if legend:
        plt.legend()

    # Ensure y-axis starts from 0 and ends at 1
    plt.ylim(0, 1)
    if measure=='MRD':
         plt.ylim(0, 2.2)

    # Add grid for better readability
    plt.grid(True, axis='y', linestyle='--', alpha=0.7)

    # for i, v in enumerate(s1):
    #     plt.text(r1[i], v + 0.02, str(round(v, 2)), ha='center', va='bottom')

    # for i, v in enumerate(s2):
    #     plt.text(r2[i], v + 0.02, str(round(v, 2)), ha='center', va='bottom')

    plt.tight_layout()
    plt.savefig(fig_name)  # Save as a file
    print('saved')

def get_statistics(arr, dir):
        groundTruth_main_dir = '/gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/ground_truths'
        groundTruth_file = { 
            "gt1": "sample1_gt.tsv",
            "gt2": "sample2_gt.tsv",
            }
        
        # label_dict = {'least':'le', 'moderately': 'md', 'most': 'mo' }

        sc = []
        pc = []
        nrmse = []
        mrd = []
        for file in arr:
            gt = file.split('.')[0].split('_')[-1]
            stats = get_stats(os.path.join(dir,file), os.path.join(groundTruth_main_dir, groundTruth_file[gt]), 
                           'transcript_name', 'transcript_name', 'tpm', 'tpm', extra=True)
            # label = file.split('_')[0]
            sc.append(stats['sc'])
            pc.append(stats['pc'])
            nrmse.append(stats['nrmse'])
            mrd.append(stats['mrd'])

        return sc, pc, nrmse, mrd

def main():
    # Day 5 is best performing

    # Exp 1
    exp1_dir_lr = '/gpfs/commons/home/sraghavendra/Results/exprmnt_2025_01_21__14_38_46'
    exp1_dir_sr = '/gpfs/commons/home/sraghavendra/Results/exprmnt_2025_01_22__14_57_21'
    exp2_dir = '/gpfs/commons/home/sraghavendra/Results/exprmnt_2025_01_22__14_58_51'
    exp4_dir = '/gpfs/commons/home/sraghavendra/Results/exprmnt_2025_01_21__14_41_09'
    exp5_dir = '/gpfs/commons/home/sraghavendra/Results/exprmnt_2025_01_22__15_00_09'

    df_dir = "out_dfs"

    exp4_file_arr_lr = []
    exp4_file_arr_sr = []

    for file in os.listdir(os.path.join(exp4_dir, df_dir)):
        if file.split('.')[0].split('_')[-1][-1]=='1':
              exp4_file_arr_lr.append(file)
        elif file.split('.')[0].split('_')[-1][-1]=='2':
              exp4_file_arr_sr.append(file)
        else:
             print('ERROR - exp 4 file name format wrong')

    exp1_file_arr_lr = os.listdir(os.path.join(exp1_dir_lr, df_dir))
    #exp4_file_arr_lr = os.listdir(os.path.join(exp4_dir_lr, df_dir))

    exp1_file_arr_sr = os.listdir(os.path.join(exp1_dir_sr, df_dir))
    #exp4_file_arr_sr = os.listdir(os.path.join(exp4_dir_sr, df_dir))

    exp2_file_arr = os.listdir(os.path.join(exp2_dir, df_dir))
    exp5_file_arr = os.listdir(os.path.join(exp5_dir, df_dir))

    #Stats
    
    exp1_stats_lr_sc, exp1_stats_lr_pc, exp1_stats_lr_nrmse, exp1_stats_lr_mrd = get_statistics(exp1_file_arr_lr, os.path.join(exp1_dir_lr, df_dir))
    exp4_stats_lr_sc, exp4_stats_lr_pc, exp4_stats_lr_nrmse, exp4_stats_lr_mrd  = get_statistics(exp4_file_arr_lr, os.path.join(exp4_dir, df_dir))

    exp1_stats_sr_sc, exp1_stats_sr_pc, exp1_stats_sr_nrmse, exp1_stats_sr_mrd = get_statistics(exp1_file_arr_sr, os.path.join(exp1_dir_sr, df_dir))
    exp4_stats_sr_sc, exp4_stats_sr_pc, exp4_stats_sr_nrmse, exp4_stats_sr_mrd = get_statistics(exp4_file_arr_sr, os.path.join(exp4_dir, df_dir))    

    exp2_stats_sc, exp2_stats_pc, exp2_stats_nrmse, exp2_stats_mrd = get_statistics(exp2_file_arr, os.path.join(exp2_dir, df_dir))
    exp5_stats_sc, exp5_stats_pc, exp5_stats_nrmse, exp5_stats_mrd = get_statistics(exp5_file_arr, os.path.join(exp5_dir, df_dir))

    models = ['Low', 'Medium', 'High']
    
    # exp1 exp4 LR bar plot SC
    fig_name = 'e1_e4_lr_div_sc.png'
    plot_bar_graph(s1=exp1_stats_lr_sc, s2=exp4_stats_lr_sc, models=models, fig_name=fig_name, 
                   e1='JOLI SS', e2='JOLI MS', measure='Spearman Correlation',
                   xlabel=False, ylabel=True, title='Simulation LR', legend=True)

    # exp1 exp4 SR bar plot SC
    fig_name = 'e1_e4_sr_div_sc.png'
    plot_bar_graph(s1=exp1_stats_sr_sc, s2=exp4_stats_sr_sc, models=models, fig_name=fig_name, 
                   e1='JOLI SS', e2='JOLI MS', measure='Spearman Correlation',
                   xlabel=False, ylabel=False, title='Simulation SR', legend=False)

    # exp 2 exp 5 bar plot SC
    fig_name = 'e2_e5_div_sc.png'
    plot_bar_graph(s1=exp2_stats_sc, s2=exp5_stats_sc, models=models, fig_name=fig_name, 
                   e1='JOLS SS', e2='JOLI MS', measure='Spearman Correlation',
                   xlabel=False, ylabel=False, title='Simulation Merged', legend=False)

    # exp1 exp4 LR bar plot PC
    fig_name = 'e1_e4_lr_div_pc.png'
    plot_bar_graph(s1=exp1_stats_lr_pc, s2=exp4_stats_lr_pc, models=models, fig_name=fig_name, 
                   e1='JOLI SS', e2='JOLI MS', measure='Pearson Correlation',
                   xlabel=False, ylabel=True, title=None, legend=False)

    # exp1 exp4 SR bar plot PC
    fig_name = 'e1_e4_sr_div_pc.png'
    plot_bar_graph(s1=exp1_stats_sr_pc, s2=exp4_stats_sr_pc, models=models, fig_name=fig_name, 
                   e1='JOLI SS', e2='JOLI MS', measure='Pearson Correlation',
                   xlabel=False, ylabel=False, title=None, legend=False)

    # exp 2 exp 5 bar plot PC
    fig_name = 'e2_e5_div_pc.png'
    plot_bar_graph(s1=exp2_stats_pc, s2=exp5_stats_pc, models=models, fig_name=fig_name, 
                   e1='JOLS SS', e2='JOLI MS', measure='Pearson Correlation',
                   xlabel=False, ylabel=False, title=None, legend=False)

    # exp1 exp4 LR bar plot MRD
    fig_name = 'e1_e4_lr_div_mrd.png'
    plot_bar_graph(s1=exp1_stats_lr_mrd, s2=exp4_stats_lr_mrd, models=models, fig_name=fig_name, 
                   e1='JOLI SS', e2='JOLI MS', measure='MRD',
                   xlabel=True, ylabel=True, title=None, legend=False)

    # exp1 exp4 SR bar plot MRD
    fig_name = 'e1_e4_sr_div_mrd.png'
    plot_bar_graph(s1=exp1_stats_sr_mrd, s2=exp4_stats_sr_mrd, models=models, fig_name=fig_name, 
                   e1='JOLI SS', e2='JOLI MS', measure='MRD',
                   xlabel=True, ylabel=False, title=None, legend=False)

    # exp 2 exp 5 bar plot MRD
    fig_name = 'e2_e5_div_mrd.png'
    plot_bar_graph(s1=exp2_stats_mrd, s2=exp5_stats_mrd, models=models, fig_name=fig_name, 
                   e1='JOLI SS', e2='JOLI MS', measure='MRD',
                   xlabel=True, ylabel=False, title=None, legend=False)


if __name__ == "__main__":
    main()