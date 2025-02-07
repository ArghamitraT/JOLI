#Code to generate SOTA comparison figures

import pandas as pd
import os
import generate_result_stat as result_process
import generate_stat as stats
from util import ExperimentFileProcessor
import matplotlib.pyplot as plt
from gen_stats import get_stats
import numpy as np


def save_fig(values, color, models, fig_name, measure, l, title, xlabel, ylabel):
    measures = {'sc': 'Spearman Correlation', 'pc': 'Pearson Correlation', 'im': 'IM', 'nrmse': 'NRMSE', 'mrd' :'MRD'}

    bar_width = 0.5 
    r1 = np.arange(len(models)) 
    bar_colors = ['salmon', 'lightseagreen', 'skyblue', 'thistle'] #lightcoral
    if l=='sr':
        r1 = [0, 1.5, 3]
        bar_colors = ['salmon', 'lightseagreen', 'skyblue']

    # plt.rcParams['xtick.labelsize'] = 8  # Sets the x-axis font size
    plt.rcParams.update({'font.size': 12})  # Sets a base font size
    plt.rcParams['axes.titlesize'] = 18  # Sets the title font size
    plt.rcParams['axes.labelsize'] = 14  # Sets the axes labels font size

    plt.figure(figsize=(4, 3))
    plt.bar(r1, values, width=bar_width, color=bar_colors) #, edgecolor='black')  # Remove fill color, add border
    if xlabel:
        plt.xlabel('Models')
    if ylabel:
        plt.ylabel(measures[measure])
    if title:
        plt.title(title)
    plt.xticks(r1, models)

    plt.ylim(0, 1)
    if measure=='im' or measure=='mrd':
        plt.ylim(0, 1.2)

    plt.grid(True, axis='y', linestyle='--', alpha=0.7)

    for i, v in enumerate(values):
        plt.text(r1[i], v - 0.15, str(round(v, 2)), ha='center', va='bottom')
    
    plt.tight_layout()
    plt.savefig(fig_name + f'_{measure}.png')
    print('saved')


def plot_bar_graph(stats, models, fig_name, title, l, ylabel, simulation=False):
    # Data preparation
    
    metrics = ['sc', 'pc', 'im']

    # Organize data into arrays for each metric
    # sc_values = [stats['kl_stats']['sc'], stats['sal_stats']['sc'], stats['e1_sstats']['sc'], stats['e4_sstats']['sc']]
    # pc_values = [stats['kl_stats']['pc'], stats['sal_stats']['pc'], stats['e1_sstats']['pc'], stats['e4_sstats']['pc']]
    # im_values = [stats['kl_stats']['im'], stats['sal_stats']['im'], stats['e1_sstats']['im'], stats['e4_sstats']['im']]

    sc_values = [stats[k]['sc'] for k in stats.keys()]
    pc_values = [stats[k]['pc'] for k in stats.keys()]
    im_values = [stats[k]['im'] for k in stats.keys()]
    if simulation:
        nrmse_values = [stats[k]['nrmse'] for k in stats.keys()]
        mrd_values = [stats[k]['mrd'] for k in stats.keys()]


    save_fig(sc_values, 'skyblue', models, fig_name, 'sc', l, title, xlabel=False, ylabel=ylabel)
    if simulation:
        save_fig(pc_values, 'lightgreen', models, fig_name, 'pc',l, None, xlabel=False, ylabel=ylabel)
    else:
        save_fig(pc_values, 'lightgreen', models, fig_name, 'pc',l, None, xlabel=False, ylabel=ylabel)
    if not simulation:
        save_fig(im_values, 'salmon', models, fig_name, 'im',l, None, xlabel=True, ylabel=ylabel)
    else:
        #save_fig(nrmse_values, 'salmon', models, fig_name, 'nrmse', l, None, xlabel=True, ylabel=ylabel)
        save_fig(mrd_values, 'salmon', models, fig_name, 'mrd', l, None, xlabel=True, ylabel=ylabel)

def main():
    # Day 5 is best performing

    # Short Read Files
    kl_r1 = '/gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/output_short_51/abundance.tsv' # kallisto
    kl_r2 = '/gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/output_short_52/abundance.tsv'
    sal_r1 = '/gpfs/commons/home/sraghavendra/SOTA/salmon/salmon_quant51/quant.sf' # salmon
    sal_r2 = '/gpfs/commons/home/sraghavendra/SOTA/salmon/salmon_quant52/quant.sf'
    e1_sr1 = '../exprmnt_2025_01_21__17_21_14/files/output_files/output_PacIllu_VIGD_token_470401_sample1_file1_ds100num1aln51short_GDlr_0.01_AlphaInitial_1.0_EMround_30_2025_1_21_06_40_34.tsv' # exp 1
    e1_sr2 = '../exprmnt_2025_01_21__17_21_14/files/output_files/output_PacIllu_VIGD_token_5426678_sample1_file1_ds100num1aln52short_GDlr_0.01_AlphaInitial_1.0_EMround_30_2025_1_21_11_33_28.tsv'
    e4_sr1 = '../exprmnt_2025_01_21__17_25_37/files/output_files/output_PacIllu_VIGD_token_24785461_sample2_file1_ds10num1aln01long_file2_ds100num1aln51short_GDlr_0.01_AlphaInitial_10000.0_EMround_30_2025_1_22_14_10_05.tsv' # exp 4
    e4_sr2 = '../exprmnt_2025_01_21__17_25_37/files/output_files/output_PacIllu_VIGD_token_1196578_sample2_file1_ds10num1aln02long_file2_ds100num1aln52short_GDlr_0.01_AlphaInitial_10000.0_EMround_30_2025_1_22_14_38_40.tsv'

    # Long Read Files
    lrkl_r1 = '/gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/output_51/abundance.tsv' #l r-kallisto
    lrkl_r2 = '/gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/output_52/abundance.tsv'
    nc_r1 = '/gpfs/commons/home/sraghavendra/SOTA/NanoCount/tx_counts_ds_51.tsv' # NanoCount
    nc_r2 = '/gpfs/commons/home/sraghavendra/SOTA/NanoCount/tx_counts_ds_52.tsv'
    oar_r1 = '/gpfs/commons/home/sraghavendra/SOTA/oarfish/day5_rep1/.quant' # Oarfish
    oar_r2 = '/gpfs/commons/home/sraghavendra/SOTA/oarfish/day5_rep2/.quant'
    # mp_r1 = '' # MPAQT
    # mp_r2 = ''
    e1_lr1 = '../exprmnt_2025_01_21__22_53_07/files/output_files/output_PacIllu_VIGD_token_10807412_sample1_file1_ds10num1aln51long_GDlr_0.01_AlphaInitial_1.0_EMround_30_2025_1_21_23_27_48.tsv' #exp 1
    e1_lr2 = '../exprmnt_2025_01_21__22_53_07/files/output_files/output_PacIllu_VIGD_token_32391888_sample1_file1_ds10num1aln52long_GDlr_0.01_AlphaInitial_1.0_EMround_30_2025_1_21_23_44_03.tsv'
    e4_lr1 = '../exprmnt_2025_01_21__17_25_37/files/output_files/output_PacIllu_VIGD_token_19932137_sample1_file1_ds10num1aln51long_file2_ds100num1aln01short_GDlr_0.01_AlphaInitial_10000.0_EMround_30_2025_1_22_08_40_46.tsv' #exp 4
    e4_lr2 = '../exprmnt_2025_01_21__17_25_37/files/output_files/output_PacIllu_VIGD_token_16541151_sample1_file1_ds10num1aln52long_file2_ds100num1aln02short_GDlr_0.01_AlphaInitial_10000.0_EMround_30_2025_1_22_09_00_04.tsv'

    # Simulation LR Files
    gt1 = '/gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/ground_truths/sample1_gt.tsv'
    gt2 = '/gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/sim_real_data/ground_truths/sample2_gt.tsv'
    kl1 = '/gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/output_sim1/abundance.tsv'
    kl2 = '/gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/output_sim2/abundance.tsv'
    oar1 = '/gpfs/commons/home/sraghavendra/SOTA/oarfish/sim2/.quant'
    oar2 = '/gpfs/commons/home/sraghavendra/SOTA/oarfish/sim2/.quant'
    kl_s1 = '/gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/output_short_sim1/abundance.tsv'
    kl_s2 = '/gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/output_short_sim2/abundance.tsv'
    nc1 = '/gpfs/commons/home/sraghavendra/SOTA/NanoCount/tx_counts_new_sim_1.tsv'
    nc2 = '/gpfs/commons/home/sraghavendra/SOTA/NanoCount/tx_counts_new_sim_2.tsv'
    e1_lr = '/gpfs/commons/home/sraghavendra/Results/exprmnt_2025_01_21__14_38_46/files/output_files/output_Simulation_VIGD_token_12497384_sample1_file1_ds100num1aln21long_GDlr_0.01_AlphaInitial_1.0_EMround_30_2025_1_21_14_42_19.tsv'
    e1_sr = '/gpfs/commons/home/sraghavendra/Results/exprmnt_2025_01_22__14_57_21/files/output_files/output_Simulation_VIGD_token_12495448_sample1_file1_ds100num1aln01short_GDlr_0.01_AlphaInitial_1.0_EMround_30_2025_1_22_16_18_38.tsv'
    e4_lr = '/gpfs/commons/home/sraghavendra/Results/exprmnt_2025_01_21__14_41_09/files/output_files/output_Simulation_VIGD_token_33423434_sample1_file1_ds100num1aln21long_file2_ds100num1aln01short_GDlr_0.01_AlphaInitial_10000.0_EMround_30_2025_1_21_16_24_36.tsv'
    e4_sr = '/gpfs/commons/home/sraghavendra/Results/exprmnt_2025_01_21__14_41_09/files/output_files/output_Simulation_VIGD_token_22767984_sample2_file1_ds100num1aln01long_file2_ds100num1aln21short_GDlr_0.01_AlphaInitial_10000.0_EMround_30_2025_1_21_16_25_37.tsv'

    ############################ REAL DATA RESULTS ################################



    print('\n###################### REAL DATA ######################\n')

    ######### SR ############

    print('\n Short Reads \n')

    kl_stats = get_stats(kl_r1, kl_r2, 'transcript_id', 'transcript_id', 'tpm', 'tpm')
    print('Real Kallisto SR: ', kl_stats)

    # sal_stats = get_stats(sal_r1, sal_r2, 'Name', 'Name', 'TPM', 'TPM')
    # # print(sal_stats)

    e1_sstats = get_stats(e1_sr1, e1_sr2, 'transcript_name', 'transcript_name', 'tpm', 'tpm')
    print('Real Exp1 SR: ', e1_sstats)

    e4_sstats = get_stats(e4_sr1, e4_sr2, 'transcript_name', 'transcript_name', 'tpm', 'tpm')
    print('Real Exp 4 SR: ', e4_sstats)

    ######### LR #########

    print('\n Long Reads \n')

    lrkl_stats = get_stats(lrkl_r1, lrkl_r2, 'transcript_id', 'transcript_id', 'tpm', 'tpm')
    print('Real Kallisto LR: ', lrkl_stats)

    # nc_stats = get_stats(nc_r1, nc_r2, 'transcript_name', 'transcript_name', 'tpm', 'tpm')
    # # print(nc_stats)

    oar_stats = get_stats(oar_r1, oar_r2, 'tname', 'tname', 'num_reads', 'num_reads')
    print('Real Oar: ', oar_stats)

    e1_lstats = get_stats(e1_lr1, e1_lr2, 'transcript_name', 'transcript_name', 'tpm', 'tpm')
    print('Real Exp 1 LR: ', e1_lstats)

    e4_lstats = get_stats(e4_lr1, e4_lr2, 'transcript_name', 'transcript_name', 'tpm', 'tpm')
    print('Real Exp 4 LR: ', e4_lstats)


    print('\n***************************************************\n')

    print('\n###################### Sim DATA ######################\n')

    ### SIMULATION STATS ###
    
    ####### SR #########

    print('\n Short Reads \n')

    # SR Kallisto STATS
    sim_kl1_sstats = get_stats(kl_s2, gt2, 'transcript_id', 'transcript_name', 'tpm', 'tpm', extra=True)
    print('Sim SR kallisto: ', sim_kl1_sstats)

    #Expt 1 SR STATS
    sim_e1sstats = get_stats(e1_sr, gt1, 'transcript_name', 'transcript_name', 'tpm', 'tpm', extra=True)
    print('Exp 1 SR: ', sim_e1sstats)

    #Expt 4 SR STATS
    sim_e4sstats = get_stats(e4_sr, gt2, 'transcript_name', 'transcript_name', 'tpm', 'tpm', extra=True)
    print('Exp 4 SR: ', sim_e4sstats)

    ######### LR #########

    print('\n Long Reads \n')

    #NC stats
    # nc_stats = get_stats(nc1, gt1, 'transcript_name', 'transcript_name', 'tpm', 'tpm')
    # print('Sim NC results: ', nc_stats)

    # LR Kallisto Sim Stats
    sim_kl1_lstats = get_stats(kl1, gt1, 'transcript_id', 'transcript_name', 'tpm', 'tpm', extra=True)
    print('Sim lr-kallisto: ', sim_kl1_lstats)

    # Oarfish Sim LR Stats
    sim_oar1_stats = get_stats(oar1, gt1, 'tname', 'transcript_name', 'num_reads', 'tpm', extra=True)
    print('Sim Oarfish: ', sim_oar1_stats)

    #Expt 1 LR STATS
    sim_el1sstats = get_stats(e1_lr, gt2, 'transcript_name', 'transcript_name', 'tpm', 'tpm', extra=True)
    print('Exp 1 LR: ', sim_el1sstats)

    #Expt 4 LR STATS
    sim_el4sstats = get_stats(e4_lr, gt2, 'transcript_name', 'transcript_name', 'tpm', 'tpm', extra=True)
    print('Exp 4 LR: ', sim_el4sstats)

    print('\n')

    ##### REAL PLOTS ####

    # SR Bar Plot
    models = ['JOLI SS', 'JOLI MS', 'kallisto']
    stats = {'e1_sstats': e1_sstats, 'e4_sstats': e4_sstats, 'kl_stats': kl_stats}
    fig_name = 'sr_real_res'
    plot_bar_graph(stats, models, fig_name, 'Illumina SR', l='sr', ylabel=False)

    # LR Bar Plot
    models = ['JOLI SS', 'JOLI MS', 'lr-kallisto', 'Oarfish']
    stats = {'e1_lstats': e1_lstats, 'e4_lstats': e4_lstats, 'lrkl_stats': lrkl_stats, 'oar_stats': oar_stats}
    fig_name = 'lr_real_res'
    plot_bar_graph(stats, models, fig_name, 'PacBio LR', l='lr', ylabel=True)

    #### Sim Bar Plots #######

    # SR Bar Plot
    models = ['JOLI SS', 'JOLI MS', 'kallisto']
    stats = {'sim_e1sstats': sim_e1sstats, 'sim_e4sstats': sim_e4sstats, 'sim_kl1_sstats': sim_kl1_sstats}
    fig_name = 'sr_sim_res'
    plot_bar_graph(stats, models, fig_name, 'Simulation SR', l='sr', ylabel=False, simulation=True)

    # LR Bar Plot
    models = ['JOLI SS', 'JOLI MS', 'lr-kallisto', 'Oarfish']
    stats = {'sim_el1sstats': sim_el1sstats, 'sim_el4sstats': sim_el4sstats, 'sim_kl1_lstats': sim_kl1_lstats, 'sim_oar1_stats': sim_oar1_stats}
    fig_name = 'lr_sim_res'
    plot_bar_graph(stats, models, fig_name, 'Simulation LR', l='lr', ylabel=True, simulation=True)

if __name__ == "__main__":
    main()