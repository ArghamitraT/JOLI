# File to get all required statistics - SC, PC, IM, MRD

import numpy as np
from scipy.stats import spearmanr, pearsonr
import pandas as pd
import datetime
import time
import os
import re
from collections import defaultdict
from scipy.interpolate import UnivariateSpline

def fraction_to_float_gen(value):
    try:
        return float(value)
    except ValueError:
        return None

def split_transcript_id(transcript):
    return transcript.split('.')[0]

def csv_tpm_processing(file_path1, file_path2, tr1, tr2, suffixes=('_quant', '_truth')):
    # Load the datasets
    our_quant = pd.read_csv(file_path1, sep="\t")
    ground_truth = pd.read_csv(file_path2, sep="\t")

    # Apply the function to the 'transcript_id' column of your DataFrame
    our_quant[tr1] = our_quant[tr1].apply(split_transcript_id)
    ground_truth[tr1] = ground_truth[tr2].apply(split_transcript_id)

    # Find common isoforms
    common_isoforms = pd.merge(our_quant, ground_truth, on=tr1, suffixes=suffixes)

    return common_isoforms

def spearman_corr_generic(file_path1, file_path2, tr1, tr2, tpm1, tpm2):
    common_isoforms = csv_tpm_processing(file_path1, file_path2, tr1, tr2)

    if f'{tpm1}_quant' in common_isoforms and f'{tpm2}_truth' in common_isoforms:
        common_isoforms[f'{tpm1}_truth'] = common_isoforms[f'{tpm2}_truth'].apply(fraction_to_float_gen)
        correlation, p_value = spearmanr(np.log(common_isoforms[f'{tpm1}_quant']+1), np.log(common_isoforms[f'{tpm2}_truth']+1))
        pcorrelation, p_value = pearsonr(np.log(common_isoforms[f'{tpm1}_quant']+1), np.log(common_isoforms[f'{tpm2}_truth']+1))

        return round(correlation, 5), round(pcorrelation, 5)
    
    elif tpm1!=tpm2:
        common_isoforms[f'{tpm2}'] = common_isoforms[f'{tpm2}'].apply(fraction_to_float_gen)
        correlation, p_value = spearmanr(np.log(common_isoforms[f'{tpm1}']+1), np.log(common_isoforms[f'{tpm2}']+1))
        pcorrelation, p_value = pearsonr(np.log(common_isoforms[f'{tpm1}']+1), np.log(common_isoforms[f'{tpm2}']+1))

        return round(correlation, 5), round(pcorrelation, 5)
    
def calculate_cv(data, tpm1, tpm2):

    if f'{tpm1}_Rep1' in data and f'{tpm2}_Rep2' in data:
        data_log = np.log(data[[f'{tpm1}_Rep1', f'{tpm2}_Rep2']] + 1)
        data['mean_abundance'] = data_log[[f'{tpm1}_Rep1', f'{tpm2}_Rep2']].mean(axis=1)
        data['std_abundance'] = data_log[[f'{tpm1}_Rep1', f'{tpm2}_Rep2']].std(axis=1)
    
    elif tpm1!=tpm2:
        data_log = np.log(data[[f'{tpm1}', f'{tpm2}']] + 1)
        data['mean_abundance'] = data_log[[f'{tpm1}', f'{tpm2}']].mean(axis=1)
        data['std_abundance'] = data_log[[f'{tpm1}', f'{tpm2}']].std(axis=1)   
    
    
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

    return data, ACVC, IM

def calculate_im_acvc(rep1, rep2, tr1, tr2, tpm1, tpm2):

    rep1_data = csv_tpm_processing(rep1, rep2, tr1, tr2, suffixes=('_Rep1', '_Rep2'))

    rep1_data, ACVC, IM = calculate_cv(rep1_data, tpm1, tpm2)
    
    return round(ACVC, 5), round(IM, 5)

def calculate_nrmse(file_path1, file_path2, tr1, tr2, tpm1, tpm2):

    common_isoforms = csv_tpm_processing(file_path1, file_path2, tr1, tr2)

    if tpm1==tpm2:
        quant = f'{tpm1}_quant'
        truth = f'{tpm2}_truth'
    else:
        quant = tpm1
        truth = tpm2

    if quant in common_isoforms and truth in common_isoforms:
        estimated = np.log(common_isoforms[quant]+1)
        ground_truth = np.log(common_isoforms[truth]+1)

        mse = np.mean((estimated - ground_truth) ** 2)
        rmse = np.sqrt(mse)
        std_dev = np.std(ground_truth)

        nrmse = rmse / std_dev
        return round(nrmse, 5)

def calculate_mrd(file_path1, file_path2, tr1, tr2, tpm1, tpm2):

    common_isoforms = csv_tpm_processing(file_path1, file_path2, tr1, tr2)

    if tpm1==tpm2:
        quant = f'{tpm1}_quant'
        truth = f'{tpm2}_truth'
    else:
        quant = tpm1
        truth = tpm2

    if quant in common_isoforms and truth in common_isoforms:
        # estimated = np.log(common_isoforms[quant]+1)
        # ground_truth = np.log(common_isoforms[truth]+1)

        estimated = common_isoforms[quant]
        ground_truth = common_isoforms[truth]

        relative_differences = np.abs(estimated - ground_truth) / np.abs(ground_truth)
        mrd = np.median(relative_differences)

        return round(mrd, 5)


def get_stats(file1, file2, tr1, tr2, tpm1, tpm2, extra=False):
    
    """ ## CALL 'spearman_corr_generic' FUNC FIRST ## """

    sc, pc = spearman_corr_generic(file1, file2, tr1, tr2, tpm1, tpm2)
    acvc, im = calculate_im_acvc(file1, file2, tr1, tr2, tpm1, tpm2)

    result_dict = {'sc': sc, 'pc': pc, 'im': im}

    if extra:
        nrmse = calculate_nrmse(file1, file2, tr1, tr2, tpm1, tpm2)
        mrd = calculate_mrd(file1, file2, tr1, tr2, tpm1, tpm2)
        result_dict['nrmse'] = nrmse
        result_dict['mrd'] = mrd

    return result_dict