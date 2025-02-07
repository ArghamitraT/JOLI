"""
This file implements EM, VI and calls GD.
The math and detail explanation is on this file: https://drive.google.com/file/d/1LGLhGvn3KRAYunf995ZVAYA4w2lgRYfr/view?usp=sharing
"""


# ~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~ #

# Third party imports
import pysam
import pandas as pd
from tqdm import tqdm
import hashlib
import pickle
import numpy as np
from scipy.special import psi, gammaln
import time
import pickle
from scipy.stats import spearmanr, pearsonr
from contextlib import redirect_stdout

from DirichletOptimizer_vector import DirichletModel 
#from dirichlet_pyro import DirichletModel


# Local imports
from NanoCount.Read import Read
from NanoCount.common import *

crnt_tm = datetime.datetime.now()


# ~~~~~~~~~~~~~~MAIN FUNCTION~~~~~~~~~~~~~~ #
class Expec_Max:

    # ~~~~~~~~~~~~~~MAGIC METHODS~~~~~~~~~~~~~~ #
    def __init__(
            self,
            file_names = [],
            alignment_file: str = "",
            count_file: str = "",
            filter_bam_out: str = "",
            min_alignment_length: int = 50,
            keep_suplementary: bool = False,
            min_query_fraction_aligned: float = 0.5,
            sec_scoring_threshold: float = 0.95,
            sec_scoring_value: str = "alignment_score",
            convergence_target: float = 0.001, 
            max_em_rounds: int = 10, 
            extra_tx_info: bool = False,
            primary_score: str = "alignment_score",
            max_dist_3_prime: int = 50,
            max_dist_5_prime: int = -1,
            verbose: bool = False,
            quiet: bool = False,
            alpha_initial: float = 1, 
            GD_lr: float = 0.01,
            process: str ='theta', #'expectation_log_theta' or 'theta'
            dirichlet_builtin: int = 1, 
            load: int = 0,
            load_filename: str = "",
            experiment_num: int = 4,
            EM_type = "MAP" #Can be VI or MAP

    ):
        """
        NOTE: File file does not use any hashing and runs EM on both ambiguous and unambiguous reads. This can be improved in future
        EXPLANATION:
        Estimate abundance of transcripts using an EM and MAP
        * alignment_file
            Sorted and indexed BAM or SAM file containing aligned ONT dRNA-Seq reads including secondary alignments
        * count_file
            Output file path where to write estimated counts (TSV format)
        * filter_bam_out
            Optional output file path where to write filtered reads selected by NanoCount to perform quantification estimation (BAM format)
        * min_alignment_length
            Minimal length of the alignment to be considered valid
        * min_query_fraction_aligned
            Minimal fraction of the primary alignment query aligned to consider the read valid
        * sec_scoring_threshold
            Fraction of the alignment score or the alignment length of secondary alignments compared to the primary alignment to be considered valid
            alignments
        * sec_scoring_value
            Value to use for score thresholding of secondary alignments either "alignment_score" or "alignment_length"
        * convergence_target
            Convergence target value of the cummulative difference between abundance values of successive EM round to trigger the end of the EM loop.
        * max_em_rounds
            Maximum number of EM rounds before triggering stop
        * extra_tx_info
            Add transcripts length and zero coverage transcripts to the output file (required valid bam/sam header)
        * primary_score
            Method to pick the best alignment for each read. By default ("alignment_score") uses the best alignment score (AS optional field), but it can be changed to
            use either the primary alignment defined by the aligner ("primary") or the longest alignment ("alignment_length"). choices = [primary, alignment_score, alignment_length]
        * keep_suplementary
            Retain any supplementary alignments and considered them like secondary alignments. Discarded by default.
        * max_dist_3_prime
            Maximum distance of alignment end to 3 prime of transcript. In ONT dRNA-Seq reads are assumed to start from the polyA tail (-1 to deactivate)
        * max_dist_5_prime
            Maximum distance of alignment start to 5 prime of transcript. In conjunction with max_dist_3_prime it can be used to select near full transcript reads
            only (-1 to deactivate).
        * verbose
            Increase verbosity for QC and debugging
        * quiet
            Reduce verbosity
        """

        # Init package
        opt_summary_dict = opt_summary(local_opt=locals())
        self.log = get_logger(name="Nanocount", verbose=verbose, quiet=quiet)

        print("Checking options and input files")
        log_dict(opt_summary_dict, self.log.debug, "Options summary")

        # Save args in self variables
        self.file_names_list = file_names
        self.alignment_file = alignment_file
        self.count_file = count_file
        self.filter_bam_out = filter_bam_out
        self.min_alignment_length = min_alignment_length
        self.min_query_fraction_aligned = min_query_fraction_aligned
        self.sec_scoring_threshold = sec_scoring_threshold
        self.sec_scoring_value = sec_scoring_value
        self.convergence_target = convergence_target
        self.max_em_rounds = max_em_rounds
        self.extra_tx_info = extra_tx_info
        self.primary_score = primary_score
        self.keep_suplementary = keep_suplementary
        self.max_dist_5_prime = max_dist_5_prime
        self.max_dist_3_prime = max_dist_3_prime
        self.all_read_dicts = {}
        self.all_ref_len_dicts = {}
        self.all_Yri = {}
        self.all_theta = {}
        self.all_Phi_ri = {}
        self.all_n = {}
        self.all_alpha = {}
        self.all_alpha_prime = {}
        self.all_isoform_indices = {}
        self.all_read_iso_prob = {}
        self.expectation_log_theta = {}
        self.elbo = {}
        self.alpha_initial = alpha_initial
        self.GD_lr = GD_lr
        self.process=process
        self.load=load
        self.load_filename = load_filename
        new_max_em_rounds = max_em_rounds
        self.all_readIsoformMap = {}
        self.all_readName = {}
        self.theta_names = {}
        self.experiment_num = experiment_num
        self.dirichlet_builtin=dirichlet_builtin
        self.EM_type = EM_type


        print("Initialise Nanocount")
        print("Parse Bam file and filter low quality alignments")
        
        start = time.time()
        start1 = time.time() # comment

        if load:
            self.load_state(max_em_rounds, count_file)
        
        else:
            self.initialize_model()
        
        self.EM_type == 'MAP'
        self.process ='theta' # 'expectation_log_theta' or 'theta'
        self.MAP()

    
        # Write out results
        print("Summarize data")
        for sample, theta_array in self.all_theta.items():
            print(f"Processing {sample}")
            
            # Retrieve the isoform names for the current sample
            theta_names = self.theta_names[sample]

            # Convert `theta_array` and `theta_names` into a dictionary
            theta_dict = {name: value for name, value in zip(theta_names, theta_array)}

            # Create DataFrame using the theta dictionary
            count_df = pd.DataFrame.from_dict(theta_dict, orient='index', columns=["raw"])
            count_df.index.name = "transcript_name"

            print("Compute estimated counts and TPM")
            # Adjusted to use the length of read_dict for the current sample
            count_df["est_count"] = count_df["raw"] * len(self.all_read_dicts[sample])
            count_df["tpm"] = count_df["raw"] * 1000000 / sum(theta_dict.values())

            # Add extra transcript info if required
            if self.extra_tx_info:
                tx_df = self._get_tx_df()
                count_df = pd.merge(count_df, tx_df, left_index=True, right_index=True, how="outer")

            # Cleanup and sort
            count_df.sort_values(by="raw", ascending=False, inplace=True)
            count_df.fillna(value=0, inplace=True)

            # Create file names and save the DataFrame to a TSV file
            if self.experiment_num == 5:
                file_name = self.create_saved_state_filename_exp5(result=self.count_file + sample)
            else:
                file_name = self.create_saved_state_filename(result=self.count_file + sample)

            
            file_name_timestamp = self.create_image_name(file_name, format="")
            count_file = f"{file_name_timestamp}.tsv" if self.count_file else None

            if count_file:
                print(f"Write file for {sample}")
                
                dir_path = os.path.dirname(count_file)
                if not os.path.exists(dir_path):
                    os.makedirs(dir_path)
                count_df.to_csv(count_file, sep="\t")
            
        end = time.time()
        interval = (end-start)/60
        print(f"time_total {interval} min")



     # ~~~~~~~~~~~~~~NEW FUNCTIONS~~~~~~~~~~~~~~ #
    def MAP(self):
        """ EM and Gradient descent """
        start1 = time.time() # comment
        
        ## Initialize the Dirichlet optimizer with the theta data
        dirichlet_optimizer = DirichletModel(self.all_alpha, self.GD_lr, self.dirichlet_builtin, self.process)
        
        # comment
        end1 = time.time()
        interval = (end1-start1)/60
        print(f"time_initialize_diri {interval} min")
        start1 = time.time() # comment
        
        # filename for saving variables
        # token = (self.count_file.split('/')[-1]).split('_')[-1]
        # model_save_path = '/'.join(self.count_file.split('/')[:-3])+'/weights/'
        token = '00000'
        model_save_path = self.count_file
        
        if self.experiment_num == 5:
            saved_state_filename = self.create_saved_state_filename_exp5(result="allWeights")
            final_save_path = self.create_image_name(model_save_path+saved_state_filename+'_token_'+token, format=".pkl")
            saved_EMstat_filename = self.create_saved_state_filename_exp5(result="allEMstats")
            final_EMstat_path = self.create_image_name(model_save_path+saved_EMstat_filename+'_token_'+token, format=".csv")
            saved_GDloss_filename = self.create_saved_state_filename_exp5(result="allGDloss")
            final_GDloss_path = self.create_image_name(model_save_path+saved_GDloss_filename+'_token_'+token, format=".pkl")
        else:
            saved_state_filename = self.create_saved_state_filename(result="allWeights")
            final_save_path = self.create_image_name(model_save_path+saved_state_filename+'_token_'+token, format=".pkl")
            saved_EMstat_filename = self.create_saved_state_filename(result="allEMstats")
            final_EMstat_path = self.create_image_name(model_save_path+saved_EMstat_filename+'_token_'+token, format=".csv")
            saved_GDloss_filename = self.create_saved_state_filename(result="allGDloss")
            final_GDloss_path = self.create_image_name(model_save_path+saved_GDloss_filename+'_token_'+token, format=".pkl")

        # Initialize an empty DataFrame to store the stats
        stats_df = pd.DataFrame(columns=[
            'EM_loop', 'Alpha_summation', 'EM_convergence', 
            'Spearman_corr_theta1_theta2', 
            'Spearman_corr_theta1_alpha', 
            'Spearman_corr_theta2_alpha',
            'pearson_corr_theta1_theta2',
            'pearson_corr_theta1_alpha', 
            'pearson_corr_theta2_alpha'
        ] + [f'ELBO_sample_{i+1}' for i in range(len(self.all_read_dicts))] 
        + [f'Convergence_sample_{i+1}' for i in range(len(self.all_read_dicts))])
        
        # Initialize an empty DataFrame to store the GD loss history
        GDloss_history_df = pd.DataFrame(columns=['GD_Loss'])
        
        # comment
        end1 = time.time()
        interval = (end1-start1)/60
        print(f"time_household_before_EM {interval} min")

        # Iterate until convergence threshold or max EM round are reached
        while self.convergence > self.convergence_target and self.em_round < self.max_em_rounds:
        # while self.em_round < self.max_em_rounds:
            self.convergence = 0
            self.em_round += 1
            elbo_values = {}
            convergence_values = {}

            start1 = time.time() # comment

            # UPDATE ALPHA
            if self.experiment_num == 4 or self.experiment_num == 5:
                self.all_alpha, GDloss_history = dirichlet_optimizer.update_alpha(self.all_unique_isoforms, self.theta_names, expectation_log_theta=self.expectation_log_theta, 
                                                              all_theta=self.all_theta)
            
            # comment
            end1 = time.time()
            interval = (end1-start1)/60
            print(f"time_GDloss {interval} min")

            start1 = time.time() # comment

            # EM
            sample_num = 0
            for sample_key in self.all_read_dicts:
                self.all_theta[sample_key], convergence = self.update_theta_vectorized(sample_key)
                self.all_Phi_ri[sample_key] = self.calculate_Z(self.all_readName[sample_key], self.all_readIsoformMap[sample_key], self.all_Yri[sample_key], self.all_read_iso_prob[sample_key], self.all_theta[sample_key])
                self.all_n[sample_key] = self.calculate_n_vectorized(sample_key)
                sample_num +=1
                self.convergence += convergence
                self.elbo[sample_key] = 0

                # Store ELBO and Convergence for this sample
                elbo_values[f'ELBO_sample_{sample_num}'] = self.elbo[sample_key]
                convergence_values[f'Convergence_sample_{sample_num}'] = convergence

                # print(f"ELBO_sample_{sample_num} {self.elbo[sample_key]}")
                print(f"Convergence_sample_{sample_num} {convergence}")

            # comment
            end1 = time.time()
            interval = (end1-start1)/60
            print(f"time_EM_calculation {interval} min")

            print("alpha_summation ", np.sum(self.all_alpha))
            self.convergence = self.convergence/sample_num
            print(f"EM_loop {self.em_round}" )
            print(f"EM_convergence {self.convergence}")

            start1 = time.time() # comment

            # generate model fit correlation
            # if self.experiment_num == 4 or self.experiment_num == 5:
            #     spearman_corr_theta1_theta2, pearson_corr_theta1_theta2 = self.spearman_pearson_corr('theta1_theta2')
            #     spearman_corr_theta1_alpha, pearson_corr_theta1_alpha = self.spearman_pearson_corr('theta1_alpha')
            #     spearman_corr_theta2_alpha, pearson_corr_theta2_alpha = self.spearman_pearson_corr('theta2_alpha')
            # else:
            #     spearman_corr_theta1_theta2, pearson_corr_theta1_theta2 = 0, 0
            #     spearman_corr_theta1_alpha, pearson_corr_theta1_alpha = self.spearman_pearson_corr('theta1_alpha')
            #     spearman_corr_theta2_alpha, pearson_corr_theta2_alpha = 0, 0


            # comment
            end1 = time.time()
            interval = (end1-start1)/60
            print(f"time_corr_calculation {interval} min")

            start1 = time.time() # comment

            # Create a DataFrame for the current iteration
            # current_stats = pd.DataFrame([{
            #     'EM_loop': self.em_round,
            #     'Alpha_summation': np.sum(self.all_alpha),
            #     'EM_convergence': self.convergence,
            #     'Spearman_corr_theta1_theta2': spearman_corr_theta1_theta2,
            #     'Spearman_corr_theta1_alpha': spearman_corr_theta1_alpha,
            #     'Spearman_corr_theta2_alpha': spearman_corr_theta2_alpha,
            #     'pearson_corr_theta1_theta2': pearson_corr_theta1_theta2,
            #     'pearson_corr_theta1_alpha': pearson_corr_theta1_alpha, 
            #     'pearson_corr_theta2_alpha': pearson_corr_theta2_alpha,
            #     **elbo_values,
            #     **convergence_values
            # }])

            # # Concatenate the current iteration's DataFrame with the main DataFrame
            # stats_df = pd.concat([stats_df, current_stats], ignore_index=True)

            # if self.experiment_num == 4 or self.experiment_num == 5:
            #     # Create a DataFrame for the loss history
            #     GDloss_history_current = pd.DataFrame({
            #         'GD_Loss': GDloss_history
            #     })
            #     # Append the loss history of the current iteration to the loss history DataFrame
            #     GDloss_history_df = pd.concat([GDloss_history_df, GDloss_history_current], ignore_index=True)
                
            #     # Check if the directory exists
            #     dir_path = os.path.dirname(final_GDloss_path)
            #     if not os.path.exists(dir_path):
            #         os.makedirs(dir_path)
                    

            #     with open(final_GDloss_path, 'wb') as f:
            #         pickle.dump(GDloss_history_df, f)

            # Save the state after each iteration
            # stats_df.to_csv(final_EMstat_path, index=False)
            # self.save_state(final_save_path)

            # comment
            end1 = time.time()
            interval = (end1-start1)/60
            print(f"time_rest {interval} min")


        print(f"Exit EM loop after {self.em_round} rounds")
        print(f"Convergence value: {self.convergence}")
        if not self.convergence <= self.convergence_target:
            print(f"Convergence target {self.convergence_target} could not be reached after {self.max_em_rounds} rounds")


    
    def create_saved_state_filename_exp5(self, result):
        
        # Iterate through the list of file names
        for index, sample_file_path in enumerate(self.file_names_list, start=1):
            result += f"_file{index}_"
            for file_path in enumerate(sample_file_path, start=1):
                # Split the file path by '/' and take the last part (the file name)
                file_name = file_path[1].split('/')[-1]
                # Extract a specific part of the file name if necessary (e.g., removing extension)
                file_identifier = ''.join(file_name.split('_')).split('.')[0]
                # Construct the string
                result += f"{file_identifier}"
        result = f"{result}_GDlr_{self.GD_lr}_AlphaInitial_{self.alpha_initial}_EMround_{self.max_em_rounds}"
        return result
    
    def create_saved_state_filename(self, result):
        
        # Iterate through the list of file names
        for fileNamelist_index, fileName_list in enumerate(self.file_names_list, start=1):
            for index, file_path in enumerate(fileName_list, start=1):
                # Split the file path by '/' and take the last part (the file name)
                file_name = file_path.split('/')[-1]
                # Extract a specific part of the file name if necessary (e.g., removing extension)
                file_identifier = ''.join(file_name.split('_')).split('.')[0]
                # Construct the string
                result += f"_file{fileNamelist_index}_{file_identifier}"
        result = f"{result}_GDlr_{self.GD_lr}_AlphaInitial_{self.alpha_initial}_EMround_{self.max_em_rounds}"
        return result


     # Add these methods to your class
    def save_state(self, filename):
        with open(filename, 'wb') as f:
            pickle.dump(self, f)

    # @staticmethod
    def load_state(self,  new_max_em_rounds, count_file):
        # with open(self.load_filename, 'rb') as f:
        #     return pickle.load(f)
        print(f"loading_weights_from {self.load_filename}")
        with open(self.load_filename, 'rb') as f:
            loaded_obj = pickle.load(f)
            # Update the current instance's __dict__ with the loaded object's __dict__
            self.__dict__.update(loaded_obj.__dict__)
        self.max_em_rounds = new_max_em_rounds
        self.count_file = count_file
        

    def create_image_name(self, name, format=".png"):
        image_name = (name+"_" + str(crnt_tm.year) + "_" + str(crnt_tm.month) + "_" + str(crnt_tm.day) + "_"
                    + time.strftime("%H_%M_%S") + format)
        return image_name

    def fraction_to_float(self, fraction):
        # Assuming fraction is in string format '1/2', '3/4', etc.
        try:
            numerator, denominator = fraction.split('/')
            return float(numerator) / float(denominator)
        except ValueError:
            # If fraction is not in expected format, return it as a float
            return float(fraction)
        
    def spearman_pearson_corr(self, comparison_criteria):
        if comparison_criteria == 'theta1_theta2':
            model_weight1 = np.array(self.all_theta['sample1'])  # Array of weights for sample1
            model_weight2 = np.array(self.all_theta['sample2'])  # Array of weights for sample2
            theta_names_1 = np.array(self.theta_names['sample1'])
            theta_names_2 = np.array(self.theta_names['sample2'])
        elif comparison_criteria == 'theta1_alpha':
            model_weight1 = np.array(self.all_theta['sample1'])
            model_weight2 = np.array(self.all_alpha)
            theta_names_1 = np.array(self.theta_names['sample1'])
            theta_names_2 = np.array(self.all_unique_isoforms)  # Assuming `all_alpha` aligns with all unique isoforms
        elif comparison_criteria == 'theta2_alpha':
            model_weight1 = np.array(self.all_theta['sample2'])
            model_weight2 = np.array(self.all_alpha)
            theta_names_1 = np.array(self.theta_names['sample2'])
            theta_names_2 = np.array(self.all_unique_isoforms)
            
        # Normalize `model_weight2` if necessary
        if comparison_criteria in ['theta1_alpha', 'theta2_alpha']:
            model_weight2 = model_weight2/np.sum(model_weight2)
        
        # Find the common isoforms and their indices in both arrays
        common_isoforms, idx1, idx2 = np.intersect1d(theta_names_1, theta_names_2, return_indices=True)

        # Extract the common values from the arrays using the indices
        weight1_values = model_weight1[idx1]
        weight2_values = model_weight2[idx2]

        # Apply logarithmic transformation
        matched_tpm = np.log(weight1_values * 1000000  + 1)
        matched_molarity = np.log(weight2_values * 1000000  + 1)

        # Calculate Spearman's and Pearson's correlation using the matched arrays
        spearman_correlation, _ = spearmanr(matched_tpm, matched_molarity)
        pearson_correlation, _ = pearsonr(matched_tpm, matched_molarity)

        # Output the results
        if comparison_criteria == 'theta1_theta2':
            print(f'Spearman_corr_theta1_theta2 {spearman_correlation}')
            print(f'Pearson_corr_theta1_theta2 {pearson_correlation}')
        elif comparison_criteria == 'theta1_alpha':
            print(f'Spearman_corr_theta1_alpha {spearman_correlation}')
            print(f'Pearson_corr_theta1_alpha {pearson_correlation}')
        elif comparison_criteria == 'theta2_alpha':
            print(f'Spearman_corr_theta2_alpha {spearman_correlation}')
            print(f'Pearson_corr_theta2_alpha {pearson_correlation}')

        return spearman_correlation, pearson_correlation

    def spearman_corr(self, comparison_criteria):
        
        if comparison_criteria == 'theta1_theta2':
            model_weight1 = self.all_theta['sample1']
            model_weight2 = self.all_theta['sample2']
        elif comparison_criteria == 'theta1_alpha':
            model_weight1 = self.all_theta['sample1']
            model_weight2 = self.all_alpha
        elif comparison_criteria == 'theta2_alpha':
            model_weight1 = self.all_theta['sample2']
            model_weight2 = self.all_alpha

        # Convert dictionary to DataFrame
        if isinstance(model_weight1, dict):
            model_weight1 = pd.DataFrame(list(model_weight1.items()), columns=['transcript_name', 'tpm'])
            model_weight2 = pd.DataFrame(list(model_weight2.items()), columns=['transcript_name', 'tpm'])
        
        # if we are correlating theta with alpha, alpha needs to be normalized
        if comparison_criteria == 'theta1_alpha' or comparison_criteria == 'theta2_alpha':
            model_weight2['tpm'] = model_weight2['tpm'] / model_weight2['tpm'].sum()

        # Initialize lists to store matched TPM and molarity values
        matched_tpm = []
        matched_molarity = []

        # Iterate over first variable to find matching molarity values in second variable
        for index, row in model_weight1.iterrows():
            # Extract the cleaned isoform name and tpm value
            cleaned_name = row['transcript_name']
            tpm_value = row['tpm']

            # Find the corresponding molarity value in second variable
            molarity_value = model_weight2.loc[model_weight2['transcript_name'] == cleaned_name, 'tpm'].values

            # If a matching isoform is found in second variable
            if len(molarity_value) > 0:
                # Append the tpm and molarity to the respective lists
                matched_tpm.append(np.log(tpm_value * 1000000  + 1))
                matched_molarity.append(np.log(molarity_value[0]* 1000000+1))  # Take the first match in case of multiple

        # Calculate Spearman's correlation using the matched lists
        correlation, p_value = spearmanr(matched_tpm, matched_molarity)

        # Output the results
        if comparison_criteria == 'theta1_theta2':
            print(f'Spearman_corr_theta1_theta2 {correlation}')
        elif comparison_criteria == 'theta1_alpha':
            print(f'Spearman_corr_theta1_alpha {correlation}')
        elif comparison_criteria == 'theta2_alpha':
            print(f'Spearman_corr_theta2_alpha {correlation}')
        
        return correlation


    def assign_alpha(self):
        # Step 1: Use `self.all_unique_isoforms` directly
        unique_isoforms = self.all_unique_isoforms

        # Step 2: Number of unique isoforms
        n = len(unique_isoforms)

        # Step 3: Initialize alpha values in an unconstrained space
        np.random.seed(5)  # Set the fixed seed for reproducibility
        alpha_raw = np.random.randn(n)

        # Step 4: Apply the softmax function to ensure they sum to 1
        alpha_softmax = np.exp(alpha_raw) / np.sum(np.exp(alpha_raw))

        # Step 5: Scale the resulting values to have the desired fixed sum
        alpha_scaled = alpha_softmax * self.alpha_initial

        # Step 6: Return the alpha values
        return alpha_scaled
    
    def assign_alpha_constant(self):
        # Step 1: Use `self.all_unique_isoforms` directly
        unique_isoforms = self.all_unique_isoforms

        # Step 2: Number of unique isoforms
        n = len(unique_isoforms)

        # Step 5: Scale the resulting values to have the desired fixed sum
        alpha_scaled = np.full(n, 3)

        # Step 6: Return the alpha values
        return alpha_scaled
    

    def calculate_theta_and_alpha_prime_0(self, all_readIsoformMap, compatibility_list):
        """
        Calculates the initial model parameter (isoform percentage) using vectors.
        """
        # Identify unique reference names and assign indices
        theta_names, inverse_indices = np.unique(all_readIsoformMap, return_inverse=True)
        
        # Calculate total abundance scores for each unique reference name
        abundance_scores = np.zeros(len(theta_names))
        abundance_scores = abundance_scores.astype(float)

        np.add.at(abundance_scores, inverse_indices, compatibility_list)

        # Normalize the abundance scores
        total_score = np.sum(abundance_scores)
        normalized_abundances = abundance_scores / total_score

        return theta_names, normalized_abundances

    def update_theta(self, sample_key):
        """ Calculates the mode of the dirichlet """

        """ EXPLANATION:
            * Eqn 5 --> {\hat{\theta}_i  = \frac{n_i + \alpha_i - 1}
            {\sum_{i=1}^{I} (\alpha_i + n_i)-I}}"""

        sample_counts = self.all_n[sample_key]
        convergence = 0
        theta_old = self.all_theta[sample_key]

        # Prepare the arrays for n_i and alpha_i
        n_i_sample = []
        alpha_i_sample = []
        isoform_order = []  # To preserve the order of isoforms

        for isoform in sample_counts:
            index = self.isoform_to_index[isoform]
            n_i_sample.append(sample_counts[isoform])
            alpha_i_sample.append(self.alpha[index])
            isoform_order.append(isoform)

        n_i_sample = np.array(n_i_sample)
        alpha_i_sample = np.array(alpha_i_sample)

        # Calculate theta_hat
        numerator = n_i_sample + alpha_i_sample - 1
        denominator = np.sum(numerator) - len(n_i_sample)
        theta_hat_sample = numerator / denominator

        # Store the results in a Counter with the same isoform order
        theta_hat = Counter({isoform: theta_hat for isoform, theta_hat in zip(isoform_order, theta_hat_sample)})

        if self.em_round > 0:
            for ref_name in theta_old.keys():
                convergence += abs(theta_old[ref_name] - theta_hat[ref_name])

        return theta_hat, convergence
    

    def update_theta_vectorized(self, sample_key):
        """
        Calculate θ̂ (theta_hat) for a given sample using the vectorized form of the equation.
        """

        # Retrieve n values and their associated names
        sample_n_values = self.all_n[sample_key]  # n_i for the sample
        sample_theta_names = self.theta_names[sample_key]  # Names corresponding to n_i
        theta_old = self.all_theta[sample_key]

        
        if self.experiment_num == 4 or self.experiment_num == 5: # multi sample
            # Retrieve alpha values and their associated names
            alpha_values = self.all_alpha  # Full list of alpha values across all samples
            alpha_names = self.all_unique_isoforms  # Names corresponding to alpha values

            # Map sample_theta_names to indices in unique_isoforms (alpha_names)
            alpha_indices = np.searchsorted(alpha_names, sample_theta_names)

            # Gather alpha values corresponding to the sample_theta_names
            sample_alpha_values = alpha_values[alpha_indices]

            # Calculate numerator: n_i + α_i - 1
            n_plus_alpha = sample_n_values + sample_alpha_values
            numerator = n_plus_alpha 

            # Calculate denominator: sum(numerator) - I
            denominator = np.sum(n_plus_alpha)

            # Compute θ̂ (theta_hat): numerator / denominator
            theta_hat = numerator / denominator
        
        else:   # single sample
            print("single sample theta calculation, alpha NOT considered")
            numerator = sample_n_values
            denominator = np.sum(numerator) 
            theta_hat = numerator / denominator
        
        if self.em_round > 0:
            # Convert theta_old to numpy array using indices for direct comparison
            convergence = np.sum(np.abs(theta_old - theta_hat))
        else:
            convergence = 0

        return theta_hat, convergence

    def calculate_n(self, compatibility_dict):
        """
        Sums up the total assignment for isoform I
        """
        """ EXPLANATION:
            * Eqn 4 -->  n_i = \sum_{r=1}^{R} z_{ri}"""

        abundance_dict = Counter()
        total = 0

        # sums over all the reads and calculates the total abundance of  each isoform
        for read_name, comp in compatibility_dict.items():
            for ref_name, score in comp.items():
                abundance_dict[ref_name] += score
                total += score
        return abundance_dict


    def calculate_Z(self, all_read_names, all_readIsoformMap, compatibility_list, p_nm, all_theta):
        """
        Vectorized EM assignment to update read-transcript compatibility based on transcript abundances.
        """
        # Find unique reads and map them to indices
        unique_reads, read_indices = np.unique(all_read_names, return_inverse=True)

        # Map isoforms (rnames) to indices in theta_names (isoforms used in theta)
        unique_rnames, isoform_indices = np.unique(all_readIsoformMap, return_inverse=True)

        # Get theta values for each isoform (isoform values are strored as sorted isoform names)
        theta_for_isoforms = all_theta[isoform_indices]
        
        # Step 1: Calculate denominator for each read
        compat_theta_product = p_nm * theta_for_isoforms
        denominators = np.zeros(len(unique_reads))
        denominators = denominators.astype(float)

        np.add.at(denominators, read_indices, compat_theta_product)

        # Step 2: Calculate Z_ri values
        Z_values = (p_nm * theta_for_isoforms) / denominators[read_indices]

        return Z_values


    def calculate_theta_0(self, compatibility_dict):
        """
        Calculates the initial model parameter, or isoform percentage
        """
        """ EXPLANATION:
            * Eqn 1 --> p_{i} = \frac{1}{R} \sum_{r=1}^{R} y_{ri}
            * Eqn 2 --> \theta_{i}^{t=0} = \frac{p_{i}}{\sum_{i=1}^{I} p_{i}} """

        abundance_dict = Counter()
        total = 0
        isoform_indices = {}
        isoform_counter = 0

        for read_name, comp in compatibility_dict.items():
            for ref_name, score in comp.items():
                abundance_dict[ref_name] += score
                #compatibility_dict[ref_name] += score
                total += score
                isoform_indices[ref_name] = isoform_counter
                isoform_counter += 1

        for ref_name in abundance_dict.keys():
            abundance_dict[ref_name] = abundance_dict[ref_name] / total

        return abundance_dict, isoform_indices
    
    def calculate_n_vectorized(self, sample_key):
        """
        Fully vectorized function to calculate the sum of Phi values for each isoform in theta_names.

        Parameters:
        - all_Phi: array of Phi values for each read.
        - all_readIsoformMap: array of isoform names for each read (same length as all_Phi).
        - theta_names: array of unique isoform names in the order you want the output.

        Returns:
        - A numpy array where each entry corresponds to the total Phi value for each isoform in theta_names.
        """
        # Convert all_readIsoformMap to numpy array for easier indexing
        all_readIsoformMap = self.all_readIsoformMap[sample_key]
        all_Phi = self.all_Phi_ri[sample_key]
        theta_names = self.theta_names[sample_key]

        # Map all_readIsoformMap to indices of theta_names
        theta_indices = np.searchsorted(theta_names, all_readIsoformMap)

        # Mask out indices where all_readIsoformMap is not in theta_names
        valid_indices = np.isin(all_readIsoformMap, theta_names)
        filtered_indices = theta_indices[valid_indices]
        filtered_Phi = all_Phi[valid_indices]

        # Sum Phi values for each index in theta_names order
        n_values = np.zeros(len(theta_names), dtype=float)  # Initialize an array of zeros
        np.add.at(n_values, filtered_indices, filtered_Phi)  # Add weights to corresponding bins

        return n_values

    
    def hash_key(self, elements):
        # Create a stable hash of the sorted elements using a cryptographic hash function
        # Convert the elements tuple to a string and encode to bytes before hashing
        return hashlib.sha256(str(sorted(elements)).encode('utf-8')).hexdigest()

    
    def is_iterable(self, obj):
        """
        Check if obj is iterable but not a string.
        """
        return isinstance(obj, (list, tuple, set, dict, frozenset)) and not isinstance(obj, (str, bytes))

    def get_compatibility_modified(self, sample_key):
        """
            This function computes compatibility scores and normalized probabilities for alignments of sequencing reads to reference isoforms.

            ### Functionality:
            1. **Input Details:**
            - Takes a sample key to retrieve read and reference length information.
            - Processes reads and alignments to calculate alignment lengths and associated isoforms.

            2. **Outputs:**
            - Normalized compatibility scores based on the number of alignments per read.
            - Probabilities for each alignment, normalized by the reference isoform length adjusted by the alignment length.

            3. **Steps:**
            - Vectorizes the handling of reads, alignments, and reference lengths for efficient computation.
            - Sorts alignments by read names and isoform names for consistent ordering.
            - Computes normalized probabilities while handling edge cases where the denominator is invalid.

            4. **Output Values:**
            - `all_read_names`: Names of reads in sorted order.
            - `all_readIsoformMap`: Associated isoforms for each read alignment.
            - `compatibility_list`: Compatibility scores inversely proportional to the number of alignments per read.
            - `p_nm`: Normalized probabilities for each alignment, adjusted by reference length and alignment length.
            """

        # Fetch reference lengths and names
        read_dict = self.all_read_dicts[sample_key]
        ref_len_dict = self.all_ref_len_dicts[sample_key]

    
        # Create a mapping of reference names to indices for quick lookup
        ref_len_dict_keys = np.array(list(ref_len_dict.keys()))
        
        ref_lengths = np.fromiter((ref_len_dict[key] for key in ref_len_dict_keys), dtype=float)
        ref_name_to_index = dict(zip(ref_len_dict_keys, np.arange(len(ref_len_dict_keys))))

        all_read_names = []
        all_alignments = []
        all_alignment_lens = []
        all_readIsoformMap = []

        # Collect all data into lists for vectorized processing
        for read_name, read in read_dict.items():
            if not self.is_iterable(read):
                read = [read]  # Wrap non-iterable read in a list

            for alignment in read:
                alignments = alignment.alignment_list
                alignment_lens = np.array([string.align_len for string in alignments])
                rnames = [string.rname for string in alignments]

                # Collect the read name and alignment details for each alignment
                all_read_names.extend([read_name] * len(alignments))
                all_alignments.extend(alignments)
                all_alignment_lens.extend(alignment_lens)
                all_readIsoformMap.extend(rnames)

        # Convert all collected data into numpy arrays for vectorized processing
        all_alignment_lens = np.array(all_alignment_lens)
        all_readIsoformMap = np.array(all_readIsoformMap)
        all_read_names = np.array(all_read_names)

        # Sort all arrays based on `all_read_names` and `all_readIsoformMap` to maintain consistent ordering
        sort_indices = np.lexsort((all_readIsoformMap, all_read_names))# Sort by `all_read_names`, then by `all_readIsoformMap`
        all_read_names = all_read_names[sort_indices]
        all_readIsoformMap = all_readIsoformMap[sort_indices]
        all_alignment_lens = all_alignment_lens[sort_indices]

        # Find the indices of reference names in ref_len_dict_keys
        ref_indices = np.vectorize(ref_name_to_index.get)(all_readIsoformMap)
        ref_lengths_for_alignments = ref_lengths[ref_indices]

        # Vectorized compatibility scores
        unique_read_names, counts = np.unique(all_read_names, return_counts=True)
        compatibility_list = np.repeat(1 / counts, counts)  # Compatibility score is the inverse of the count of alignments for each read

        # Vectorized probability calculation
        denominators = ref_lengths_for_alignments - all_alignment_lens + 1
        p_nm = np.where(denominators <= 0, np.nan, 1 / denominators)

        # Handle problematic probabilities
        problematic_indices = np.isnan(p_nm) | (p_nm < 0)
        p_nm[problematic_indices] = 1  # Replace problematic probabilities with 1
        
        return all_read_names, all_readIsoformMap, compatibility_list, p_nm
    

    def initialize_model(self):
        start = time.time()

        if self.experiment_num == 4 or self.experiment_num == 1:
            for fileName_list_index,fileName_list in enumerate(self.file_names_list, start=1):
                for index, file_name in enumerate(fileName_list, start=1):
                # for index, file_name in enumerate(self.file_names_list, start=1):
                    print(f"sample_{index} {file_name}")
                    
                    pkl_file_name = file_name+'_read_dicts.pkl'
                    with open(pkl_file_name, 'rb') as file:
                        read_dict = pickle.load(file)
                        print(f"opened {pkl_file_name}")

                    pkl_file_name = file_name+'_ref_len_dicts.pkl'
                    with open(pkl_file_name, 'rb') as file:
                        ref_len_dict = pickle.load(file)
                        print(f"opened {pkl_file_name}")


                    # Store the dictionaries in the universal dictionary with a sample key
                    sample_key = f'sample{fileName_list_index}'
                    self.all_read_dicts[sample_key] = read_dict
                    self.all_ref_len_dicts[sample_key] = ref_len_dict
        
        elif self.experiment_num == 2:
            for fileName_list_index, fileName_list in enumerate(self.file_names_list, start=1):
                for index, file_name in enumerate(fileName_list, start=1):
                # for index, file_name in enumerate(self.file_names_list, start=1):
                    print(f"sample_{index} {file_name}")
                    
                    pkl_file_name = file_name+'_read_dicts.pkl'
                    with open(pkl_file_name, 'rb') as file:
                        read_dict = pickle.load(file)
                        print(f"opened {pkl_file_name}")

                    pkl_file_name = file_name+'_ref_len_dicts.pkl'
                    with open(pkl_file_name, 'rb') as file:
                        ref_len_dict = pickle.load(file)
                        print(f"opened {pkl_file_name}")


                    # Merge read_dict
                    if fileName_list_index ==1:
                        merged_read_dict = read_dict
                        merged_ref_len_dict = ref_len_dict
                    else:
                        for key, value in read_dict.items():
                            if key in merged_read_dict:
                                for subkey, subvalue in value.items():
                                    merged_read_dict[key][subkey]+=subvalue
                            else:
                                merged_read_dict[key]=value
                
                        # Merge ref_len_dict
                        for key, value in ref_len_dict.items():
                            if key in merged_ref_len_dict:
                                continue
                            else:
                                merged_ref_len_dict[key] = value
            # Store the dictionaries in the universal dictionary with a sample key
            sample_key = f'sample1'
            self.all_read_dicts[sample_key] = merged_read_dict
            self.all_ref_len_dicts[sample_key] = merged_ref_len_dict
        
        elif self.experiment_num == 5:
            for fileName_list_index, fileName_list in enumerate(self.file_names_list, start=1):
                for index, file_name in enumerate(fileName_list, start=1):
                    print(f"sample_{index} {file_name}")
                    
                    pkl_file_name = file_name+'_read_dicts.pkl'
                    with open(pkl_file_name, 'rb') as file:
                        read_dict = pickle.load(file)
                        print(f"opened {pkl_file_name}")

                    pkl_file_name = file_name+'_ref_len_dicts.pkl'
                    with open(pkl_file_name, 'rb') as file:
                        ref_len_dict = pickle.load(file)
                        print(f"opened {pkl_file_name}")


                    # Merge read_dict
                    if index ==1:
                        merged_read_dict = read_dict
                        merged_ref_len_dict = ref_len_dict
                    else:
                        for key, value in read_dict.items():
                            if key in merged_read_dict:
                                for subkey, subvalue in value.items():
                                    merged_read_dict[key][subkey]+=subvalue
                            else:
                                merged_read_dict[key]=value
                
                        # Merge ref_len_dict
                        for key, value in ref_len_dict.items():
                            if key in merged_ref_len_dict:
                                continue
                            else:
                                merged_ref_len_dict[key] = value
            # Store the dictionaries in the universal dictionary with a sample key
                sample_key = f'sample{fileName_list_index}'
                self.all_read_dicts[sample_key] = merged_read_dict
                self.all_ref_len_dicts[sample_key] = merged_ref_len_dict

        
        end = time.time()
        interval = (end-start)/60
        print(f"time_parse {interval} min")


        # COMMENT
        print("parsing done")

        if self.filter_bam_out:
            print("Write selected alignments to BAM file")
            self._write_bam()

        """ EXPLANATION
            * Yri --> Binary Compatibility Matrix
            * theta --> isoform percentage/abundance/quantity
            * Phi_ri --> Expectation of A (true asignment) matrix
            * n --> # of reads for each isoform
        """

        
        all_unique_isoforms_set = set()

        start1 = time.time() # comment
        
        # All the initial calculation
        for sample_key in self.all_read_dicts:
            all_read_names, all_readIsoformMap, compatibility_list, p_nm = self.get_compatibility_modified(sample_key)
            print(f"read_num {len(all_read_names)}")
            self.all_Yri[sample_key], self.all_read_iso_prob[sample_key] = compatibility_list, p_nm
            self.all_readIsoformMap[sample_key], self.all_readName[sample_key] = all_readIsoformMap, all_read_names

            theta_names, all_theta = self.calculate_theta_and_alpha_prime_0(all_readIsoformMap, compatibility_list)
            self.all_theta[sample_key], self.theta_names[sample_key] = all_theta, theta_names
            
            # Update the set of unique isoforms
            all_unique_isoforms_set.update(theta_names)


            self.all_Phi_ri[sample_key] = self.calculate_Z(all_read_names, all_readIsoformMap, compatibility_list, p_nm, all_theta)
            if self.EM_type == 'MAP':
                self.all_n[sample_key] = self.calculate_n_vectorized(sample_key)

            

        # Convert the set of unique isoforms to a NumPy array
        all_unique_isoforms_set = sorted(all_unique_isoforms_set)
        self.all_unique_isoforms = np.array(list(all_unique_isoforms_set))

        # self.all_alpha = self.assign_alpha_constant()
        self.all_alpha = self.assign_alpha()

        # comment
        end1 = time.time()
        interval = (end1-start1)/60
        print(f"time_initialize {interval} min")

        # COMMENT
        print("Initiation")

        
        # EM loop to calculate abundance and update read-transcript compatibility
        print("Start EM abundance estimate")

        self.em_round = 0
        self.convergence = 1

    # ~~~~~~~~~~~~~~PRIVATE METHODS~~~~~~~~~~~~~~ #

    def _parse_bam(self, file_name):
        """
        Parse Bam/Sam file, group alignments per reads, filter reads based on
        selection criteria and return a dict of valid read/alignments
        """
        # Parse bam files
        read_dict = defaultdict(Read)
        ref_len_dict = OrderedDict()
        c = Counter()

        aligned_read = file_name

        with pysam.AlignmentFile(aligned_read) as bam:

            # Collect reference lengths in dict
            for name, length in zip(bam.references, bam.lengths):
                ref_len_dict[name] = length

            for idx, alignment in enumerate(bam):
                if alignment.is_unmapped:
                    c["Discarded unmapped alignments"] += 1
                elif alignment.is_reverse:
                    c["Discarded negative strand alignments"] += 1
                # elif not self.keep_suplementary and alignment.is_supplementary:
                #     c["Discarded supplementary alignments"] += 1
                # elif self.min_alignment_length > 0 and alignment.query_alignment_length < self.min_alignment_length:
                #     c["Discarded short alignments"] += 1
                # elif self.max_dist_3_prime >= 0 and alignment.reference_end <= ref_len_dict[
                #     alignment.reference_name] - self.max_dist_3_prime:
                #     c["Discarded alignment with invalid 3 prime end"] += 1
                # elif self.max_dist_5_prime >= 0 and alignment.reference_start >= self.max_dist_5_prime:
                #     c["Discarded alignment with invalid 5 prime end"] += 1
                else:
                    c["Valid alignments"] += 1
                    read_dict[alignment.query_name].add_pysam_alignment(pysam_aligned_segment=alignment, read_idx=idx)

        # Write filtered reads counters
        log_dict(
            d=c,
            logger=self.log.info,
            header="Summary of alignments parsed in input bam file",
        )

        # Filter alignments
        filtered_read_dict = defaultdict(Read)
        c = Counter()

        for query_name, read in read_dict.items():
            # Check if best alignment is valid
            best_alignment = read.get_best_alignment(primary_score=self.primary_score)

            # In case the primary alignment was removed by filters
            if best_alignment:
                if best_alignment.align_score == 0:
                    c["Reads with zero score"] += 1
                elif best_alignment.align_len == 0:
                    c["Reads with zero len"] += 1
                elif best_alignment.query_fraction_aligned < self.min_query_fraction_aligned:
                    c["Reads with low query fraction aligned"] += 1
                else:
                    filtered_read_dict[query_name].add_alignment(best_alignment)
                    c["Reads with valid best alignment"] += 1
                    for alignment in read.get_secondary_alignments_list(primary_score=self.primary_score):

                        # Filter out secondary alignments based on minimap alignment score
                        if self.sec_scoring_value == "alignment_score" and alignment.align_score / best_alignment.align_score < self.sec_scoring_threshold:
                            c["Invalid secondary alignments"] += 1

                        # Filter out secondary alignments based on minimap alignment length
                        elif self.sec_scoring_value == "alignment_length" and alignment.align_len / best_alignment.align_len < self.sec_scoring_threshold:
                            c["Invalid secondary alignments"] += 1

                        # Select valid secondary alignments
                        else:
                            c["Valid secondary alignments"] += 1
                            filtered_read_dict[query_name].add_alignment(alignment)
            else:
                c["Reads without best alignment"] += 1

        if not "Valid secondary alignments" in c:
            self.log.error(
                "No valid secondary alignments found in bam file. Were the reads aligned with minimap `-p 0 -N 10` options ?")

        # Write filtered reads counters
        log_dict(d=c, logger=self.log.info, header="Summary of reads filtered")

        return filtered_read_dict, ref_len_dict

    def _write_bam(self):
        """"""
        c = Counter()

        # Make list of alignments idx to select
        selected_read_idx = set()
        for read in self.read_dict.values():
            for alignment in read.alignment_list:
                selected_read_idx.add(alignment.read_idx)
                c["Alignments to select"] += 1

        # Select from original bam file and write to output bam file
        with pysam.AlignmentFile(self.alignment_file) as bam_in:
            with pysam.AlignmentFile(self.filter_bam_out, "wb", template=bam_in) as bam_out:
                for read_idx, alignment in enumerate(bam_in):
                    if read_idx in selected_read_idx:
                        bam_out.write(alignment)
                        c["Alignments written"] += 1
                    else:
                        c["Alignments skipped"] += 1

        log_dict(d=c, logger=self.log.info, header="Summary of alignments written to bam")

    def _get_compatibility(self):
        """"""
        compatibility_dict = defaultdict(dict)

        for read_name, read in self.read_dict.items():
            # Check if read is iterable and not a string; if not, make it a list
            if not self.is_iterable(read):
                read = [read]  # Wrap non-iterable read in a list

            for alignment in read:
                for string in alignment.alignment_list:
                    compatibility_dict[read_name][string.rname] = \
                        score = 1.0 / alignment.n_alignment

        return compatibility_dict

    def _calculate_abundance(self):
        """
        Calculate the abundance of the transcript set based on read-transcript compatibilities
        """
        abundance_dict = Counter()
        total = 0
        convergence = 0

        for read_name, comp in self.compatibility_dict.items():
            for ref_name, score in comp.items():
                abundance_dict[ref_name] += score
                total += score

        for ref_name in abundance_dict.keys():
            abundance_dict[ref_name] = abundance_dict[ref_name] / total

            if self.em_round > 1:
                convergence += abs(self.abundance_dict[ref_name] - abundance_dict[ref_name])

        if self.em_round == 1:
            self.convergence = 1
        else:
            self.convergence = convergence

        return abundance_dict

    def _update_compatibility(self):
        """
        Update read-transcript compatibility based on transcript abundances
        """
        compatibility_dict = defaultdict(dict)

        for read_name, comp in self.compatibility_dict.items():
            total = 0
            for ref_name in comp.keys():
                total += self.abundance_dict[ref_name]

            for ref_name in comp.keys():
                compatibility_dict[read_name][ref_name] = self.abundance_dict[ref_name] / total

        return compatibility_dict

    def _get_tx_df(self):
        """
        Extract transcript info from bam file header
        """
        try:
            with pysam.AlignmentFile(self.alignment_file) as bam:
                references = bam.references
                lengths = bam.lengths

            return pd.DataFrame(index=references, data=lengths, columns=["transcript_length"])
        # If any error return empty DataFrame silently
        except Exception:
            return pd.DataFrame()
        