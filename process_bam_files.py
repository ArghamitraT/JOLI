
# ~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~ #

# Third party imports
import pysam
import time
import pickle

# Local imports
from NanoCount.Read import Read
from NanoCount.common import *

crnt_tm = datetime.datetime.now()


# ~~~~~~~~~~~~~~MAIN FUNCTION~~~~~~~~~~~~~~ #
class process_bam:

    # ~~~~~~~~~~~~~~MAGIC METHODS~~~~~~~~~~~~~~ #
    def __init__(
            self,
            file_names = [],
            alignment_file: str = "",
            count_file: str = "",
            filter_bam_out: str = "",
            min_alignment_length: int = 50,
            keep_supplementary: bool = False,
            keep_neg_strand: bool = False,
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
            process: str ='theta', #(AT) 'expectation_log_theta' or 'theta'
            dirichlet_builtin: int = 1, 
            load: int = 0,
            load_filename: str = "",
            experiment_num: int = 4,
            parse_original: int = 0,
            process_dir: str = "",
            bam_dir: str = ""
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
        self.keep_supplementary = keep_supplementary
        self.max_dist_5_prime = max_dist_5_prime
        self.max_dist_3_prime = max_dist_3_prime
        self.keep_neg_strand = keep_neg_strand
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
        self.all_rnames = {}
        self.all_readName = {}
        self.theta_names = {}
        self.experiment_num = experiment_num
        self.dirichlet_builtin=dirichlet_builtin
        self.parse_original = parse_original
        self.process_dir = process_dir
        self.bam_dir = bam_dir


        print("Initialise Nanocount")
        print("Parse Bam file and filter low quality alignments")
        
        start = time.time()
        start1 = time.time() # comment

        if load:
            self.load_state(max_em_rounds, count_file)
        
        else:
            self.initialize_model()
        
        
    

    def initialize_model(self):
        
        start = time.time()
        for index, file_name in enumerate(self.file_names_list, start=1):
            main_dir = self.bam_dir
            file_name = main_dir+'/'+file_name
            print(f"sample_{index} {file_name}")
            
            if self.parse_original:
                read_dict, ref_len_dict = self._parse_bam_original(file_name=file_name)
            else:
                read_dict, ref_len_dict = self._parse_bam(file_name=file_name)

            main_dir =  self.process_dir 
            pkl_file_name_read_dict = file_name.split('/')[-1].split('.')[0] + '_read_dicts.pkl'
            pkl_file_name_ref_len_dict = file_name.split('/')[-1].split('.')[0] + '_ref_len_dicts.pkl'

            read_dict_path = os.path.join(main_dir, pkl_file_name_read_dict)
            with open(read_dict_path, 'wb') as file:
                    pickle.dump(read_dict, file)
                    print(f"Saved {read_dict_path}")


            # Check if the ref_len_dict file already exists
            ref_len_dict_path = os.path.join(main_dir, pkl_file_name_ref_len_dict)
            with open(ref_len_dict_path, 'wb') as file:
                    pickle.dump(ref_len_dict, file)
                    print(f"Saved {ref_len_dict_path}")

            if not os.path.exists(ref_len_dict_path):
                with open(ref_len_dict_path, 'wb') as file:
                    pickle.dump(ref_len_dict, file)
                    print(f"Saved {ref_len_dict_path}")
            else:
                print(f"File {ref_len_dict_path} already exists. Skipping save.")


            # Store the dictionaries in the universal dictionary with a sample key
            sample_key = f'sample{index}'
            self.all_read_dicts[sample_key] = read_dict
            self.all_ref_len_dicts[sample_key] = ref_len_dict

        

    # ~~~~~~~~~~~~~~PRIVATE METHODS~~~~~~~~~~~~~~ #

    def _parse_bam_original(self, file_name):
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
        # with pysam.AlignmentFile(self.alignment_file) as bam:
            # Collect reference lengths in dict
            for name, length in zip(bam.references, bam.lengths):
                ref_len_dict[name] = length

            for idx, alignment in enumerate(bam):
                if alignment.is_unmapped:
                    c["Discarded unmapped alignments"] += 1
                elif not self.keep_neg_strand and alignment.is_reverse:
                    c["Discarded negative strand alignments"] += 1
                elif not self.keep_supplementary and alignment.is_supplementary:
                    c["Discarded supplementary alignments"] += 1
                elif self.min_alignment_length > 0 and alignment.query_alignment_length < self.min_alignment_length:
                    c["Discarded short alignments"] += 1
                elif self.max_dist_3_prime >= 0 and alignment.reference_end <= ref_len_dict[alignment.reference_name] - self.max_dist_3_prime:
                    c["Discarded alignment with invalid 3 prime end"] += 1
                elif self.max_dist_5_prime >= 0 and alignment.reference_start >= self.max_dist_5_prime:
                    c["Discarded alignment with invalid 5 prime end"] += 1
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
            self.log.error("No valid secondary alignments found in bam file. Were the reads aligned with minimap `-p 0 -N 10` options ?")

        # Write filtered reads counters
        log_dict(d=c, logger=self.log.info, header="Summary of reads filtered")
        return filtered_read_dict, ref_len_dict

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

