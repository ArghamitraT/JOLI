# This file contains helper functions

import os
import re

class ExperimentFileProcessor:
    def __init__(self, directory):
        self.directory = directory
        self.regex_patterns = {
    "exp1": {
        "simulation": {
            "pattern": r'output_Simulation_VIGD_token_(\d+)_sample(\d+)_file(\d+)_ds(\d+)num(\d+)aln(\d+)(long|short)_GDlr_(0\.\d+)_AlphaInitial_(\d+(?:\.\d+)?)_EMround_(\d+)_',
            "labels": ["token", "sample", "file_num1", "ds_percentage1", "num1", "aln_replica1", "length1", "GDlr", "AlphaInitial", "EMround"]
        },
        "real": {
            "pattern": r'output_PacIllu_VIGD_token_(\d+)_sample(\d+)_file(\d+)_ds(\d+)num(\d+)aln(\d+)(long|short)_GDlr_(0\.\d+)_AlphaInitial_(\d+(?:\.\d+)?)_EMround_(\d+)_',
            "labels": ["token", "sample", "file_num1", "ds_percentage1", "num1", "aln_replica1", "length1", "GDlr", "AlphaInitial", "EMround"]
        }
    },
    "exp2": {
        "simulation": {
            "pattern": r'output_Simulation_VIGD_token_(\d+)_sample(\d+)_file(\d+)_ds(\d+)num(\d+)aln(\d+)(long|short)_file(\d+)_ds(\d+)num(\d+)aln(\d+)(long|short)_GDlr_(0\.\d+)_AlphaInitial_(\d+(?:\.\d+)?)_EMround_(\d+)_',
            "labels": ["token", "sample", "file_num1", "ds_percentage1", "num1", "aln_replica1", "length1", "file_num2", "ds_percentage2", "num2", "aln_replica2", "length2", "GDlr", "AlphaInitial", "EMround"]
        },
        "real": {
            "pattern": r'output_PacIllu_VIGD_token_(\d+)_sample(\d+)_file(\d+)_ds(\d+)num(\d+)aln(\d+)(long|short)_file(\d+)_ds(\d+)num(\d+)aln(\d+)(long|short)_GDlr_(0\.\d+)_AlphaInitial_(\d+(?:\.\d+)?)_EMround_(\d+)_',
            "labels": ["token", "sample", "file_num1", "ds_percentage1", "num1", "aln_replica1", "length1", "file_num2", "ds_percentage2", "num2", "aln_replica2", "length2", "GDlr", "AlphaInitial", "EMround"]
        }
    },
    "exp4": {
        "simulation": {
            "pattern": r'output_Simulation_VIGD_token_(\d+)_sample(\d+)_file(\d+)_ds(\d+)num(\d+)aln(\d+)(long|short)_file(\d+)_ds(\d+)num(\d+)aln(\d+)(long|short)_GDlr_(0\.\d+)_AlphaInitial_(\d+(?:\.\d+)?)_EMround_(\d+)_',
            "labels": ["token", "sample", "file_num1", "ds_percentage1", "num1", "aln_replica1", "length1", "file_num2", "ds_percentage2", "num2", "aln_replica2", "length2", "GDlr", "AlphaInitial", "EMround"]
        },
        "real": {
            "pattern": r'output_PacIllu_VIGD_token_(\d+)_sample(\d+)_file(\d+)_ds(\d+)num(\d+)aln(\d+)(long|short)_file(\d+)_ds(\d+)num(\d+)aln(\d+)(long|short)_GDlr_(0\.\d+)_AlphaInitial_(\d+(?:\.\d+)?)_EMround_(\d+)_',
            "labels": ["token", "sample", "file_num1", "ds_percentage1", "num1", "aln_replica1", "length1", "file_num2", "ds_percentage2", "num2", "aln_replica2", "length2", "GDlr", "AlphaInitial", "EMround"]
        }
    },
    "exp5": {
        "simulation": {
            # "pattern": r'output_Simulation_VIGD_token_(\d+)_sample(\d+)_file(\d+)_.*_file(\d+)_.*_GDlr_(0\.\d+)_AlphaInitial_(\d+(?:\.\d+)?)_EMround_(\d+)_',
            # "labels": ["token", "sample", "file_num1", "file_num2", "GDlr", "AlphaInitial", "EMround"]
            "pattern": r'output_Simulation_VIGD_token_(\d+)_sample(\d+)_file(\d+)_(.*?)_file(\d+)_(.*?)_GDlr_(0\.\d+)_AlphaInitial_(\d+(?:\.\d+)?)_EMround_(\d+)_',
            "labels": ["token", "sample", "file_num1", "filename1", "file_num2", "filename2", "GDlr", "AlphaInitial", "EMround"]
        },
        "real": {
            "pattern": r'output_PacIllu_VIGD_token_(\d+)_sample(\d+)_file(\d+)_.*_file(\d+)_.*_GDlr_(0\.\d+)_AlphaInitial_(\d+(?:\.\d+)?)_EMround_(\d+)_',
            "labels": ["token", "sample", "file_num1", "file_num2", "GDlr", "AlphaInitial", "EMround"]
        }
    }
}

    def process_files(self, experiment, data_type, file):
        """
        Process file names for a specific experiment and data type.

        :param experiment: str, experiment name (e.g., 'exp1', 'exp4')
        :param data_type: str, 'simulation' or 'real'
        :return: Tuple (labels, extracted_data)
        """
        if experiment not in self.regex_patterns or data_type not in self.regex_patterns[experiment]:
            raise ValueError(f"Invalid experiment '{experiment}' or data type '{data_type}'.")

        # Get the pattern and labels
        pattern_info = self.regex_patterns[experiment][data_type]
        file_pattern = re.compile(pattern_info["pattern"])
        labels = pattern_info["labels"]

        # Extract data from matching files
        match = file_pattern.search(file)
        if match:
            extracted_data  = dict(zip(labels, match.groups()))
        
            # Add file_name1 based on experiment logic
            if experiment in ["exp2", "exp4"]:
                sample = extracted_data["sample"]
                if sample == '1':
                    file_num1 = extracted_data["num1"]
                    replica1 = extracted_data["aln_replica1"][1]
                    dayF1 = extracted_data["aln_replica1"][0]
                    length = extracted_data["length1"]
                    ds_prctF1 = extracted_data["ds_percentage1"]
                elif sample == '2':
                    file_num1 = extracted_data["num2"]
                    replica1 = extracted_data["aln_replica2"][1]
                    dayF1 = extracted_data["aln_replica2"][0]
                    length = extracted_data["length2"]
                    ds_prctF1 = extracted_data["ds_percentage2"]

                file_name1 = f'ds{ds_prctF1}num{file_num1}aln{dayF1}{replica1}{length}'

            elif experiment == "exp1":
                file_num1 = extracted_data["num1"]
                replica1 = extracted_data["aln_replica1"][1]
                dayF1 = extracted_data["aln_replica1"][0]
                ds_prctF1 = extracted_data["ds_percentage1"]
                length = extracted_data["length1"]

                file_name1 = f'ds{ds_prctF1}num{file_num1}aln{dayF1}{replica1}{length}'

            elif experiment == "exp5":
                sample = extracted_data["sample"]
                if sample == '1':
                    file_name1 = extracted_data["filename1"]
                elif sample == '2':
                    file_name1 = extracted_data["filename2"]
                pattern = re.compile(r'ds(\d+)num(\d+)aln(\d+)(long|short)ds(\d+)num(\d+)aln(\d+)(long|short)')
                match = pattern.search(file_name1)
                labels = ["ds_percentage1", "num1", "aln_replica1", "length1", "ds_percentage2", "num2", "aln_replica2", "length2"]
                extracted_data_filename  = dict(zip(labels, match.groups()))
                replica1 = extracted_data_filename["aln_replica1"][1]
                dayF1 = extracted_data_filename["aln_replica1"][0]
                
            else:
                raise ValueError(f"Unsupported experiment type '{experiment}'.")


            if experiment == "exp5":
                extracted_data["file_name"] = file_name1
                extracted_data["result_aln_replica"] = f"{dayF1}{replica1}"
            
            else:
                extracted_data["file_name"] = file_name1
                extracted_data["result_file_num"] = file_num1
                extracted_data["result_ds_percentage"] = ds_prctF1
                extracted_data["result_aln_replica"] = f"{dayF1}{replica1}"
                extracted_data["result_length"] = length
            
            return extracted_data
        else:
            return

        # for file in os.listdir(self.directory):
        #     match = file_pattern.search(file)
        #     if match:
        #         extracted_data.append(match.groups())

        
