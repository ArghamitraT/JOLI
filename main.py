from JOLI import Expec_Max
import os
import argparse
from process_bam_files import process_bam


def create_load_file_path():
    default_load_filepath = '/gpfs/commons/home/atalukder/RNA_Splicing/files/results/exprmnt_2024_08_10__02_05_36/weights'
    # Initialize an empty string to store the result
    fileName = "allWeights"

    # Iterate through the list of file names
    for index, file_path in enumerate(file_names_list, start=1):
        # Split the file path by '/' and take the last part (the file name)
        file_name = file_path.split('/')[-1]
        # Extract a specific part of the file name if necessary (e.g., removing extension)
        file_identifier = ''.join(file_name.split('_')).split('.')[0]
        # Construct the string
        fileName += f"_file{index}_{file_identifier}"
    fileName = f"{fileName}_GDlr_{GD_lr}_AlphaInitial_{alpha_initial}.0_EMround_{last_EM_round}"
    file_name = fileName.strip()
    
    # Search for files that match this pattern in the specified directory
    file_list = []
    for file in os.listdir(default_load_filepath):
        if file.startswith(file_name) and file.endswith('.pkl'):
            file_list.append(os.path.join(default_load_filepath, file))
    
    if not file_list:
        print(f"No file found matching the format {fileName}'")
        return
    elif len(file_list)>1:
        print(f"More than 1 file found matching the format {fileName}'")
        return
    else:
        return file_list[0]

######### parameters #############

# hardcoded
EM_default = 'MAP'
load = 0
load_filename = "generic"
EM_type = EM_default
dirichlet_builtin = 0
process = 'theta'

# process mode default paths
process_bam_required_default = 0
parse_original = 0 #if 1 then, runs original nanocount parsing
process_saved_dir_default = '/gpfs/commons/home/atalukder/RNA_Splicing/data/Argha/RNA_Splicing/data/trial/pklfiles/' 
bam_dir_default = '/gpfs/commons/home/atalukder/RNA_Splicing/data/Argha/RNA_Splicing/data/trial/'

# quantification mode default paths
default_sample1 =  ['aln_E00']
default_sample2 =  ['aln_E22']
data_folder_default = '/gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/PacBio_data_Liz/transcriptome_aln_pklfiles/'
# output_file_default = os.path.join(os.getcwd(), 'results/exprmntSingleRun_2024_00_00__00_00_00/files/output_files/outputTRIAL_PacIllu_VIGD_token_00000')
experiment_num_default = 4
output_file_default = '/gpfs/commons/home/atalukder/RNA_Splicing/data/Argha/RNA_Splicing/data/trial/pklfiles/'

parser = argparse.ArgumentParser(description="Process BAM files and output results.")
parser.add_argument("--data_folder", type=str, default=data_folder_default,
                    help="Path for the processed.")
parser.add_argument("--output_path", type=str, default=output_file_default,
                    help="Path for the output file.")
parser.add_argument("--sample1", type=str, nargs='+', default=default_sample1, help="Sample1 (LR) file name(s). as a list []")
parser.add_argument("--sample2", type=str, nargs='+', default=default_sample2, help="Sample2 (SR) file name(s). as a list []")
parser.add_argument("--GD_lr", type=float, default=0.01, help="Learning rate for dirichlet gradient descent.")
parser.add_argument("--alpha_initial", type=float, default=1, help="The fixed sum value of alpha")
parser.add_argument("--max_em_rounds", type=int, default=30, help="The maximum EM iterations")
parser.add_argument("--experiment_num", type=int, default=experiment_num_default, help="Different experiment setup, 1: for 1 sample, 2 for merged, 4 for multisample, 5 for merged multisample")

parser.add_argument("--process_bam_required", type=int, default=process_bam_required_default, help="if 1 just process the files and save as pkl files. do not run the EM")
parser.add_argument("--process_saved_dir",  type=str, default=process_saved_dir_default,
                    help="Path to save the processed bam files")
parser.add_argument("--bam_dir",  type=str, default=bam_dir_default,
                    help="Path of the original bam files")

######### parameters ##############

# Parse the arguments
args = parser.parse_args()
input_folder = args.data_folder
output_file = args.output_path
experiment_num = args.experiment_num
sample1 = [input_folder +file for file in args.sample1]
sample2 = [input_folder +file for file in args.sample2]
GD_lr = args.GD_lr
alpha_initial = args.alpha_initial
max_em_rounds = args.max_em_rounds


for file in sample1:
    print(file)


# Print all the parameters
last_EM_round = 25
process_bam_required = args.process_bam_required
if process_bam_required:
    experiment_num = 1
else:
    experiment_num = experiment_num

if experiment_num == 1:
    file_names_list = [sample1]
else:
    file_names_list = [sample1, sample2]

print("Output_file_path", output_file)
print("GD_lr", GD_lr)
print("alpha_initial", alpha_initial)
print("max_em_rounds", max_em_rounds)
print("dirichlet_builtin", dirichlet_builtin)
print("Inference_type", EM_type)
print("dirichlet_process", process)

print("experiment_num", experiment_num)
if experiment_num == 1:
    print("Single sample, no gradient decent")
elif experiment_num == 2:
    print("2 samples merged as one, no gradient decent")
elif experiment_num == 4:
    print(print("Multi-sample with gradient decent"))
elif experiment_num == 5:
    print(print("Merged multi-sample with gradient decent"))


if load:
    print("load_filename", load_filename)
    if load_filename == "generic":
        load_filename= create_load_file_path()

process_dir = args.process_saved_dir
bam_dir = args.bam_dir

if process_bam_required:
    print ("data processing mode.......................")
    file_names_list = [item for sublist in file_names_list for item in sublist]
    file_names_list = [file.split('/')[-1]+'.bam' for file in file_names_list]
    process_bam (file_names=file_names_list, 
           count_file=output_file, 
           GD_lr=GD_lr, 
           alpha_initial=alpha_initial, 
           max_em_rounds=max_em_rounds,
           load=load,
           load_filename=load_filename,
           experiment_num = 4,
           dirichlet_builtin = dirichlet_builtin,
           parse_original = parse_original,
           process_dir = process_dir,
           bam_dir = bam_dir)

else:
    Expec_Max (file_names=file_names_list, 
            count_file=output_file, 
            GD_lr=GD_lr, 
            alpha_initial=alpha_initial, 
            max_em_rounds=max_em_rounds,
            load=load,
            load_filename=load_filename,
            experiment_num = experiment_num,
            dirichlet_builtin = dirichlet_builtin,
            EM_type=EM_type,
            process=process)

print("#########END###########")
