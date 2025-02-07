# JOLI: JOint Long and short Isoform quantification

JOLI builds upon the standard **Expectation-Maximization (EM) framework** by incorporating two key extensions:  
1️⃣ **Joint analysis of short-read (SR) and long-read (LR) data**  
2️⃣ **Sharing information between multiple similar samples ("multisample" analysis)**  

Our proposed method simultaneously leverages read information from multiple short and long-read samples, exploiting their **complementarity** to improve accuracy, consistency, and reproducibility.  

JOLI operates on multiple related samples (e.g., from the same tissue type), under the assumption that the **true isoform abundances for these samples are correlated**.  

Using an **empirical Bayes framework**, we learn a common prior across the samples that allows the model to capture shared variability and facilitates **robust and consistent inference**.  

To estimate the **Maximum A Posteriori (MAP) parameters**, we employ:  
- The **EM algorithm** within each sample  
- **Gradient descent** for the shared prior parameters  

JOLI builds upon the foundations of  **[NanoCount](https://github.com/a-slide/NanoCount)** as part of its development.

## Technical files, folders and their functionality

### **main.py**
This file is the entry point for launching the main program.

### **process_bam_files.py**
This file is used to process minimap aligned bam files

### **JOLI.py**
This file runs the EM algorithm of JOLI

### **DirichletOptimizer_vector.py**
This file updates the shared dirichlet (empirical Bayes framework) prior with Gradient descent

### **Results**
Holds the scripts for result analysis

### **Sample Data**
To test the script with **sample data**, you can download the dataset from the following link:

🔗 **[Download Sample Data](https://drive.google.com/drive/folders/1aWmm-ZAsBTUhnqNyg3uNMNIMKvaV3QPN?usp=drive_link)**

Once downloaded, place the files in the appropriate directory before running the script.


### **Environments**
Has .yml and .txt files required to install and run JOLI

## Running the Program

To set up and run **JOLI**, follow these steps:

### Clone and Environment setup 
1. clone the repository to your local machine: `git clone https://github.com/ArghamitraT/JOLI.git`
2. Inside the environments folder, you will find the following files:

```bash
    JOLI_conda_requirements.txt
    JOLI_pip_requirements.txt
    JOLI.yml
```
    
You can create the Conda environment using either of the following methods:
🔹 Option 1: Using JOLI.yml (Recommended)

This method installs all dependencies from a single YML file.

```bash
conda env create -f Environments/JOLI.yml
conda activate JOLI
```

🔹 Option 2: Using Conda & Pip Requirements Files

Alternatively, you can manually install dependencies step by step:

```bash
conda create --name JOLI --file Environments/JOLI_conda_requirements.txt
conda activate JOLI
pip install -r Environments/JOLI_pip_requirements.txt
```

### Command Line Arguments

JOLI provides several command-line flags to configure the execution. Below are the available options:

- **`--data_folder`** (`str`): Path for the processed data (output of `process_bam_files.py`, pkl files).  
- **`--output_path`** (`str`): Path for the output abundance file.  
- **`--sample1`** (`str (list)`): Sample1 file name(s).  
- **`--sample2`** (`str (list)`): Sample2 file name(s).  
- **`--GD_lr`** (`float`, default=`0.01`): Learning rate for Dirichlet gradient descent.  
- **`--alpha_initial`** (`float`, default=`1`): The fixed sum value of alpha. For single sample (SS), use `1`; for multi-sample (MS), depends on number of isoforms, for 150k isoforms we used `10e5`.  
- **`--max_em_rounds`** (`int`, default=`30`): The maximum number of EM iterations.  
- **`--experiment_num`** (`int`): Defines the experiment setup:
  - `1` for single sample (JOLI SS)
  - `4` for multi-sample (JOLI MS)
- **`--process_bam_required`** (`int`): If set to `1`, only processes the BAM files and saves them as `.pkl` files (EM will not run).  
- **`--process_saved_dir`** (`str`): Path to save the processed BAM files.  
- **`--bam_dir`** (`str`): Path to the original BAM files.  

---

### Alignment of fasta reads to bam files
JOLI needs Minimap2 or STAR aligned .bam files. We used Minimap2 to align our Long reads. Our used command to align the fasta read files with transcriptome was `minimap2 -ax splice -uf -C5 transcriptome.fna reads.fq > aln.sam`; then use samtools to convert the .sam files into .bamf files: `samtools view -bS aln.sam > aln.bam`. To align the short reads we used STAR. Commands are:

```bash
STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir <output dir of STAR> \
--limitGenomeGenerateRAM 150000000000 \
--genomeFastaFiles <transcriptome .fna file>
```


```bash
STAR --genomeDir <output dir of STAR from previous step> \
--runThreadN 6 \
--readFilesIn <pair_ended_read_fastaFile_1> <pair_ended_read_fastaFile_2> \
--outFileNamePrefix <bam file dir> \
--outSAMtype BAM SortedByCoordinate
```

### Processing bam files

To process .bam files you can run following command. This script should give 2 pkl files as output: filename_read_dicts.pkl and filename_ref_len_dicts.pkl

```bash
cd <directory of the repo>
python main.py --process_bam_required 1 \
               --process_saved_dir '<folder to save processed pkl files>' \
               --bam_dir '<folder containing BAM files>' \
               --sample1 '<BAM file 1>' '<BAM file 2>'  # (without .bam extension)

```
Example:
```bash
cd /Users/jondoe/Documents/JOLI
python main.py --process_bam_required 1 \
               --process_saved_dir '/Users/jondoe/Documents/pklfiles/' \
               --bam_dir '/Users/jondoe/Documents/bamfiles/' \
               --sample1 'aln_E0' 'aln_E2'
```

### Run JOLI for Isoform Quantification
The output should be .tsv files with isoform abundance. The .tsv file has 4 columns: transcript_name, raw, est_count, and tpm. 

Example: Run JOLI SS
```bash
cd /Users/jondoe/Documents/JOLI
python main.py --data_folder '/Users/jondoe/Documents/pklfiles/' \
               --output_path '/Users/jondoe/Documents/results/' \
               --sample1 'aln_E0' \
               --experiment_num 1\
               --alpha_initial 1\
               --max_em_rounds 30
```
The naming of the JOLI SS output would be "sample1_file1_<sample1_file_name>_alnE0short_GDlr_001_AlphaInitial_1_EMround_<date&time>.tsv" which would contain the abundanc of sample 1 file


Example: Run JOLI MS
```bash
cd /Users/jondoe/Documents/JOLI
python main.py --data_folder '/Users/jondoe/Documents/pklfiles/' \
               --output_path '/Users/jondoe/Documents/results/' \
               --sample1 'aln_E0' \
               --sample2 'aln_E2' \
               --experiment_num 4\
               --alpha_initial 10000\
               --max_em_rounds 30
```
JOLI MS would give 2 tsv files as output. The naming would be 
"sample1_file1_<sample1_file_name>_file2_<sample2_file_name>_GDlr_001_AlphaInitial_1_EMround_30_<date&time>.tsv"
"sample2_file1_<sample1_file_name>_file2_<sample2_file_name>_GDlr_001_AlphaInitial_1_EMround_30_<date&time>.tsv"

Here sample1 means the file contains the abundance of file1_<sample1_file_name> and sample2 means it has the abundance of file2_<sample2_file_name>
