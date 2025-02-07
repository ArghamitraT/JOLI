#!/bin/bash
#SBATCH --job-name=oar_quant                        # Job name
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=sraghavendra@nygenome.org       # Where to send mail  
#SBATCH --mem=50G                                   # Job memory request. Different units can be specified using the suffix [K|M|G|T]
#SBATCH --output=stdout_%j.log                      # Standard output and error log

~/miniconda3/bin/conda init bash
source /gpfs/commons/home/sraghavendra/.bashrc
conda activate mpaqt

# oarfish -j 6 -a /gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/PacBio_data_Liz/transcriptome_aln/subfolder/ds_10_num1_aln_01_long.bam -o day0_rep1 --filter-group no-filters --model-coverage

# oarfish -j 16 --reads /gpfs/commons/home/sraghavendra/PacBio/reads/long/downsampled/ds_01.fastq --reference /gpfs/commons/home/sraghavendra/PacBio/reference/sota/genomic.fna --seq-tech pac-bio -o /gpfs/commons/home/sraghavendra/SOTA/oarfish/day0_rep1/ --filter-group no-filters --model-coverage

# oarfish -j 16 --reads /gpfs/commons/home/sraghavendra/PacBio/reads/long/downsampled/ds_02.fastq --reference /gpfs/commons/home/sraghavendra/PacBio/reference/sota/genomic.fna --seq-tech pac-bio -o /gpfs/commons/home/sraghavendra/SOTA/oarfish/day0_rep2/ --filter-group no-filters --model-coverage

# oarfish -j 16 --reads /gpfs/commons/home/sraghavendra/PacBio/reads/long/downsampled/ds_51.fastq --reference /gpfs/commons/home/sraghavendra/PacBio/reference/sota/genomic.fna --seq-tech pac-bio -o /gpfs/commons/home/sraghavendra/SOTA/oarfish/day5_rep1/ --filter-group no-filters --model-coverage

# oarfish -j 16 --reads /gpfs/commons/home/sraghavendra/PacBio/reads/long/downsampled/ds_52.fastq --reference /gpfs/commons/home/sraghavendra/PacBio/reference/sota/genomic.fna --seq-tech pac-bio -o /gpfs/commons/home/sraghavendra/SOTA/oarfish/day5_rep2/ --filter-group no-filters --model-coverage

# oarfish -j 16 --reads /gpfs/commons/home/sraghavendra/Simulation/lrgasp-simulation/sim_result/human_simulated_job_correct/PacBio.simulated.fasta --reference /gpfs/commons/home/sraghavendra/Simulation/lrgasp-simulation/sim_result/ref_data/human.transcripts.fasta --seq-tech pac-bio -o /gpfs/commons/home/sraghavendra/SOTA/oarfish/sim1/ --filter-group no-filters --model-coverage

# oarfish -j 16 --reads /gpfs/commons/home/sraghavendra/Simulation/lrgasp-simulation/sim_result_2/human_simulated_job_correct/PacBio.simulated.fasta --reference /gpfs/commons/home/sraghavendra/Simulation/lrgasp-simulation/sim_result_2/ref_data/human.transcripts.fasta --seq-tech pac-bio -o /gpfs/commons/home/sraghavendra/SOTA/oarfish/sim2/ --filter-group no-filters --model-coverage

# oarfish -j 16 -a /gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/PacBio_data_Liz/transcriptome_aln/subfolder/ds_10_num1_aln_02_long.bam -o day0_rep2 --filter-group no-filters --model-coverage

# oarfish -j 16 -a /gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/PacBio_data_Liz/transcriptome_aln/subfolder/ds_10_num1_aln_51_long.bam -o day5_rep1 --filter-group no-filters --model-coverage

# oarfish -j 16 -a /gpfs/commons/home/spark/knowles_lab/Argha/RNA_Splicing/data/PacBio_data_Liz/transcriptome_aln/subfolder/ds_10_num1_aln_52_long.bam -o day5_rep2 --filter-group no-filters --model-coverage

# oarfish -j 16 -a /gpfs/commons/home/sraghavendra/SOTA/oarfish/alignments/sim1_aln_pb.bam -o /gpfs/commons/home/sraghavendra/SOTA/oarfish/sim1/ --filter-group no-filters --model-coverage

# oarfish -j 16 -a /gpfs/commons/home/sraghavendra/SOTA/oarfish/alignments/sim2_aln_pb.bam -o /gpfs/commons/home/sraghavendra/SOTA/oarfish/sim2/ --filter-group no-filters --model-coverage

mkdir /gpfs/commons/home/sraghavendra/SOTA/oarfish/day5_rep1_new /gpfs/commons/home/sraghavendra/SOTA/oarfish/day5_rep2_new 
mkdir /gpfs/commons/home/sraghavendra/SOTA/oarfish/day5_rep1_hifi /gpfs/commons/home/sraghavendra/SOTA/oarfish/day5_rep2_hifi

oarfish -j 16 -a /gpfs/commons/home/sraghavendra/SOTA/oarfish/alignments/aln_pb_51_new.bam -o /gpfs/commons/home/sraghavendra/SOTA/oarfish/day5_rep1_new/ --filter-group no-filters --model-coverage

oarfish -j 16 -a /gpfs/commons/home/sraghavendra/SOTA/oarfish/alignments/aln_pb_52_new.bam -o /gpfs/commons/home/sraghavendra/SOTA/oarfish/day5_rep2_new/ --filter-group no-filters --model-coverage

# oarfish -j 16 -a /gpfs/commons/home/sraghavendra/Simulation/alignments/aln_pb_2.bam -o sim2 --filter-group no-filters --model-coverage

oarfish -j 16 -a /gpfs/commons/home/sraghavendra/SOTA/oarfish/alignments/aln_pb_51_hifi.bam -o /gpfs/commons/home/sraghavendra/SOTA/oarfish/day5_rep1_hifi/ --filter-group no-filters --model-coverage

oarfish -j 16 -a /gpfs/commons/home/sraghavendra/SOTA/oarfish/alignments/aln_pb_52_hifi.bam -o /gpfs/commons/home/sraghavendra/SOTA/oarfish/day5_rep2_hifi/ --filter-group no-filters --model-coverage