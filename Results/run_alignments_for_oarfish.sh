#!/bin/bash
#SBATCH --job-name=pb_aln_sim                       # Job name
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=sraghavendra@nygenome.org       # Where to send mail  
#SBATCH --mem=20G                                   # Job memory request
#SBATCH --output=stdout_%j.log                      # Standard output and error log

module load SAMtools


/gpfs/commons/home/sraghavendra/SOTA/oarfish/minimap2/minimap2 --eqx -N 100 -ax map-pb /gpfs/commons/home/sraghavendra/Simulation/lrgasp-simulation/sim_result/ref_data/human.transcripts.fasta\
 /gpfs/commons/home/sraghavendra/Simulation/lrgasp-simulation/sim_result/human_simulated_job_correct/PacBio.simulated.fasta \
 | samtools view -@4 -b -o /gpfs/commons/home/sraghavendra/SOTA/oarfish/alignments/sim1_aln_pb.bam

/gpfs/commons/home/sraghavendra/SOTA/oarfish/minimap2/minimap2 --eqx -N 100 -ax map-pb /gpfs/commons/home/sraghavendra/Simulation/lrgasp-simulation/sim_result_2/ref_data/human.transcripts.fasta\
 /gpfs/commons/home/sraghavendra/Simulation/lrgasp-simulation/sim_result_2/human_simulated_job_correct/PacBio.simulated.fasta \
 | samtools view -@4 -b -o /gpfs/commons/home/sraghavendra/SOTA/oarfish/alignments/sim2_aln_pb.bam

/gpfs/commons/home/sraghavendra/SOTA/oarfish/minimap2/minimap2 --eqx -N 100 -ax map-pb /gpfs/commons/home/sraghavendra/PacBio/reference/sota/transcriptome.fna \
 /gpfs/commons/home/sraghavendra/PacBio/reads/long/downsampled/ds_51.fastq \
 | samtools view -@4 -b -o /gpfs/commons/home/sraghavendra/SOTA/oarfish/alignments/aln_pb_51.bam

/gpfs/commons/home/sraghavendra/SOTA/oarfish/minimap2/minimap2 --eqx -N 100 -ax map-pb /gpfs/commons/home/sraghavendra/PacBio/reference/sota/transcriptome.fna \
 /gpfs/commons/home/sraghavendra/PacBio/reads/long/downsampled/ds_52.fastq \
 | samtools view -@4 -b -o /gpfs/commons/home/sraghavendra/SOTA/oarfish/alignments/aln_pb_52.bam