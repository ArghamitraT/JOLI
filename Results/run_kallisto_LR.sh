#!/bin/bash
#SBATCH --job-name=lr_kallisto_pip                  # Job name
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=sraghavendra@nygenome.org       # Where to send mail
#SBATCH --mem=100G                                  # Job memory request. Different units can be specified using the suffix [K|M|G|T]
#SBATCH --output=stdout_%j.log                      # Standard output and error log

~/miniconda3/bin/conda init bash
source /gpfs/commons/home/sraghavendra/.bashrc
conda activate NanoCount_5

###### INDEX ###########
# kb ref -k 63 -i new_index.idx -g /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/t2g.txt\
#  -f1 /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/transcriptome.fasta \
#   /gpfs/commons/home/sraghavendra/PacBio/reference/sota/genomic.fna\
#    /gpfs/commons/home/sraghavendra/PacBio/reference/sota/genomic.gtf

#######################################

output_dir='output_sim2'
reads_file='PacBio.simulated.fasta'
reads_dir='/gpfs/commons/home/sraghavendra/Simulation/lrgasp-simulation/sim_result_2/human_simulated_job_correct/'
index_file='new_index.idx'

/gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/kallisto/build/src/kallisto bus -x bulk --threshold 0.8 -t 32 --long --unmapped\
 -i /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$index_file\
  $reads_dir$reads_file\
  -o /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$output_dir/

/gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/bustools/build/src/bustools sort -t 30 /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$output_dir/output.bus\
 -o /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$output_dir/sorted.bus

/gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/bustools/build/src/bustools count /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$output_dir/sorted.bus\
					-t /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$output_dir/transcripts.txt\
					-e /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$output_dir/matrix.ec\
					-o /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$output_dir/count --cm -m \
					-g /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/t2g.txt #/gpfs/commons/home/sraghavendra/PacBio/reference/sota/genomic.gtf ## TO CHANGE 

/gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/kallisto/build/src/kallisto quant-tcc -t 32 --long -P PacBio\
 /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$output_dir/count.mtx\
					-i /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$index_file\
					-e /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$output_dir/count.ec.txt\
					-o /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$output_dir/

output_dir='output_52'
reads_file='ds_52.fastq'
reads_dir='/gpfs/commons/home/sraghavendra/PacBio/reads/long/downsampled/'
index_file='new_index.idx'

/gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/kallisto/build/src/kallisto bus -x bulk --threshold 0.8 -t 32 --long --unmapped\
 -i /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$index_file\
  $reads_dir$reads_file\
  -o /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$output_dir/

/gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/bustools/build/src/bustools sort -t 30 /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$output_dir/output.bus\
 -o /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$output_dir/sorted.bus

/gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/bustools/build/src/bustools count /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$output_dir/sorted.bus\
					-t /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$output_dir/transcripts.txt\
					-e /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$output_dir/matrix.ec\
					-o /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$output_dir/count --cm -m \
					-g /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/t2g.txt #/gpfs/commons/home/sraghavendra/PacBio/reference/sota/genomic.gtf ## TO CHANGE 

/gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/kallisto/build/src/kallisto quant-tcc -t 32 --long -P PacBio\
 /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$output_dir/count.mtx\
					-i /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$index_file\
					-e /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$output_dir/count.ec.txt\
					-o /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$output_dir/