#!/bin/bash
#SBATCH --job-name=sr_kallisto_pip                  # Job name
#SBATCH --mail-type=END,FAIL                        # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=sraghavendra@nygenome.org       # Where to send mail
#SBATCH --mem=80G                                  # Job memory request. Different units can be specified using the suffix [K|M|G|T]
#SBATCH --output=stdout_%j.log                      # Standard output and error log


output_dir='output_short_sim1'
reads_file1='Illumina.simulated_1.fq'
reads_file2='Illumina.simulated_2.fq'
reads_dir='/gpfs/commons/home/sraghavendra/Simulation/lrgasp-simulation/sim_result/human_simulated_job_correct/'
index_file='new_index.idx'

# /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/kallisto_old/kallisto bus -x bulk --threshold 0.8 -t 32 --paired\
#  -i /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/human_k-31.idx\
#   /gpfs/commons/home/sraghavendra/Simulation/lrgasp-simulation/sim_result_2/human_simulated_job_correct/Illumina.simulated_1.fq \
#   /gpfs/commons/home/sraghavendra/Simulation/lrgasp-simulation/sim_result_2/human_simulated_job_correct/Illumina.simulated_2.fq \
#   -o /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/output_short_sim2/

/gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/kallisto/build/src/kallisto bus -x bulk --threshold 0.8 -t 32 --paired\
 -i /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$index_file\
  $reads_dir$reads_file1 \
  $reads_dir$reads_file2 \
  -o /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$output_dir/

/gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/bustools/build/src/bustools sort -t 30 /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$output_dir/output.bus\
 -o /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$output_dir/sorted.bus

/gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/bustools/build/src/bustools count /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$output_dir/sorted.bus\
					-t /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$output_dir/transcripts.txt\
					-e /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$output_dir/matrix.ec\
					-o /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$output_dir/count --cm -m \
					-g /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/t2g.txt #/gpfs/commons/home/sraghavendra/PacBio/reference/sota/genomic.gtf

/gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/kallisto/build/src/kallisto quant-tcc -t 32\
 /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$output_dir/count.mtx\
					-i /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$index_file\
					-e /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$output_dir/count.ec.txt\
					-o /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$output_dir/

# output_dir='output_short_51'
# reads_file1='ill_51_R1.fastq'
# reads_file2='ill_51_R2.fastq'
# reads_dir='/gpfs/commons/home/sraghavendra/PacBio/reads/short/'
# index_file='new_index.idx'

# /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/kallisto/build/src/kallisto bus -x bulk --threshold 0.8 -t 32 --paired\
#  -i /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$index_file\
#   $reads_dir$reads_file1 \
#   $reads_dir$reads_file2 \
#   -o /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$output_dir/

# /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/bustools/build/src/bustools sort -t 30 /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$output_dir/output.bus\
#  -o /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$output_dir/sorted.bus

# /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/bustools/build/src/bustools count /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$output_dir/sorted.bus\
# 					-t /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$output_dir/transcripts.txt\
# 					-e /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$output_dir/matrix.ec\
# 					-o /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$output_dir/count --cm -m \
# 					-g /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/t2g.txt #/gpfs/commons/home/sraghavendra/PacBio/reference/sota/genomic.gtf

# /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/kallisto/build/src/kallisto quant-tcc -t 32\
#  /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$output_dir/count.mtx\
# 					-i /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$index_file\
# 					-e /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$output_dir/count.ec.txt\
# 					-o /gpfs/commons/home/sraghavendra/SOTA/lr-kallisto/$output_dir/