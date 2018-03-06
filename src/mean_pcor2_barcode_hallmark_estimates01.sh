#!/bin/sh
#SBATCH -J mean_pcor2_hallmark # A single job name for the array
#SBATCH -c 8 # Number of cores
#SBATCH -n 1 # Number of cores
#SBATCH -N 1 # All cores on one machine
#SBATCH -p serial_requeue,irizarry # Partition
#SBATCH --mem 24000 # Memory request
#SBATCH -t 0-06:35 # (D-HH:MM)
#SBATCH -o /net/irizarryfs01/srv/export/irizarryfs01/share_root/ypitajuarez/PCxN/log/mean_pcor2_hallmark_%a.out # Standard output
#SBATCH -e /net/irizarryfs01/srv/export/irizarryfs01/share_root/ypitajuarez/PCxN/log/mean_pcor2_hallmark_%a.err # Standard error
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=yered.h@gmail.com # Email to which notifications will be sent


source new-modules.sh
module load R/3.2.2-fasrc01
export R_LIBS_USER=$HOME/apps/R/3.2.2-fasrc01:$R_LIBS_USE

R CMD BATCH "--args ${SLURM_ARRAY_TASK_ID} ${SLURM_CPUS_PER_TASK}" /net/irizarryfs01/srv/export/irizarryfs01/share_root/ypitajuarez/PCxN/src/mean_pcor2_barcode/mean_pcor2_barcode_hallmark_estimates01.R /net/irizarryfs01/srv/export/irizarryfs01/share_root/ypitajuarez/PCxN/log/mean_pcor2_hallmark_"${SLURM_ARRAY_TASK_ID}".Rout
