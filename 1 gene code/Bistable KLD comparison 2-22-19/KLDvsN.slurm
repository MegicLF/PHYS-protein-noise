#!/bin/bash

#SBATCH --mail-user=john.j.vastola@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --mem=4G
#SBATCH --output=KLDvsN_job.out

module load Anaconda3


#########################

python KLDvsN_best.py 

