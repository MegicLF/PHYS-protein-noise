#!/bin/bash
#SBATCH --mail-user=john.j.vastola@vanderbilt.edu,linghui.feng@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=3G
#SBATCH --time=4:00:00
#SBATCH --output=kld_out.txt

module load Anaconda3
python cleanKLD_main.py
