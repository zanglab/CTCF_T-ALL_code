#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DPH7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DPH7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DPM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DPM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DPM3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DPP10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DPP3 
