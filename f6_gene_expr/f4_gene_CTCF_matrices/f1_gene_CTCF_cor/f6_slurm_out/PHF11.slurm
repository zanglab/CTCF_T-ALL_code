#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PHF11.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PHF11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PHF12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PHF13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PHF14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PHF19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PHF2 
