#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GAL.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GAL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GAL3ST1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GAL3ST2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GAL3ST3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GAL3ST4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GALC 
