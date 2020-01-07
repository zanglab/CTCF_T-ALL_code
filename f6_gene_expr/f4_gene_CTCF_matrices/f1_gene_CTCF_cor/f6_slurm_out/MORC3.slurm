#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MORC3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MORC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MORC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MORF4L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MORF4L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MORN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MORN2 
