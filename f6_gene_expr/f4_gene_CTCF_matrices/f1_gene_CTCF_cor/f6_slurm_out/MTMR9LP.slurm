#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MTMR9LP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MTMR9LP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTNR1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTNR1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTO1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTOR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTPAP 
