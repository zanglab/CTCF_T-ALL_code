#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ELP2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ELP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ELP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ELP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ELP5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ELP6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ELSPBP1 
