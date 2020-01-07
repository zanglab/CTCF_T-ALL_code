#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ECSCR.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ECSCR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ECSIT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ECT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ECT2L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EDA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EDA2R 
