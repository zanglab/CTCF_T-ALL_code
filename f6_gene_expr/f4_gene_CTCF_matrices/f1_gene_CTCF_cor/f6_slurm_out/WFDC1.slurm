#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o WFDC1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i WFDC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WFDC10A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WFDC10B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WFDC11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WFDC12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WFDC13 
