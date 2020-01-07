#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SCRT2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SCRT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCTR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCUBE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCUBE2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCUBE3 
