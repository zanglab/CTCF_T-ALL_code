#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PLCG1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PLCG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLCG2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLCH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLCH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLCL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLCL2 
