#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DACT2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DACT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DACT3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DAD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DAG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DAGLA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DAGLB 
