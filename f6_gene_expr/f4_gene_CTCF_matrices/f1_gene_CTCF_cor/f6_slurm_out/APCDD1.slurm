#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o APCDD1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i APCDD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APCDD1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APCS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APEH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APELA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i APEX1 
