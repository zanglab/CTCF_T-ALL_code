#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HIST1H2AE.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H2AE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H2AG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H2AH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H2AI 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H2AJ 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H2AK 
