#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HIST1H2BD.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H2BD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H2BE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H2BF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H2BG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H2BH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H2BI 
