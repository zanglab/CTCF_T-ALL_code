#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MAF1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MAF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAFA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAFB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAFF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAFG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAFIP 
