#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PAPOLG.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PAPOLG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAPPA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAPPA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAPSS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAPSS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAQR3 
