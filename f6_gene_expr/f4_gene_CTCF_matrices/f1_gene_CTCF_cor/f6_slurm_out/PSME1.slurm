#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PSME1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PSME1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSME2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSME3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSME4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSMF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSMG1 
