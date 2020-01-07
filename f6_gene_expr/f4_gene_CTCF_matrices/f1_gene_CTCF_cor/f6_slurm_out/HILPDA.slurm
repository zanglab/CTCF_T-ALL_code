#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HILPDA.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HILPDA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HILS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HINFP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HINT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HINT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HINT3 
