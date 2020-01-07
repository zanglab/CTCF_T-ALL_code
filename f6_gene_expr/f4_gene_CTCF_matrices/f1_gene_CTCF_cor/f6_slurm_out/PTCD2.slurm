#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PTCD2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PTCD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTCD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTCH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTCH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTCHD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTCHD3 
