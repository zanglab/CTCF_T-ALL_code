#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ALKBH5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ALKBH5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALKBH6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALKBH7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALKBH8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALLC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALMS1 
