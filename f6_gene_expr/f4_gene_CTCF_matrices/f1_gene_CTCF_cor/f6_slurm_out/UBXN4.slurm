#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o UBXN4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i UBXN4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBXN6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBXN7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBXN8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UCHL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UCHL3 
