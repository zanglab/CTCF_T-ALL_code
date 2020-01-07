#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ASCL5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ASCL5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASF1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASF1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASGR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASGR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ASH1L 
