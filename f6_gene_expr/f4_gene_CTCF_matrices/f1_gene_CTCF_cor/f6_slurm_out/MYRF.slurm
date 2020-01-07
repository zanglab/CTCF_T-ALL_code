#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MYRF.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MYRF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYRFL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYRIP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYSM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYT1L 
