#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RADIL.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RADIL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAET1E 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAET1G 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAET1K 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAET1L 
