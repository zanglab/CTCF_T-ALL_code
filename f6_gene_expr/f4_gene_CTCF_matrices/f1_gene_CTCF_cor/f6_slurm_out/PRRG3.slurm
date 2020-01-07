#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PRRG3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PRRG3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRRG4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRRT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRRT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRRT3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRRT4 
