#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PGD.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PGD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PGF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PGGT1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PGK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PGK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PGLS 
