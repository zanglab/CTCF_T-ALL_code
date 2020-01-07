#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PGM2L1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PGM2L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PGM3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PGM5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PGM5P2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PGP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PGPEP1 
