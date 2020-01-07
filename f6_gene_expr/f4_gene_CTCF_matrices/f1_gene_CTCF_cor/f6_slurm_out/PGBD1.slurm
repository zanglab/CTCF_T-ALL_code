#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PGBD1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PGBD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PGBD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PGBD4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PGBD5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PGC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PGCP1 
