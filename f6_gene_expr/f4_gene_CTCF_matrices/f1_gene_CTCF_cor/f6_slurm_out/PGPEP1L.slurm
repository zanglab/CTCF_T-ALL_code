#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PGPEP1L.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PGPEP1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PGR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PGRMC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PGRMC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PGS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PHACTR1 
