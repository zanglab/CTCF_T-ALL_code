#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PGAM2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PGAM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PGAM4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PGAM5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PGAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PGAP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PGAP3 
