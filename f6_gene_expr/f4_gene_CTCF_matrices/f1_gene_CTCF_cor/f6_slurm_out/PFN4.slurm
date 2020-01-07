#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PFN4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PFN4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PGA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PGA4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PGA5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PGAM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PGAM1P5 
