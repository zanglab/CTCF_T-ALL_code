#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PHGR1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PHGR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PHIP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PHKA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PHKA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PHKB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PHKG1 
