#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NHLH2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NHLH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NHLRC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NHLRC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NHLRC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NHLRC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NHP2 
