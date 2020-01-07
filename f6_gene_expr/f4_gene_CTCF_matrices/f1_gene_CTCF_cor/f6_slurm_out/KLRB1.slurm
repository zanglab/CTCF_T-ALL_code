#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KLRB1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KLRB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLRC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLRC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLRC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLRC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KLRC4-KLRK1 
