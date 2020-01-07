#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TPGS2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TPGS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TPH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TPH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TPI1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TPI1P2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TPI1P3 
