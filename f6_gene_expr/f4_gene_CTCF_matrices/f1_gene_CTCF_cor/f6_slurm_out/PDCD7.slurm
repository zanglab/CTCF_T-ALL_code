#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PDCD7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PDCD7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDCL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDCL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDCL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDCL3P4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDE10A 
