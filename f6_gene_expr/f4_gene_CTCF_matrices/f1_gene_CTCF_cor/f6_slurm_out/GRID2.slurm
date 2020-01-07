#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GRID2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GRID2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRID2IP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRIFIN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRIK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRIK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRIK3 
