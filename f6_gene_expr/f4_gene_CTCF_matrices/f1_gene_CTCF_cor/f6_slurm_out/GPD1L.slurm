#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GPD1L.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GPD1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPER1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPHA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPHB5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPHN 
