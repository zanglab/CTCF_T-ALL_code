#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LRRTM1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRTM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRTM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRTM3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRTM4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRSAM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRTM1 
