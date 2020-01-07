#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PRLH.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PRLH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRLHR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRLR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRM3 
