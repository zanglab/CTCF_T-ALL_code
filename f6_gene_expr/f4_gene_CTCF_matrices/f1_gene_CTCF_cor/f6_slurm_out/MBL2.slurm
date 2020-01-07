#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MBL2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MBL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MBLAC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MBLAC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MBNL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MBNL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MBNL3 
