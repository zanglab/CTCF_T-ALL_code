#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PXMP4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PXMP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PXN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PXT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PXYLP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PYCARD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PYCR1 
