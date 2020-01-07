#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FDFT1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FDFT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FDPS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FDPSP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FDX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FDXACB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FDXR 
