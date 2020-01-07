#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLF2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLFN11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLFN12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLFN12L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLFN13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLFN14 
