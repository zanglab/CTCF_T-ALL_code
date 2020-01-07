#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MFSD5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MFSD5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MFSD6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MFSD6L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MFSD8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MFSD9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MGA 
