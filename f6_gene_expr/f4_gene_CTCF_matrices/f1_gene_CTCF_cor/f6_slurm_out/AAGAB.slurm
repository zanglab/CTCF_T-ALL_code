#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o AAGAB.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i AAGAB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AAK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AAMDC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AAMP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AANAT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AAR2 
