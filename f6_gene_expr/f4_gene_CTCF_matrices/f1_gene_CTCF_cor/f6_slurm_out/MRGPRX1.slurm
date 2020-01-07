#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MRGPRX1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MRGPRX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRGPRX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRGPRX3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRGPRX4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRI1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRLN 
