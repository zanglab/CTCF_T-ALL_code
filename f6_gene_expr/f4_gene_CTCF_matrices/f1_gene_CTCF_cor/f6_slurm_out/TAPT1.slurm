#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TAPT1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TAPT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TARBP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TARBP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TARDBP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TARM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TARS 
