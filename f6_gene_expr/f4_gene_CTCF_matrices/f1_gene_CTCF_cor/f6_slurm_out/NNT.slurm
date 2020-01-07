#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NNT.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NNT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NOA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NOB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NOBOX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NOC2L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NOC3L 
