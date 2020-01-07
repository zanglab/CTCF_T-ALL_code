#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NOC4L.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NOC4L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NOCT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NOD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NOD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NODAL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NOG 
