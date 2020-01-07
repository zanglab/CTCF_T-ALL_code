#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SH3BP5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SH3BP5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SH3BP5L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SH3D19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SH3D21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SH3GL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SH3GL1P1 
