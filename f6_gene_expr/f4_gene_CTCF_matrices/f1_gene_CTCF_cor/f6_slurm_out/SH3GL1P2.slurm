#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SH3GL1P2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SH3GL1P2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SH3GL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SH3GL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SH3GLB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SH3GLB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SH3KBP1 
