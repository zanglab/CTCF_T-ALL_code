#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SH2D4A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SH2D4A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SH2D4B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SH2D5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SH2D6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SH2D7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SH3BGR 
