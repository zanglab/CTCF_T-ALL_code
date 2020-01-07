#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SYT1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SYT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYT10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYT11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYT12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYT13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SYT14 
