#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o EXPH5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i EXPH5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EXT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EXT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EXTL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EXTL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EXTL3 
