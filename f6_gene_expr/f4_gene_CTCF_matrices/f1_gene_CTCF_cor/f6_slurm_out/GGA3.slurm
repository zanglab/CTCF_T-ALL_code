#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GGA3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GGA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GGACT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GGCT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GGCX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GGH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GGN 
