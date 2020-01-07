#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PANK4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PANK4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PANO1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PANX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PANX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PANX3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAOX 
