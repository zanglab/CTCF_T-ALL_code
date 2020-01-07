#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DMRTC2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DMRTC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DMTF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DMTN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DMWD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DMXL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DMXL2 
