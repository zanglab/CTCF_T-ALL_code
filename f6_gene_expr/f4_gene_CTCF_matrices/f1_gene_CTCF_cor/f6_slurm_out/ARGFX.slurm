#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ARGFX.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ARGFX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARGFXP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARGLU1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGAP10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGAP11A 
