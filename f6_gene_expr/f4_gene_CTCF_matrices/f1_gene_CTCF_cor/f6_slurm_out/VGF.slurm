#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o VGF.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i VGF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VGLL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VGLL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VGLL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VGLL4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VHL 
