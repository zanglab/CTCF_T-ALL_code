#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GLYATL1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GLYATL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GLYATL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GLYATL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GLYCAM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GLYCTK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GLYR1 
