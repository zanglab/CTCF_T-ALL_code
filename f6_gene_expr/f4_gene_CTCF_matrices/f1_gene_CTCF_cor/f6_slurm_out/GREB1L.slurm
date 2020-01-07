#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GREB1L.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GREB1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GREM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GREM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRHL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRHL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRHL3 
