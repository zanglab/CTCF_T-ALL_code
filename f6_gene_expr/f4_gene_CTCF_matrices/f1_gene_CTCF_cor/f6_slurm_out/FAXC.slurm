#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FAXC.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FAXC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAXDC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBLIM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FBLL1 
