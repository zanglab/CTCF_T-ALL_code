#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GH1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GHDC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GHITM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GHR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GHRH 
