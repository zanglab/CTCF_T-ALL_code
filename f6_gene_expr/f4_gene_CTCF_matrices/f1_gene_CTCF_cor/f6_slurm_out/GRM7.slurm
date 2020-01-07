#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GRM7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GRM7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRM8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRPEL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRPEL2 
