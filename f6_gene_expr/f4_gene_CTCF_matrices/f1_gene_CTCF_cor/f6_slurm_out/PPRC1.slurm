#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PPRC1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PPRC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPT2-EGFL8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPTC7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPWD1 
