#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CRTAP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CRTAP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRTC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRTC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRTC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRY1 
