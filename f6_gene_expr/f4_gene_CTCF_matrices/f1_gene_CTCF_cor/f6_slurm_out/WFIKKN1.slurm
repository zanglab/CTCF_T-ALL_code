#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o WFIKKN1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i WFIKKN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WFIKKN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WFS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WHAMM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WHAMMP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WHAMMP3 
