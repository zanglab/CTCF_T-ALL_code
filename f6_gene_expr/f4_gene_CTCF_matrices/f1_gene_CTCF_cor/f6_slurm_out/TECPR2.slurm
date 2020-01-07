#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TECPR2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TECPR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TECR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TECRL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TECTA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TECTB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TEDDM1 
