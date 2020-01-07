#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GYG2P1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GYG2P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GYPA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GYPB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GYPC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GYPE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GYS1 
