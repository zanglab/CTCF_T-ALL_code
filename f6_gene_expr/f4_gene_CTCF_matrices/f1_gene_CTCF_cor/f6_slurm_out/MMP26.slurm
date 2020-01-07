#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MMP26.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MMP26 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MMP27 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MMP28 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MMP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MMP7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MMP8 
