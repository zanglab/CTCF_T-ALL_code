#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MMP20.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MMP20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MMP21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MMP23A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MMP23B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MMP24 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MMP25 
