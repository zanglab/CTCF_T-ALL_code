#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o AUNIP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i AUNIP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AUP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AURKA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AURKAIP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AURKAPS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AURKB 
