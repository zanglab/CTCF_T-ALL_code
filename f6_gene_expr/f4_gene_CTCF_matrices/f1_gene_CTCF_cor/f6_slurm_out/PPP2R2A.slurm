#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PPP2R2A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP2R2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP2R2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP2R2C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP2R2D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP2R3A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP2R3B 
