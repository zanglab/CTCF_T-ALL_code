#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PPP1R9A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP1R9A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP1R9B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP2CA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP2CB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP2R1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP2R1B 
