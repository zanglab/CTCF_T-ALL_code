#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PPP1R12A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP1R12A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP1R12B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP1R12C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP1R13B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP1R13L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP1R14A 
