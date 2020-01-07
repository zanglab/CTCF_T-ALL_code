#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PPP2R3C.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP2R3C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP2R5A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP2R5B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP2R5C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP2R5D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP2R5E 
