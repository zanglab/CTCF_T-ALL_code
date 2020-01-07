#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PPP5C.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP5C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP5D1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP6C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP6R1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP6R2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP6R3 
