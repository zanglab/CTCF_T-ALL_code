#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PPP1R2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP1R2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP1R21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP1R26 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP1R27 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP1R32 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP1R35 
