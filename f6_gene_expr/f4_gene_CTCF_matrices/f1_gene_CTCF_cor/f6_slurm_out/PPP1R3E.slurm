#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PPP1R3E.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP1R3E 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP1R3F 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP1R3G 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP1R42 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP1R7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP1R8 
