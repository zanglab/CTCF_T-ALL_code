#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PPP1R16B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP1R16B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP1R17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP1R18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP1R1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP1R1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP1R1C 
