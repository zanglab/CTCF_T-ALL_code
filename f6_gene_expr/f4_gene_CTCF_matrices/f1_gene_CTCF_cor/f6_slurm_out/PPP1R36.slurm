#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PPP1R36.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP1R36 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP1R37 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP1R3A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP1R3B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP1R3C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP1R3D 
