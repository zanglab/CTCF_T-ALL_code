#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PPP3CA.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP3CA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP3CB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP3CC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP3R1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP3R2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP4C 
