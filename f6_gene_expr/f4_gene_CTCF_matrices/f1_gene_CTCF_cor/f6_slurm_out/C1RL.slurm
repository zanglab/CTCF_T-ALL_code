#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C1RL.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C1RL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1S 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1orf100 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1orf105 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1orf109 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1orf112 
