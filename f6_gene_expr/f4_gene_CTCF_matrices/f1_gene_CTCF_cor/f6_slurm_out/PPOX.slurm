#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PPOX.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PPOX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP1CA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP1CB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP1CC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP1R10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP1R11 
