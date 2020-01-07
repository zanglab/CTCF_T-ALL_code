#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PPP1R14B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP1R14B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP1R14C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP1R14D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP1R15A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP1R15B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPP1R16A 
