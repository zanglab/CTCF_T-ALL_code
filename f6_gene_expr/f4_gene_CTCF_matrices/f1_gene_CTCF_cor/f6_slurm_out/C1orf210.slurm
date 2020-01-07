#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C1orf210.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C1orf210 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1orf216 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1orf226 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1orf228 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1orf35 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C1orf43 
