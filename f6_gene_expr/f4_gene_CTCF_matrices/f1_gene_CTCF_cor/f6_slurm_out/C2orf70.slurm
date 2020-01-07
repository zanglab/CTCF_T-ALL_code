#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C2orf70.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C2orf70 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C2orf71 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C2orf72 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C2orf73 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C2orf74 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C2orf76 
