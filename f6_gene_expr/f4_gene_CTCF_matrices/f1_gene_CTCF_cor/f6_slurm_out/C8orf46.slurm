#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C8orf46.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C8orf46 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C8orf48 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C8orf58 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C8orf59 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C8orf74 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C8orf76 
