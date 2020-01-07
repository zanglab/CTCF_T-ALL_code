#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C4orf36.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C4orf36 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C4orf45 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C4orf46 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C4orf47 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C4orf48 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C4orf51 
