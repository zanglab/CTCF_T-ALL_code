#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HERPUD2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HERPUD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HES1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HES2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HES3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HES4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HES5 
