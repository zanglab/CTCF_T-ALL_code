#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C21orf33.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C21orf33 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C21orf58 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C21orf59 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C21orf62 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C21orf91 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C22orf15 
