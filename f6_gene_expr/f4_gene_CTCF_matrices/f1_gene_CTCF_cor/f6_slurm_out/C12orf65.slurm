#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C12orf65.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C12orf65 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C12orf66 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C12orf71 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C12orf73 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C12orf74 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C12orf75 
