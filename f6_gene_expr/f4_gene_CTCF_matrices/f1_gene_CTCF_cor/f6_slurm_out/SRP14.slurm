#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SRP14.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SRP14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SRP19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SRP54 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SRP68 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SRP72 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SRP9 
