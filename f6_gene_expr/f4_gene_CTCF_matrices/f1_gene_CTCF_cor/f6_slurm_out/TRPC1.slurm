#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TRPC1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TRPC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRPC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRPC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRPC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRPC4AP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRPC5 
