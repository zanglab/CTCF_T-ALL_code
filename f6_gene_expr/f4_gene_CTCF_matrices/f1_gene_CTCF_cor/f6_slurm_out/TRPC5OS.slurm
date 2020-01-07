#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TRPC5OS.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TRPC5OS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRPC6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRPC7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRPM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRPM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TRPM3 
