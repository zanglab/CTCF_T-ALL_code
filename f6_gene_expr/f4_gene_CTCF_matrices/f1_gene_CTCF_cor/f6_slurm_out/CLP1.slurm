#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CLP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CLP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLPB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLPP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLPS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLPSL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLPSL2 
