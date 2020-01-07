#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GCNT4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GCNT4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GCNT7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GCOM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GCSAM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GCSAML 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GCSAML-AS1 
