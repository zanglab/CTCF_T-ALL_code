#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IGFL1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IGFL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IGFL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IGFL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IGFL4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IGFLR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IGFN1 
