#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o B3GALNT2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i B3GALNT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i B3GALT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i B3GALT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i B3GALT4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i B3GALT5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i B3GALT6 
