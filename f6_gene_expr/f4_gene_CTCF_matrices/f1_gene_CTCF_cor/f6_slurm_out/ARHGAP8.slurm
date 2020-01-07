#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ARHGAP8.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGAP8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGAP9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGDIA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGDIB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGDIG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGEF1 
