#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GDF3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GDF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GDF5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GDF6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GDF7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GDF9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GDI1 
