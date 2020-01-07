#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CD38.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CD38 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD3D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD3E 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD3EAP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD3G 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD4 
