#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PDE6G.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PDE6G 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDE6H 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDE7A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDE7B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDE8A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDE8B 
