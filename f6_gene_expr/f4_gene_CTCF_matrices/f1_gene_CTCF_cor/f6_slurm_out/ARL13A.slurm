#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ARL13A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ARL13A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARL13B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARL14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARL14EP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARL14EPL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARL15 
