#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FAM83E.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM83E 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM83F 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM83G 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM83H 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM84A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM84B 
