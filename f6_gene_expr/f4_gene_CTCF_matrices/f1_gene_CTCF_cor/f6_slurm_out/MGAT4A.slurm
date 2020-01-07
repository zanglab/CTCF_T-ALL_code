#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MGAT4A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MGAT4A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MGAT4B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MGAT4C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MGAT4D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MGAT4EP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MGAT5 
