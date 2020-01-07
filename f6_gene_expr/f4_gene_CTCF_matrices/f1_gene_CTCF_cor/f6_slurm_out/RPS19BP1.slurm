#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RPS19BP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS19BP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS23 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS24 
