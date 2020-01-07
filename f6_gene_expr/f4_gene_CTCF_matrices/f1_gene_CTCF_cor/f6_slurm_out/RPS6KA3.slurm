#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RPS6KA3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS6KA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS6KA4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS6KA5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS6KA6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS6KB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPS6KB2 
