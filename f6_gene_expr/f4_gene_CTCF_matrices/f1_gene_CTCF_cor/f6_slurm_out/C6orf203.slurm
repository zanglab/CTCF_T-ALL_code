#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C6orf203.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C6orf203 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C6orf222 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C6orf223 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C6orf226 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C6orf229 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C6orf47 
