#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C8B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C8B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C8G 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C8orf33 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C8orf34 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C8orf37 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C8orf44 
