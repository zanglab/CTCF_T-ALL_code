#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C5orf30.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C5orf30 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C5orf34 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C5orf38 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C5orf42 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C5orf46 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C5orf47 
