#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LAGE3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LAGE3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LAIR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LAIR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LAMA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LAMA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LAMA3 
