#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TOP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TOP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TOP1MT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TOP1P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TOP2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TOP2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TOP3A 
