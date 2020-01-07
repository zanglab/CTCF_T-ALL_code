#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KIF18B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KIF18B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIF19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIF1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIF1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIF1C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIF20A 
