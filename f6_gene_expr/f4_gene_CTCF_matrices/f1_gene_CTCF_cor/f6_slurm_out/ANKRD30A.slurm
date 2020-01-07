#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ANKRD30A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD30A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD30B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD30BL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD30BP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD30BP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANKRD31 
