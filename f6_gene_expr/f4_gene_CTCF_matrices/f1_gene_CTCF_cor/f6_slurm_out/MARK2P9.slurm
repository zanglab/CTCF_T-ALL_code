#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MARK2P9.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MARK2P9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MARK3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MARK4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MARS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MARS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MARVELD1 
