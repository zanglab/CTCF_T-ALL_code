#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MORN3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MORN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MORN4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MORN5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MOS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MOSPD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MOSPD2 
