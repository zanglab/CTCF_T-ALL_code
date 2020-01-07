#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o BAP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i BAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BARD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BARHL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BARHL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BARX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BARX2 
