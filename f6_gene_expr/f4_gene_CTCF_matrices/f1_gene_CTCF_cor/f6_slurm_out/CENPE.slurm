#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CENPE.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CENPE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CENPF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CENPH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CENPI 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CENPJ 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CENPK 
