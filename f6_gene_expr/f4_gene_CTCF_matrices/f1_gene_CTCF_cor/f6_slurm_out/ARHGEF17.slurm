#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ARHGEF17.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGEF17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGEF18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGEF19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGEF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGEF25 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGEF26 
