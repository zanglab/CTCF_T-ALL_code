#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ARHGEF38.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGEF38 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGEF39 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGEF4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGEF40 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGEF5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGEF6 
