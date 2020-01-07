#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ARHGEF28.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGEF28 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGEF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGEF33 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGEF34P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGEF35 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGEF37 
