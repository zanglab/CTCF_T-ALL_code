#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ARHGAP11B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGAP11B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGAP12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGAP15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGAP17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGAP18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ARHGAP19 
