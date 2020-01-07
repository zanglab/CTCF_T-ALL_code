#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RCVRN.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RCVRN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RD3L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RDH10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RDH11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RDH12 
