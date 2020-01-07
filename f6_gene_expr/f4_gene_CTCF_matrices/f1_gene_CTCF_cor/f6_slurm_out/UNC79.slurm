#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o UNC79.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i UNC79 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UNC80 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UNC93A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UNC93B1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UNCX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UNG 
