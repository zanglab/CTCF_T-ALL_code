#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GDE1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GDE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GDF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GDF10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GDF11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GDF15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GDF2 
