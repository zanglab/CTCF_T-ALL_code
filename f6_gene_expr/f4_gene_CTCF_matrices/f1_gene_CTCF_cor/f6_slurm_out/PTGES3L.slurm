#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PTGES3L.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PTGES3L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTGES3L-AARSD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTGFR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTGFRN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTGIR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTGIS 
