#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CDC14B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CDC14B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDC14C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDC16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDC20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDC20B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDC23 
