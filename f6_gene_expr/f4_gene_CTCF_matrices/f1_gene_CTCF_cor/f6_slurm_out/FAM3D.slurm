#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FAM3D.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM3D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM43A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM43B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM45A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM46A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM46B 
