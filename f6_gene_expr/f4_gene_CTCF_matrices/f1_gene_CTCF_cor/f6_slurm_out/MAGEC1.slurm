#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MAGEC1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MAGEC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAGEC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAGEC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAGED1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAGED2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAGED4 
