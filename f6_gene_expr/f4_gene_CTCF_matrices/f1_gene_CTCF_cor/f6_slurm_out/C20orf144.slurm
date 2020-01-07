#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C20orf144.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C20orf144 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C20orf173 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C20orf194 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C20orf196 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C20orf202 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C20orf203 
