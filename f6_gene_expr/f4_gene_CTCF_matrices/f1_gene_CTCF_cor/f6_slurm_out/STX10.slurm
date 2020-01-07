#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o STX10.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i STX10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STX11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STX12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STX16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STX16-NPEPL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i STX17 
