#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PEX1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PEX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PEX10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PEX11A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PEX11B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PEX11G 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PEX12 
