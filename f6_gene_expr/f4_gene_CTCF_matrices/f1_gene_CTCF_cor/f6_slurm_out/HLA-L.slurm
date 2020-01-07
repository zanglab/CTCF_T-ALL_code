#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HLA-L.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HLA-L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HLCS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HLF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HLTF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HLX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HM13 
