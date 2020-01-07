#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CRYGC.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CRYGC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRYGD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRYGN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRYGS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRYL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRYM 
