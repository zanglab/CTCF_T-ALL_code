#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CRYBB1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CRYBB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRYBB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRYBB2P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRYBB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRYBG3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CRYGB 
