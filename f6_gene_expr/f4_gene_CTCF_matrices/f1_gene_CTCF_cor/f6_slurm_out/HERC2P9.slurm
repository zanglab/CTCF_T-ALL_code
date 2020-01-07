#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HERC2P9.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HERC2P9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HERC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HERC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HERC5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HERC6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HERPUD1 
