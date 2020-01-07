#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HSP90AB4P.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HSP90AB4P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSP90B1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSP90B2P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSPA12A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSPA12B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HSPA13 
