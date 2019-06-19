#!/bin/sh
#SBATCH -p RG-Cristobal-Uauy # partition (queue)
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=5000
#SBATCH -J Circos
#SBATCH -o Circos.%N.%j.out
#SBATCH -e Circos.%N.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sophie.harrington@jic.ac.uk

source perl-5.22.1
source circos-0.69.3

work='/nbi/Research-Groups/NBI/Cristobal-Uauy/Sophie'

circos -conf $work/Scripts/Circos/circos_K2282.conf