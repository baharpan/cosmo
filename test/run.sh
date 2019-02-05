#!/bin/bash
#count the k-mers
mkdir -p kmc_temp
./kmc -b -ci12 -fq -k32 -cs250 new.fastq list.kmc kmc_temp
./kmc_tools sort list.kmc list.kmc.sorted

echo "list.kmc.sorted" > list

#construct the de Bruijn graph
../cosmo-pack -k list

#construct the read-colored matrix with reduced number of colors:
../reduce_color new.fastq 20000
