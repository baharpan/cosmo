# Recoloring  the  Colored  de  Bruijn  Graph
 A heuristic algorithm to recolor a colored DBG, with orders of magnitude less number of colors. In this project colors are the sequencing reads.
 ### Prerequisites

VARI and its all Prerequisites

### Installing
```
#Fetch software and setup directories
git clone https://github.com/cosmo-team/cosmo
cd cosmo
git checkout Recoloring
Follow all steps listed in Building notes of VARI (https://github.com/cosmo-team/cosmo/tree/VARI#building-notes)
```
# Input
fastq files.
# Output
Color matrix constructed with orders of magnitude less number of colors and compressed using Elias-Fano coding 
## Running 
```
#count the k-mers
mkdir -p kmc_temp
ls -1 --color=no *.fastq |xargs -l -i echo "./kmc -b -fq -k32 -ci12 -cs250 {} {}.kmc kmc_temp" >kmercount.sh
source kmercount.sh
ls -1 --color=no *.fastq |xargs -l -i echo "./kmc_tools sort {}.kmc {}.kmc.sorted " >kmercountsort.sh
source kmercountsort.sh
ls -1 --color=no *.fastq |xargs -l -i echo "{}.kmc.sorted" > filtered_kmc2_list

#construct the de Bruijn graph
./cosmo-pack -k filtered_kmc2_list

#construct the read-colored matrix with reduced number of colors:
./reduce_color <fastq file> <number of reads>


```
## Paper
https://link.springer.com/chapter/10.1007/978-3-030-00479-8_1

## Authors
Bahar Alipanahi, Alan Kuhnle and Christina Boucher

