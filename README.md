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
Color matrix constructed with orders of magnitude less number of colors and compressed using Elias-Fano coding. And the file labels.txt which is a list of colors with same labels. 
## Running 
```
#count the k-mers
mkdir -p kmc_temp
./kmc -b -ci0 -fq -k32 -cs250 <fastq file> <output.kmc> kmc_temp
./kmc_tools sort <output.kmc> <output.kmc.sorted>

#construct the de Bruijn graph
./cosmo-pack -k <output.kmc.sorted>

#construct the read-colored matrix with reduced number of colors:
./reduce_color <fastq file> <number of reads>

#example: (please check the paths before running)
cd test
bash run.sh



```
## Paper
https://link.springer.com/chapter/10.1007/978-3-030-00479-8_1

## Authors
Bahar Alipanahi, Alan Kuhnle and Christina Boucher

