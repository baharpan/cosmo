# LueVari
LueVari is a reference free resistom SNP caller, based on the read-colored de Bruijn graph. LuVari is an extension of VARI (https://github.com/cosmo-team/cosmo/tree/VARI), in which the reads are stored as colors to allow the read-coherent traversal of de Bruijn graph. 
### Prerequisites

VARI and its all Prerequisites



### Installing
```
#Fetch software and setup directories
git clone https://github.com/cosmo-team/cosmo
cd cosmo
git checkout LueVari
Follow all steps listed in Building notes of VARI (https://github.com/cosmo-team/cosmo/tree/VARI#building-notes)
```
# Input
.fastq files 
# Output
.fasta file of gene-sized sequences spanning the SNPs, with specifying the SNP index and varying nucleotides.  
## Running LueVari
```
#count the k-mers
mkdir -p kmc_temp
ls -1 --color=no *.fasta |xargs -l -i echo "~/kmc -b -fq -k32 -ci0 -cs250 {} {}.kmc kmc_temp" >kmercount.sh
source kmercount.sh
ls -1 --color=no *.fasta |xargs -l -i echo "~/kmc_tools sort {}.kmc {}.kmc.sorted " >kmercountsort.sh
source kmercountsort.sh
ls -1 --color=no *.fasta |xargs -l -i echo "{}.kmc.sorted" > filtered_kmc2_list

#build the succinct de Bruijn graph
./cosmo-pack -k filtered_kmc2_list 


#construct the read-colored matrix
./bubbles_matrix -n <number of reads> -k <k value> -i <input file>

#SNP calling
./bubbles -c <minimum coverage> -u <number of output sequences> -b <max length of output sequences> -n <number of reads> -k <k value>

```
## Authors
Bahar Alipanahi, Martin D. Muggli, Musa Jundi, Noelle Noyes, and Christina Boucher
