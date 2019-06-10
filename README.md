# LueVari
LueVari is a reference free metagenome SNP caller, based on the read-colored de Bruijn graph. LuVari is an extension of VARI (https://github.com/cosmo-team/cosmo/tree/VARI), in which the reads are stored as colors to allow the read-coherent traversal of de Bruijn graph. 
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
.fastq file (read set) 
# Output
.fasta file of gene-sized sequences spanning the SNPs, with specifying the SNP index and varying nucleotides.  
## Running LueVari
```
#count the k-mers (-k should be same in KMC, bubble_matrix and bubbles. Note that this version of BOSS does not support k values > 32. This issue will be fixed).
mkdir -p kmc_temp
ls -1 --color=no *.fastq |xargs -l -i echo "~/kmc -b -fq -k32 -ci0 -cs250 {} {}.kmc kmc_temp" >kmercount.sh
source kmercount.sh
ls -1 --color=no *.fastq |xargs -l -i echo "~/kmc_tools sort {}.kmc {}.kmc.sorted " >kmercountsort.sh
source kmercountsort.sh
ls -1 --color=no *.fastq |xargs -l -i echo "{}.kmc.sorted" > kmers_list

#build the succinct de Bruijn graph
./cosmo-pack -k kmers_list 


#construct the read-colored matrix
./bubbles_matrix -n <number of reads> -k <k value> -i <input file>

#SNP calling
./bubbles -c <minimum coverage> -u <number of output sequences> -b <max length of output sequences> -n <number of reads> -k <k value>

```
## Authors
Bahar Alipanahi, Martin D. Muggli, Musa Jundi, Noelle Noyes, and Christina Boucher
