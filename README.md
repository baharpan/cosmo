# LueVari
LueVari is a reference free resistom SNP caller, based on the read-colored de Bruijn graph. LuVari is an extension of VARI (https://github.com/cosmo-team/cosmo/tree/VARI), in which the reads are stored as colors to allow the read-coherent traversal of de Bruijn graph. 
### Prerequisites

VARI and its all Prerequisites

Bowtie (https://github.com/BenLangmead/bowtie)

### Installing
```
#Fetch software and setup directories
git clone https://github.com/cosmo-team/cosmo
cd cosmo
git checkout LueVari
Follow all steps listed in Building notes of VARI (https://github.com/cosmo-team/cosmo/tree/VARI#building-notes)
```
# Input
Either .fasta or .fastq files.
# Output
.fasta file of gene-sized sequences spanning the SNPs, with specifying the SNP index and varying Nucleotides.  
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
cosmo-pack -k filtered_kmc2_list 

#align the k-mers to reads
bowtie-build  <reads.fasta> flat_bowtie
bowtie -f  -a -v0 flat_bowtie <kmers.fasta> outputFile

#prepare the kmers and the reads to fill the read-colored matrix
matrix_prepare.py outputFile kmers.txt

#construct the read-colored matrix
fast_matrix

#SNP calling
bubbles

```
## Authors
Bahar Alipanahi, Martin D. Muggli, Musa Jundi, Noelle Noyes, and Christina Boucher
