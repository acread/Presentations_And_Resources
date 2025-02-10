## Nanopore long read sequencing lab

In this lab we will dive deeper into how Nanopore sequencing works and the DNA requirements to \
maximize data generation, learn some of the \
newest technologies that have come to market, and do a hands-on base-calling and read-mapping \
workflow

This markdown file will be accompanied by a powerpoint presentation

## Let's do some computational stuff!

Log onto your MSI account and navigate to the Nanopore folder in the agro5431/shared directory \
**NOTE: you need to sign onto mangi!**
`````
ssh 'username'@mangi.msi.umn.edu
`````

Each group has a total of 3 fast5 raw read files output by a MinIon run \
Each fast5 contains 4,000 raw read 'squiggles'

## Basecalling
We will be using **Dorado** on MSI for basecalling \


Create a .sh file to submit to MSI - note that this script makes use of the GPU on the v100 partition
You need to update everything that is in quotes
`````
#!/bin/bash -l
#SBATCH -p v100                                             
#SBATCH --gres=gpu:v100:1
#SBATCH --time=0:20:00
#SBATCH --ntasks=5
#SBATCH --mem=80g
#SBATCH --tmp=64g
#SBATCH --job-name=Basecall_class
#SBATCH --mail-type=ALL
#SBATCH --mail-user='x500'@umn.edu

mkdir /home/agro5431/shared/Nanopore/'x500'/'x500'_guppy_basecalled

module load dorado/0.7.1

#usage - call dorado, point to a basecall model, point to a folder of raw reads, tell it we want fastq output (default is bam), direct it to our new output folder
#note that by default dorado will use the GPU - it can be run using CPU by adding --device cpu

dorado basecaller /home/agro5431/shared/LongReadLab/dna_r9.4.1_e8_hac@v3.3 \
 /home/agro5431/shared/LongReadLab/Table1_fast5s \
 --emit-fastq > 'x500'_guppy_basecalled/calls.fastq

`````

### How long do you think it would take to basecall using a CPU instead of a GPU?

As the .fastqs begin to populate, take a few minutes to look at the files 
### what is the anatomy of a fastq file?
Some good stuff here: https://en.wikipedia.org/wiki/FASTQ_format#FAST5_and_HDF5_evolutions

### What do the read lengths look like?
save this little script with the .sh extension \
This will output the length of each read in the fastqs in the current directory\
**Note need to run 'chmod +x "script.sh" prior to running "./script.sh"**
`````
for i in *.fastq; do
cat $i | awk '{if(NR%4==2) print length($1)}' >  ${i}.readslength.txt
done
`````
### What is the sum of the sequenced bases?


## Mapping the reads to the Setaria reference genome
We will use Minimap2 for this - https://lh3.github.io/minimap2/minimap2.html \
I saved a copy of the Setaria reference genome in the /home/agro5431/shared/LongReadLab folder - the below script points to it

`````
#!/bin/bash -l
#SBATCH --time=0:15:00
#SBATCH --ntasks=1
#SBATCH --mem=32g
#SBATCH --tmp=10g
#SBATCH --job-name=minimap

#Andy adding samtools module load
module load samtools/1.14
module load minimap2/2.17

#usage - point to a reference sequence, give basecalled reads to map, tell it we're using Nanopore (ont) output as a sam
minimap2 /home/agro5431/shared/LongReadLab/Sviridis_500_v2.0.fa  *.fastq -ax map-ont > 'x500'.sam

`````

## Did all the reads map to the provided genome?
We can output the unmapped reads using samtools - https://www.htslib.org/doc/ \
we can run this on the command line pointing to our 'x500'.sam from the above step
`````
module load samtools
samtools view 'x500'.sam -f 4 > 'x500'_unmapped.txt
`````
### Let's BLAST a few of the unmapped reads
https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&BLAST_SPEC=&LINK_LOC=blasttab&LAST_PAGE=blastn \

You can view your unmapped.txt file however you would like (head/tail/less) \
Copy one of the DNA sequences and paste it into the NCBI browser 
### What does the sequence hit?  Try with 2-3 additional sequences.
