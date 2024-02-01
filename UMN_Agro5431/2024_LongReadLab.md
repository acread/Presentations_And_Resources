## Nanopore long read sequencing lab

In this lab we will dive deeper into how Nanopore sequencing works and the DNA requirements to \
maximize data generation, demo a flow-cell check using the Nanopore software, learn some of the \
newest technologies that have come to market, and do a hands-on base-calling and read-mapping \
workflow

This markdown file will be accompanied by a powerpoint presentation

## Let's do some computational stuff!

Log onto your MSI account and navigate to the Nanopore folder in the agro5431/shared directory \
**NOTE: you need to sign onto mangi!**
`````
ssh 'username'@mangi.msi.umn.edu
`````

I did some Nanopore sequencing of the model monocot Setaria viridis

Each table will work with the same set of 11 fast5 raw read files output by a MinIon run, but every person should do each step\
Each fast5 contains 4,000 raw read 'squiggles'

Make a new directory called "'x500'_fast5s" in your folder and copy the fast5s from /home/agro5431/shared/LongReadLab/'yourtablenumber'_fast5s
`````
mkdir 'x500_fast5s'
cd 'x500_fast5s'
cp /home/agro5431/shared/LongReadLab/'yourtablenumber'_fast5s/*.fast5 .
`````

## Basecalling
We will be using **guppy** on MSI for basecalling \

Create a .sh file to submit to MSI - note that this script makes use of the GPU on the v100 partition
You need to update everything that is in quotes\
**NOTE: I use VIM to create new files, you can use whatever you're comfortable with**

`````
vim 'x500'_basecall.sh
`````
copy and paste in the script \
Be sure to change the parts in quotes to match your files
`````
#!/bin/bash -l
#SBATCH -p v100                                             
#SBATCH --gres=gpu:v100:1
#SBATCH --time=1:00:00
#SBATCH --ntasks=5
#SBATCH --mem=40g
#SBATCH --tmp=32g
#SBATCH --job-name=Basecall_class

mkdir "x500"_basecalled

module load guppy/3.2.4

#the command to call the Nanopore basecaller
/common/software/install/migrated/guppy/4.2.2-gpu/bin/guppy_basecaller -i "x500"_fast5s  -r -s "x500"_basecalled \
  --config dna_r9.4.1_450bps_hac_prom.cfg  --device CUDA:0
`````

Run the script -- it make take awhile becasue so many people are submitting at once,\
When I tested it took about 10 minutes
`````
sbatch 'x500'_basecall.sh
`````

you can check the status of your job in a few ways
`````
squeue -u 'x500'

less slurm-'xxxxxxxxxx'.out

#check the contents of your 'x500'_basecalled folder
`````

### How long do you think it would take to basecall using a CPU instead of a GPU?

As the .fastqs begin to populate, take a few minutes to look at the files 
### what is the anatomy of a fastq file?
Some good stuff here: https://en.wikipedia.org/wiki/FASTQ_format#FAST5_and_HDF5_evolutions

### What do the read lengths look like?
copy and paste this little script as 'x500'_lens.sh in your basecalled folder \
This will output the length of each read in the fastqs in the current directory\
**Note need to run 'chmod +x 'x500'_lens.sh prior to running ./'x500'_lens.sh**
`````
for i in *.fastq; do
cat $i | awk '{if(NR%4==2) print length($1)}' >  ${i}.readslength.txt
done
`````
### What is the sum of the sequenced bases?

### Combine the 11 .fastqs into a single file
`````
cd /home/agro5431/'x500'/'x500'_basecalled
cat *.fastq > combined.fastq
`````


## Mapping the reads to the Setaria reference genome using minimap2
`````
vim /home/agro5431/'x500'/minimap_align.sh
`````
I saved a copy of the Setaria reference genome in the /home/agro5431/LongReadLab folder - the below script points to it

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

#arguments for minimap: ref genome, fastqs to align - the asterix means all files that end in .fastq, -ax tell it we're using nanopore data
minimap2 /home/agro5431/shared/LongReadLab/Sviridis_500_v2.0.fa  /home/agro5431/'x500'/'x500'_basecalled/combined.fastq -ax map-ont > "x500"_SetariaAlignment.sam

#I think we can get rid of this step for the exercise
#convert the sam to a bam and sort it
#samtools view -bS "x500"_SetariaAlignment.sam > "x500"_SetariaAlignment.bam
#samtools sort "x500"_SetariaAlignment.bam -o "x500"_SetariaAlignment.bam
`````

## Did all the reads map to the provided genome?
We can output the unmapped reads using samtools - we can run this right on the command line \
We will point to our "x500"_SetariaAlignment.sam from the above step
`````
module load samtools/1.14
samtools view "x500"_SetariaAlignment.sam -f 4 > "x500"_unmapped.txt
`````
### Let's BLAST a few of the unmapped reads
https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&BLAST_SPEC=&LINK_LOC=blasttab&LAST_PAGE=blastn \

You can view your unmapped.txt file however you would like (head/tail/less) \
Copy one of the DNA sequences and paste it into the NCBI browser \
Note: This is a little more challenging if you're using a PC - copy/paste doesn't play as nicely...
### What does the sequence hit?  Try with 2-3 additional sequences.
