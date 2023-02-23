

## It starts with purified nucleic acid (DNA or RNA)
Isolation of high molecular weight nucleic acid can be challenging - largely depends on tissue type \
Plants seem particularly challenging, probably because of polysaccharides and secondary metabolites \

### How can you tell if you have high molecular weight DNA?

<details>
  <summary>Agilent TapeStation/Bioanalyzer Output</summary>
<img width="276" alt="image" src="https://user-images.githubusercontent.com/43852873/218581783-36197ece-baff-4a4e-afe9-5a1cf123ca2d.png">
</details>
  
  
Nanodrop should be used for absorbance - targetting 260/230 2.0-2.2 and 260/280 ~2.0
### What do low 260/230 and 260/280 indicate?

Qubit or other flourescent method for quantification

### Choosing a library prep and flowcell
Go check out some of the options https://nanoporetech.com/products/kits

Flow-cell demo - qc and/or loading a library

### Beyond the basics, what is Nanopore capable of?


### What is the fast5 format and what data is stored here?

<details>
  <summary>Zooming in on fast5 squiggles</summary>
<img width="631" alt="image" src="https://user-images.githubusercontent.com/43852873/218881708-58308251-bd3a-461b-879b-42b98f65f459.png">
</details>

## Let's do some computational stuff!

Log onto your MSI account and navigate to the Nanopore folder in the agro5431/shared directory
`````
cd /home/agro5431/shared/Nanopore_Lab
`````

Copy the folder that corresponds to your group to your scratch.global folder - each group \
has a total of 11 fast5 raw read files output by a MinIon run \

````
cp -r /scratch.global/read0094/Nanopore_Class/Group'X' "your_scratch_folder"
````

### What is the range of file sizes for your group's .fast5s - why are they different?
`````
ls -lh
`````

## Basecalling
We will be using **guppy** on MSI for basecalling \
The version of guppy that is available is not the most recent - I'm working on getting this updated\
If you are interested in downloading guppy for your own basecalling - https://community.nanoporetech.com/downloads \
**Note: you need a Nanopore Community log-in to access the downloads page**

Here are the available configs based on flow-cell + kit used to generate the data 
`````
dna_r10.3_450bps_hac
dna_r10.3_450bps_hac_prom
dna_r10_450bps_hac
dna_r9.4.1_450bps_hac
dna_r9.4.1_450bps_hac_prom
dna_r9.5_450bps
rna_r9.4.1_70bps_hac
rna_r9.4.1_70bps_hac_prom
`````

Create a .sb file to submit to MSI - note that this script makes use of the GPU on the v100 partition
`````
#!/bin/bash -l
#SBATCH -p v100                                             
#SBATCH --gres=gpu:v100:1
#SBATCH --time=0:10:00
#SBATCH --ntasks=5
#SBATCH --mem=40g
#SBATCH --tmp=32g
#SBATCH --job-name=Basecall_class

mkdir "your_initials"basecalled

module load guppy/3.2.4

/panfs/roc/msisoft/guppy/4.2.2-gpu/bin/guppy_basecaller -i "fast5folder"  -r -s "your_initials"basecalled \
  --config dna_r9.4.1_450bps_hac_prom.cfg  --device CUDA:0
`````

### How long do you think it would take to basecall using a CPU instead of a GPU?

As the .fastqs begin to populate, take a few minutes to look at the files 
### what is the anatomy of a fastq file?
Some good stuff here: https://en.wikipedia.org/wiki/FASTQ_format#FAST5_and_HDF5_evolutions

### Rank these three from highest to lowest quality based on the fastq quality code:
1.  {|}{|}{|}{|}{|}{|}{|}{|}
2.  CDECDECDECDECDECDECDECDE
3.  cdecdecdecdecdecdecdecde
4.  #$%#$%#$%#$%#$%#$%#$%#$%

### What do the read lengths look like?
save this little script with the .sh extension \
This will output the length of each read in the fastqs in the current directory\
**Note need to run 'chmod +x "script.sh" prior to running "./script.sh"**
`````
for i in *.fastq; do
cat $i | awk '{if(NR%4==2) print length($1)}' >  ${i}.readslength.txt
done
`````
### What is the sum of the sequenced bases? (you can do this before or after merging the fastqs)

### Combine the 11 .fastqs into a single file BUT please keep the 11 individual files as well (we'll use these later)
`````
cat *.fastq > combined.fastq
`````


## Mapping the reads to the Setaria reference genome
I have a copy of the Setaria viridis reference genome saved in my scratch folder, the below script points to it

`````
#!/bin/bash -l
#SBATCH --time=0:15:00
#SBATCH --ntasks=1
#SBATCH --mem=32g
#SBATCH --tmp=10g
#SBATCH --job-name=minimap

#Andy added this samtools module load Feb 23
module load samtools/1.14
module load minimap2/2.17

minimap2 /scratch.global/read0094/Nanopore_Class/Sviridis_500_v2.0.fa  *.fastq -ax map-ont > "your_output".sam

samtools view -bS "your_output".sam > "your_output".bam
samtools sort <input.bam> -o <output.bam>
`````

Copy your <output.bam> file to your group's folder in my scratch
`````
cp <output.bam> /scratch.global/read0094/Nanopore_Class/Group"X"
`````

I'll visualize the alignments using IGV (Integrated Genomics Viewer) \
https://software.broadinstitute.org/software/igv/download

### Did all the reads map to the provided genome?
We can output the unmapped reads using samtools - we can run this right on the command line \
We will point to our "your_output".sam from the above step
`````
module load samtools
samtools view "your_output".sam -f 4 > unmapped.txt
`````
### Let's BLAST a few of the unmapped reads
https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&BLAST_SPEC=&LINK_LOC=blasttab&LAST_PAGE=blastn \

You can view your unmapped.txt file however you would like (head/tail/less) \
Copy one of the DNA sequences and paste it into the NCBI browser 
### What does the sequence hit?  Try with 2-3 additional sequences.

### Is there a higher throughput way to see what all of our unmapped reads align to?



## Homework Questions
1.  Our version of Guppy did not include the plant specific base-calling model. What is the \
    name of this config?
2.  The fast5 format is not efficient. Recently an alternative format has been proposed. \
      -What is the name of this file format? \
      -What is the name of the binary version of the new format
4.  How long is your longest read? \
      -DNA passes through the pore at 450 bases/second, how long did it take this read to traverse the pore?
5.  Perform an online nucleotide BLAST with one of your reads \
      -What organism does it match? \
      -What is the length and identity of the match? \
      -What can you infer about the base-calling accuracy? 
6.  You calculated the sum of the bases called - what coverage of the Setaria genome is this? \
      **Hint: it is possible that this is less than one**
7.  How many reads were able to map to the Setaria genome?
8.  We saw that some of our unmapped reads mapped to common bacteria -- name three potential sources of bacterial \
    DNA in our sequencing library.
9.  Each group got a subset of 11 fast5 files -- in reality, I generated 170 fast5s from my run. \
    If the read length for all fast5s is about equal to the mean of your subset, what depth coverage \
    of Setaria did I get from this single MinIon flowcell?
10.  Using R, generate boxplots that show the distribution of fastq lengths from each of your 11 fast5s \
    **Note: this means you will have 11 bars.** \
    Display a text label indicating the mean read length for each fastq \
    Use a MetBrewer color palette for your bars https://github.com/BlakeRMills/MetBrewer \
    **Note: there are a few palettes with 11 colors, if you choose one with fewer you will need to figure out** \
    **how to deal with this (either adding colors or creating a continuous palette)**


