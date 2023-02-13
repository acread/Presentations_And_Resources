## Today's Agenda
1.  


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

Copy the folder that corresponds to your group to your scratch.global folder - each group \
has a total of 11 fast5 raw read files output by a MinIon run \

````
cp -r /scratch.global/read0094/Nanopore_Class/Group'X' "your_scratch_folder"
````

### What is the range of file sizes for your group's .fast5s - why are they different?

## Basecalling
We will be using **guppy** on MSI for basecalling \
The version of guppy that is available is not the most recent - I'm working on getting this updated

Here are the available configs based on flow-cell + kit used to generate the data \
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
#SBATCH --time=1:00:00
#SBATCH --ntasks=5
#SBATCH --mem=40g
#SBATCH --tmp=32g
#SBATCH --job-name=Basecall_class

mkdir basecalled

module load guppy

/panfs/roc/msisoft/guppy/4.2.2-gpu/bin/guppy_basecaller -i "fast5folder"  -r -s basecalled \
  --config dna_r9.4.1_450bps_hac_prom.cfg  --device CUDA:0
`````

As the .fastqs begin to populate, take a few minutes to look at the files \
## what is the anatomy of a fastq file?
Some good stuff here: https://en.wikipedia.org/wiki/FASTQ_format#FAST5_and_HDF5_evolutions


## What do the read lengths look like?
save this little script with the .sh extension \
This will output the length of each read in the fastqs in the current directory\
**Note need to run 'chmod +x "script.sh" prior to running "./script.sh"**
`````
for i in *.fastq; do
cat $i | awk '{if(NR%4==2) print length($1)}' >  ${i}.readslength.txt
done
`````

## Combine the 11 .fastqs into a single file BUT please keep the 11 individual files as well (we'll use these later)
`````
cat *.fastq > combined.fastq
`````


## Mapping the reads to the Setaria reference genome
I have a copy of the Setaria viridis reference genome saved in my scratch folder, you can point to it 

`````
#!/bin/bash -l
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --mem=32g
#SBATCH --tmp=10g
#SBATCH --job-name=minimap


module load minimap2

minimap2 /scratch.global/read0094/Nanopore_Class/Sviridis_500_v2.0.fa  *.fastq -ax map-ont > "your_output".sam

samtools view -bS "your_output".sam > "your_output".bam
samtools sort <input.bam> -o <output.bam>
samtools rmdup <input.bam> <output.bam>
`````

