# Evaluating Allele-Specific Expression based on global transcription using ASElux
## Alfredo Lebron - Civelek Lab - August 10th, 2019
### Overall Goal
Allele-specific expression (ASE) is a genetic phenomenon in which one allele in
a heterozygous site is expressed more than another. ASE sites across the genome
indicate the biased expression of a transcript and are valuable clues of
genetic regulation. For example, ASE in enhancer regions imply a difference in
enhancer activity depending on the allele. This document explains how to
evaluate ASE from PROseq data using ASElux. To find genome-wide ASE sites, this
analysis requires genome coordinates based on global transcription, as defined
in the GTF format.

### Folder Organization
The image below displays the folder structure for the ASE analysis,
along with all of the files that you should have when done:
file:///home/ace08/Dropbox/SRIP2019/Progress_Reports/folder_pic.png

### Installing ASElux
Create a folder that will contain the ASElux file and analysis files. I named it
"newASE" since I was previously running an older version of ASElux. From here,
enter the following commands to get the ASElux directory:
```
#get from github repository
cd ~/newASE
module load gcc/7.1.0
module load openmpi/2.1.5
git clone https://github.com/abl0719/ASElux.git
cd ASElux
make -f Makefile
```
These commands will run ASE on demo files. Note that the path specifying the
directory of the executable will be different depending on the system and folder
locations.
```
#test executable on demo
ASElux=/nv/vol192/civeleklab/alfredo/newASE/ASElux/bin/ASElux
export PATH=/nv/vol192/civeleklab/alfredo/newASE/ASElux/bin:$PATH
ASElux build --gtf annotation.gtf --ref genome.fa --out demo
ASElux align --fa --pe --readLen 50 --index demo --vcf snps.vcf --seqFiles read1.fa read2.fa --out demo
```
### Setting up necessary files
Create a directory named "ASE_Analysis". In it, run these commands to create the
other necessary directories:
```
mkdir reference
mkdir run_ASE
mkdir vcf
```
In the "vcf" directory, copy over the VCF files from Hannah Zeidler's
"ASE Analysis" (2019). In the reference directory, copy over the hg38.fa file
referenced in Hannah's documentation, as well as the fixed_best_HMM.gtf file
from "Creating a custom GTF file using groHMM for ASE analysis" (Lebron, 2019).

### Creating indices to run ASE by chromosome
In order to assist in troubleshooting and speed up processing time, ASElux can
be run in parallel by chromosome. This requires indices to be built for every
chromosome by splitting up the GTF file. In this case, the GTF coordinates
were generated with groHMM based on global transcription of the smooth muscle
cell data I worked on. More information on how to get these coordinates can be
found in the previously mentioned documentation. To run this in parallel, a
SLURM script is used.

First, create a "chromlist.txt" file with the following content in the reference
directory:
```
chr1
chr2
chr3
chr4
chr5
chr6
chr7
chr8
chr9
chr10
chr11
chr12
chr13
chr14
chr15
chr16
chr17
chr18
chr19
chr20
chr21
chr22
chrX
chrY
chrM
```
Create arrays based on the chromosomes with these commands:
```
chr=$(cat chromlist.txt)
index=1
for ii in ${chr}
do
echo ${ii} > array${index}
index=$(expr ${index} + 1)
done
```
Next, create a SLURM script. I named it "org_main_20190810.sh". Note that it initializes a script
that is mentioned later. Also, not that the "dir" variable will be different
depending on the system used.
```
#!/bin/bash
#SBATCH -A CivelekLab
#SBATCH --ntasks=1
#SBATCH --time=72:00:00
#SBATCH --partition=standard
#SBATCH --mem 256000

dir=/nv/vol192/civeleklab/alfredo/newASE/ASE_Analysis
refdir=$dir/reference
cd ${refdir}
module load gcc/7.1.0
module load openmpi/2.1.5

cnt=${SLURM_ARRAY_TASK_ID}
${refdir}/index_build_20190810.sh ${cnt}
```
Next create a script that will split the "fixed_best_HMM.gtf" file by chromosome
and create indexes for ASElux to use. Again, note the difference in some
directories. I named this file "index_build_20190810.sh".
```
#!/bin/bash
dir=/nv/vol192/civeleklab/alfredo/newASE/ASE_Analysis
refdir=$dir/reference
vcfdir=$dir/vcf

cnt="$1"

cd ${refdir}

chrname=$(cat array${cnt})
gtf=${refdir}/fixed_best_HMM.gtf
cat ${gtf} | grep -w ${chrname} > ${chrname}.gtf

ASElux=/nv/vol192/civeleklab/alfredo/newASE/ASElux/bin/ASElux
export PATH=/nv/vol192/civeleklab/alfredo/newASE/ASElux/bin:$PATH

ASElux build --gtf ${refdir}/${chrname}.gtf --ref hg38.fa --out ${chrname}_index
```
To run these scripts, enter the commands below:
```
chmod u+x *.sh
sbatch --array=1-25 org_main_20190810.sh
```
The result of running this will be a reference directory with index files that look
like this for every chromosome number:
```
chr#.gtf
chr#_index.annotation
chr#_index_gene.sa
chr#_index_genome.sa
```
To evaluate ASE, create a folder describing a donor and experimental condition.
For example, I named one folder "Z99P"; Z99 is the donor and "P" stands for
PGDF treatment. In this folder, copy in the 25 array files from the reference.
You could also create them with the code that uses the "chromlist.txt" file, but the
end result will be the same: 25 files named array# with the corresponding
chromosome.

Similar to what was done before, create a SLURM script with this content (I
  named it "org_main_Z99P_20190813.sh"):
```
#!/bin/bash
#SBATCH -A CivelekLab
#SBATCH --ntasks=1
#SBATCH --time=72:00:00
#SBATCH --partition=standard
#SBATCH --mem 256000

dir=/nv/vol192/civeleklab/alfredo/newASE/ASE_Analysis
rundir=$dir/run_ASE/Z99P
cd ${rundir}
module load gcc/7.1.0
module load openmpi/2.1.5

cnt=${SLURM_ARRAY_TASK_ID}
${rundir}/Z99P_20190813.sh ${cnt}
```
The script above references "Z99P_20190813.sh", which contains this:
```
dir=/nv/vol192/civeleklab/alfredo/newASE/ASE_Analysis
rundir=$dir/run_ASE/Z99P
refdir=$dir/reference
fastqdir=/nv/vol192/civeleklab/alfredo/ASE1hr
vcfdir=$dir/vcf

cd $rundir

cnt="$1"

donor=Z99
chrname=$(cat $refdir/array${cnt})

index=$refdir/${chrname}_index
readLen=31
fastq=$fastqdir/${donor}_t1h_PDGF_cat.fastq
vcf=$vcfdir/${donor}_chr.vcf
output=${donor}_PDGF_${chrname}_results.txt
export PATH=/nv/vol192/civeleklab/alfredo/newASE/ASElux/bin:$PATH
ASElux=/nv/vol192/civeleklab/alfredo/newASE/ASElux/bin/ASElux
$ASElux align --fq --se --readLen $readLen --index $index --vcf $vcf --seqFiles $fastq --out $output
```
Run this by entering this command:
```
chmod u+x *.sh
sbatch --array=1-25 org_main_Z99P_20190813.sh
```
The end result will be text files containing the ASE results by chromosome.
Here is an example of a few of the files:
```
head *results*

Z99_PDGF_chr10_results.txt <==
10:120674475	1	1	0.5
10:131310207	2	0	1
10:120678856	2	2	0.5
10:122541384	0	41	0
10:125492105	5	1	0.833333
10:121222466	0	2	0
10:117283987	4	1	0.8
10:104778812	9	0	1
10:132301149	3	0	1
10:103163607	4	0	1

==> Z99_PDGF_chr11_results.txt <==
11:129723725	22	0	1
11:134019675	3	0	1
11:106058586	0	13	0
11:119934081	3	25	0.107143
11:128412297	1	1	0.5
11:109981606	0	2	0
11:130462511	10	16	0.384615
11:121047241	2	0	1
11:119938625	1	4	0.2
11:118937014	0	1	0

==> Z99_PDGF_chr12_results.txt <==
12:104354728	2	0	1
12:107331989	0	1	0
12:123058706	4	0	1
12:122479308	0	1	0
12:122389765	0	1	0
12:122395777	0	1	0
12:122396395	0	1	0
12:122397476	0	1	0
12:128357168	36	32	0.529412
12:128875948	0	6	0

##more results from other chromosomes will be displayed
##some chromosomes may not have results (X,Y,M,etc.)
```
To run ASElux on other donors and conditions, make the corresponding directory
and create an org_main SLURM script that references the script that runs ASElux in
that directory. Make sure to change the FASTQ name, donor name, directories, and
the corresponding VCF if necessary.
