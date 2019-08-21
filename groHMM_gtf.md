# Creating a custom GTF file using groHMM for ASE analysis
## Alfredo Lebron - Civelek Lab - July 5th, 2019
### Overall Goal
ASElux uses a prepared index to detect allele specific expression. Generally, an index is prepared using an hg38 FASTA file along with an hg38 GTF file. For more information, refer to Hannah Zeidler's "ASE Analysis" vignette (April 23rd, 2019). Sometimes, there is an issue with restricting an ASE analysis to the chromosome coordinates set by an hg38 GTF annotation. Since ASE often occurs outside of genes, limiting an analysis to set gene coordinates will not fully is=llustrate the regulatory effects of ASE in a given condition. Thus, I present a method of creating a customized GTF file that covers global transcription regions using the groHMM package. The final output GTF can be used in the same way as other GTFs in an ASE analysis, and can be used to broadly detect ASE across the entire transcriptome. My code and workflow is adapted from Warren Anderson's and Michael Guertin's "Identification of transcriptional units" vignette (November 7th, 2018).

### Setting up groHMM
This part of the workflow requires R, which in my case was run on UVA's Rivanna cluster. Make sure to look up UVA's high-performance computing documentation for more information on how to log on to the Rivanna server and how to set up a personal lab directory. This process requires significant CPU and RAM usage, so requesting a job using this command before any part of this analysis is basically required:
```
ijob -A civeleklab -p standard --mem 128000 --time=10:00:00
```
I've found that Rivanna will sometimes take very long to grant this job request. If this happens, try lowering the memory request to less than 128 gigabytes. If that still doesn't work, it's probably best to just wait and try again later.

During this analysis, you will need these commands to run the code that I used. Just be aware that Rivanna is subject to updates and changes that will likely change the specific commands:
```
module load gcc/7.1.0
module load openmpi/2.1.5
module load R/3.6.0
module load bedtools/2.26.0
```
Type in "R" in the command line and R should load up. For this to work, you will need to load these packages:
```
library(dplyr)
library(bigWig)
library(groHMM)
library(GenomicFeatures)
```
For the most part, instructions for package installation are pretty easy to find online, but people tend to have an issue with installing bigWig, especially on Rivanna. Just use these commands if you have trouble:
```
install.packages('devtools", repos='[1]http://cran.cnr.berkeley.edu')

devtools::install_github("andrelmartins/bigWig/bigWig")
```
A BAM file is needed. This is obtained from experimental results. Make sure your working directory contains this file, and enter these commands into R ("expr" is the main R object you will want to keep from this block of code):
```
## get merged bam
bam.file="smc_merged.bam"
data <- readGAlignments(bam.file)
data_gr <- granges(data)
data_gr <- sort(data_gr)
expr <- keepStandardChromosomes(data_gr, pruning.mode="coarse")
```
The next block of code requires an hg38 BED annotation file called "hg38.wide.annotation.bed". I generated this file by first sourcing the hg38.gtf file described in Hannah Zeilder's previously mentioned documentation. Then, I used the following shell script to convert the GTF into a BED file. Note that this script includes a way to obtain an hg38 GTF file if Hannah's is not available:
```
###########################################################################
## get comprehensive annotations
## this will be used to estimate TTS
###########################################################################

# get gencode comprehensive annotations
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_21/gencode.v21.annotation.gtf.gz
gunzip gencode.v21.annotation.gtf
mv gencode.v21.annotation.gtf hg38.gtf

# get 'comprehensive' protein coding genes
grep 'gene_type "protein_coding"' hg38.gtf | cut -f1,4,5,7,9 | tr ";" "\t" > hg38_unfilt1.bed
grep 'gene_id' hg38_unfilt1.bed > hg38_unfilt2.bed
grep 'gene_name' hg38_unfilt2.bed > hg38_protcode.bed
awk '{for(i=5;i<=NF;i++){if($i~/^gene_name/){a=$(i+1)}} print $1,$2,$3,a,"na",$4}' hg38_protcode.bed | tr " " "\t" > hg38.annotation.bed
rm *unfilt*
```
This will generate a file with a head that looks like this:
```
chr1	69091	70008	"OR4F5"	na	+
chr1	69091	70008	"OR4F5"	na	+
chr1	69091	70008	"OR4F5"	na	+
chr1	69091	70005	"OR4F5"	na	+
chr1	69091	69093	"OR4F5"	na	+
chr1	70006	70008	"OR4F5"	na	+
chr1	70006	70008	"OR4F5"	na	+
chr1	182393	184158	"FO538757.3"	na	+
chr1	182393	184158	"FO538757.3"	na	+
chr1	182393	182746	"FO538757.3"	na	+
chr1	182709	182746	"FO538757.3"	na	+
chr1	182709	182711	"FO538757.3"	na	+
chr1	183114	183240	"FO538757.3"	na	+
chr1	183114	183240	"FO538757.3"	na	+
chr1	183922	184158	"FO538757.3"	na	+
chr1	183922	184155	"FO538757.3"	na	+
chr1	184156	184158	"FO538757.3"	na	+
chr1	182393	182708	"FO538757.3"	na	+
chr1	184156	184158	"FO538757.3"	na	+
chr1	184923	200322	"FO538757.2"	na	-
```
As you'll notice, there are repeating rows with different start and end sites. For this analysis, the eventual transcription units won't need to be exact; they just need to broadly cover all possible regions of transcription across the genome. So, groHMM will need to use a GTF with the widest possible coverage as a baseline to predict transcription units. This R code will read this GTF file and parse out the unique genes with the earliest start site and farthest end site:
```
#Set the working directory where the .bed file is located
setwd("~/working/directory/")

#read the .bed file as a data frame in R
hg38.annotation.bed <- read.table("hg38.annotation.bed",header=FALSE, stringsAsFactors=F)

#Only include unique rows (remove duplicates)
unique.hg38.annotation.bed <- unique(hg38.annotation.bed)

##create new data frame by taking the original and extracting the rows with the earliest start site and
##farthest end site for each gene

#naming
bed0 = unique.hg38.annotation.bed
names(bed0) = c("chr","start","end","gene","zz","strand")

#create an empty vector that will hold the results of the loop
widebed0 <- c()

#run loop
for (ii in unique(bed0$gene)){
  ind <- which(bed0$gene==ii)
  start <- min(bed0$start[ind])
  end <- max(bed0$end[ind])
  out <- c(bed0$chr[ind][1],start,end,ii,"na",bed0$strand[ind][1])
  widebed0 <- rbind(widebed0,out)
}

#convert output to data frame since we have different data types and they are being coerced into a matrix of characters
widebed0 <- as.data.frame(widebed0)

#give columns descriptive names
names(widebed0) <- c("chr","start","end","gene","zz","strand")

#remove rownames
rownames(widebed0)=NULL

#save as a .bed file
write.table(widebed0,file="hg38.wide.annotation.bed",quote=F,sep="\t",row.names=F)
```
### Running groHMM
After "hg38.wide.annotation.bed" is created, the groHMM package can run. Two main parameters are used to run the groHMM prediction algorithm. They are the UTS and LtProbB, or the "log-transformed transition probability of switching from transcribed state to non-transcribed state and the variance of the emission probability for reads in the non-transcribed state, respectively" (refer to groHMM documentation for more detailed information on what these probabilities mean and how the package works). After this code is run, multiple GRanges objects will be created based on combinations of various UTS and LtProbB parameters. Ultimately, you must select the GRanges object with the lowest density of reads in untranscribed regions and make a GTF. Furthermore, this GRanges object should have merge and disassociation errors in the first quartile. Here is the code for the analysis:
```
##get gene annotations based on TSS id

#loading the data from the .bed file
gene.ann0 = read.table("hg38.wide.annotation.bed",sep="\t",header=T,stringsAsFactors=F)

#create a GRanges object from the data
gene.ann = makeGRangesFromDataFrame(gene.ann0[,-5], seqnames.field="chr",
                                    start.field="start",end.field="end",
                                    strand.field="strand",
                                    starts.in.df.are.0based=TRUE,
                                    keep.extra.columns=TRUE)

#organize
gene.ann <- sort(gene.ann)
gene.ann <- keepStandardChromosomes(gene.ann, pruning.mode="coarse")
unique(seqnames(gene.ann))
unique(seqnames(expr))

# window parameters
Fp <- windowAnalysis(expr, strand="+", windowSize=50)
Fm <- windowAnalysis(expr, strand="-", windowSize=50)
# set the number of cores for grohmm
options(mc.cores=getCores(40))

# specify parameters for variation
vars = c(1,2,3,4,5,6,8,10,12,15,20,25,30,40,50,60,80,100,150,200,300)
ltpr = c(-10,-20,-50,-100,-150,-200,-300,-400,-500)
LtProbB = sapply(ltpr,function(x){rep(x,length(vars))}) %>% as.vector
UTS = rep(vars,length(ltpr))
tune <- data.frame(LtProbB=LtProbB, UTS=UTS)

#run groHMM
hmm.vars <- mclapply(seq_len(nrow(tune)), function(x) {
  hmm <- detectTranscripts(Fp=Fp, Fm=Fm, LtProbB=tune$LtProbB[x],
                           UTS=tune$UTS[x], threshold=1)
  e <- evaluateHMMInAnnotations(hmm$transcripts, gene.ann)
  return(list(hmm=hmm, eval=e$eval))
},mc.cores=getOption("mc.cores"), mc.silent=FALSE)

#give names based on parameters and save as .RData
names(hmm.vars) = apply(tune,1,function(x)paste0(x[2],"_",x[1]))
save(hmm.vars, file="grohmm_vars.RData")
```
Next up is creating an R object that contains sensitivity information. Merged bigWig files are required for this.
```
# load hmm data
load("grohmm_vars.RData")

# load bigwigs
bw.plus = load.bigWig("smc_merged_plus.bigWig")
bw.minus = load.bigWig("smc_merged_minus.bigWig")

# loop through hmms, check reads in complement, document metrics
eval.dat = c()

for(ii in 1:length(hmm.vars)){

  # get complement to hmm regions
  hmm.test = hmm.vars[[ii]]$hmm
  transcripts = hmm.test[["transcripts"]]
  if(is.null(transcripts)==TRUE){next}
  hmm.bed0 <- data.frame(chr=seqnames(transcripts),
                         start=start(transcripts)-1,
                         end=end(transcripts),
                         names=c(rep(".", length(transcripts))),
                         scores=c(rep(".", length(transcripts))),
                         strand=strand(transcripts))
  hmmp = hmm.bed0 %>% filter(strand=="+")
  hmmm = hmm.bed0 %>% filter(strand=="-")
  write.table(hmmp,"hmmp.bed",sep="\t",quote=F,col.names=F,row.names=F)
  write.table(hmmm,"hmmm.bed",sep="\t",quote=F,col.names=F,row.names=F)
  command1=paste('sort -k1,1 -k2,2n', 'hmmp.bed', '> hmmp.sorted.bed')
  command2=paste('sort -k1,1 -k2,2n', 'hmmm.bed', '> hmmm.sorted.bed')
  system(command1); system(command2)
  comm1 = "bedtools complement -i hmmp.sorted.bed -g hg38.sorted.bed > compp.bed"
  comm2 = "bedtools complement -i hmmm.sorted.bed -g hg38.sorted.bed > compm.bed"
  system(comm1); system(comm2)
  complplus = read.table("compp.bed",stringsAsFactors=F)
  complminus = read.table("compm.bed",stringsAsFactors=F)
  system(paste0("rm hmmp.bed hmmm.bed hmmp.sorted.bed hmmm.sorted.bed compp.bed compm.bed"))
  # map reads to the complement regions
  plus.bed = complplus %>% mutate(gene=".",xy=".",strand="+")
  minus.bed = complplus %>% mutate(gene=".",xy=".",strand="-")
  names(plus.bed)[1:3] = names(minus.bed)[1:3] = c("chr","start","end")
  cnt.plus = bed.region.bpQuery.bigWig(bw.plus, plus.bed)
  cnt.minus = bed.region.bpQuery.bigWig(bw.minus, minus.bed)
  # get count density metrics
  len.p = plus.bed$end - plus.bed$start
  len.m = minus.bed$end - minus.bed$start
  density.plus = cnt.plus / len.p
  density.minus = cnt.minus / len.m
  counts = c(cnt.plus, cnt.minus)
  densities = c(density.plus, density.minus)
  max_cnt = max(counts);
  max_den = max(densities)
  mean_cnt = mean(counts); mean_den = mean(densities)
  med_cnt = median(counts); med_den = median(densities)
  sum_cnt = sum(counts);
  sum_den = sum(densities)
  # combine evaluation metrics
  pars0 = strsplit(names(hmm.vars)[ii],"_")[[1]]
  LtProbB = pars0[2]
  UTS = pars0[1]
  hmm.eval = hmm.vars[[ii]]$eval
  new = data.frame(LtProbB, UTS, hmm.eval,
                   max_cnt, mean_cnt, med_cnt, sum_cnt,
                   max_den, mean_den, med_den, sum_den)
  eval.dat = rbind(eval.dat, new)
} # ii, hmm loop

save(eval.dat, file="hmmEval.RData")
```
### Selecting the best set of TU coordinates
With the RData files now created, there are some simple commands in R that can be used to select the GRanges object with the greatest transcription unit coverage. First load the RData file and sort the GRanges objects:
```
#load the data
load("hmmEval.RData")
load("grohmm_vars.RData")

#sort data based on merge errors, dissiciation errors, and sums of reads outside of TU
eval.dat.sorted1 = eval.dat[with(eval.dat, order(merged, dissociated, sum_cnt)),] %>%
  select(merged, dissociated, sum_cnt)
eval.dat.sorted2 = eval.dat[with(eval.dat, order(dissociated, merged, sum_cnt)),] %>%
  select(merged, dissociated, sum_cnt)
eval.dat.sorted3 = eval.dat[with(eval.dat, order(sum_cnt, dissociated, merged)),] %>%
  select(merged, dissociated, sum_cnt)

#preview each sorted data frame

head(eval.dat.sorted1)
    merged dissociated   sum_cnt
169   3383         432  96244978
170   3454         536  55574440
148   3642         548 100991615
171   3709         612  34938064
172   3909         633  25821547
149   3959         703  60909182

head(eval.dat.sorted2)
merged dissociated   sum_cnt
169   3383         432  96244978
170   3454         536  55574440
148   3642         548 100991615
171   3709         612  34938064
172   3909         633  25821547
173   4023         636  20683627

head(eval.dat.sorted3)
merged dissociated  sum_cnt
179   4704         797 12384006
178   4371         729 12544452
180   4955         849 12649531
177   4296         697 12854009
181   5137         885 13244494
176   4229         672 13533299
```
This is just an example based on what I did, but I'll explain the general process. First, look at the quartiles of the sum count density outside of the predicted transcription units:
```
quantile(eval.dat$sum_cnt)
       0%       25%       50%       75%      100%
 12384006  38254115  91629053 155365885 193436291
```
In this case, we learn that the lowest value for read counts outside transcription units is 12384006. We can look at which GRanges object this value refers to with this code:
```
#find out position of GRanges object with the lowest read counts outside of TUs
which(eval.dat$sum_cnt==12384006)

[1] 179

#preview information on the object
eval.dat[179,]

    LtProbB UTS merged dissociated total  errorRate txSize max_cnt mean_cnt
179    -500  20   4704         797  5501 0.08514557  44934 2573280  276.651
    med_cnt  sum_cnt max_den   mean_den      med_den  sum_den
179       1 12384006 155.777 0.02654926 6.535948e-05 1188.451
```
Based on this information, the GRanges object with the lowest density of read counts outside of transcription units is the one generated with the LtProbB and UTS parameters at -500 and 20, respectively. This file can be converted into a data frame and written as a .bed file, then as a GTF.
```
#make R object from list of GRanges objects in grohmm_vars
hmm.best = hmm.vars[["20_-500"]]$hmm

#extract GRanges object
transcripts = hmm.best[["transcripts"]]

#preview GRanges object
head(transcripts)

GRanges object with 6 ranges and 2 metadata columns:
      seqnames        ranges strand |  type           ID
         <Rle>     <IRanges>  <Rle> | <Rle>  <character>
  [1]     chr1   26950-27099      + |    tx  chr1_26950+
  [2]     chr1   38900-38999      + |    tx  chr1_38900+
  [3]     chr1   55800-56899      + |    tx  chr1_55800+
  [4]     chr1   79000-82249      + |    tx  chr1_79000+
  [5]     chr1   92650-92749      + |    tx  chr1_92650+
  [6]     chr1 128750-130299      + |    tx chr1_128750+
  -------
  seqinfo: 25 sequences from an unspecified genome; no seqlengths

#make an R object in the .bed format
hmm.bed0 <- data.frame(chr=seqnames(transcripts),
                       start=start(transcripts)-1,
                       end=end(transcripts),
                       names=c(rep(".", length(transcripts))),
                       scores=c(rep(".", length(transcripts))),
                       strand=strand(transcripts))

#make an R object in gtf format
emptyvals=rep(".",length(hmm.bed0[,1]))
hmm.gtf0 <- data.frame(chr=hmm.bed0$chr,src=rep("ASElux",length(hmm.bed0[,1])),feature=rep("transcript",
                       length(hmm.bed0[,1])),start=hmm.bed0$start,end=hmm.bed0$end,score=emptyvals,
                       strand=hmm.bed0$strand,frame=emptyvals,attr=emptyvals)

#write gtf file
write.table(hmm.gtf0,file="best_HMM.gtf",quote=F,sep="\t",row.names=F,col.names=F)
```
### Correctly formatting the GTF
ASElux requires that chromosome coordinates include "gene", "transcript", and "exon" in the feature column. I found this out when the GTF file created above didn't work while building the index for ASElux. Instead of properly building the index, ASElux would put out a fatal "segmentation fault" error. [This](https://github.com/abl0719/ASElux/issues/5) thread started by Nick Garza on the ASElux Github page gives some insight on the specifics for the GTF format and helped me figure out what was wrong with the above file. Use this code in the Unix command line to sort the GTF by chromosome number, although this might not be necessary:
```
sort -k1,1 -k2,2n best_HMM1.gtf > sorted_best_HMM.gtf
```
Use this code in R to duplicate each set of coordinates and add "gene", "transcript", and "exon" to the feature column.
```
nwork <- read.table("sorted_best_HMM.gtf",stringsAsFactors=F,header=F)

#adding gene, transcript, exon
fixhmm <- c()
for (i in 1:length(nwork[,1])){
  fixhmm <- rbind(fixhmm,nwork[i,],nwork[i,],nwork[i,])
  fixhmm$V3 <- c("gene","transcript","exon")
}

#create new GTF file
write.table(fixhmm,"fixed_best_HMM.gtf",row.names=F,col.names=F,sep="\t",quote=F)
```
The head of the newly created GTF should look something like this:
```
chr1	ASElux	gene	100783049	101419549	.	+	.	.
chr1	ASElux	transcript	100783049	101419549	.	+	.	.
chr1	ASElux	exon	100783049	101419549	.	+	.	.
chr1	ASElux	gene	101416349	101507649	.	-	.	.
chr1	ASElux	transcript	101416349	101507649	.	-	.	.
chr1	ASElux	exon	101416349	101507649	.	-	.	.
chr1	ASElux	gene	101430249	101431149	.	+	.	.
chr1	ASElux	transcript	101430249	101431149	.	+	.	.
chr1	ASElux	exon	101430249	101431149	.	+	.	.
chr1	ASElux	gene	101442349	101642399	.	+	.	.
```
With that done, "fixed_best_HMM.gtf" should be ready to be used for an ASElux index.
