#Running DADA2 on MPC communities from the beginning!
#Glade Dlott
#July 22 2019
#Pipeline for 16S amplicon processing

#Glade Dlott, 2019_June_17

#From: https://benjjneb.github.io/dada2/ITS_workflow.html accessed 2019 June 2019

library(dada2)
packageVersion("dada2")

library(ShortRead)
packageVersion("ShortRead")

library(Biostrings)
packageVersion("Biostrings")

path <- "~/Desktop/MPCREVISION2019/MPC_Illumina_Data/"
list.files(path)

fnFs <- sort(list.files(path, pattern = ".R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = ".R2.fastq", full.names = TRUE))

#designate primers

FWD <- "GGMTTAGATACCC"
REV <- "CCGYCAATTYMTTTRAGTTT"

#Test primer removal

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

#Remove Ns

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))

#Next step completely crashed the computer last time; disabling multithread... ran in about 25 minutes

filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = FALSE, compress = FALSE)

#Count primers

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

#OUTPUT:
#Forward Complement Reverse RevComp
#FWD.ForwardReads 1661953          0       0       1
#FWD.ReverseReads       1          0       0      30
#REV.ForwardReads       0          0       0      21
#REV.ReverseReads 1638942          0       0       0

#Assume that reverse primers were truncated already. 

#remove primers

cutadapt <- "/Users/gadlott/miniconda3/bin/cutadapt"
system2(cutadapt, args = "--version")

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, "-m", 50,
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

#Look for primers:

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

#1 'incorrect' forward/reverse reads remain. Not sure why but its literally one, so, ignoring?

#Proceeding: only concern: 

#Bases preceding removed adapters:
#A: 0.7%
#C: 0.6%
#G: 1.5%
#T: 97.2%
#none/other: 0.0%
#WARNING:
#  The adapter is preceded by "T" extremely often.
#The provided adapter sequence could be incomplete at its 3' end.

#This happened pretty often, may be a source of trouble.


# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = ".R1.fastq", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = ".R2.fastq", full.names = TRUE))

#trying to remove future 0-samples:

cutFs <- cutFs[-104]
cutFs <- cutFs[-187]

cutRs <- cutRs[-104]
cutRs <- cutRs[-187]

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)
sample.names

#read quality profiles
plotQualityProfile(cutFs[1:2])
plotQualityProfile(cutRs[1:2])

#filter and trim!
filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

#Add in truncation; expect a 142-length amplicon; use truncation lengths 77F, 77R (12 bp overlap).

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), truncLen=c(77,77),
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, verbose = TRUE, multithread = FALSE)
head(out)

#lost quite a few reads and also all of: oct4 blank, 29r3. Attempting to rename as I did before. I really think this shoukd work without having to childishly remove files are repeat cutadapt. 

#Did it. I missed one without an error before because I removed two, sequentially. When you remove one, the position of the other one changes.

#Learn Error rates

errF <- learnErrors(filtFs, multithread = FALSE)
errR <- learnErrors(filtRs, multithread = FALSE)

#note - weird divercence from expected printout here. Must have changed verbose code.

plotErrors(errF, nominalQ = TRUE)

#This looks even worse than in the fungal dataset. Why is the error rate the same regardless of read quality?

#dereplicate reads:

#NOTE here is the source of error for 0-samples.

derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

names(derepFs) <- sample.names
names(derepRs) <- sample.names

#sample inference

dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

#merge paired reads, minimum 12

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

#construct sequence table

seqtab <- makeSequenceTable(mergers)
dim(seqtab)

#remove chimeras

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

#sequence length distribution

table(nchar(getSequences(seqtab.nochim)))

#..Why are all the sequences 253 bp? -primers? 

#track reads

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, 
                                                                       getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", 
                     "nonchim")
rownames(track) <- sample.names
head(track)

#This is it. Good.

#assign taxonomy

###

#RDP
#taxa <- assignTaxonomy(seqtab.nochim, "~/Desktop/MPC_Illumina_Data/rdp_train_set_16.fa", multithread=FALSE)
#taxa <- addSpecies(taxa, "~/Desktop/MPC_Illumina_Data/rdp_species_assignment_16.fa")

#GreenGenes
#taxaGG <- assignTaxonomy(seqtab.nochim, "~/Desktop/MPC_Illumina_Data/GreenGenes/gg_13_8_train_set_97.fa", multithread=FALSE)
#taxaGG <- addSpecies(taxa, "~/Desktop/MPC_Illumina_Data/GreenGenes/rdp_species_assignment_14.fa")

#Silva
taxaSilva <- assignTaxonomy(seqtab.nochim, "~/Desktop/MPC_Illumina_Data/Silva/silva_nr_v132_train_set.fa", multithread=FALSE)
taxaSilva <- addSpecies(taxa, "~/Desktop/MPC_Illumina_Data/Silva/silva_species_assignment_v132.fa")

#RDP
#taxa.print <- taxa  # Removing sequence rownames for display only
#rownames(taxa.print) <- NULL
#head(taxa.print)

#GreenGenes
#taxa.printGG <- taxaGG  # Removing sequence rownames for display only
#rownames(taxa.printGG) <- NULL
#head(taxa.printGG)

#Silva
taxa.printSilva <- taxaSilva  # Removing sequence rownames for display only
rownames(taxa.printSilva) <- NULL
head(taxa.printSilva)

###

library(phyloseq); packageVersion("phyloseq")

library(Biostrings); packageVersion("Biostrings")

library(ggplot2); packageVersion("ggplot2")

theme_set(theme_bw())

#create phyloseq object; already have otu table (sectab.nochim), taxonomy (taxa), only need sample data

#create a dataframe in phyloseq style. If you’re skilled with R it may be faster for you to put this together ad-hoc; i will import mine from a .csv.

samdf <- read.csv("~/Desktop/MPC_Illumina_Data/MPC_Data_2019_07_18.csv")
rownames(samdf) <- sample.names

#create phyloseq object

#RDP:
#psMPC <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                  #sample_data(samdf), 
                  #tax_table(taxa))
#psMPC

#GreenGenes:
#psMPC_GG <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                  #sample_data(samdf), 
                  #tax_table(taxaGG))
#psMPC_GG

#Silva:
psMPC_Silva <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                  sample_data(samdf), 
                  tax_table(taxaSilva))
psMPC_Silva

#export phyloseq object

#RDP
#saveRDS(psMPC, "~/Desktop/MPC_Illumina_Data/MPC_16S_phyloseq_2019_07_23.rds")

#Greengenes
#saveRDS(psMPC_GG, "~/Desktop/MPC_Illumina_Data/GreenGenes/MPC_GreenGenes_16S_phyloseq_2019_07_23.rds")

#Silva
saveRDS(psMPC_Silva, "~/Desktop/MPC_Illumina_Data/Silva/MPC_Silva_16S_phyloseq_2019_07_23.rds")

### END ###

#Glade Dlott, 2019 July 23, Verify earlier findings...

#Look at MPC phyloseq OTU table

#OTU(ASV) Table
write.csv(otu_table(psMPC), "~/Desktop/MPC_Illumina_Data/MPC_otu_table_2019_07_23.csv")

#RDP Taxonomy
#write.csv(tax_table(psMPC), "~/Desktop/MPC_Illumina_Data/MPC_taxonomy_table_2019_07_23.csv")

#Greengenes Taxonomy
#write.csv(tax_table(psMPC_GG), "~/Desktop/MPC_Illumina_Data/GreenGenes/MPC_GreenGenes_taxonomy_table_2019_07_23.csv")

#Silva Taxonomy
write.csv(tax_table(psMPC_Silva), "~/Desktop/MPC_Illumina_Data/Silva/MPC_Silva_taxonomy_table_2019_07_23.csv")

