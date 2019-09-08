# Loading libraries
library(dada2); packageVersion("dada2")
library(data.table)
library(phyloseq)
library(openssl)

### 16S Pipeline
# Step 1: Directories and filenames
###############################################################################################################

path_fastq <- "Demultiplexed_fastq_files/Run1_16S/fastq_16S/" # Directory containing the fastq files after unzipping.
path_main <- getwd() # Working directory
dir.create(paste0(path_main,"/16S"))
path_results <- paste0(path_main,"/16S")

# Assuming forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path_fastq, pattern="R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path_fastq, pattern="R2_001.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_SXXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_S\\d+"), `[`, 1)

# Step 2: Sequence quality filtering
###############################################################################################################

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path_results, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path_results, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
# Filter and trim reads
# Forward and reverse reads are truncated at 226 bp and 220 bp for 16S, and 225 and 216 for 18S, respectively
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=c(0,0), 
                     truncLen=c(226,220),maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,compress=TRUE, multithread=T)

#Step 4: Learn error rates
###############################################################################################################

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

#Step 5: Dereplication
###############################################################################################################

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Step 6: Removing sequencing erros from the learned error rates
###############################################################################################################

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

# Step 7: Merging forward and reverse reads
###############################################################################################################

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, minOverlap = 10, maxMismatch = 0, verbose=TRUE)

# Step 8: Construct feature table
###############################################################################################################

seqtab <- makeSequenceTable(mergers)
# Inspect distribution oflengths
table(nchar(getSequences(seqtab)))

# Step 9: Remove chimeras
###############################################################################################################
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

# Step 10: Track read number throughout the pipeline
###############################################################################################################

getN <- function(x) sum(getUniques(x))
track <- as.data.frame(cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim)))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
Sample_names <- as.data.frame(sample.names)
head(track)
track = cbind(Sample_names,track)
write.table(as.data.frame(track), paste0(path_results,"/Read_counts.csv"), sep=",", quote = F, row.names = F)

# Step 11: Exporting the features' table and features' fasta file for taxonomic assignment in Qiime2
###############################################################################################################

otab <- otu_table(seqtab.nochim, taxa_are_rows=FALSE)
colnames(otab) <- md5(colnames(seqtab.nochim))
Names = as.character(colnames(otab))
otab = t(otab)
otab = cbind(ASVs = rownames(otab), otab)
write.table(otab, paste0(path_results,"/dada_table.txt"),quote=FALSE,sep="\t", row.names = FALSE)
uniquesToFasta(getUniques(seqtab.nochim), fout= paste0(path_results,"/rep_ASVs.fasta"),ids=as.character(as.list(md5(names(getUniques(seqtab.nochim))))))














### 18S Pipeline
# Step 1: Directories and filenames
###############################################################################################################

path_fastq <- "Demultiplexed_fastq_files/Run1_18S/fastq_18S/" # Directory containing the fastq files after unzipping.
path_main <- getwd() # Working directory
dir.create(paste0(path_main,"/18S"))
path_results <- paste0(path_main,"/18S")

# Assuming forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path_fastq, pattern="R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path_fastq, pattern="R2_001.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_SXXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_S\\d+"), `[`, 1)

# Step 2: Sequence quality filtering
###############################################################################################################

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path_results, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path_results, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
# Filter and trim reads
# Forward and reverse reads are truncated at 226 bp and 220 bp for 16S, and 225 and 216 for 18S, respectively
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=c(0,0), 
                     truncLen=c(225,216),maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,compress=TRUE, multithread=T)

#Step 4: Learn error rates
###############################################################################################################

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

#Step 5: Dereplication
###############################################################################################################

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Step 6: Removing sequencing erros from the learned error rates
###############################################################################################################

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

# Step 7: Merging forward and reverse reads
###############################################################################################################

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, minOverlap = 10, maxMismatch = 0, verbose=TRUE)

# Step 8: Construct feature table
###############################################################################################################

seqtab <- makeSequenceTable(mergers)
# Inspect distribution oflengths
table(nchar(getSequences(seqtab)))

# Step 9: Remove chimeras
###############################################################################################################
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

# Step 10: Track read number throughout the pipeline
###############################################################################################################

getN <- function(x) sum(getUniques(x))
track <- as.data.frame(cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim)))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
Sample_names <- as.data.frame(sample.names)
head(track)
track = cbind(Sample_names,track)
write.table(as.data.frame(track), paste0(path_results,"/Read_counts.csv"), sep=",", quote = F, row.names = F)

# Step 11: Exporting the features' table and features' fasta file for taxonomic assignment in Qiime2
###############################################################################################################

otab <- otu_table(seqtab.nochim, taxa_are_rows=FALSE)
colnames(otab) <- md5(colnames(seqtab.nochim))
Names = as.character(colnames(otab))
otab = t(otab)
otab = cbind(ASVs = rownames(otab), otab)
write.table(otab, paste0(path_results,"/dada_table.txt"),quote=FALSE,sep="\t", row.names = FALSE)
uniquesToFasta(getUniques(seqtab.nochim), fout= paste0(path_results,"/rep_ASVs.fasta"),ids=as.character(as.list(md5(names(getUniques(seqtab.nochim))))))
