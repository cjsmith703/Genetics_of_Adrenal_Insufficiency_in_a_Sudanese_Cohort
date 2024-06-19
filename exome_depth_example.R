#title "Exome_depth_example_with build38)"
#author: Chris J Smith

####Install Dependencies##################################################
#if (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")

#BiocManager::install(c("Biostrings", "IRanges", "Rsamtools", "GenomicRanges", "GenomicAlignments"))
#install.packages("aod")
#install.packages("VGAM")

#Install exome depth from tar.gz file if not added back to CRAN

############################################################################

#Load packages
library(tidyverse)
library(ExomeDepth)
library(data.table)

#set bam files to character vector
bamFiles <- c("/path/to/WES/bam/proband.recal.bam", 
              "/path/to/WES/bam/unrelated_panel_1.recal.bam", #panel samples should be from the same sequencing batch
              "/path/to/WES/bam/unrelated_panel_2.recal.bam", 
              "/path/to/WES/bam/unrelated_panel_3.recal.bam", 
              "/path/to/WES/bam/unrelated_panel_4.recal.bam",
              "/path/to/WES/bam/unrelated_panel_5.recal.bam",
              "/path/to/WES/bam/unrelated_panel_6.recal.bam",
              "/path/to/WES/bam/unrelated_panel_7.recal.bam",
              "/path/to/WES/bam/unrelated_panel_8recal.bam",
              "/path/to/WES/bam/unrelated_panel_9.recal.bam",
              "/path/to/WES/bam/unrelated_panel_10.recal.bam",)

fasta <- c("/path/to/Homo_sapiens_assembly38.fasta")

bed <- read.csv("path/to/hg38_bed.csv")

bed <- bed %>% #manipulate bed file to match bam file
  mutate(chrom = gsub('chr', '', chrom)) %>%
  rename(chromosome = chrom) %>%
  rename(name = gene_name)

#Generate read count data
my.counts <- getBamCounts(bed.frame = bed,
                          bam.files = bamFiles,
                          include.chr = TRUE,
                          referenceFasta = fasta)

#set as dataframe and remove 'chr'
my.counts.dafr <- as(my.counts, 'data.frame')
my.counts.dafr$chromosome <- gsub(as.character(my.counts.dafr$chromosome), 
                                   pattern = 'chr', 
                                   replacement = '')
head(my.counts.dafr)

#build reference set
my.test<- my.counts$proband.recal.bam

my.ref.samples <- c("unrelated.panel.1.recal.bam",
                    "unrelated.panel.2.recal.bam",
                    "unrelated.panel.3.recal.bam",
                    "unrelated.panel.4.recal.bam",
                    "unrelated.panel.5.recal.bam",
                    "unrelated.panel.6.recal.bam",
                    "unrelated.panel.7.recal.bam",
                    "unrelated.panel.8.recal.bam",
                    "unrelated.panel.9.recal.bam",
                    "unrelated.panel.10.recal.bam")

my.reference.set <- as.matrix(my.counts.dafr[, my.ref.samples])

my.choice <- select.reference.set(test.counts = my.test,
                                  reference.counts = my.reference.set,
                                  bin.length = (my.counts.dafr$end - my.counts.dafr$start)/1000,
                                  n.bins.reduced = 10000)

print(my.choice[[1]])

my.matrix <- as.matrix(my.counts.dafr[, my.choice$reference.choice, drop = FALSE])

#selects most similar bams form panel
my.reference.selected <- apply(X = my.matrix,
                               MARGIN = 1,
                               FUN = sum)
#call cnvs
all.exons <- new('ExomeDepth',
                 test = my.test,
                 reference = my.reference.selected,
                 formula = 'cbind(test, reference) ~ 1')


all.exons <- CallCNVs(x = all.exons,
                      transition.probability = 10^-4,
                      chromosome = my.counts.dafr$chromosome,
                      start = my.counts.dafr$start,
                      end = my.counts.dafr$end,
                      my.counts.dafr$exon)

#filter out common cnvs
data(Conrad.hg19)
head(Conrad.hg19.common.CNVs)

all.exons <- AnnotateExtra(x = all.exons,
                           reference.annotation = Conrad.hg19.common.CNVs,
                           min.overlap = 0.5,
                           column.name = 'Conrad.hg19')


#convert to dataframe
cnv_analysis <- as.data.frame(all.exons@CNV.calls)

head(cnv_analysis)

#re-read in bed file
bed <- read.csv("hg38.csv")

bed <- bed %>%
  mutate(chrom = gsub('chr', '', chrom))

#function to add gene name for build 38 to cnv analysis 

add_gene_name <- function(data, bed) {
  gene_names <- character(nrow(data))  
  for (i in seq_len(nrow(data))) {
      filtered_bed <- subset(bed, chrom == data$chromosome[i] & start >= data$start[i] & end <= data$end[i])
      if (nrow(filtered_bed) > 0) {
        gene_names[i] <- paste(filtered_bed$gene_name, collapse = ",")
      } else if (nrow(filtered_bed) == 0) {
        filtered_bed <- subset(bed, chrom == data$chromosome[i] & start <= data$start[i] & end >= data$end[i])
        gene_names[i] <- paste(filtered_bed$gene_name, collapse = ",")
      }
        else {
          gene_names[i] <- NA
      }
  }
  data$gene_name <- gene_names
  return(data)
}

#run function on cnv analysis
annotated_cnv <- add_gene_name(cnv_analysis, bed)

head(annotated_cnv)

#write csv to output folder
write.csv(annotated_cnv, "output/proband_exome_depth.csv")


#plot gene of interest
plot (all.exons,
      sequence = '1',
      xlim = c(44280000, 44400000),
      count.threshold = 20,
      main = 'GENE',
      cex.lab = 0.8,
      with.gene = TRUE)


