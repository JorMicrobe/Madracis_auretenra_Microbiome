#-------------------------------------------------------------------------------
#Microbiome analysis based on amplicon sequencing 16S rRNA
#ANALISIS DEL MICROBIOMA DEL CORAL Madracis auretenra BASADO EN SECUENCIACION AMPLICON DEL GEN 16S ARNr REGION HIPERVARIABLE V4
#PIPELINE DADA2
#Instalar los paquetes requeridos y la ruta a las secuencias .fastq

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.18")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")

library(dada2); packageVersion("dada2")
library(BiocManager)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(tidyr)
library(stringr)
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")
library(vegan)
library(lessR)
library(ape)
library(data.table)

#Adicionar la ruta donde están las sequencias en bruto (raw sequences)

path <- "/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/16S_Seq_Ma_2024"

list.files(path)

#-------------------------------------------------------------------------------
#ANALISIS PRELIMINAR DE LAS SECUENCIAS
#Organizar las lecturas de acuerdo a su nombre

fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
head(sample.names)

#Inspeccionar la calidad de las sequencias Forward y Reverse

plotQualityProfile(fnFs[1:2:3:4:5:6])
plotQualityProfile(fnFs[7:8:9:10:11:12])
plotQualityProfile(fnFs[13:14:15:16:17:18])
plotQualityProfile(fnFs[19:20:21:22:23:24])
plotQualityProfile(fnFs[25:26:27:28:29:30])
plotQualityProfile(fnFs[31:32:33:34:35:36])
plotQualityProfile(fnFs[37:38:39:40:41:42])
plotQualityProfile(fnFs[43:44:45:46:47:48])
plotQualityProfile(fnFs[49:50:51:52:53:54])
plotQualityProfile(fnFs[55:56:57:58:59:60])
plotQualityProfile(fnFs[61:62:63:64:65:66])
plotQualityProfile(fnFs[67:68:69:70:71:72])

plotQualityProfile(fnRs[1:2:3:4:5:6])
plotQualityProfile(fnRs[7:8:9:10:11:12])
plotQualityProfile(fnRs[13:14:15:16:17:18])
plotQualityProfile(fnRs[19:20:21:22:23:24])
plotQualityProfile(fnRs[25:26:27:28:29:30])
plotQualityProfile(fnRs[31:32:33:34:35:36])
plotQualityProfile(fnRs[37:38:39:40:41:42])
plotQualityProfile(fnRs[43:44:45:46:47:48])
plotQualityProfile(fnRs[49:50:51:52:53:54])
plotQualityProfile(fnRs[55:56:57:58:59:60])
plotQualityProfile(fnRs[61:62:63:64:65:66])
plotQualityProfile(fnRs[67:68:69:70:71:72])

#Filtrar y cortar las secuencias con base en el análisis de calidad
#La región V4 del gen 16S ARNr tienen un tamaño de 254 pb por ello no cortar tanto para no perder información

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(220,200),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)

head(out)

#Verificar la tasa de errores

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#Hacer la dereplicación, es decir, combinar las secuencias similares en una única

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

names(derepFs) <- sample.names
names(derepRs) <- sample.names

#Aplicar el algoritmo de inferencia de la muestra

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
dadaFs[[1]]

#Fusionar las lecturas pareadas Forward y Reverse (el resultado son lecturas sin ruido)

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
head(mergers[[1]])

#-------------------------------------------------------------------------------
#CONSTRUCCION DE LA TABLA DE ASVs, ASIGNACION TAXONOMICA Y REMOCION DE SECUENCIAS NO DESEADAS
#Construir la tabla de secuencias variantes amplicon (ASVs)

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

#Remover las quimeras

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

#Hacer seguimiento de lecturas a lo largo del pipeline

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)
write.csv(track, "/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/track_Ma.csv")

#Asignar la taxonomia con la base de datos SILVA hasta nivel de genero (>97% de identidad)

taxa <- assignTaxonomy(seqtab.nochim, "/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/16S_Seq_Ma_2024/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/16S_Seq_Ma_2024/silva_species_assignment_v138.1.fa.gz")
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)

#Remover las lecturas asignadas a Cloroplastos

is.chloroplast <- taxa[,"Order"] %in% "Chloroplast"
seqtab.nochloro <- seqtab.nochim[,!is.chloroplast]
dim(seqtab.nochloro)
taxa.nochloro <- taxa[!is.chloroplast,]
dim(taxa.nochloro)

#Remover las lecturas asignadas a Mitocondrias

is.mitochondria <- taxa.nochloro[,"Family"] %in% "Mitochondria"
seqtab.nomito <- seqtab.nochloro[,!is.mitochondria]
taxa.nomito <- taxa.nochloro[!is.mitochondria,]
dim(seqtab.nomito)
dim(seqtab.nochloro)
dim(seqtab.nochim)

#Organizar la tabla y crear el archivo de excel que contiene los taxa

seqtab_nochim <- seqtab.nomito
dim(seqtab_nochim)
taxa <- taxa.nomito
head(taxa)
write.csv(taxa, "/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/Ma_taxa_2024.csv")

library(lessR)
library(S4Vectors)

rows <- nrow(taxa)
asvs <- lessR::to("asv", rows)
sample.asvs <- seqtab_nochim
flipped.sample.asvs <- t(sample.asvs)
sample.asvs <- flipped.sample.asvs
seqs <- rownames(sample.asvs)
length(seqs)
seqs <- DNAStringSet(seqs)
names(seqs) <- asvs
dim(seqs)
rownames(taxa) <- asvs
rownames(sample.asvs) <- asvs

dim(sample.asvs)
length(asvs)

head(sample.asvs)

#Crear los archivos que contienen las ASVs asignadas a las taxa y las ASVs de cada muestra

write.csv(taxa, "/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/Ma_asvs_taxa_2024.csv")
write.csv(sample.asvs, "/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/Ma_sample.asvs_taxa_2024.csv")
