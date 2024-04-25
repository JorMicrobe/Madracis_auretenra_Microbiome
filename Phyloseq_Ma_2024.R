#-------------------------------------------------------------------------------
#ANALISIS DEL MICROBIOMA BASADO EN SECUENCIACION AMPLICON DEL GEN 16S ARNr
#ANALISIS BASADO EN PHYLOSEQ Y OTROS
#Instalar y cargar los paquetes necesarios

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.18")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DECIPHER")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("Biostrings", "seqLogo"))

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")

library(BiocManager)
BiocManager::install("microbiome")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ANCOMBC")

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
library(breakaway)
library(microbiome)
library(ggpubr)
library(kableExtra)
library (DESeq2)
library(ComplexHeatmap)
library(ANCOMBC)
library(DECIPHER)
library(phangorn)
#-------------------------------------------------------------------------------
#PHYLOSEQ pipeline (modificado por Jordan Ruiz, 2024)

#Hacer los anlisis usando el paquete de Phyloseq (McMurdie y Holmes, 2013)
#El obejto de phyloseq tiene 4 componentes: la tabla de otus, los datos de las muestras, la tabla de taxa, y el arbol filogenetico

#Datos de las muestras (Metadata)
#Importar los metadatos y fusionar las columnas codigo y muestra
Metadata <- read.csv("/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/Metadata_2.csv", sep = ",", header = TRUE, check.names = FALSE)
Metadata$Sample <- paste(Metadata$Code, Metadata$Sample, sep = "-")
colnames(Metadata)
head(Metadata)
Metadata2 <- select(Metadata, -Code)
head(Metadata2)
Year <- as.factor(Metadata2$Year)
Season <- as.factor(Metadata2$Season)
Site <- as.factor(Metadata2$Site)
Health_Status <- as.factor(Metadata2$Health_Status)
Algae_Touching <- as.factor(Metadata2$Algae_Touching)
Source <- as.factor(Metadata2$Source)

#Tabla de otus y tabla de taxa
#Crear las matrices que van a hacer parte del objeto de phyloseq
asv_mat<- read.csv("/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/Ma_sample.asvs_taxa_decontam_2024_2.csv", sep = ",", header = TRUE, check.names = FALSE)
tax_mat<- read.csv("/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/Ma_asvs_taxa_corrected_2024_2.csv", sep = ",", header = TRUE, check.names = FALSE)

head(asv_mat)
head(tax_mat)

asv_mat <- asv_mat %>% 
  tibble::column_to_rownames("asv")

tax_mat <- tax_mat %>%
  tibble::column_to_rownames("asv")

Metadata2 <- Metadata2 %>%
  tibble::column_to_rownames("Sample")

asv_mat <- as.matrix(asv_mat)
tax_mat <- as.matrix(tax_mat)

OTU = otu_table(asv_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
Metadata = sample_data(Metadata2)

#OBJETO DE PHYLOSEQ Y EXPLORACION DE LOS DATOS
#Crear el objeto de phyloseq
Ma <- phyloseq(OTU, TAX, Metadata)
Ma

#Arbol filogenetico:
#Hacer el alineamiento de las secuencias para crear el arbol filogenetico
random_tree = rtree(ntaxa(Ma), rooted=TRUE, tip.label=taxa_names(Ma))
plot(random_tree)

#Unir el objeto de phyloseq con el arbol
Ma = merge_phyloseq(Ma, random_tree)
Ma

sample_names(Ma)
rank_names(Ma)
sample_variables(Ma)

head(Metadata)

total = median(sample_sums(Ma))
standf = function(x, t=total) round(t * (x / sum(x)))
Ma = transform_sample_counts(Ma, standf)

plot_bar(Ma, fill = "Phylum")

plot_bar(Ma, fill = "Phylum") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

#Revisar el numero de lecturas por muestras, si hay alguna muestra con muy pocas secuencias, se debe remover
#Todas las muestras tienen entre 17000 y 17400 secuencias por lo que no se debe remover ninguna
readcount = data.table(as(sample_data(Ma), "data.frame"),
                       TotalReads = sample_sums(Ma), 
                       keep.rownames = TRUE)
setnames(readcount, "rn", "SampleID")
ggplot(readcount, aes(TotalReads)) + geom_histogram() + ggtitle("Sequencing Depth")
head(readcount[order(readcount$TotalReads), c("SampleID", "TotalReads")])

write.csv(readcount, 
          file="/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/read_count_table.csv")

#Calcular el tamano de las librerias por muestra
tab <- otu_table(asv_mat, taxa_are_rows = TRUE)
class(tab) <- "matrix"
tab <- t(tab) # transpose observations to rows
sum_seq <- rowSums(tab)
plot(sum_seq, ylim=c(0,50000), main=c("Number of counts per sample"), xlab=c("Samples"), ylab=c("Number of reads"))
sum_seq
mean(sum_seq)
sd(sum_seq)
min(sum_seq)
max(sum_seq)

#Hacer la curva de rarefraccion para determinar la profundidad de la secuenciacon
#La profundidad nos dice si se secuencio lo suficiente para representar la diversidad de la comunidad
rare <- rarecurve(tab, step=10000, lwd=2, ylab="OTU",  label=F)
par(cex = 1, las = 1)
leg.txt <- c("Urban", "Protected")
lty_vector <- c(2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1)
rare <- rarecurve(tab, 
                  step = 100, 
                  lwd=1.3, 
                  xlab = "Sample size", 
                  ylab = "Number of ASV observed", 
                  xlim = c(-50, 50000), 
                  ylim = c(-5, 1000), 
                  label = FALSE, 
                  lty = lty_vector)
legend(32000, 900, leg.txt, lty = c(2,1), lwd = 1.3, box.lwd = 0.6, cex = 1)

source("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/load-extra-functions.R")

rare1 <- ggrare(Ma, step = 100, color = "Season_Year", plot = TRUE, parallel = TRUE, se = FALSE)
figure_rare <- rare1 + theme_classic() +
  labs(x = "\nSample size", y = "Number of ASV observed\n") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        axis.ticks.length = unit(.25, "cm"), 
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "gray")) +
  scale_x_continuous(labels=comma) +
  scale_color_manual(values=c("#238A8DFF", "#FDE725FF", "#6ece58", "#440154FF")) +
  theme(plot.margin = margin(0.25,0.25,0.25,0.25, "inches"))
figure_rare

ggsave("/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/Plot_rarefraction_curve_2024_season.png", 
       figure_rare, width = 8, height = 6, dpi = 600)

#Filtrar la taxonomia, remover filos que tienen solo 1 y el NA los cuales son artifactos
table(tax_table(Ma) [, "Phylum"], exclude = NULL)
Ma <- subset_taxa(Ma, !is.na(Phylum) & !Phylum %in% c("", "Acetothermia", "Aenigmarchaeota", "Asgardarchaeota", 
                                                      "Calditrichota", "Deferrisomatota", "Elusimicrobiota", "Halobacterota", "NKB15", "RCP2-54", "Zixibacteria", "<NA>"))

Ma

plot_bar(Ma, fill = "Phylum") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

ggsave("/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/bar_plot_phylum_filtered_2024.png", 
       plot = (plot_bar(Ma, fill = "Phylum")), width = 8, height = 6, dpi = 600)

#Hacer un analisis de la prevalencia de los filos en las muestras
prevdf <- apply(X = otu_table(Ma),
                MARGIN = ifelse(taxa_are_rows(Ma), yes = 1, no = 2),
                FUN = function(x){sum(x > 0)})
prevdf <- data.frame(Prevalence = prevdf,
                     TotalAbundance = taxa_sums(Ma),
                     tax_table(Ma))
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

prevdf1 <- subset(prevdf, Phylum %in% get_taxa_unique(Ma, "Phylum"))
fig_prev <- ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(Ma),color=Phylum)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

fig_prev

ggsave("/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/plot_phylum_prevalence_2024.png", 
       fig_prev, width = 8, height = 6, dpi = 600)

#Remover las ASVs que estan por debajo de una prevalencia del 5%
prevalenceThreshold <- 0.05 * nsamples(Ma)
keepTaxa <- rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
Ma1 <- prune_taxa(keepTaxa, Ma)
saveRDS(Ma1, "Ma1.rds")

Ma1

total = median(sample_sums(Ma))
standf = function(x, t=total) round(t * (x / sum(x)))
Ma1 = transform_sample_counts(Ma1, standf)

plot_bar(Ma1, fill = "Phylum") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

plot_bar(Ma1, fill = "Family") + 
  geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack")

write.csv(cbind(data.frame(otu_table(Ma1)),
                tax_table(Ma1)), 
          file="/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/filtered_otu_table.csv")

#Hacer la rarefraccion (normalizar) de los datos para los analisis de diversidad
set.seed(111); .Random.seed
Ma_rarefied <- rarefy_even_depth(Ma1, 
                               sample.size = min(colSums(otu_table(Ma1))), 
                               rngseed = 63)
Ma_rarefied

head(phyloseq::sample_sums(Ma_rarefied))

write.csv(cbind(data.frame(otu_table(Ma_rarefied)),
                tax_table(Ma_rarefied)), 
          file="/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/filtered_rarefied_otu_table.csv")


#-------------------------------------------------------------------------------
#ANALISIS FILOGENETICO Y ARBOLES

#Aglomerar las taxa para la construccion del arbol
length(get_taxa_unique(Ma1, taxonomic.rank = "Genus"))
Matree = tax_glom(Ma1, "Genus", NArm = TRUE)
h1 = 0.4
Matree1 = tip_glom(Matree, h = h1)
Matree1 = plot_tree(Matree1, method = "treeonly",
                    ladderize = "left", title = "By Height") +
  theme(plot.title = element_text(size = 12))
Matree1

(subset_taxa(Ma1, Family=="Endozoicomonadaceae") -> Ma1.en)
prune_samples(sample_sums(Ma1.en)>=5, Ma1.en) -> Ma1.en
tree1 <- plot_tree(Ma1.en, color="Source", label.tips="Genus", shape="Site")
tree1 <- tree1 + theme(legend.text = element_text(size = 12),
                       legend.title = element_text(size = 14),
                       plot.margin = margin(0.25, 0.25, 0.25, 0.25, "inches"),
                       legend.box.background = element_rect(colour = "gray")) +
  scale_color_manual(values=c("#440154FF", "#6ece58", "#FDE725FF")) +
  guides(color = guide_legend(override.aes = list(size = 5)))
tree1

ggsave("/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/plot_tree1_2024.png", 
       tree1, width = 8, height = 8, dpi = 600)

(subset_taxa(Ma1, Family=="Vibrionaceae") -> Ma1.vib)
prune_samples(sample_sums(Ma1.vib)>=5, Ma1.vib) -> Ma1.vib
tree2 <- plot_tree(Ma1.vib, color="Source", label.tips="Genus", shape="Site")
tree2 <- tree2 + theme(legend.text = element_text(size = 12),
                       legend.title = element_text(size = 14),
                       plot.margin = margin(0.25, 0.25, 0.25, 0.25, "inches"),
                       legend.box.background = element_rect(colour = "gray")) +
  scale_color_manual(values=c("#440154FF", "#6ece58", "#FDE725FF")) +
  guides(color = guide_legend(override.aes = list(size = 6)))
tree2

ggsave("/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/plot_tree2_2024.png", 
       tree2, width = 8, height = 8, dpi = 600)

(subset_taxa(Ma1, Family=="Rhodobacteraceae") -> Ma1.rod)
prune_samples(sample_sums(Ma1.rod)>=5, Ma1.rod) -> Ma1.rod
tree3 <- plot_tree(Ma1.rod, color="Source", label.tips="Genus", shape="Site")
tree3 <- tree3 + theme(legend.text = element_text(size = 12),
                       legend.title = element_text(size = 14),
                       plot.margin = margin(0.25, 0.25, 0.25, 0.25, "inches"),
                       legend.box.background = element_rect(colour = "gray")) +
  scale_color_manual(values=c("#440154FF", "#6ece58", "#FDE725FF")) +
  guides(color = guide_legend(override.aes = list(size = 6)))
tree3

ggsave("/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/plot_tree3_2024.png", 
       tree3, width = 8, height = 8, dpi = 600)

figure_tree<- ggarrange(tree1, tree2, tree3,
                        labels=c("A","B","C"),
                        ncol=3,nrow=1)

figure_tree <- figure_tree + theme(plot.margin = margin(0.25, 0.25, 0.25, 0.25, "inches"))
figure_tree

ggsave("/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/plot_tree_2024.png", 
       figure_tree, width = 16, height = 10, dpi = 600)


#-------------------------------------------------------------------------------
#ANALISIS DE LA DIVERSIDAD DE LA COMUNIDAD MICROBIANA parte 1
#DIVERSIDAD ALFA
#Hacer el analisis de diversidad alfa (diversidad dentro de la muestra)
#Usar el objeto filtrado y normalizado (rarefraccion) para el analisis
#Este es el analisis univariado de la diversidad

#Revisar la correlacion entre las lecturas totales y la riqueza observada
#La diversidad aumenta a medida que se obtienen mas lecturas como es de esperarse
#hacerlo sobre el objeto original de phyloseq (Ma)
figure_obs_rich <- ggplot(data = data.frame("total_reads" =  phyloseq::sample_sums(Ma),
                         "observed" = phyloseq::estimate_richness(Ma, measures = "Observed")[, 1]),
       aes(x = total_reads, y = observed)) +
  geom_point() +
  geom_smooth(method="lm", se = TRUE) +
  labs(x = "\nTotal Reads", y = "Number of ASV observed\n") +
  theme_classic() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        axis.ticks.length = unit(.25, "cm")) +
  theme(plot.margin = margin(0.25,0.25,0.25,0.25, "inches")) +
  scale_x_continuous(labels=comma)
figure_obs_rich

ggsave("/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/plot_observed_richness_2024.png", 
       figure_obs_rich, width = 8, height = 6, dpi = 600)

figure_rare_rich <- ggarrange(figure_rare,figure_obs_rich,
                   labels=c("A","B"),
                   ncol=2,nrow=1)

figure_rare_rich2 <- figure_rare_rich + theme(plot.margin = margin(0.25,0.25,0.25,0.25, "inches"))
figure_rare_rich2

ggsave("/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/plot_rare_richness_2024_2.png", 
       figure_rare_rich2, width = 14, height = 6, dpi = 600)


#Determinar los indices de diversidad, riqueza y dominancia para cada factor y exportar la tabla
plot_richness(Ma_rarefied) 

adiv <- data.frame(
  "Observed" = phyloseq::estimate_richness(Ma_rarefied, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(Ma_rarefied, measures = "Shannon"),
  "Simpson" = phyloseq::estimate_richness(Ma_rarefied, measures = "Simpson"),
  "Chao1" = phyloseq::estimate_richness(Ma_rarefied, measures = "Chao1"),
  "InvSimpson" = phyloseq::estimate_richness(Ma_rarefied, measures = "InvSimpson"),
  "Source" = phyloseq::sample_data(Ma_rarefied)$Source, 
  "Site" = phyloseq::sample_data(Ma_rarefied)$Site,
  "Season" = phyloseq::sample_data(Ma_rarefied)$Season,
  "Season_Year" = phyloseq::sample_data(Ma_rarefied)$Season_Year,
  "Algae_Touching" = phyloseq::sample_data(Ma_rarefied)$Algae_Touching,
  "Health_Status" = phyloseq::sample_data(Ma_rarefied)$Health_Status)

head(adiv)

write.csv(adiv, file="/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/alpha_div.csv")

#Graficar los indices observado, shannon y dominancia por sitio y fuente

adiv %>%
  gather(key = metric, value = value, c("Observed", "Shannon", "Simpson")) %>%
  mutate(metric = factor(metric, levels = c("Observed", "Shannon", "Simpson"))) %>%
  ggplot(aes(x = Source, y = value)) +
  geom_boxplot(aes(color = Season)) +
  geom_jitter(aes(color = Season), height = 0, width = .2) +
  labs(x = "", y = "") +
  facet_wrap(~ metric, scales = "free")

rich1 <- adiv %>%
  gather(key = metric, value = value, c("Observed", "Shannon", "Simpson")) %>%
  mutate(metric = factor(metric, levels = c("Observed", "Shannon", "Simpson"))) %>%
  ggplot(aes(x = Source, y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = Site), height = 0, width = .2) +
  labs(x = "", y = "Alpha Diversity Measure") +
  facet_wrap(~ metric, scales = "free") +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        axis.ticks.length = unit(.25, "cm"),
        strip.text = element_text(size = 14)) +
  theme(plot.margin = margin(2, 2, 2, 2, "cm")) +
  scale_color_manual(values=c("#238A8DFF", "#FDE725FF")) +
  theme(plot.margin = margin(2, 2, 2, 2, "cm"))
rich1

rich2 <- adiv %>%
  gather(key = metric, value = value, c("Observed", "Shannon", "Simpson")) %>%
  mutate(metric = factor(metric, levels = c("Observed", "Shannon", "Simpson"))) %>%
  ggplot(aes(x = Source, y = value)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(aes(color = Season_Year), height = 0, width = .2) +
  labs(x = "", y = "Alpha Diversity Measure") +
  facet_wrap(~ metric, scales = "free") +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        axis.ticks.length = unit(.25, "cm"),
        strip.text = element_text(size = 14)) +
  theme(plot.margin = margin(2, 2, 2, 2, "cm")) +
  scale_color_manual(values=c("#238A8DFF", "#FDE725FF", "#404788FF", "#55C667FF")) +
  theme(plot.margin = margin(2, 2, 2, 2, "cm"))
rich2

ggsave("/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/boxplot_alpha_diversity_site_2024.png", 
       rich1, width = 10, height = 8, dpi = 600)

adiv %>%
  group_by(Source) %>%
  dplyr::summarise(median_observed = median(Observed),
                   median_shannon = median(Shannon),
                   median_simpson = median(Simpson))

#Hacer analisis estadistico (univariado) no parametrico para la diversidad, obvervada y shannon con base en cada factor
#Homegeneidad de varianzas con test de Levene y diferencias con test de Wilcox
library(car)

#Source
leveneTest(y = adiv$Shannon, group = adiv$Source, center = "median")
kruskal.test(Shannon ~ Source, data = adiv)
kruskal.test(Observed ~ Source, data = adiv)
TukeyHSD_Shannon <- TukeyHSD(aov(Shannon ~ Source, data =  adiv))
TukeyHSD_Shannon_df <- data.frame(TukeyHSD_Shannon$Source)
TukeyHSD_Shannon_df$measure = "Shannon"
TukeyHSD_Shannon_df$shapiro_test_pval = (shapiro.test(residuals(aov(Shannon ~ Source, data =  adiv))))$p.value
TukeyHSD_Shannon_df
TukeyHSD_Observed <- TukeyHSD(aov(Observed ~ Source, data =  adiv))
TukeyHSD_Observed_df <- data.frame(TukeyHSD_Observed$Source)
TukeyHSD_Observed_df$measure = "Observed"
TukeyHSD_Observed_df$shapiro_test_pval = (shapiro.test(residuals(aov(Observed ~ Source, data =  adiv))))$p.value
TukeyHSD_Observed_df

#Sitio
leveneTest(y = adiv$Shannon, group = adiv$Site, center = "median")
kruskal.test(Shannon ~ Site, data = adiv)
kruskal.test(Observed ~ Site, data = adiv)
TukeyHSD_Shannon <- TukeyHSD(aov(Shannon ~ Site, data =  adiv))
TukeyHSD_Shannon_df <- data.frame(TukeyHSD_Shannon$Site)
TukeyHSD_Shannon_df$measure = "Shannon"
TukeyHSD_Shannon_df$shapiro_test_pval = (shapiro.test(residuals(aov(Shannon ~ Site, data =  adiv))))$p.value
TukeyHSD_Shannon_df
TukeyHSD_Observed <- TukeyHSD(aov(Observed ~ Site, data =  adiv))
TukeyHSD_Observed_df <- data.frame(TukeyHSD_Observed$Site)
TukeyHSD_Observed_df$measure = "Observed"
TukeyHSD_Observed_df$shapiro_test_pval = (shapiro.test(residuals(aov(Observed ~ Site, data =  adiv))))$p.value
TukeyHSD_Observed_df

#Estacion climatica
leveneTest(y = adiv$Shannon, group = adiv$Season_Year, center = "median")
kruskal.test(Shannon ~ Season_Year, data = adiv)
kruskal.test(Observed ~ Season_Year, data = adiv)
TukeyHSD_Shannon <- TukeyHSD(aov(Shannon ~ Season_Year, data =  adiv))
TukeyHSD_Shannon_df <- data.frame(TukeyHSD_Shannon$Season_Year)
TukeyHSD_Shannon_df$measure = "Shannon"
TukeyHSD_Shannon_df$shapiro_test_pval = (shapiro.test(residuals(aov(Shannon ~ Season_Year, data =  adiv))))$p.value
TukeyHSD_Shannon_df
TukeyHSD_Observed <- TukeyHSD(aov(Observed ~ Season_Year, data =  adiv))
TukeyHSD_Observed_df <- data.frame(TukeyHSD_Observed$Season_Year)
TukeyHSD_Observed_df$measure = "Observed"
TukeyHSD_Observed_df$shapiro_test_pval = (shapiro.test(residuals(aov(Observed ~ Season_Year, data =  adiv))))$p.value
TukeyHSD_Observed_df

#Contacto con alga
leveneTest(y = adiv$Shannon, group = adiv$Algae_Touching, center = "median")
kruskal.test(Shannon ~ Algae_Touching, data = adiv)
kruskal.test(Observed ~ Algae_Touching, data = adiv)
TukeyHSD_Shannon <- TukeyHSD(aov(Shannon ~ Algae_Touching, data =  adiv))
TukeyHSD_Shannon_df <- data.frame(TukeyHSD_Shannon$Algae_Touching)
TukeyHSD_Shannon_df$measure = "Shannon"
TukeyHSD_Shannon_df$shapiro_test_pval = (shapiro.test(residuals(aov(Shannon ~ Algae_Touching, data =  adiv))))$p.value
TukeyHSD_Shannon_df
TukeyHSD_Observed <- TukeyHSD(aov(Observed ~ Algae_Touching, data =  adiv))
TukeyHSD_Observed_df <- data.frame(TukeyHSD_Observed$Algae_Touching)
TukeyHSD_Observed_df$measure = "Observed"
TukeyHSD_Observed_df$shapiro_test_pval = (shapiro.test(residuals(aov(Observed ~ Algae_Touching, data =  adiv))))$p.value
TukeyHSD_Observed_df

#Estado de salud
leveneTest(y = adiv$Shannon, group = adiv$Health_Status, center = "median")
kruskal.test(Shannon ~ Health_Status, data = adiv)
kruskal.test(Observed ~ Health_Status, data = adiv)
TukeyHSD_Shannon <- TukeyHSD(aov(Shannon ~ Health_Status, data =  adiv))
TukeyHSD_Shannon_df <- data.frame(TukeyHSD_Shannon$Health_Status)
TukeyHSD_Shannon_df$measure = "Shannon"
TukeyHSD_Shannon_df$shapiro_test_pval = (shapiro.test(residuals(aov(Shannon ~ Health_Status, data =  adiv))))$p.value
TukeyHSD_Shannon_df
TukeyHSD_Observed <- TukeyHSD(aov(Observed ~ Health_Status, data =  adiv))
TukeyHSD_Observed_df <- data.frame(TukeyHSD_Observed$Health_Status)
TukeyHSD_Observed_df$measure = "Observed"
TukeyHSD_Observed_df$shapiro_test_pval = (shapiro.test(residuals(aov(Observed ~ Health_Status, data =  adiv))))$p.value
TukeyHSD_Observed_df

#Otros plots y analisis estadisticos de la diversidad alfa
rich = estimate_richness(Ma_rarefied, measures = c("Observed", "Shannon"))
wilcox.shannon <- pairwise.wilcox.test(rich$Shannon, 
                                        sample_data(Ma_rarefied)$Site, 
                                        p.adjust.method = "BH")
tab.observed <- wilcox.shannon$p.value %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "group1") %>%
  gather(key="group2", value="p.adj", -group1) %>%
  na.omit()
tab.observed

wilcox.shannon <- pairwise.wilcox.test(rich$Shannon, 
                                       sample_data(Ma_rarefied)$Season, 
                                       p.adjust.method = "BH")
tab.observed <- wilcox.shannon$p.value %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "group1") %>%
  gather(key="group2", value="p.adj", -group1) %>%
  na.omit()
tab.observed

p <- plot_richness(Ma_rarefied, 
                   x="sample", 
                   color="Site", 
                   measures=c("Observed","Shannon", "Chao1"), 
                   nrow = 1)
print(p)

p$data %>% head()

ggplot(p$data,aes(Site,value,colour=Site)) +
  facet_grid(variable ~ Site, drop=T,scale="free",space="fixed") +
  geom_boxplot(outlier.colour = NA,alpha=1)

diver1 <- ggplot(p$data,aes(Source,value, fill=Site)) +
  facet_grid(variable ~ Source, drop=T,scale="free",space="fixed") +
  geom_boxplot(outlier.colour = NA,alpha=0.4, 
               position = position_dodge(width=0.9)) + 
  geom_point(size=2,position=position_jitterdodge(dodge.width=0.9)) +
  ylab("Diversity index") + xlab(NULL) + 
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        axis.ticks.length = unit(.25, "cm"),
        strip.text = element_text(size = 14)) +
  theme(plot.margin = margin(2, 2, 2, 2, "cm")) +
  scale_color_manual(values=c("#238A8DFF", "#FDE725FF")) +
  scale_fill_manual(values=c("#238A8DFF", "#FDE725FF")) +
  theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), 
                           "inches"))
diver1

diver2 <- ggplot(p$data,aes(Source,value,fill=Season)) +
  facet_grid(variable ~ Source, drop=T,scale="free",space="fixed") +
  geom_boxplot(outlier.colour = NA,alpha=0.4, 
               position = position_dodge(width=0.9)) + 
  geom_point(size=2,position=position_jitterdodge(dodge.width=0.9)) +
  ylab("Diversity index")  + xlab(NULL) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        axis.ticks.length = unit(.25, "cm"),
        strip.text = element_text(size = 14)) +
  theme(plot.margin = margin(2, 2, 2, 2, "cm")) +
  scale_color_manual(values=c("#238A8DFF", "#FDE725FF")) +
  scale_fill_manual(values=c("#238A8DFF", "#FDE725FF")) +
  theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), 
                           "inches"))
diver2

figure_diversity <- ggarrange(diver1,diver2,
                   labels=c("A","B"),
                   ncol=1,nrow=2)

figure_diversity_index <- figure_diversity + theme(plot.margin = margin(1.5, 1.5, 1.5, 1.5, "cm"))
figure_diversity_index 

#Hacer otros analisis y graficas con la significancia
# Generamos un objeto `phyloseq` sin taxa que sume 0 reads
Ma2 <- prune_taxa(taxa_sums(Ma_rarefied) > 0, Ma_rarefied)
Ma2
# Calculamos los índices de diversidad
#(Asegurarse de cargar los paquetes global para que funcione)
tab <- diversity(Ma2, index = "all")
# Y finalmente visualizamos la tabla de resultados
head(tab) %>%
  kable(format = "html", col.names = colnames(tab), digits = 2) %>%
  kable_styling() %>%
  kableExtra::scroll_box(width = "100%", height = "310px")

Ma2_meta <- meta(Ma2)

head(Ma2_meta) %>%
  kable(format = "html", col.names = colnames(Ma2_meta), digits = 2) %>%
  kable_styling() %>%
  kableExtra::scroll_box(width = "100%", height = "500px")

Ma2_meta$Shannon <- tab$diversity_shannon
Ma2_meta$Simpson <- tab$diversity_gini_simpson

Shannon <- tab$shannon
Simpson <- tab$gini_simpson

metadf <- Ma2_meta
metadf$Source.fact <- as.factor(metadf$Source)
metadf$Health_Status.fact <- as.factor(metadf$Health_Status)
metadf$Algae_Touching.fact <- as.factor(metadf$Algae_Touching)
metadf$Season.fact <- as.factor(metadf$Season)
metadf$Site.fact <- as.factor(metadf$Site)
metadf$Season_Year.fact <- as.factor(metadf$Season_Year)

# Obtenemos las variables desde nuestro objeto `phyloseq`
sources <- levels(metadf$Source.fact)
status <- levels(metadf$Health_Status.fact)
algae <- levels(metadf$Algae_Touching.fact)
season <- levels(metadf$Season.fact)
site <- levels(metadf$Site.fact)
season.year <- levels(metadf$Season_Year.fact)
# Creamos una lista de lo que queremos comparar
pares.sources <- combn(seq_along(sources), 2, simplify = FALSE, FUN = function(i)sources[i])
pares.status <- combn(seq_along(status), 2, simplify = FALSE, FUN = function(i)status[i])
pares.algae <- combn(seq_along(algae), 2, simplify = FALSE, FUN = function(i)algae[i])
pares.season <- combn(seq_along(season), 2, simplify = FALSE, FUN = function(i)season[i])
pares.site <- combn(seq_along(site), 2, simplify = FALSE, FUN = function(i)site[i])
pares.season.year <- combn(seq_along(season.year), 2, simplify = FALSE, FUN = function(i)season.year[i])
# Imprimimos en pantalla el resultado
print(pares.sources)
print(pares.status)
print(pares.algae)
print(pares.season)
print(pares.site)
print(pares.season.year)

boxplot1 <- ggboxplot(Ma2_meta, x = "Source", y = "Shannon",
               add = "jitter", fill = "Source", palette = c("#404788FF", "#FDE725FF", "#238A8DFF"), alpha = 0.5) +
  labs(x = "\nSource", y = "Shannon\n") +
  theme(axis.title = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(.25, "cm"),
        axis.text.x = element_text(angle = 90)) +
  theme(legend.position = "none") +
  theme(plot.margin = margin(0.25,0.25,0.25,0.25, "inches"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  ylim(3.2, 5.0)

boxplot2 <- ggboxplot(Ma2_meta, x = "Health_Status", y = "Shannon",
                      add = "jitter", fill = "Health_Status", palette = c("#404788FF", "#238A8DFF", "#55C667FF"), alpha = 0.5) +
  labs(x = "\nHealth Status", y = "Shannon\n") +
  theme(axis.title = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(.25, "cm"),
        axis.text.x = element_text(angle = 90)) +
  theme(legend.position = "none") +
  theme(plot.margin = margin(0.25,0.25,0.25,0.25, "inches"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  ylim(3.2, 5.0)

boxplot3 <- ggboxplot(Ma2_meta, x = "Algae_Touching", y = "Shannon", 
                    add = "jitter", fill = "Algae_Touching", palette = c("#404788FF", "#238A8DFF", "#55C667FF"), alpha = 0.5) +
  labs(x = "\nAlgae Touching", y = "Shannon\n") +
  theme(axis.title = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(.25, "cm"),
        axis.text.x = element_text(angle = 90)) +
  theme(legend.position = "none") +
  theme(plot.margin = margin(0.25,0.25,0.25,0.25, "inches"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  ylim(3.2, 5.0)

boxplot4 <- ggboxplot(Ma2_meta, x = "Season_Year", y = "Shannon", 
                    add = "jitter", fill = "Season_Year", palette = c("#404788FF", "#238A8DFF", "#55C667FF", "#FDE725FF"), alpha = 0.5) +
  labs(x = "\nSeason", y = "Shannon\n") +
  theme(axis.title = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(.25, "cm"),
        axis.text.x = element_text(angle = 90)) +
  theme(legend.position = "none") +
  theme(plot.margin = margin(0.25,0.25,0.25,0.25, "inches"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  ylim(3.2, 5.5)

boxplot5 <- ggboxplot(Ma2_meta, x = "Site", y = "Shannon", 
                    add = "jitter", fill = "Site", palette = c("#404788FF", "#238A8DFF", "#55C667FF"), alpha = 0.5) +
  labs(x = "\nSite", y = "Shannon\n") +
  theme(axis.title = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(.25, "cm"),
        axis.text.x = element_text(angle = 90)) +
  theme(legend.position = "none") +
  theme(plot.margin = margin(0.25,0.25,0.25,0.25, "inches"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  ylim(3.2, 5.0)

boxplot6 <- ggboxplot(Ma2_meta, x = "Source", y = "Simpson",
                      add = "jitter", fill = "Source", palette = c("#404788FF", "#FDE725FF", "#238A8DFF"), alpha = 0.5) +
  labs(x = "\nSource", y = "Simpson\n") +
  theme(axis.title = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(.25, "cm"),
        axis.text.x = element_text(angle = 90)) +
  theme(legend.position = "none") +
  theme(plot.margin = margin(0.25,0.25,0.25,0.25, "inches"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  ylim(0.95, 1)

boxplot7 <- ggboxplot(Ma2_meta, x = "Season_Year", y = "Simpson", 
                      add = "jitter", fill = "Season_Year", palette = c("#404788FF", "#238A8DFF", "#55C667FF", "#FDE725FF"), alpha = 0.5) +
  labs(x = "\nSeason", y = "Simpson\n") +
  theme(axis.title = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(.25, "cm"),
        axis.text.x = element_text(angle = 90)) +
  theme(legend.position = "none") +
  theme(plot.margin = margin(0.25,0.25,0.25,0.25, "inches"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  ylim(0.95, 1.01)

boxplot8 <- ggboxplot(Ma2_meta, x = "Site", y = "Simpson", 
                      add = "jitter", fill = "Site", palette = c("#404788FF", "#238A8DFF", "#55C667FF"), alpha = 0.5) +
  labs(x = "\nSite", y = "Simpson\n") +
  theme(axis.title = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(.25, "cm"),
        axis.text.x = element_text(angle = 90)) +
  theme(legend.position = "none") +
  theme(plot.margin = margin(0.25,0.25,0.25,0.25, "inches"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        axis.line = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  ylim(0.95, 1)

b1 <- boxplot1 + stat_compare_means(comparisons = pares.sources)
print(b1)

b2 <- boxplot2 + stat_compare_means(comparisons = pares.status)
print(b2)

b3 <- boxplot3 + stat_compare_means(comparisons = pares.algae)
print(b3)

b4 <- boxplot4 + stat_compare_means(comparisons = pares.season.year)
print(b4)

b5 <- boxplot5 + stat_compare_means(comparisons = pares.site)
print(b5)

b6 <- boxplot6 + stat_compare_means(comparisons = pares.sources)
print(b6)

b7 <- boxplot7 + stat_compare_means(comparisons = pares.season.year)
print(b7)

b8 <- boxplot8 + stat_compare_means(comparisons = pares.site)
print(b8)

figure_diversities <- ggarrange(b5,b4,b1,b8,b7,b6,
                   labels=c("A","B","C","D","E","F"),
                   ncol=3,nrow=2)

figure_diversities_index <- figure_diversities + theme(plot.margin = margin(0.25,0.25,0.25,0.25, "inches"))
figure_diversities_index

ggsave("/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/plot_diversities_2024.png", 
       figure_diversities_index, width = 18, height = 12, dpi = 600)


#-------------------------------------------------------------------------------
#ANALISIS DE LA DIVERSIDAD DE LA COMUNIDAD MICROBIANA parte 2
#DIVERSIDAD BETA
#Hacer el analisis de diversidad beta (diversidad entre muestras)
#Usar el objeto filtrado y normalizado (rarefraccion) para el analisis
#Este es el analisis multivariado de la diversidad
library(vegan)
library(phyloseq)
library(ggplot2)
library(dplyr)

#DENDOGRAMA
library(dendextend)
#Clustering y dendogramas de las muestras
ps_rel_otu <- data.frame(phyloseq::otu_table(Ma_rarefied))
ps_rel_otu <- t(ps_rel_otu)
bc_dist <- vegan::vegdist(ps_rel_otu, method = "bray")
#Save as dendrogram
ward <- as.dendrogram(hclust(bc_dist, method = "ward.D2"))
#Provide color codes
meta <- data.frame(phyloseq::sample_data(Ma_rarefied))
colorCode <- c(Seawater = "#FDE725FF", Mucus = "#6ece58", Tissue = "#440154FF")
labels_colors(ward) <- colorCode[meta$Source][order.dendrogram(ward)]
#Plot
dendogram <- plot(ward)

png(filename = "Dendogram.png", width = 15, height = 11, res = 1200, units = "in")
dendogram <- plot(ward)
print(dendogram)
dev.off()

#ANALISIS Y PCoA
#Determinar la disimilutud como medidad de distancia
Ma_rarefied %>% transform_sample_counts(function(x) x/sum(x)) %>%
  otu_table() %>%
  t() %>%
  sqrt() %>%
  as.data.frame() %>%
  vegdist(binary=F, method = "bray") -> dist

#Correr analisis de ordenacion PCoA con la distancia
#Bray Curtis
ord <- ordinate(Ma_rarefied,"PCoA",dist)
ord$vectors

pcoa <- plot_ordination(Ma_rarefied, 
                ord,
                color = "Site", 
                shape="Site", 
                label= "Sample_ID") + 
  stat_ellipse(geom = "polygon", aes(group = Source, label = Source, fill = Source), 
               alpha = 0.5, size = 1, linetype = 2) +
  stat_ellipse(aes(group = Site), linetype = 2, label = Year) +
  geom_vline(xintercept = c(0), color = "grey", linetype = 1) +
  geom_hline(yintercept = c(0), color = "grey", linetype = 1) +
  labs(x = "PCo1 [62.0%]", y = "PCo2 [11.1%]") +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        plot.margin = margin(0.25, 0.25, 0.25, 0.25, "inches"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  geom_point(aes(size=Shannon)) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_color_manual(values=c("#238A8DFF", "#FDE725FF")) +
  scale_fill_manual(values=c("#404788FF", "#FDE725FF", "#238A8DFF"))

pcoa1 <- pcoa +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        axis.ticks.length = unit(.25, "cm"),
        strip.text = element_text(size = 14),
        legend.key.size = unit(1, "cm"),
        plot.title = element_text(size = 12)) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  ggtitle("ADONIS: p=0.001, R2=0.5865 \nNPMANOVA: p=0.001") 

pcoa1

ggsave("/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/plot_PCoA_beta_2024.png", 
       pcoa1, width = 10, height = 10, dpi = 600)

#Hacer prueba estadistica de significancia con PERMANOVA (Adonis)
#(actualizar o recargar los paquetes vegan y phyloseq por si no funciona)
adonis(dist ~ get_variable(Ma_rarefied, "Site"), permutations = 1000)$aov.tab
adonis(dist ~ get_variable(Ma_rarefied, "Source"), permutations = 1000)$aov.tab
adonis(dist ~ get_variable(Ma_rarefied, "Season"), permutations = 1000)$aov.tab

#Pair-wise PERMANOVA para compartimento (fuente)
cbn <- combn(x=unique(Metadata2$Source), m = 2)
p <- c()

for(i in 1:ncol(cbn)){
  Ma_subs <- subset_samples(Ma_rarefied, Source %in% cbn[,i])
  metadata_sub <- data.frame(sample_data(Ma_subs))
  permanova_pairwise <- adonis(phyloseq::distance(Ma_subs, method = "bray") ~ Source, 
                               data = metadata_sub)
  p <- c(p, permanova_pairwise$aov.tab$`Pr(>F)`[1])
}

p.adj <- p.adjust(p, method = "BH")
p.table <- cbind.data.frame(t(cbn), p=p, p.adj=p.adj)
p.table

#Pair-wise PERMANOVA para sitio
cbn <- combn(x=unique(Metadata2$Site), m = 2)
p <- c()

for(i in 1:ncol(cbn)){
  Ma_subs <- subset_samples(Ma_rarefied, Site %in% cbn[,i])
  metadata_sub <- data.frame(sample_data(Ma_subs))
  permanova_pairwise <- adonis(phyloseq::distance(Ma_subs, method = "bray") ~ Site, 
                               data = metadata_sub)
  p <- c(p, permanova_pairwise$aov.tab$`Pr(>F)`[1])
}

p.adj <- p.adjust(p, method = "BH")
p.table <- cbind.data.frame(t(cbn), p=p, p.adj=p.adj)
p.table

#Pair-wise PERMANOVA para estacion climatica
cbn <- combn(x=unique(Metadata2$Season), m = 2)
p <- c()

for(i in 1:ncol(cbn)){
  Ma_subs <- subset_samples(Ma_rarefied, Season %in% cbn[,i])
  metadata_sub <- data.frame(sample_data(Ma_subs))
  permanova_pairwise <- adonis(phyloseq::distance(Ma_subs, method = "bray") ~ Season, 
                               data = metadata_sub)
  p <- c(p, permanova_pairwise$aov.tab$`Pr(>F)`[1])
}

p.adj <- p.adjust(p, method = "BH")
p.table <- cbind.data.frame(t(cbn), p=p, p.adj=p.adj)
p.table

#Dispersion
disper1 <- boxplot(betadisper(dist, 
                   get_variable(Ma_rarefied, "Site")),las=2, 
        main=paste0("Multivariate Dispersion Test Bray-Curtis "," pvalue = ", 
                    permutest(betadisper(dist, get_variable(Ma_rarefied, "Site")))$tab$`Pr(>F)`[1]))
disper1

disper2 <- boxplot(betadisper(dist, 
                   get_variable(Ma_rarefied, "Season_Year")),las=2, 
        main=paste0("Multivariate Dispersion Test Bray-Curtis "," pvalue = ", 
                    permutest(betadisper(dist, get_variable(Ma_rarefied, "Season_Year")))$tab$`Pr(>F)`[1]))
disper2

disper3 <- boxplot(betadisper(dist, 
                   get_variable(Ma_rarefied, "Source")),las=2, 
        main=paste0("Multivariate Dispersion Test Bray-Curtis "," pvalue = ", 
                    permutest(betadisper(dist, get_variable(Ma_rarefied, "Source")))$tab$`Pr(>F)`[1]))
disper3

#ANOSIM
plot(anosim(dist, get_variable(Ma_rarefied, "Site"))
     ,main="ANOSIM Bray-Curtis "
     ,las=2)
plot(anosim(dist, get_variable(Ma_rarefied, "Season_Year"))
     ,main="ANOSIM Bray-Curtis "
     ,las=2)
plot(anosim(dist, get_variable(Ma_rarefied, "Source"))
     ,main="ANOSIM Bray-Curtis "
     ,las=2)

#Graficar el NMDS con base en la distancia de jaccard
dist = phyloseq::distance(Ma_rarefied, method="jaccard", binary = TRUE)
ordination = ordinate(Ma_rarefied, method="NMDS", distance=dist)
 
pcoa2 <- plot_ordination(Ma_rarefied, 
                        ordination,
                        color = "Source", 
                        shape="Season_Year", 
                        label= "Sample_ID") + 
  stat_ellipse(geom = "polygon", aes(group = Season_Year, label = Season_Year, fill = Season_Year), 
               alpha = 0.5, size = 1, linetype = 2) +
  stat_ellipse(aes(group = Source), linetype = 2) +
  labs(x = "NMDS1", y = "NMDS2") +
  geom_vline(xintercept = c(0), color = "grey", linetype = 1) +
  geom_hline(yintercept = c(0), color = "grey", linetype = 1) +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        plot.margin = margin(0.25, 0.25, 0.25, 0.25, "inches"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  geom_point(aes(size=Shannon)) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_color_manual(values=c("#404788FF", "#FDE725FF", "#238A8DFF")) +
  scale_fill_manual(values=c("#440154FF", "#6ece58", "#B8DE29FF", "#482878")) +
  ggtitle("ANOSIM: p=0.001, R=0.331 \nSTRESS=0.065") +
  ylim(-1.6, 1.2)

pcoa3 <- pcoa2 +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        axis.ticks.length = unit(.25, "cm"),
        strip.text = element_text(size = 14),
        legend.key.size = unit(1, "cm"),
        plot.title = element_text(size = 12)) +
  guides(color = guide_legend(override.aes = list(size = 5)))

pcoa3

ggsave("/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/plot_NMDS_beta_2024.png", 
       pcoa3, width = 10, height = 10, dpi = 600)

figure_ord<- ggarrange(pcoa1, pcoa3,
                   labels=c("A","B"),
                   ncol=1,nrow=2)

figure_beta_ord <- figure_ord + theme(plot.margin = margin(0.25, 0.25, 0.25, 0.25, "inches"))
figure_beta_ord

ggsave("/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/plot_PCOA+NMDS_beta_2024.png", 
       figure_beta_ord, width = 8, height = 14, dpi = 600)

## Simper analysis, extract abundance matrix from a phyloseq object
Ma1_OTUs = as(otu_table(Ma1), "matrix")

# transpose so we have the OTUs as columns
if(taxa_are_rows(Ma1)){Ma1_OTUs <- t(Ma1_OTUs)}

# Coerce the object to a data.frame
Ma1_OTUs = as.data.frame(Ma1_OTUs)

# running the simper analysis on the dataframe and the variable of interest "time"
Ma1_simper <- simper(Ma1_OTUs, Metadata2$Source, permutations = 100)

# printing the top OTUs
print(Ma1_simper)


#-------------------------------------------------------------------------------
#ABUNDANCIA RELATIVA Y VISUALIZACION DE LA COMPOSICION DE LA COMUNIDAD MICROBIANA

#Visualizar y hacer analisis de la diversidad con base en la abundancia relativa
#Combinar todo en el nivel filo y hallar la abundancia relativa
#Hacer bar plot y box plot con la abundancia relativa para el filo y el genero

Ma1_phylum <- tax_glom(Ma1, "Phylum", NArm = TRUE)
Ma1_phylum_relabun <- transform_sample_counts(Ma1_phylum, function(OTU) OTU/sum(OTU) * 100)
taxa_abundance_table_phylum <- psmelt(Ma1_phylum_relabun)
StackedBarPlot_phylum <- taxa_abundance_table_phylum %>% 
  ggplot(aes(x =Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  labs(x = "",
       y = "Relative Abundance",
       title = "Phylum Relative Abundance") +
  facet_grid(~ Site + Season, scales = "free") +
  theme(
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12))

StackedBarPlot_phylum

ggsave("/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/plot_box_abundance_phylum_site_2024.png", 
       BoxPlot_phylum, width = 17, height = 5, dpi = 600)

Ma1_genus <- tax_glom(Ma1, "Genus", NArm = TRUE)
Ma1_genus_relabun <- transform_sample_counts(Ma1_genus, function(OTU) OTU/sum(OTU) * 100)
taxa_abundance_table_genus <- psmelt(Ma1_genus_relabun)
StackedBarPlot_genus <- taxa_abundance_table_genus %>% 
  ggplot(aes(x =Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity") +
  labs(x = "",
       y = "Relative Abundance",
       title = "Genus Relative Abundance") +
  facet_grid(~Source, scales = "free") +
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 10),
    strip.text = element_text(size = 12))

StackedBarPlot_genus

#BARRAS APILADAS, BOXPLOTS Y % DE ABUNDANCIA RELATIVA (TAXA RARA Y POR GRUPOS)
#Otra forma de hacer el analisis de la abundancia relativa basada en barras apiladas
#GENUS
Ma1.rel = transform_sample_counts(Ma1, function(x) x/sum(x)*100)
# agglomerate taxa
glom <- tax_glom(Ma1.rel, taxrank = 'Genus', NArm = FALSE)
Ma1.melt <- psmelt(glom)
# change to character for easy-adjusted level
Ma1.melt$Genus <- as.character(Ma1.melt$Genus)

Ma1.melt <- Ma1.melt %>%
  group_by(Source, Site, Season, Algae_Touching, Health_Status, Genus) %>%
  mutate(median=median(Abundance))
# select group median > 1
keep <- unique(Ma1.melt$Genus[Ma1.melt$median > 2.5])
Ma1.melt$Genus[!(Ma1.melt$Genus %in% keep)] <- "< 2.5%"
#to get the same rows together
Ma1.melt_sum <- Ma1.melt %>%
  group_by(Sample, Source, Site, Season, Algae_Touching, Health_Status, Genus) %>%
  summarise(Abundance=sum(Abundance))



bar1 <- ggplot(Ma1.melt_sum, aes(x = Sample, y = Abundance, fill = Genus)) + 
  geom_bar(stat = "identity", aes(fill=Genus), alpha = 0.9) + 
  labs(x="", y="Relative abundance (%)") +
  facet_wrap(~Site + Source + Season, scales= "free_x", nrow=1) +
  theme_classic() + 
  theme(strip.background = element_blank(), 
        axis.text.x.bottom = element_text(angle = -90)) +
  theme(axis.title=element_text(size=12), 
        axis.text=element_text(size=10), 
        legend.text=element_text(size=10),
        legend.title = element_text(size = 12),
        axis.text.x.bottom = element_blank(),
        strip.text = element_text(size = 12),
        axis.ticks.length = unit(.25, "cm")) +
  scale_fill_viridis_d()
bar1

ggsave("/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/plot_bar_relative_genus_2024.png", 
       bar1, width = 17, height = 5, dpi = 600)

write.csv(Ma1.melt_sum, 
          file="/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/relative_abundance_genus_table.csv")

BoxPlot_genus1 <- Ma1.melt_sum %>% 
  ggplot(aes(x = Genus, y = Abundance, fill = Genus)) +
  geom_boxplot() +
  labs(x = "",
       y = "Relative Abundance (%)") +
  facet_grid(~Source, scales = "free") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.text = element_text(color = "black"),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text = element_text(size = 12),
        axis.ticks.length = unit(.25, "cm")) +
  scale_fill_viridis_d()
BoxPlot_genus1

BoxPlot_genus2 <- Ma1.melt_sum %>% 
  ggplot(aes(x = Genus, y = Abundance, fill = Genus)) +
  geom_boxplot() +
  labs(x = "",
       y = "Relative Abundance (%)") +
  facet_grid(~Health_Status + Algae_Touching, scales = "free") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.text = element_text(color = "black"),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        strip.text = element_text(size = 12),
        axis.ticks.length = unit(.25, "cm")) +
  scale_fill_viridis_d()
BoxPlot_genus2

ggsave("/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/plot_box_abundance_genus_source_2024.png", 
       BoxPlot_genus1, width = 12, height = 6, dpi = 600)

ggsave("/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/plot_box_abundance_genus_health_2024.png", 
       BoxPlot_genus2, width = 14, height = 6, dpi = 600)


#PHYLUM
glom <- tax_glom(Ma1.rel, taxrank = 'Phylum', NArm = FALSE)
Ma1.melt <- psmelt(glom)
# change to character for easy-adjusted level
Ma1.melt$Phylum <- as.character(Ma1.melt$Phylum)

Ma1.melt <- Ma1.melt %>%
  group_by(Source, Site, Season, Algae_Touching, Health_Status, Phylum) %>%
  mutate(median=median(Abundance))

#to get the same rows together
Ma1.melt_sum <- Ma1.melt %>%
  group_by(Sample, Source, Site, Season, Algae_Touching, Health_Status, Phylum) %>%
  summarise(Abundance=sum(Abundance))

bar2 <- ggplot(Ma1.melt_sum, aes(x = Sample, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity", aes(fill=Phylum), alpha = 0.9) + 
  labs(x="", y="Relative abundance (%)") +
  facet_wrap(~Site + Source + Season, scales= "free_x", nrow=1) +
  theme_classic() + 
  theme(strip.background = element_blank(), 
        axis.text.x.bottom = element_text(angle = -90)) +
  theme(axis.title=element_text(size=14), 
        axis.text=element_text(size=12, color = "black"), 
        legend.text=element_text(size=12),
        legend.title = element_text(size = 14),
        axis.text.x.bottom = element_blank(),
        strip.text = element_text(size = 12),
        axis.ticks.length = unit(.25, "cm")) +
  scale_fill_viridis_d()
bar2

ggsave("/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/plot_bar_relative_phylum_2024.png", 
       bar2, width = 17, height = 5, dpi = 600)

write.csv(Ma1.melt_sum, 
          file="/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/relative_abundance_phylum_table.csv")

BoxPlot_phylum <- Ma1.melt_sum %>% 
  ggplot(aes(x =Phylum, y = Abundance, fill = Phylum)) +
  geom_boxplot() +
  labs(x = "",
       y = "Relative Abundance (%)") +
  facet_grid(~ Source + Health_Status, scales = "free") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 10),
        strip.text = element_text(size = 12)) +
  scale_fill_viridis_d()
BoxPlot_phylum

ggsave("/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/plot_box_abundance_phylum_health_2024.png", 
       BoxPlot_phylum, width = 17, height = 5, dpi = 600)

#FAMILY
glom <- tax_glom(Ma1.rel, taxrank = 'Family', NArm = FALSE)
Ma1.melt <- psmelt(glom)
# change to character for easy-adjusted level
Ma1.melt$Family <- as.character(Ma1.melt$Family)

Ma1.melt <- Ma1.melt %>%
  group_by(Source, Site, Season, Algae_Touching, Health_Status, Family) %>%
  mutate(median=median(Abundance))

# select group median > 2.5%
keep <- unique(Ma1.melt$Family[Ma1.melt$median > 2.5])
Ma1.melt$Family[!(Ma1.melt$Family %in% keep)] <- "< 2.5%"

#to get the same rows together
Ma1.melt_sum <- Ma1.melt %>%
  group_by(Sample, Source, Site, Season, Algae_Touching, Health_Status, Family) %>%
  summarise(Abundance=sum(Abundance))

bar3 <- ggplot(Ma1.melt_sum, aes(x = Sample, y = Abundance, fill = Family)) + 
  geom_bar(stat = "identity", aes(fill=Family), alpha = 0.9) + 
  labs(x="", y="Relative abundance (%)") +
  facet_wrap(~Site + Source + Season, scales= "free_x", nrow=1) +
  theme_classic() + 
  theme(strip.background = element_blank(), 
        axis.text.x.bottom = element_text(angle = -90)) +
  theme(axis.title=element_text(size=14), 
        axis.text=element_text(size=12, color = "black"), 
        legend.text=element_text(size=12),
        legend.title = element_text(size = 14),
        axis.text.x.bottom = element_blank(),
        strip.text = element_text(size = 12),
        axis.ticks.length = unit(.25, "cm")) +
  scale_fill_viridis_d()
bar3

ggsave("/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/plot_bar_relative_family_2024.png", 
       bar3, width = 17, height = 5, dpi = 600)

write.csv(Ma1.melt_sum, 
          file="/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/relative_abundance_family_table.csv")

BoxPlot_family <- Ma1.melt_sum %>% 
  ggplot(aes(x =Family, y = Abundance, fill = Family)) +
  geom_boxplot() +
  labs(x = "",
       y = "Relative Abundance (%)") +
  facet_grid(~Site + Season + Source, scales = "free") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 10),
        strip.text = element_text(size = 12)) +
  scale_fill_viridis_d()
BoxPlot_family

ggsave("/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/plot_box_abundance_family_season+site_2024.png", 
       BoxPlot_family, width = 17, height = 5, dpi = 600)

#Unir los dos barplots de la abundancia relativa (phylum + genus)
figure_bar_abundance <- ggarrange(bar2, bar3,
                                  labels = c("A", "B"),
                                  ncol=1, nrow=2)

figure_bar_abundance <- figure_bar_abundance + theme(plot.margin = margin(0.25, 0.25, 0.25, 0.5, "inches"))
figure_bar_abundance

ggsave("/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/plot_bar_relative_abundance_2024.png", 
       figure_bar_abundance, width = 13, height = 9, dpi = 600)

#Otra forma de graficar usando burbujas con la abundancia relativa (%)
#FAMILY
Ma1.rel = transform_sample_counts(Ma1, function(x) x/sum(x)*100)
# agglomerate taxa
glom <- tax_glom(Ma1.rel, taxrank = 'Genus', NArm = FALSE)
Ma1.melt <- psmelt(glom)
# change to character for easy-adjusted level
Ma1.melt$Genus <- as.character(Ma1.melt$Genus)

Ma1.melt <- Ma1.melt %>%
  group_by(Source, Site, Season, Genus) %>%
  mutate(median=median(Abundance))
# select group median > 1
keep <- unique(Ma1.melt$Genus[Ma1.melt$median > 1.0])
Ma1.melt$Genus[!(Ma1.melt$Genus %in% keep)] <- "< 1.0%"
#to get the same rows together
Ma1.melt_sum <- Ma1.melt %>%
  group_by(Sample, Source, Site, Season, Genus) %>%
  summarise(Abundance=sum(Abundance))

bubble1 <- ggplot(Ma1.melt_sum, aes(x = Sample, y = Abundance, fill = Genus)) + 
  geom_point(aes(size = Abundance, fill = Genus), alpha = 0.75, shape = 21) + 
  scale_size_continuous(limits = c(0.000001, 100), range = c(1,17), breaks = c(1,10,50,75)) + 
  labs( x= "", y = "", size = "Relative Abundance (%)", fill = "Genus") +
  facet_grid(~Source, scales = "free") +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme_classic() +
  theme(
    axis.text = element_blank(),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    strip.text = element_text(size = 12),
    axis.text.x.bottom = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
    legend.position = "right") +
  scale_fill_viridis_d()
bubble1

ggsave("/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/plot_bubble_genus_abundance_2024.png", 
       bubble1, width = 17, height = 5, dpi = 600)


#-------------------------------------------------------------------------------
#ANALISIS DIFERENCIAL DE LA ABUNDANCIA (DESeq2)

#Hacer un analisis diferencial de la abundancia con DEseq2 para detectar los grupos bacterianos por condicion
#Factorizar para DESeq2
sample_data(Ma1)$Source <- as.factor(sample_data(Ma1)$Source)
Ma1.taxa <- tax_glom(Ma1, taxrank = 'Family', NArm = FALSE)
#Comparacion pairwise
Ma1.taxa.sub <- subset_samples(Ma1.taxa, Source %in% c("Mucus", "Tissue"))
#filtrar los caracteres con >90% de ceros
Ma1.taxa.pse.sub <- prune_taxa(rowSums(otu_table(Ma1.taxa.sub) == 0) < ncol(otu_table(Ma1.taxa.sub)) * 0.9, Ma1.taxa.sub)
Ma1_ds = phyloseq_to_deseq2(Ma1.taxa.pse.sub, ~ Source)
#Estimador alternativo "cada gen contiene una muestra con un cero"
ds <- estimateSizeFactors(Ma1_ds, type="poscounts")
ds = DESeq(ds, test="Wald", fitType="parametric")
alpha = 0.05 
res = results(ds, alpha=alpha)
res = res[order(res$padj, na.last=NA), ]
res
taxa_sig = rownames(res[1:13, ]) # seleccionar 13 con los p.adj valores mas bajos
Ma1.taxa.rel <- transform_sample_counts(Ma1, function(x) x/sum(x)*100)
Ma1.taxa.rel.sig <- prune_taxa(taxa_sig, Ma1.taxa.rel)
#Hacer la matriz
Ma1.taxa.rel.sig <- prune_samples(colnames(otu_table(Ma1.taxa.pse.sub)), Ma1.taxa.rel.sig)
matrix <- as.matrix(data.frame(otu_table(Ma1.taxa.rel.sig)))
rownames(matrix) <- as.character(tax_table(Ma1.taxa.rel.sig)[, "Family"])
metadata_sub <- data.frame(sample_data(Ma1.taxa.rel.sig))
colnames(matrix) = NULL

#Definir la anotacion
annotation_col = data.frame(
  Site = as.factor(metadata_sub$Site),
  Season = as.factor(metadata_sub$Season),
  `Algae Touching` =as.factor(metadata_sub$Algae_Touching),
  `Health Status` =as.factor(metadata_sub$Health_Status),
  Source = as.factor(metadata_sub$Source), 
  check.names = FALSE)
rownames(annotation_col) = rownames(metadata_sub)

annotation_row = data.frame(Phylum = as.factor(tax_table(Ma1.taxa.rel.sig)[, "Phylum"]))
rownames(annotation_row) = rownames(matrix)

#Definir la anotacion de los colores
phylum_col = RColorBrewer::brewer.pal(length(levels(annotation_row$Phylum)), "Paired")
names(phylum_col) = levels(annotation_row$Phylum)
ann_colors = list(
  Source = c(`Mucus` = "#440154FF", `Tissue` = "#FDE725FF"),
  `Site` = c(`Protected` = "#440154FF", `Urban` = "#FDE725FF"),
  `Season` = c(`Dry` = "#440154FF", `Rainy` = "#FDE725FF"),
  `Algae Touching` = c(`Algae_No` = "#440154FF", `Algae_Yes` = "#FDE725FF"),
  `Health Status` = c(`Healthy` = "#440154FF", `Stressed` = "#FDE725FF"),
  `Phylum` = c(`Bacteroidota` = "#440154FF", `Cyanobacteria` = "#482677ff", `Firmicutes` = "#404788ff",
               `Proteobacteria` = "#FDE725FF", `Spirochaetota` = "#dce319ff"))
library(circlize)
col_fun = colorRamp2(c(-4, 0, 4), c("yellow", "blue", "red"))
col_fun(seq(-3, 3))

Deseq_heatmap <- ComplexHeatmap::pheatmap(matrix, scale= "row", 
                         annotation_col = annotation_col, 
                         annotation_row = annotation_row, 
                         annotation_colors = ann_colors, 
                         col = col_fun, 
                         heatmap_legend_param = list( 
                           title = "Abundance",
                           legend_height = unit(4, "cm"),
                           title_position = "leftcenter-rot"))
Deseq_heatmap
                         
png(filename = "Deseq_Heatmap_family2.png", width = 10, height = 10, res = 1200, units = "in")
Heatmap <- draw(Deseq_heatmap, heatmap_legend_side = "left", annotation_legend_side = "right")
dev.off()

#Hacer un ANCOM-CB Analisis basado en la abundancia relativa para datos de microbioma
#el ANCOM-CB detecta que taxones cambiaron su abundancia relativa significativamente
#ancombc
library(tidyverse)
library(DT)
Ma1.taxa.sub <- subset_samples(Ma1.taxa, Source %in% c("Mucus", "Tissue", 
                                                       Site %in% c("Protected", "Urban", 
                                                                   Season %in% c("Dry", "Rainy", 
                                                                                 Algae_Touching %in% c("Yes", "No",
                                                                                                       Health_Status %in% c("Healthy", "Stressed"))))))
out = ancombc(data = Ma1.taxa.sub, assay_name = "counts", 
              tax_level = "Family", phyloseq = NULL, 
              formula = "Source + Site + Season + Algae_Touching + Health_Status", 
              p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, 
              group = "Source", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
              max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE,
              n_cl = 1, verbose = TRUE)
res <- out$res
res
res_global = out$res_global

tab_lfc = res$lfc
col_name = c("Taxon", "Intercept", "Mucus-Tissue", "Protected-Urban", "Dry-Rainy", "Yes-No", "Healthy-Stressed")
colnames(tab_lfc) = col_name
tab_lfc %>% 
  datatable(caption = "Log Fold Changes from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

tab_se = res$se
colnames(tab_se) = col_name
tab_se %>% 
  datatable(caption = "SEs from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

tab_w = res$W
colnames(tab_w) = col_name
tab_w %>% 
  datatable(caption = "Test Statistics from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

tab_p = res$p_val
colnames(tab_p) = col_name
tab_p %>% 
  datatable(caption = "P-values from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

tab_q = res$q
colnames(tab_q) = col_name
tab_q %>% 
  datatable(caption = "Adjusted p-values from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

tab_diff = res$diff_abn
colnames(tab_diff) = col_name
tab_diff %>% 
  datatable(caption = "Differentially Abundant Taxa from the Primary Result")

#que bacterias contribuyen mas a las diferencias en la composicion de las comunidades entre los grupos (beta-diversidad)
#hacer un analisis PERMANOVA
pseq <- Ma1

# Pick relative abundances (compositional) and sample metadata 
pseq.rel <- microbiome::transform(pseq, "compositional")
otu <- abundances(pseq.rel)
meta <- meta(pseq.rel)

# samples x species as input
library(vegan)
permanova <- adonis(t(otu) ~ Source,
                    data = meta, permutations=999, method = "bray")

# P-value
print(as.data.frame(permanova$aov.tab)["Source", "Pr(>F)"])
#revisar la homogeneidad de condicion
dist <- vegdist(t(otu))
anova(betadisper(dist, meta$Source))
permutest(betadisper(dist, meta$Source), pairwise = TRUE)
#taxa que contribuyen mas a las diferencias entre compartimentos o fuente (source)
coef <- coefficients(permanova)["Source1",]
top.coef <- coef[rev(order(abs(coef)))[1:40]]
par(mar = c(3, 6, 3, 6))
barplot(sort(top.coef), horiz = T, las = 1, main = "Top taxa")


#-------------------------------------------------------------------------------
#REDES DE CO-OCURRENCIA

#Redes de co-ocurrencia para inferir sobre potenciales interacciones ecológicas
#Estas interacciones pueden ser directas o indirectas y nos permiten comenzar a descifrar mecanismos ecológicos detrás de la composición de una comunidad microbiana
#Hacer red de distancia
plot_net(Ma1, type = "taxa", point_label = "Family", point_size = 10, point_alpha = 0.5, maxdist = 0.5, color = "Phylum", distance = "bray", laymeth = "auto") +
  scale_fill_viridis_d()

library(SpiecEasi)
se.mb.Ma1 <- spiec.easi(Ma1, method='mb', lambda.min.ratio=1e-2,
                         nlambda=20, icov.select.params=list(rep.num=50))
se.mb.Ma1
ig2.mb <- adj2igraph(getRefit(se.mb.Ma1),  vertex.attr=list(name=taxa_names(Ma1)))
ig2.mb
red_phylum <- plot_network(ig2.mb, Ma1, type='taxa', color="Phylum") + scale_color_viridis_d() +
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.margin = margin(0.25, 0.25, 0.25, 0.25, "inches"),
        legend.box.background = element_rect(colour = "gray"))
red_phylum
red_family <- plot_network(ig2.mb, Ma1, type='taxa', color="Family") + scale_color_viridis_d() +
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.margin = margin(0.25, 0.25, 0.25, 0.25, "inches"),
        legend.box.background = element_rect(colour = "gray"))
red_family
red_genus <- plot_network(ig2.mb, Ma1, type='taxa', color="Genus") + scale_color_viridis_d() +
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.margin = margin(0.25, 0.25, 0.25, 0.25, "inches"),
        legend.box.background = element_rect(colour = "gray"))
red_genus

ggsave("/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/plot_network_family_2024.png", 
       red_family, width = 14, height = 8, dpi = 600)

figure_net<- ggarrange(red_phylum, red_family,
                   labels=c("A","B"),
                   ncol=1,nrow=2)

figure_network <- figure_net + theme(plot.margin = margin(0.25, 0.25, 0.25, 0.25, "inches"))
figure_network

ggsave("/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/plot_network_coocurrence_2024.png", 
       figure_network, width = 14, height = 14, dpi = 600)


#-------------------------------------------------------------------------------
#ANALISIS DE CORRELACION DE METADATOS AMBIENTALES Y ESTRUCTURA DE LA COMUNIDAD MICROBIANA
#ANALISIS DE CORRESPONDENCIA CANONICA ENTRE VARIABLES AMBIENTALES VS. VARIABLES BIOLOGICAS

#Ligar los datos ambientales con los datos de la comunidad microbiana
#Usar los OTU con una abundancia relativa sobre 0.5 para el analisis e importar datos
Ma1_filter <- filter_taxa(Ma1, function(x) mean(x) > 0.05, TRUE)
otutable <- t(otu_table(Ma1_filter))
sampledf <- data.frame(sample_data(Ma1_filter))
#Revisar que analisis aplica para los datos, si CCA o RDA
dca <- decorana(otutable)
dca 
summary(dca)
#si los valores son < 3 usar RDA pero si son > 4 usar CCA
#Envfit + CCA
#escalar a unidad de varianza
df <- data.frame(scale(sampledf[,-c(1:8)])) %>% 
  na.omit()
df
#probar todos los factores ambientales y graficar
ord_all <- cca(otutable ~ ., data=df)
bray <- vegdist(otutable, "bray")
ord_all <- capscale(bray~., df)
plot(ord_all, type = "p", scaling = "sites")
#analisis de inflacion de la varianza para detectar colinealidad (dependencias lineales)
temp <- vif.cca(ord_all)
temp
#dejar las variables que tengan un valor < 10 de VIF
select_para <- names(temp[temp < 10])
select_para
ord <- cca(otutable ~ ., df[,select_para]) 
#encajar los vectores de las variables ambientales en la ordenacion
fit <- envfit(ord, df[,select_para], perm = 999, display = "lc", scaling = "sites")
fit$vectors #revisar significancia
#extraer las variables con mejor significancia (p<0.05)
spp.scrs <- as.data.frame(scores(fit, display = "vectors"))
spp.scrs
pval <-  fit$vectors$pvals
pval 
#flechas de vectores
fdat <- cbind(spp.scrs, Vector = rownames(spp.scrs), pval)
#seleccionar solo las variables significativas
bestEnvVariables<-rownames(fdat)[fdat$pval<=0.05]
#rehacer el analisis con las variables significativas
eval(parse(text=paste("ord1 <- cca(otutable ~ ",do.call(paste,c(as.list(bestEnvVariables),sep=" + ")),",data=sampledf)",sep="")))
summary(ord1)
#hacer el annova de la ordenacion
anova.cca(ord1, perm=9999)
anova.cca(ord1, by="margin", perm=9999) # marginal effects of the terms 
anova.cca(ord1, by="terms", perm=9999) # sequential
anova.cca(ord1, by="axis") # axis, slow
drop1(ord1, test="perm")
#reencajar las variables ambientales en la ordenacion
fit1 <- envfit(ord1,sampledf[,bestEnvVariables], perm = 999, display = "lc", scaling = "sites")
plot(ord1)
plot(ord1, type="n")
points(ord1, display = "sites", 
       col = as.numeric(Source),
       pch=16)
plot(fit1, col = "red", cex=1.2, axis=TRUE, p.max = 0.05)
#graficar usando ggplot2
spp.scrs <- data.frame(scores(fit1, display = "vectors"))
pval <-  fit1$vectors$pvals
#datos para envfit flechas
spp.scrs <- cbind(spp.scrs, Vector = rownames(spp.scrs), pval) # vector table
scrs <- as.data.frame(scores(ord1, display = "sites")) # sample table
scrs <- cbind(scrs, sampledf)
spp.scrs1<- subset(spp.scrs, pval<=0.05) #extracts relevant environment vectors from envifit
spp.scrs1
#biplot
library(digest)
CCA_ord <- ggplot(scrs) +
  geom_point(mapping = aes(x = CCA1, y = CCA2, colour = Source), alpha = 0.8, size = 5) +
  geom_segment(data = spp.scrs1,
               aes(x = 0, xend = CCA1*2.5, y = 0, yend = CCA2*2.5),
               arrow = arrow(length = unit(0.3, "cm")), colour = "#440154FF", size=1.0) +
  geom_text(data = spp.scrs1, aes(x = CCA1*3, y = CCA2*3, label = Vector),
            size = 6) +
  stat_ellipse(aes(x = CCA1, y = CCA2, group = Source), linetype = 2) +
  stat_ellipse(geom = "polygon", aes(x = CCA1, y = CCA2, group = Season_Year, fill = Season_Year), 
               alpha = 0.3, size = 2, linetype = 2) +
  scale_fill_manual(values=c("#440154FF", "#6ece58", "#B8DE29FF", "#FDE725FF")) +
  theme_classic() +
  scale_color_manual(values =c("#404788FF", "#FDE725FF", "#238A8DFF"), guide = guide_legend(ncol=1)) +
  theme(legend.position="right") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14, color = "black"),
        axis.ticks.length = unit(.25, "cm")) +
  theme(legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        legend.key.size = unit(1, "cm"),
        plot.title = element_text(size = 14),
        legend.background = element_blank()) +
  geom_vline(xintercept = c(0), color = "gray", linetype = 1) +
  geom_hline(yintercept = c(0), color = "gray", linetype = 1) +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        plot.margin = margin(0.25, 0.25, 0.25, 0.25, "inches"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  ggtitle("ANOVA: p=0.001, F=1.566")
  
CCA_ord

ggsave("/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/plot_CCA_ord_stats_2024_2.png", 
       CCA_ord, width = 10, height = 10, dpi = 600)

#hacer el cca para la relacion de las variables fisicoquimicas y la abundancia de los taxa
spe <- as.data.frame(scores(ord1, display = "species"))
tax <- tax_table(Ma1_filter) 
otu <- otu_table(Ma1_filter) 
abundance <- rowMeans(x=otu)# calculate mean abundance of each OTU 
df2 <- data.frame(tax, Abundance = abundance) 
df2 <- subset(df2, abundance >=0.05)
df2
spe.scrs <- data.frame(scores(fit1, display = "vectors"))
pval <-  fit1$vectors$pvals
#data for the envfit arrows
spe.scrs <- cbind(spe.scrs, Vector = rownames(spe.scrs), pval) # vector table
spe <- as.data.frame(scores(ord1, display = "species")) # sample table
spe <- cbind(spe, df2)
spe.scrs1<- subset(spe.scrs, pval<=0.05) #extracts relevant environment vectors from envifit
spe.scrs1
#ggplot
CCA_ord2 <- ggplot(spe) +
  geom_point(mapping = aes(x = CCA1, y = CCA2, colour = Phylum), alpha = 0.8, size = 4) +
  geom_segment(data = spe.scrs1,
               aes(x = 0, xend = CCA1*2.5, y = 0, yend = CCA2*2.5),
               arrow = arrow(length = unit(0.25, "cm")), colour = "#440154FF", size=1) +
  geom_text(data = spe.scrs1, aes(x = CCA1*3, y = CCA2*3, label = Vector),
            size = 6) +
  xlim(-2, 4) +
  theme_classic() +
  theme(axis.title=element_text(size=14), axis.text=element_text(size=14), legend.text=element_text(size=14)) +
  scale_color_viridis_d(guide = guide_legend(ncol=1)) +
  theme(legend.position="right") +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14, color = "black"),
        axis.ticks.length = unit(.25, "cm")) +
  theme(legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        plot.title = element_text(size = 14),
        legend.background = element_blank()) +
  geom_vline(xintercept = c(0), color = "gray", linetype = 1) +
  geom_hline(yintercept = c(0), color = "gray", linetype = 1) +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        plot.margin = margin(0.25, 0.25, 0.25, 0.25, "inches"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
  ggtitle("ANOVA: p=0.001, F=1.766")

CCA_ord2

ggsave("/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/plot_CCA_ord2_taxa_2024_2.png", 
       CCA_ord2, width = 10, height = 10, dpi = 600)

figure_CCA<- ggarrange(CCA_ord, CCA_ord2,
                   labels=c("A","B"),
                   ncol=1,nrow=2)

figure_CCA_ord <- figure_CCA + theme(plot.margin = margin(0.25, 0.25, 0.25, 0.25, "inches"))
figure_CCA_ord

ggsave("/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/plot_CCA_ord_taxa+samples_2024_2.png", 
       figure_CCA_ord, width = 8, height = 14, dpi = 600)


#-------------------------------------------------------------------------------
#ANALISIS DEL CORE MICROBIOME

#Definir el CORE microbioma (core-microbiome) y hacer analisis de taxones compartidos (Venn)
#El "core" se define como un grupo de taxones que se detectan en una marcada fraccion de la comunidad a un umbral de abundancia dada
library(microbiomeutilities)
library(RColorBrewer)
Ma.rel <- microbiome::transform(Ma1, "compositional")
Ma.rel.f <- format_to_besthit(Ma.rel)
#Set different detection levels and prevalence
prevalences <- seq(.05, 1, .05) #0.5 = 95% prevalence
detections <- 10^seq(log10(1e-2), log10(.2), length = 10)
#(1e-3) = 0.001% abundance; change "-3" to -2 to increase to 0.01%
core_plot <- plot_core(Ma.rel.f, plot.type = "heatmap", 
               colours = rev(brewer.pal(10, "Spectral")),
               min.prevalence = 0.6, 
               prevalences = prevalences, 
               detections = detections) +
  xlab("Detection Threshold (Relative Abundance (%))") +
  theme_classic()

print(core_plot)

#Familia
Ma.rel.f.fam <- aggregate_taxa(Ma.rel.f, "Family")
prevalences <- seq(.05, 1, .05)
detections <- round(10^seq(log10(1e-5), log10(.2), length = 10), 3)

core_plot1 <- plot_core(Ma.rel.f.fam, 
                plot.type = "heatmap",
                colours = rev(brewer.pal(10, "Spectral")),
                prevalences = prevalences, 
                detections = detections, min.prevalence = .5) +
  xlab("Relative Abundance (%)")
core_plot1 <- core_plot1 + theme_classic() + ylab("Family") +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=14, color = "black"), 
        legend.text=element_text(size=14),
        legend.title = element_text(size = 16),
        legend.key.height = unit(1, "cm"),
        axis.line = element_blank(),
        axis.ticks.length = unit(.25, "cm")) +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        plot.margin = margin(0.25, 0.25, 0.25, 0.25, "inches"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))
core_plot1
#Genero
Ma.rel.f.gen <- aggregate_taxa(Ma.rel.f, "Genus")
prevalences <- seq(.05, 1, .05)
detections <- round(10^seq(log10(1e-5), log10(.2), length = 10), 3)

core_plot2 <- plot_core(Ma.rel.f.gen, 
                plot.type = "heatmap", 
                colours = rev(brewer.pal(10, "Spectral")),
                prevalences = prevalences, 
                detections = detections, min.prevalence = .5) +
  xlab("Relative Abundance (%)")
core_plot2 <- core_plot2 + theme_classic() + ylab("Genus") +
  theme(axis.title=element_text(size=16), 
        axis.text=element_text(size=14, colour = "black"), 
        legend.text=element_text(size=14),
        legend.title = element_text(size = 16),
        legend.key.height = unit(1, "cm"),
        axis.line = element_blank(),
        axis.ticks.length = unit(.25, "cm")) +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        plot.margin = margin(0.25, 0.25, 0.25, 0.25, "inches"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5))
core_plot2

figure_core <- ggarrange(core_plot1, core_plot2,
                         labels = c("A", "B"),
                         ncol = 2, nrow = 1)

figure_core <- figure_core + theme(plot.margin = margin(0.25, 0.25, 0.25, 0.25, "inches"))
figure_core 
  
ggsave("/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/plot_core_prevalence_2024.png", 
       figure_core, width = 16, height = 6, dpi = 600)

#Otra forma de determinar el "core"
Ma1.core <- core(Ma1, detection = .2/100, prevalence = 50/100)
Ma1.core

core.taxa.standard <- core_members(Ma.rel, detection = 0, prevalence = 50/100)
pseq.core <- core(Ma.rel, detection = 0, prevalence = .5)
pseq.core
pseq.core2 <- aggregate_rare(Ma.rel, "Genus", detection = 0, prevalence = .5)
pseq.core2
core.taxa <- taxa(pseq.core)
core.taxa
core.abundance <- sample_sums(core(Ma.rel, detection = .01, prevalence = .95))
det <- c(0, 0.1, 0.5, 2, 5, 20)/100
prevalences <- seq(.05, 1, .05)

plot_core(Ma.rel, 
          prevalences = prevalences, 
          detections = det, 
          plot.type = "lineplot") + 
  xlab("Relative Abundance (%)")


#-------------------------------------------------------------------------------
#ANALISIS DE CONJUNTOS Y PHYLOGRUPOS COMPARTIDOS ENTRE MUESTRAS (DIAGRAMA DE VENN)

#Analisis de conjuntos para determinar taxones compartidos y unicos con base en los grupos
library(MicEco)

#Para fuente (mucus, tejido y agua)
ps_venn(Ma1, group = "Source", quantities = list(type=c("counts")))
#teniendo en cuenta el peso de la abundancia relativa
ps_venn(Ma1, group = "Source", quantities = list(type=c("percent")), weight = TRUE)
#incluyendo taxones presentes en al menos el 5% de las muestras en cada grupo
ps_venn(Ma1, group = "Source", quantities = list(type=c("counts")), fraction = 0.05)
#con color
v1 <- ps_venn(Ma1, group = "Source", fill = c("#404788FF", "#238A8DFF", "#FDE725FF"), 
        alpha = 0.5, quantities = list(type=c("percent", "counts"), fraction = 0.05, weight = TRUE, col = "grey16"))
#Para Sitio (protegido, urbano)
v2 <- ps_venn(Ma1, group = "Site", fill = c("#404788FF", "#238A8DFF", "#55C667FF"), 
                       alpha = 0.5, quantities = list(type=c("percent", "counts"), fraction = 0.05, weight = TRUE, col = "grey16"))
#Para Epoca climatica (seca, lluviosa)
v3 <- ps_venn(Ma1, group = "Season_Year", fill = c("#404788FF", "#238A8DFF", "#55C667FF", "#FDE725FF"), 
                       alpha = 0.5, quantities = list(type=c("percent", "counts"), fraction = 0.05, weight = TRUE, col = "grey16"))
#Para Contacto con alga (sin y con)
v4 <- ps_venn(Ma1, group = "Algae_Touching", fill = c("#73D055FF", "#FDE725FF", "#31688e"), 
                       alpha = 0.5, quantities = list(type=c("percent", "counts"), fraction = 0.05, weight = TRUE, col = "grey16"))
#Para Estado de salud (sano y estresado)
v5 <- ps_venn(Ma1, group = "Health_Status", fill = c("#73D055FF", "#FDE725FF", "#31688e"), 
                       alpha = 0.5, quantities = list(type=c("percent", "counts"), fraction = 0.05, weight = TRUE, col = "grey16"))


#figura de grupos (venn) + diversidad alfa
#SOURCE
figure_venn_alfa1 <- ggarrange(v1, b1, b6,
                              labels=c("A","B", "C"),
                              ncol=3,nrow=1)
figure_venn_alfa1 <- figure_venn_alfa1 + theme(plot.margin = margin(0.25,0.25,0.25,0.25, "inches")) +
  theme(panel.background = element_rect(fill = "white", colour = "white"))
figure_venn_alfa1

#SEASON
figure_venn_alfa2 <- ggarrange(v3, b4, b7,
                              labels=c("A","B", "C"),
                              ncol=3,nrow=1)
figure_venn_alfa2 <- figure_venn_alfa2 + theme(plot.margin = margin(0.25,0.25,0.25,0.25, "inches")) +
  theme(panel.background = element_rect(fill = "white", colour = "white"))
figure_venn_alfa2

#SITE
figure_venn_alfa3 <- ggarrange(v2, b5, b8,
                              labels=c("A","B", "C"),
                              ncol=3,nrow=1)
figure_venn_alfa3 <- figure_venn_alfa3 + theme(plot.margin = margin(0.25,0.25,0.25,0.25, "inches")) +
  theme(panel.background = element_rect(fill = "white", colour = "white"))
figure_venn_alfa3


ggsave("/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/plot_venn_alfa_source_2024.png", 
       figure_venn_alfa1, width = 12, height = 6, dpi = 600)
ggsave("/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/plot_venn_alfa_season_2024.png", 
       figure_venn_alfa2, width = 12, height = 6, dpi = 600)
ggsave("/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/plot_venn_alfa_site_2024.png", 
       figure_venn_alfa3, width = 12, height = 6, dpi = 600)


#-------------------------------------------------------------------------------
#ANALISIS DE CORRELACION DE LAS VARIABLES AMBIENTALES
#ANALISIS DE COMPONENTES PRINCIPALES PCA
#cargar los metadatos
env_metadata<-read.csv("/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/Metadata_Environmental_2.csv", sep = ",", header = TRUE, check.names = FALSE)
head(env_metadata)
library(caret)
library(textshape)
library(ggfortify)
library(radiant.data)

df_alpha <- env_metadata %>%  
  column_to_rownames("Sample_ID") %>% 
  dplyr::select(where(is.numeric)) %>%
  na.omit() %>% 
  scale(center = T, scale = T) %>% 
  as.data.frame()
cor_matrix <- cor(df_alpha)
highly_correlated <- caret::findCorrelation(cor_matrix, cutoff = 0.9) 
# none
# df_alpha <- df_alpha[, - highly_correlated]

# lets see if any env vars can be packed neatly in PCs
pca_env <- prcomp(df_alpha, scale = FALSE)
biplot(pca_env)
#plot
pca_plot <- ggplot2::autoplot(pca_env,
                  data = 
                    (df_alpha %>% 
                        rownames_to_column("Sample_ID") %>% 
                        left_join(., env_metadata[,c("Season_Year", "Site", "Sample_ID")], 
                                  by = "Sample_ID")), 
                  colour = 'Season_Year',
                  shape = 'Site',
                  size = 5, 
                  alpha = 0.6,
                  loadings = TRUE, 
                  loadings.colour = "black", 
                  loadings.label = TRUE, 
                  loadings.label.size = 6, 
                  loadings.label.colour = "black",
                  loadings.label.repel = TRUE) +
  stat_ellipse(aes(group = Season_Year, label = Season_Year), linetype = 2) +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 1) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 1) +
  scale_colour_viridis_d() +
  theme_classic() +
  theme(legend.background = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14, color = "black"),
        axis.ticks.length = unit(.25, "cm"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        plot.title = element_text(size = 14))

pca_plot1 <- pca_plot + theme(panel.background = element_rect(fill = "white", colour = "black"),
                 plot.margin = margin(0.25, 0.25, 0.25, 0.25, "inches"),
                 panel.border = element_rect(colour = "black", fill=NA, size=0.5))
pca_plot1

ggsave("/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/plot_PCA_enviromental_2024.png", 
       pca_plot1, width = 7, height = 7, dpi = 600)

summary(pca_env)$importance %>% 
  as.data.frame() %>% 
  dplyr::select("PC1", "PC2", "PC3", "PC4")

# Extract the loadings from the PCA result and convert to tidy format 
loadings_df <-   tidy(pca_env, matrix = "rotation")

# PC1
# Create the ggplot2 object for plotting the loadings
loadings_df %>% filter(PC == "1") %>% 
  #x=reorder(class,-amount,sum)
  ggplot( aes(x = reorder(column, -value), y = value, group = factor(column))) +
  geom_bar(stat = "identity", position = "dodge", 
           fill= "gray90", color = "gray20") +
  geom_text(aes(label = column), position = position_dodge(width = 0.9), 
            #vjust = "inward",  
            angle = 90 , hjust = "inward") +
  labs(x = "Variables of PC1", y = "Loadings") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        legend.position = "none")  -> plot_pca1

# PC2 
loadings_df %>% filter(PC == "2") %>% 
  #x=reorder(class,-amount,sum)
  ggplot( aes(x = reorder(column, -value), y = value, group = factor(column))) +
  geom_bar(stat = "identity", position = "dodge", 
           fill= "gray90", color = "gray20") +
  geom_text(aes(label = column), position = position_dodge(width = 0.9), 
            #vjust = "inward",  
            angle = 90 , hjust = "inward") +
  labs(x = "Variables of PC2", y = "Loadings") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        legend.position = "none") -> plot_pca2

# PC3
loadings_df %>% filter(PC == "3") %>% 
  #x=reorder(class,-amount,sum)
  ggplot( aes(x = reorder(column, -value), y = value, group = factor(column))) +
  geom_bar(stat = "identity", position = "dodge", 
           fill= "gray90", color = "gray20") +
  geom_text(aes(label = column), position = position_dodge(width = 0.9), 
            #vjust = "inward",  
            angle = 90 , hjust = "inward") +
  labs(x = "Variables of PC3", y = "Loadings") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        legend.position = "none") -> plot_pca3

pca_variables <- plot_pca1 + plot_pca2 + plot_pca3
pca_variables

ggsave("/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/plot_loadings_PCA_2024.png", 
       pca_variables, width = 14, height = 6, dpi = 600)

#Figura de todos las ordenaciones (PCA + CCAs) de las variables ambientales
figure_ord<- ggarrange(pca_plot1, CCA_ord, CCA_ord2,
                       labels=c("A","B","C"),
                       ncol=3,nrow=1)

figure_ords <- figure_ord + theme(plot.margin = margin(0.25, 0.25, 0.25, 0.25, "inches"))
figure_ords

ggsave("/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/plot_odinations_PCA+CCA_2024_2.png", 
       figure_ords, width = 22, height = 12, dpi = 600)


#-------------------------------------------------------------------------------
#ANALISIS ECOLOGICO DE LA RELACION DE LAS VARIABLES AMBIENTALES CON LA ESTRUCTURA DE LA COMUNIDAD MICROBIANA
#cargar los paquetes y datos

library(vegan)
library(ggplot2)

env_data <- read.csv("/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/Env_Metadata_Samples_2024.csv", sep = ",", header = TRUE, check.names = FALSE, row.names = 1)
otu_data <- read.csv("/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/RarefiedTaxa_Samples_Abundance_2024.csv", sep = ",", header = TRUE, check.names = FALSE, row.names = 1)
otu_data<-t(otu_data)
otu_data<-otu_data[rowSums(otu_data)>200,]
env_data<-env_data[rownames(otu_data),]

sel_env<-c("Temperature","Salinity","pH","Oxygen","Ammonia","Phosphates","Nitrates","Nitrites","TSS")
sel_env_label <- list(
  'Temperature'="Temperature",
  'Salinity'="Salinity",
  'pH'="pH",
  'Oxygen'="Oxygen",
  'Ammonia'="Ammonia",
  'Phosphates'="Phosphates",
  'Nitrates'="Nitrates",
  'Nitrites'="Nitrites",
  'TSS'="TSS")
sel_env_label<-t(as.data.frame(sel_env_label))
sel_env_label<-as.data.frame(sel_env_label)
colnames(sel_env_label)<-c("Trans")
sel_env_label$Trans<-as.character(sel_env_label$Trans)

env_data_filtered<-env_data[,sel_env]
otu_data_filtered<-otu_data[rownames(env_data_filtered),]
env_data_filtered
x<-log((otu_data+1)/(rowSums(otu_data_filtered)+dim(otu_data_filtered)[2]))
x<-x[,order(colSums(x),decreasing=TRUE)]
N<-80
taxa_list<-colnames(x)[1:N]
taxa_list<-taxa_list[!grepl("Unknown",taxa_list)]
N<-length(taxa_list)
x<-data.frame(x[,colnames(x) %in% taxa_list])
y<-env_data_filtered

grouping_info<-data.frame(row.names=rownames(otu_data),
                          t(as.data.frame(strsplit(rownames(otu_data),"-"))))
head(grouping_info)
grouping_info
groups<-grouping_info[,1]
groups
method<-"kendall"

corr<-NULL
for(g in colnames(x)){
  for(h in colnames(y)){
    for(i in unique(groups)){
      a<-x[groups==i,g,drop=F]
      b<-y[groups==i,h,drop=F]
      tmp<-c(g,h,cor(a[complete.cases(b),],b[complete.cases(b),],use="everything",method=method),cor.test(a[complete.cases(b),],b[complete.cases(b),],method=method)$p.value,i)
      if(is.null(corr)){
        corr<-tmp  
      }
      else{
        corr<-rbind(corr,tmp)
      }    
    }
  }
}

corr<-data.frame(row.names=NULL,corr)
corr
colnames(corr)<-c("Taxa","Env","Correlation","Pvalue","Type")
corr
corr$Pvalue<-as.numeric(as.character(corr$Pvalue))
corr$AdjPvalue<-rep(0,dim(corr)[1])
corr$Correlation<-as.numeric(as.character(corr$Correlation))
corr

adjustment_label<-c("NoAdj","AdjEnvAndType","AdjTaxaAndType","AdjTaxa","AdjEnv")
adjustment<-1
if(adjustment==1){
  corr$AdjPvalue<-corr$Pvalue
} else if (adjustment==2){
  for(g in unique(corr$Env)){
    for(h in unique(corr$Type)){
      sel<-corr$Env==g & corr$Type==h
      corr$AdjPvalue[sel]<-p.adjust(corr$Pvalue[sel],method="BH")
    }
  }
} else if (adjustment==3){
  for(g in unique(corr$Taxa)){
    for(h in unique(corr$Type)){
      sel<-corr$Taxa==g & df$Type==h
      corr$AdjPvalue[sel]<-p.adjust(corr$Pvalue[sel],method="BH")
    }
  }
} else if (adjustment==4){
  for(g in unique(corr$Taxa)){
    sel<-corr$Taxa==g
    corr$AdjPvalue[sel]<-p.adjust(corr$Pvalue[sel],method="BH")
  }
} else if (adjustment==5){
  for(g in unique(corr$Env)){
    sel<-corr$Env==g
    corr$AdjPvalue[sel]<-p.adjust(corr$Pvalue[sel],method="BH")
  }
}

#Now we generate the labels for signifant values
corr$Significance<-cut(corr$AdjPvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))

corr<-corr[complete.cases(corr),]
corr

write.csv(corr, file="/Users/jordanruiz/Desktop/Doctoral_Thesis_Analysis/Objective_3_Microbiome_Madracis/Kendall_correlation_env_taxa.csv")

#We want to reorganize the Env data based on they appear
corr$Env<-factor(corr$Env,as.character(corr$Env))

#We use the function to change the labels for facet_grid in ggplot2
Env_labeller <- function(variable,value){
  return(sel_env_label[as.character(value),"Trans"])}

heat_map_env <- ggplot(aes(x=Type, y=Taxa, fill=Correlation), data=corr)
heat_map_env <- heat_map_env + geom_tile() + scale_fill_gradient2(low="blue", mid="white", high="red") + theme_classic()
heat_map_env <- heat_map_env + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5))
heat_map_env <- heat_map_env + geom_text(aes(label=Significance), color="black", size=3)+labs(y=NULL, x=NULL, fill=method)
heat_map_env <- heat_map_env + facet_grid(. ~ Env, drop=TRUE,scale="free",space="free_x",labeller=Env_labeller)
heat_map_env

png(filename = "Heatmap_env2.png", width = 13, height = 11, res = 1200, units = "in")
print(heat_map_env)
dev.off()


#-------------------------------------------------------------------------------
#OTROS...

Ma_Site <- merge_samples(Ma1, "Site")
plot_bar(Ma_Site, fill = "Family") + 
  geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack")

Ma_Season <- merge_samples(Ma1, "Season")
plot_bar(Ma_Season, fill = "Family") + 
  geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack")

Ma_Health <- merge_samples(Ma1, "Health_Status")
plot_bar(Ma_Health, fill = "Family") + 
  geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack")

Ma_Algae <- merge_samples(Ma1, "Algae_Touching")
plot_bar(Ma_Algae, fill = "Family") + 
  geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack")

Ma_Source <- merge_samples(Ma1, "Source")
plot_bar(Ma_Source, fill = "Family") + 
  geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack")



Ma_Vibrio <- subset_taxa(Ma1, Family %in% c("Vibrionaceae"))

Ma_Firmi <- subset_taxa(Ma, Phylum %in% c("Firmicutes"))

Ma_Actino <- subset_taxa(Ma, Phylum %in% c("Actinobacteriota"))

Ma_Proteo <- subset_taxa(Ma, Phylum %in% c("Proteobacteria"))

Ma_Gamma <- subset_taxa(Ma1, Class %in% c("Gammaproteobacteria"))

Ma_Health <- subset_samples(Ma, Health_Status %in% c("Healthy", "Stressed"))

plot_bar(Ma_Vibrio, x="Genus", fill = "Genus", facet_grid = Site~Season) +
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")

plot_bar(Ma_Health, x="Phylum", fill = "Phylum", facet_grid = Site~Source) +
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

plot_bar(Ma_Firmi, x="Family", fill = "Family", facet_grid = Site~Source) +
  geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack")

plot_bar(Ma_Actino, x="Family", fill = "Family", facet_grid = Site~Source) +
  geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack")

plot_bar(Ma_Proteo, x="Class", fill = "Class", facet_grid = Site~Source) +
  geom_bar(aes(color=Class, fill=Class), stat="identity", position="stack")

plot_bar(Ma_Gamma, x="Family", fill = "Family", facet_grid = Site~Source) +
  geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack")

plot_bar(Ma_Gamma, x="Family", fill = "Family", facet_grid = Health_Status~Source) +
  geom_bar(aes(color=Family, fill=Family), stat="identity", position="stack")

Ma_Psedoal <- subset_taxa(Ma1, Family %in% c("Pseudoalteromonadaceae"))

plot_bar(Ma_Psedoal, x="Genus", fill = "Genus", facet_grid = Site~Source~Season) +
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")


plot_heatmap(Ma1, method = "NMDS", distance = "bray")

Ma_abund <- filter_taxa(Ma1, function(x) sum(x > total*0.009) > 0, TRUE)
Ma_abund

plot_heatmap(Ma_abund, method = "NMDS", distance = "bray")

plot_heatmap(Ma_abund, method = "MDS", distance = "bray", 
             taxa.label = "Genus", taxa.order = "Genus", 
             trans=NULL, low="beige", high="cyan4", na.value="beige")


plot_heatmap(Ma_Gamma, method = "NMDS", distance = "bray", 
             taxa.label = "Family", taxa.order = "Family", 
             low="beige", high="red", na.value="beige")


plot_richness(Ma1, measures=c("Chao1", "Shannon"))
plot_richness(Ma1, measures=c("Chao1", "Shannon"), x="Site", color="Season")
plot_richness(Ma, measures=c("Chao1", "Shannon"), x="Source", color="Health_Status")

Ma.ord <- ordinate(Ma1, "NMDS", "bray")
plot_ordination(Ma1, Ma.ord, type="taxa", color="Phylum", shape= "Kingdom", 
                title="ASVs")

plot_ordination(Ma, Ma.ord, type="taxa", color="Kingdom", 
                title="ASVs", label="Phylum") + 
  facet_wrap(~Kingdom, 2)

plot_ordination(Ma1, Ma.ord, type="samples", color="Site", 
                shape="Season", title="Samples") + geom_point(size=3)

plot_ordination(Ma, Ma.ord, type="split", color="Class", 
                shape="Source", title="biplot", label = "Site") +  
  geom_point(size=3)

plot_net(Ma_rarefied, distance = "(A+B-2*J)/(A+B)", type = "taxa", 
         maxdist = 0.7, color="Phylum", point_label="Genus")

plot_net(Ma_rarefied, distance = "(A+B-2*J)/(A+B)", type = "taxa", 
         maxdist = 0.8, color="Family", point_label="Phylum") 







plot_richness(Ma_rarefied, x="Site", measures=c("Observed", "Shannon")) +
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90))

plot_richness(Ma_rarefied, x="Season", measures=c("Observed", "Shannon")) +
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90))

plot_richness(Ma_rarefied, x="Health_Status", measures=c("Observed", "Shannon")) +
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90))

plot_richness(Ma_rarefied, x="Algae_Touching", measures=c("Observed", "Shannon")) +
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90))

plot_richness(Ma_rarefied, x="Source", measures=c("Observed", "Shannon")) +
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90))


















