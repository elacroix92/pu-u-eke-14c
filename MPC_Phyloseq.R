#Glade Dlott 2019_07_26

library(phyloseq); packageVersion("phyloseq")

library(Biostrings); packageVersion("Biostrings")

library(ggplot2); packageVersion("ggplot2")

theme_set(theme_bw())

psMPC <- readRDS("~/Desktop/MPC_Illumina_Data/MPC_16S_phyloseq_2019_07_23.rds")

psMPC

#Make ordinations in phyloseq

#remove 0-samples

Bacteria_No0s <- subset_samples(psMPC, sample_sums(psMPC) > 0)

Bacteria_No0s

#ordinate bacteria

Bacteria_No0s_Ordination <- ordinate(Bacteria_No0s, "NMDS", "bray")

#plot ordination

head(sample_data(Bacteria_No0s))

plot_ordination(Bacteria_No0s, Bacteria_No0s_Ordination, label = "horizon", color = "distance_mm", title = "Bacteria")+theme_classic()

#Trim unrealistic samples

sample_sums(Bacteria_No0s)

Bacteria_No0s_Trim <- subset_samples(Bacteria_No0s,
                                  sample_sums(Bacteria_No0s) > 200)

Bacteria_No0s_Trim

Bacteria_No0s_Trim_Ordination <- ordinate(Bacteria_No0s_Trim, "NMDS", "bray")
plot_ordination(Bacteria_No0s_Trim, Bacteria_No0s_Trim_Ordination, label = "horizon", color = "distance_mm", title = "Bacteria")+theme_classic()

#color by per_t_delta_14C

plot_ordination(Bacteria_No0s_Trim, Bacteria_No0s_Trim_Ordination, label = "horizon", color = "per_t_delta_14C", title = "Bacteria")+theme_classic()

#color by age_14C_age_bp

plot_ordination(Bacteria_No0s_Trim, Bacteria_No0s_Trim_Ordination, label = "horizon", color = "age_14C_age_bp", title = "Bacteria")+theme_classic()

#No seeiming influence of anything but horizon and distance, horizons 4 and 5 are weird.

#Heatmaps!
  
#Make heatmaps in phyloseq

#Subset 100 most abundant ASVs

Top100BacteriaASV_names <- names(sort(taxa_sums(Bacteria_No0s_Trim), decreasing = TRUE))[1:100]

Top100BacteriaASV <- prune_taxa(Top100BacteriaASV_names, Bacteria_No0s_Trim)

#Make heatmap, ordered by sampling process

plot_heatmap(Top100BacteriaASV, "NMDS", "bray", "horizon", "Class", sample.order = "sample_id", low="#FFFFCC", high="#000033", na.value="white")

#Same heatmap, by tree

plot_heatmap(Top100BacteriaASV_No0, "NMDS", "bray", "proc", "Genus", sample.order = "tree", low="#FFFFCC", high="#000033", na.value="white")

#Top 10?

#Subset 10 most abundant ASVs

Top10BacteriaASV_names <- names(sort(taxa_sums(Bacteria_No0s_Trim), decreasing = TRUE))[1:10]

Top10BacteriaASV <- prune_taxa(Top10BacteriaASV_names, Bacteria_No0s_Trim)

#Make heatmap, ordered by sampling process

plot_heatmap(Top10BacteriaASV, "NMDS", "bray", "horizon", "Class", sample.order = "sample_id", low="#FFFFCC", high="#000033", na.value="white")

#Not much to say or speculate.

