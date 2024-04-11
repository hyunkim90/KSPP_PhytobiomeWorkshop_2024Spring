######## 2024_Spring_KSPP workshop part.2
### 0. Install and load packages and set the working directory
packages <- c("phyloseq", "vegan", "ggplot2", "lawstat", "FSA", "rcompanion")
install_packages <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
if (length(install_packages) > 0) {
  install.packages(install_packages)
}

library(phyloseq)
library(vegan)
library(ggplot2)
library(lawstat)
library(FSA)
library(rcompanion)

# Set the working directory
setwd("C:/Users/WIN/Desktop/KSPP_PhytobiomeWorkshop_2024Spring-main")
getwd()
wd <- "./Session_1/Output/"
setwd(wd)
getwd()

#### 1. Data preparation
## 1-1. Import biom, metadata, tree file
a <- import_biom("otu_table_final.biom")
b <- import_qiime_sample_data("sample_metadata.tsv")
c <- read_tree("tree.nwk")

# merge phyloseq object
run <- merge_phyloseq(a,b,c)
run

## 1-2. Change taxonomy rank names
colnames(tax_table(run))
colnames(tax_table(run)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
colnames(tax_table(run))

## 1-3. Filter taxa and rarefy
# "subset_taxa()" can subset for your purpose
# For example, if you want to subset Chloroplast order, use subset_taxa(phyloseq, Order=="Chloroplast")
# If you want to filter the Chloroplast order, use subset_taxa(phyloseq, Order!="Chloroplast")
# "==" means check left side and right side are same.
# "!=" means check left side and right side are different.
# Using this function, you can remain or delete taxa you want.
# We already filtered chloroplast, mitochondira and others in Qiime2.
run <- subset_taxa(run, Kingdom != "Unassigned") # remove "Unassigned" kingdom
run

# rarefaction
rare.run <- rarefy_even_depth(run, rngseed = 1024)
rare.run

## 1-4. Save the phyloseq object as .Rdata
save.wd <- "./../../Session_2"
setwd(save.wd)
getwd()

phyloseq_object <- rare.run
save(phyloseq_object, file="phyloseq_object.RData")

#### 2. Alpha diversity
## 2-1. Make a boxplot using phyloseq basic function (embeded ggplot2)
plot_richness(rare.run, x="Branch_number", color="Branch_number")
# measures=c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "invSimpson", "Fisher"), default is all index

## 2-2. Make a boxplot using phyloseq and ggplot2 advanced
plot_richness(rare.run, x="Branch_number", color="Branch_number", measures=c("Observed","Chao1","Simpson","Shannon")) +
  geom_boxplot() +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.y = element_text(size=10)) + theme(axis.text.x = element_text(size=10, angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.title = element_text(size=15, hjust=0.5)) +
  theme(legend.position="none") +
  theme(strip.background=element_rect(fill="lightgrey")) +
  theme(strip.text = element_text(size=12)) +
  labs(x="\n Branch", y="Alpha diversity \n") +
  scale_color_manual(values=c("Branch_2"="#FF5765", "Branch_3"="#01DEE6", "Branch_4"="#8A6FDF", "Branch_5"="#7CB53A"))

## 2-3. Calculate alpha diversity
alpha <- estimate_richness(
  rare.run,
  split=TRUE # Logical (set TRUE -> separate set of richness estimates be performed for each sample,
  # set FALSE -> pool all samples and estimate richness of the entire set)
  # measures: c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher"))
)
head(alpha)

## 2-4. Modify dataframe of "alpha" object for plotting
alpha$SampleID <- row.names(alpha) # add the "SampleID" column
alpha$Branch_number <- sample_data(rare.run)$Branch_number # add the "Branch_number" column of sample_metadata
rownames(alpha) <- NULL # remove the rowname
head(alpha)

## 2-5. Observed OTU
# Statistical tests
shapiro.test(alpha$Observed) # normality test
kruskal.test(Observed~Branch_number, data=alpha) # Kruskal Wallis Test
DT <- dunnTest(Observed~Branch_number, data=alpha)
PT <- DT$res
cld <- cldList(P.adj ~ Comparison,data = PT, threshold = 0.05)
cld
names(cld) <- c("Branch_number", "Letter", "MonoLetter")
cld

# plotting
ggplot(data = alpha, aes(x=Branch_number, y=Observed, fill=Branch_number)) +
  geom_boxplot(lwd=0.8, width = 0.8, outlier.shape = 16) +
  theme_bw() + 
  scale_fill_manual(values=c("Branch_2"="#FF5765", "Branch_3"="#01DEE6",
                             "Branch_4"="#8A6FDF", "Branch_5"="#7CB53A")) +
  #geom_point(position='jitter',shape=1, alpha=0.8, colour="black", size=1.5) + 
  xlab("\n Branch") + ylab("Observed OTU \n") +
  theme(axis.title.x = element_text(size = 15, hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15, hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(size=12, angle = 0, hjust = 0.5, vjust=0.5, face='bold',color='black'))+
  theme(axis.text.y = element_text(size = 15, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none") + geom_text(data=cld, aes(label=Letter, x=Branch_number, y=31))

## 2-6. Shannon diversity
shapiro.test(alpha$Shannon) # normality test
kruskal.test(Shannon~Branch_number, data=alpha) # Kruskal Wallis Test (nonparametric test)
DT <- dunnTest(Shannon~Branch_number, data=alpha)
PT <- DT$res
cld <- cldList(P.adj ~ Comparison,data = PT, threshold = 0.05)
cld
names(cld) <- c("Branch_number", "Letter", "MonoLetter")
cld

# plotting
ggplot(data = alpha, aes(x=Branch_number, y=Shannon, fill=Branch_number)) +
  geom_boxplot(lwd=0.8, width = 0.8, outlier.shape = 16) +
  theme_bw() + scale_fill_manual(values=c("Branch_2"="#FF5765", "Branch_3"="#01DEE6",
                                          "Branch_4"="#8A6FDF", "Branch_5"="#7CB53A")) +
  #geom_point(position='jitter',shape=1, alpha=0.8, colour="black", size=1.5) + 
  xlab("\n Branch") + ylab("Shannon \n") +
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none") +
  geom_text(data=cld, aes(label=Letter, x=Branch_number, y=2.5))

## 2-7. Simpson
shapiro.test(alpha$Simpson) # normality test
kruskal.test(Simpson~Branch_number, data=alpha) # Kruskal Wallis Test (nonparametric test)
DT <- dunnTest(Simpson~Branch_number, data=alpha)
PT <- DT$res
cld <- cldList(P.adj ~ Comparison,data = PT, threshold = 0.05)
cld
names(cld) <- c("Branch_number", "Letter", "MonoLetter")
cld

# plotting
ggplot(data = alpha, aes(x=Branch_number, y=Simpson, fill=Branch_number)) +
  geom_boxplot(lwd=0.8, width = 0.8, outlier.shape=16) +
  theme_bw() + scale_fill_manual(values=c("Branch_2"="#FF5765", "Branch_3"="#01DEE6",
                                          "Branch_4"="#8A6FDF", "Branch_5"="#7CB53A")) +
  #geom_point(position='jitter',shape=1, alpha=0.8, colour="black", size=1.5) + 
  xlab("\n Branch") + ylab("Simpson \n") +
  theme(axis.title.x = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15,hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust=0.4,size=12, face='bold',color='black'))+
  theme(axis.text.y = element_text(size=15, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none") +
  geom_text(data=cld, aes(label=Letter, x=Branch_number, y=1.0))

## 2-8. all index plotting
# change the dataframe for plotting all index
alpha.melt <- reshape2::melt(alpha)
head(alpha.melt)
names(alpha.melt) <- c("SampleID", "Branch_number", "Index", "Value")
head(alpha.melt)

# plotting
alpha.melt.wanted <- subset(alpha.melt, Index %in% c("Observed", "Shannon", "Simpson", "InvSimpson"))
ggplot(data = alpha.melt.wanted, aes(x=Branch_number, y=Value, fill=Branch_number)) +
  geom_boxplot(lwd=0.8, width = 0.8, outlier.shape=16) +
  facet_wrap(vars(Index), scales = "free") +
  theme_bw() + 
  scale_fill_manual(values=c("Branch_2"="#FF5765", "Branch_3"="#01DEE6",
                             "Branch_4"="#8A6FDF", "Branch_5"="#7CB53A")) +
  #geom_point(position='jitter',shape=1, alpha=0.8, colour="black", size=1.5) + 
  xlab("\n Branch") + ylab("Alpha diversity \n") +
  theme(axis.title.x = element_text(size = 15, hjust = 0.5, face='bold')) + 
  theme(axis.title.y = element_text(size = 15, hjust = 0.5, face='bold')) + 
  theme(axis.text.x = element_text(size=12, angle = 0, hjust = 0.5, vjust=0.5, face='bold',color='black'))+
  theme(axis.text.y = element_text(size = 15, face='bold',color='black'))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none") +
  theme(aspect.ratio = 1) +
  theme(panel.spacing=unit(1, "cm"))

#### 3. Beta diversity
## 3-1. Bray-curtis dissimilarity and PCoA ordination
ra = transform_sample_counts(run, function(x) 100 * x/sum(x))
bray_distance <- phyloseq::distance(ra, "bray")
ordu1 <- ordinate(run, "PCoA", "bray")

# PERMANOVA
ra_df <- data.frame(sample_data(ra))
p <- adonis2(bray_distance ~ Branch_number, data=ra_df, permutations = 999)
p
p <- p$`Pr(>F)`
p <- p[1]
p

# plotting
bray_pcoa <- plot_ordination(ra, ordu1, color="Branch_number")+
  geom_point(size=3) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15)) +
  theme(legend.title = element_text(size = 15)) +
  theme(legend.text = element_text(size = 12)) +
  theme(aspect.ratio = 1) +
  annotate("text", label=paste("p-value =",round(p, 3)), x=Inf, y=-Inf,
           hjust=1.05, vjust=-0.6) +
  scale_color_manual(values = c("Branch_2"="#FF5765", "Branch_3"="#01DEE6",
                                "Branch_4"="#8A6FDF", "Branch_5"="#7CB53A"))
bray_pcoa

# change the axis title
explained <- ordu1$values$Relative_eig
bray_pcoa <- plot_ordination(ra, ordu1, color="Branch_number")+
  geom_point(size=3) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15)) +
  theme(legend.title = element_text(size = 15)) +
  theme(legend.text = element_text(size = 12)) +
  coord_fixed(ratio=1) +
  scale_color_manual(values = c("Branch_2"="#FF5765", "Branch_3"="#01DEE6",
                                "Branch_4"="#8A6FDF", "Branch_5"="#7CB53A")) +
  labs(x = paste("PCoA 1 [", round(explained[1] * 100, 2), "%]", sep=""), 
       y = paste("PCoA 2 [", round(explained[2] * 100, 2), "%]", sep="")) +
  annotate("text", label=paste("p-value =",round(p, 3)), x=Inf, y=-Inf,
           hjust=1.05, vjust=-0.6)
bray_pcoa

## 3-2. unweighted Unifrac and PCoA ordination
unweighted_unifrac_distance <- phyloseq::distance(ra, method="unifrac", weighted=FALSE)
ordu2 <- ordinate(run, "PCoA", "unifrac", weighted=FALSE)

# PERMANOVA
ra_df <- data.frame(sample_data(ra))
p <- adonis2(unweighted_unifrac_distance ~ Branch_number, data=ra_df, permutations = 999)
p
p <- p$`Pr(>F)`
p <- p[1]
p

# plotting
uw_unifrac_pcoa <- plot_ordination(ra, ordu2, color="Branch_number")+
  geom_point(size=3) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15)) +
  theme(legend.title = element_text(size = 15)) +
  theme(legend.text = element_text(size = 12)) +
  theme(aspect.ratio = 1) +
  annotate("text", label=paste("p-value =",round(p, 3)), x=Inf, y=-Inf,
           hjust=1.05, vjust=-0.6) +
  scale_color_manual(values = c("Branch_2"="#FF5765", "Branch_3"="#01DEE6",
                                "Branch_4"="#8A6FDF", "Branch_5"="#7CB53A"))
uw_unifrac_pcoa

# change the axis title
explained <- ordu2$values$Relative_eig
uw_unifrac_pcoa <- plot_ordination(ra, ordu2, color="Branch_number")+
  geom_point(size=3) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15)) +
  theme(legend.title = element_text(size = 15)) +
  theme(legend.text = element_text(size = 12)) +
  theme(aspect.ratio = 1) +
  scale_color_manual(values = c("Branch_2"="#FF5765", "Branch_3"="#01DEE6",
                                "Branch_4"="#8A6FDF", "Branch_5"="#7CB53A")) +
  labs(x = paste("PCoA 1 [", round(explained[1] * 100, 2), "%]", sep=""), 
       y = paste("PCoA 2 [", round(explained[2] * 100, 2), "%]", sep="")) +
  annotate("text", label=paste("p-value =",round(p, 3)), x=Inf, y=-Inf,
           hjust=1.05, vjust=-0.6)
uw_unifrac_pcoa

## 3-3. Bray-curtis dissimilarity and NMDS ordination
ordu3 <- ordinate(ra, "NMDS", "bray")
stress <- ordu3$stress # The stress value is indicators of the quality of the NMDS. 
# The stress value represents the degree to which data has failed to project in to a low-dimensional space 
# while maintaining the similarity of the original data, and usually has a value between 0 and 1.
# In general, the smaller the stress value, the better the NMDS result is for the original data.
# stress < 0.2 - good
# stress < 0.1 - very good
# stress < 0.05 - perfect

p <- adonis2(bray_distance ~ Branch_number, data=ra_df, permutations = 999)
p <- p$`Pr(>F)`
p <- p[1]
p

# plotting
bray_NMDS <- plot_ordination(ra, ordu3, color="Branch_number") + 
  geom_point(size=3) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15)) +
  theme(legend.title = element_text(size = 15)) +
  theme(legend.text = element_text(size = 12)) +
  annotate("text", label=paste("stress =",round(stress, 2),"\n p-value =",round(p, 3)), x=Inf, y=-Inf,
           hjust=1.05, vjust=-0.6) +
  theme(aspect.ratio = 1) +
  scale_color_manual(values = c("Branch_2"="#FF5765", "Branch_3"="#01DEE6",
                                "Branch_4"="#8A6FDF", "Branch_5"="#7CB53A"))
bray_NMDS

## 3-4. Unweighted Unifrac NMDS
ordu4 <- ordinate(run, "NMDS", "unifrac", weighted=FALSE)
stress <- ordu4$stress

p <- adonis2(unweighted_unifrac_distance ~ Branch_number, data=ra_df, permutations = 999)
p <- p$`Pr(>F)`
p <- p[1]
p

# plotting
uw_unifrac_NMDS <- plot_ordination(ra, ordu4, color="Branch_number") + 
  geom_point(size=3) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(size = 12)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 15)) +
  theme(axis.title.y = element_text(size = 15)) +
  theme(legend.title = element_text(size = 15)) +
  theme(legend.text = element_text(size = 12)) +
  annotate("text", label=paste("stress =",round(stress, 2),"\n p-value =",round(p, 3)), x=Inf, y=-Inf,
           hjust=1.05, vjust=-0.6) +
  theme(aspect.ratio = 1) +
  scale_color_manual(values = c("Branch_2"="#FF5765", "Branch_3"="#01DEE6",
                                "Branch_4"="#8A6FDF", "Branch_5"="#7CB53A"))
uw_unifrac_NMDS