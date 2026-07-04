library(ape)
library(phytools)
library(ggtree)
library(khroma)
library(stringr)

tree <- read.tree(file = "gtdb_alignment.faa.treefile")
tree

tree.rooted <- midpoint_root(tree = tree)
tree.rooted

plot(tree.rooted)

df <- read.delim(file = "gtdbtk.bac120.summary.tsv",
                 header = T, sep = "\t",
                 stringsAsFactors = F,
                 check.names = F)
head(df)

df$phylum <- unlist(lapply(strsplit(x = df$classification, split = ";"), "[[", 2))
sort(table(df$phylum))

phylum.count <- as.data.frame(table(df$phylum), stringsAsFactors = F)
head(phylum.count)
class(phylum.count$Var1)
phylum.count$other <- phylum.count$Var1
phylum.count$other[which(phylum.count$Freq<60)] <- "Other"
table(phylum.count$other)

phylum_other <- c()
for (i in 1:length(df$phylum)){
  test <- df$phylum[i]
  phylum_other <- c(phylum_other, phylum.count[which(phylum.count$Var1==test), "other"])
}
phylum_other
df$phylum_other <- phylum_other
head(df)

df$user_genome <- gsub(pattern = "_in", replacement = "", x = df$user_genome)
head(df)

table(df$user_genome%in%tree.rooted$tip.label)
table(tree.rooted$tip.label%in%df$user_genome)
table(df$user_genome==tree.rooted$tip.label)
table(tree.rooted$tip.label==df$user_genome)
df <- df[match(tree.rooted$tip.label, df$user_genome),]
table(df$user_genome==tree.rooted$tip.label)
table(tree.rooted$tip.label==df$user_genome)

df$phylum_other <- gsub(pattern = "p__", replacement = "", x = df$phylum_other)

rownames(df) <- df$user_genome

muted <- color("muted")
p <- ggtree(tree.rooted, layout="fan", open.angle=10) %<+% df
p + geom_tippoint(mapping=aes(color=phylum_other),
                size=0.5,
                show.legend=TRUE) + 
  scale_color_manual(values = c("Actinomycetota" = "#CC6677",
                                "Bacillota" = "#332288",
                                "Bacillota_A" =  "#882255",
                                "Bacteroidota" = "#117733",
                                "Campylobacterota" = "#DDCC77",
                                "Cyanobacteriota" = "#88CCEE",
                                "Desulfobacterota" = "#999933",
                                "Other" = "#DDDDDD",
                                "Pseudomonadota" = "#44AA99",
                                "Spirochaetota" = "#AA4499")) +
  geom_treescale(offset = 50)

ml <- read.delim(file = "model_list.tsv",
                 header = T, sep = "\t",
                 stringsAsFactors = F, check.names = F)
head(ml)

table(ml$assembly_accession%in%tree.rooted$tip.label)
table(tree.rooted$tip.label%in%ml$assembly_accession)
table(ml$assembly_accession==tree.rooted$tip.label)
table(tree.rooted$tip.label==ml$assembly_accession)
ml <- ml[match(tree.rooted$tip.label, ml$assembly_accession),]
table(ml$assembly_accession==tree.rooted$tip.label)
table(tree.rooted$tip.label==ml$assembly_accession)

table(ml$assembly_accession==df$user_genome)
table(df$user_genome==ml$assembly_accession)

ml$matchname <- unlist(lapply(strsplit(x = ml$file_path, split = "\\/"), "[[", 4))
ml$matchname <- gsub(pattern = ".xml.gz", replacement = "", x = ml$matchname)
ml$matchname2 <- str_replace_all(ml$organism_name, "[^[:alnum:]]", "_")
head(ml)

carveme <- read.delim(file = "CarveMe_organism_names.csv",
                      header = F, sep = ",",
                      stringsAsFactors = F,
                      check.names = F)
head(carveme)

agora <- read.delim(file = "AGORA_organism_names.csv",
                      header = F, sep = ",",
                      stringsAsFactors = F,
                      check.names = F) 
head(agora)

table(carveme$V1%in%agora$V1)
table(agora$V1%in%carveme$V1)

table(agora$V1%in%ml$matchname)
table(carveme$V1%in%ml$matchname)

table(agora$V1%in%ml$matchname2)
sanity <- c(ml$matchname, ml$matchname2)
table(agora$V1%in%sanity)

table(carveme$V1%in%ml$matchname)
table(agora$V1%in%sanity)
table(ml$matchname%in%carveme$V1)
table(ml$matchname%in%agora$V1)
table(ml$matchname2%in%agora$V1)

ml$dataset <- ifelse(test = ml$matchname%in%agora$V1, yes = "AGORA + CarveMe",
                     no = ifelse(test = ml$matchname2%in%agora$V1, yes = "AGORA + CarveMe",
                                 no = "CarveMe"))
table(ml$dataset)

head(ml)
head(df)
table(ml$assembly_accession%in%df$user_genome)
table(df$user_genome%in%ml$assembly_accession)
table(ml$assembly_accession==df$user_genome)
table(df$user_genome==ml$assembly_accession)
df <- cbind(df, ml)

p2 <- ggtree(tree.rooted, layout="fan", open.angle=10) %<+% df
p2 <- 
p2 + geom_tippoint(mapping=aes(color=phylum_other),
                  size=0.5,
                  show.legend=TRUE) + 
  scale_color_manual(values = c("Actinomycetota" = "#CC6677",
                                "Bacillota" = "#332288",
                                "Bacillota_A" =  "#882255",
                                "Bacteroidota" = "#117733",
                                "Campylobacterota" = "#DDCC77",
                                "Cyanobacteriota" = "#88CCEE",
                                "Desulfobacterota" = "#999933",
                                "Other" = "#DDDDDD",
                                "Pseudomonadota" = "#44AA99",
                                "Spirochaetota" = "#AA4499")) +
  geom_treescale(offset = 50) 

df.datasets <- data.frame(df$dataset)
head(df.datasets)
colnames(df.datasets) <- c("dataset")
rownames(df.datasets) <- df$user_genome
head(df.datasets)
table(df.datasets$dataset)

table(rownames(df.datasets)==tree.rooted$tip.label)
table(tree.rooted$tip.label==rownames(df.datasets))

gheatmap(p = p2, data = df.datasets, offset=0, width=0.1, color = NA, colnames = F) +
  scale_fill_manual(values = c("deeppink3", "white"))


p3 <- ggtree(tree.rooted, layout="fan", open.angle=10) %<+% df

p3 <- p3 + 
  geom_point2(aes(subset = dataset %in% c("AGORA + CarveMe")),
              color = "dodgerblue", size = 0.5)

df.phyla <- data.frame(df$phylum_other)
head(df.phyla)
colnames(df.phyla) <- c("Phylum")
rownames(df.phyla) <- df$user_genome
head(df.phyla)
table(df.phyla$Phylum)

table(rownames(df.phyla)==tree.rooted$tip.label)
table(tree.rooted$tip.label==rownames(df.phyla))

gheatmap(p = p3, data = df.phyla, offset=0, width=0.1, color = NA, colnames = F) +
  scale_fill_manual(values = c("Actinomycetota" = "#CC6677",
                                "Bacillota" = "#332288",
                                "Bacillota_A" =  "#882255",
                                "Bacteroidota" = "#117733",
                                "Campylobacterota" = "#DDCC77",
                                "Cyanobacteriota" = "#88CCEE",
                                "Desulfobacterota" = "#999933",
                                "Other" = "#DDDDDD",
                                "Pseudomonadota" = "#44AA99",
                                "Spirochaetota" = "#AA4499")) +
  geom_treescale(offset = 50) 
