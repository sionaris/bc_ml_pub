# This script is used to compare our signature with other published signatures

# Libraries #####
library(dplyr)
library(ggplot2)
library(openxlsx)
library(genefu)
library(org.Hs.eg.db)

# Annotation: files and function ######
# ID_Map
official = org.Hs.egSYMBOL
mapped_genes_official = mappedkeys(official)
official_df = as.data.frame(official[mapped_genes_official])
official_df = official_df %>% dplyr::rename(EntrezGene.ID = gene_id, Gene.Symbol = symbol)
official_df$HGNC_Official = "Yes"
official_df = official_df[-which(duplicated(official_df$Gene.Symbol)==T),]

alias = org.Hs.egALIAS2EG
mapped_genes_alias = mappedkeys(alias)
alias_df = as.data.frame(alias[mapped_genes_alias])
alias_df = alias_df %>% dplyr::rename(EntrezGene.ID = gene_id, Gene.Symbol = alias_symbol)
alias_df = alias_df[-which(alias_df$Gene.Symbol %in% official_df$Gene.Symbol),]
alias_df$HGNC_Official = "No"

ID_Map = rbind(official_df, alias_df) %>% distinct()
ID_Map$EntrezGene.ID = as.numeric(ID_Map$EntrezGene.ID)
ID_Map = ID_Map[order(ID_Map$EntrezGene.ID),] %>%
  dplyr::rename(probe=Gene.Symbol) %>%
  dplyr::select(probe, EntrezGene.ID, HGNC_Official)

# Aliases
aliases_for_join = alias_df %>% dplyr::rename(Alias = Gene.Symbol)
Aliases = official_df %>% inner_join(aliases_for_join,
                                     by = "EntrezGene.ID", multiple = "all") %>%
  dplyr::select(Alias, Gene.Symbol, EntrezGene.ID) %>%
  dplyr::rename(probe = Alias, HGNC_Symbol = Gene.Symbol,
                Entrez = EntrezGene.ID) %>%
  distinct()

rm(alias, alias_df, aliases_for_join, official, mapped_genes_alias, mapped_genes_official)

# Import our gene expression signature #####
NRS = read.xlsx("data/Output sets/Downstream/DGEA.xlsx") %>%
  dplyr::filter(adj.P.Val < 0.05) %>%
  dplyr::select(Gene.Symbol, EntrezGene.ID)

# Genefu signatures #####

# ClaudinLow
# Entrez ID's are matched to official gene symbols. Unmapped Entrez ID's are discarded
data("claudinLowData"); claudinLow = as.data.frame(list(EntrezGene.ID = claudinLowData$fnames)) %>%
  inner_join(official_df, by = "EntrezGene.ID") %>%
  dplyr::select(Gene.Symbol, EntrezGene.ID)
rm(claudinLowData)

# pam50
data("pam50"); pam50 = pam50$centroids.map %>% 
  dplyr::select(Gene.Symbol = probe.centroids, EntrezGene.ID)

# scmgene is easy to generate manually
scmgene = as.data.frame(rbind(c("ESR1", 2099), c("ERBB2", 2064), c("AURKA", 6790)))
colnames(scmgene) = c("Gene.Symbol", "EntrezGene.ID")

# scmod1
data("scmod1.robust")
scmod1_ESR1_module = as.data.frame(scmod1.robust$mod$ESR1) %>%
  inner_join(official_df, by = "EntrezGene.ID") %>%
  dplyr::select(Gene.Symbol, EntrezGene.ID)

scmod1_ERBB2_module = as.data.frame(scmod1.robust$mod$ERBB2) %>%
  inner_join(official_df, by = "EntrezGene.ID") %>%
  dplyr::select(Gene.Symbol, EntrezGene.ID)

scmod1_AURKA_module = as.data.frame(scmod1.robust$mod$AURKA) %>%
  inner_join(official_df, by = "EntrezGene.ID") %>%
  dplyr::select(Gene.Symbol, EntrezGene.ID)
rm(scmod1.robust)

# scmod2
data("scmod2.robust")
scmod2_ESR1_module = as.data.frame(scmod2.robust$mod$ESR1) %>%
  inner_join(official_df, by = "EntrezGene.ID") %>%
  dplyr::select(Gene.Symbol, EntrezGene.ID)

scmod2_ERBB2_module = as.data.frame(scmod2.robust$mod$ERBB2) %>%
  inner_join(official_df, by = "EntrezGene.ID") %>%
  dplyr::select(Gene.Symbol, EntrezGene.ID)

scmod2_AURKA_module = as.data.frame(scmod2.robust$mod$AURKA) %>%
  inner_join(official_df, by = "EntrezGene.ID") %>%
  dplyr::select(Gene.Symbol, EntrezGene.ID)
rm(scmod2.robust)

# endoPredict
data("sig.endoPredict"); endoPredict = sig.endoPredict %>% 
  dplyr::select(Gene.Symbol = symbol, EntrezGene.ID)
rm(sig.endoPredict)

# OncotypeDX
data("sig.oncotypedx"); OncotypeDX = sig.oncotypedx %>% 
  dplyr::select(Gene.Symbol = symbol, EntrezGene.ID)
rm(sig.oncotypedx)

# Mammaprint
data("sig.gene70"); Mammaprint = sig.gene70 %>%
  dplyr::select(Gene.Symbol = NCBI.gene.symbol, EntrezGene.ID) %>%
  dplyr::filter(!is.na(Gene.Symbol))
rm(sig.gene70)

# GGI
data("sig.ggi"); GGI = sig.ggi%>%
  dplyr::select(Gene.Symbol = NCBI.gene.symbol, EntrezGene.ID) %>%
  dplyr::filter(!is.na(Gene.Symbol))
rm(sig.ggi)

# PIK3CA
data("sig.pik3cags"); sig.pik3cags$EntrezGene.ID = as.character(sig.pik3cags$EntrezGene.ID)
PIK3CA = sig.pik3cags %>%
  inner_join(official_df, by = "EntrezGene.ID") %>%
  dplyr::select(Gene.Symbol, EntrezGene.ID)
rm(sig.pik3cags)

# METABRIC Integrative Clusters
library(iC10)
data("Map.Exp"); iC10 = Map.Exp %>% dplyr::rename(Gene.Symbol = Gene_symbol) %>%
  left_join(official_df, by = "Gene.Symbol") %>%
  dplyr::select(Gene.Symbol, EntrezGene.ID)
rm(Map.Exp)

# IHC4 (ESR1, ERBB2, PGR, MKI67)
IHC4 = as.data.frame(rbind(c("ESR1", 2099), c("ERBB2", 2064), 
                           c("PGR", 5241), c("MKI67", 4288)))
colnames(IHC4) = c("Gene.Symbol", "EntrezGene.ID")

# MetaGX (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6584731/)
MetaGX = read.xlsx("Signatures/MetaGX.xlsx") %>%
  dplyr::select(Gene.Symbol, EntrezGene.ID)

# Danaher (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5319024/)
# Supplementary Material
Danaher = read.xlsx("Signatures/Danaher.xlsx", sheet = 4) %>%
  dplyr::rename(Gene.Symbol = Gene) %>%
  left_join(official_df, by = "Gene.Symbol") %>%
  dplyr::select(Gene.Symbol, EntrezGene.ID)

# MCP (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5073889/)
# Genes retrieved from GitHub repo 
# (https://github.com/ebecht/MCPcounter/blob/master/Signatures/genes.txt)
MCP = read.table("Signatures/MCP.txt")
colnames(MCP) = MCP[1,]
MCP = MCP[-1,] %>% 
  dplyr::select(Gene.Symbol = `HUGO symbols`, EntrezGene.ID = ENTREZID)

# Immunophenoscore (https://pubmed.ncbi.nlm.nih.gov/28052254/)
# Supplementary Material
Immunophenoscore = read.xlsx("Signatures/Immunophenoscore.xlsx")
colnames(Immunophenoscore) = Immunophenoscore[1,]
Immunophenoscore = Immunophenoscore[-1,] %>%
  dplyr::rename(Gene.Symbol = Metagene) %>%
  left_join(official_df, by = "Gene.Symbol") %>%
  dplyr::select(Gene.Symbol, EntrezGene.ID)

# Create a signature list #####
Signatures = list(pam50 = pam50,
                  scmod1_AURKA = scmod1_AURKA_module,
                  scmod1_ERBB2_module = scmod1_ERBB2_module,
                  scmod1_ESR1_module = scmod1_ESR1_module,
                  scmod2_AURKA = scmod2_AURKA_module,
                  scmod2_ERBB2_module = scmod2_ERBB2_module,
                  scmod2_ESR1_module = scmod2_ESR1_module,
                  scmgene = scmgene,
                  IHC4 = IHC4,
                  claudinLow = claudinLow,
                  endoPredict = endoPredict,
                  OncotypeDX = OncotypeDX,
                  Mammaprint = Mammaprint,
                  GGI = GGI,
                  PIK3CA = PIK3CA,
                  iC10 = iC10,
                  MetaGX = MetaGX,
                  Danaher = Danaher,
                  MCP = MCP,
                  Immunophenoscore = Immunophenoscore)

# Cleaning up
rm(pam50, scmod1_AURKA_module, scmod1_ERBB2_module, scmod1_ESR1_module,
   scmod2_AURKA_module, scmod2_ERBB2_module, scmod2_ESR1_module, scmgene,
   IHC4, claudinLow, endoPredict, OncotypeDX, Mammaprint, GGI, PIK3CA, iC10,
   MetaGX, Danaher, MCP, Immunophenoscore)

# Overlaps with our signature (NRS) #####

# Preparing output
GS_overlap = as.data.frame(list(NRSGS = NRS$Gene.Symbol))
Entrez_overlap = as.data.frame(list(NRSE = NRS$EntrezGene.ID))

# Overlaps
for (i in 1:length(Signatures)){
  GS_overlap[,i+1] = "No"
  Entrez_overlap[,i+1] = "No"
  yes_ind_gs = which(NRS$Gene.Symbol %in% Signatures[[i]]$Gene.Symbol)
  yes_ind_entrez = which(NRS$EntrezGene.ID %in% Signatures[[i]]$EntrezGene.ID)
  GS_overlap[yes_ind_gs, i+1] = "Yes"
  Entrez_overlap[yes_ind_entrez, i+1] = "Yes"
  rm(yes_ind_entrez, yes_ind_gs)
}

rm(i); gc()
colnames(GS_overlap) = colnames(Entrez_overlap) = c("NRS", names(Signatures))

# Add total overlaps in the last row of both output sets
GS_overlap[nrow(GS_overlap) + 1,] = c("Total", 
                                     lapply(GS_overlap[,2:ncol(GS_overlap)], function(x) sum(x == "Yes")))
Entrez_overlap[nrow(Entrez_overlap) + 1,] = c("Total", 
                                     lapply(Entrez_overlap[,2:ncol(GS_overlap)], function(x) sum(x == "Yes")))

# Heatmap of NRS vs signatures
# Load required libraries
library(ComplexHeatmap)
library(circlize)
library(rcartocolor)

# Convert your dataset's relevant columns to a matrix
overlap_mat <- GS_overlap[1:166,2:ncol(GS_overlap)]

# Convert 'Yes' and 'No' to numeric values
overlap_mat[overlap_mat == 'Yes'] <- 1
overlap_mat[overlap_mat == 'No'] <- 0
overlap_mat[, colnames(overlap_mat)] = lapply(overlap_mat[, colnames(overlap_mat)], as.numeric)
overlap_mat = as.matrix(overlap_mat)

# Convert the matrix back to numeric
overlap_mat <- apply(overlap_mat, 2, as.numeric)
rownames(overlap_mat) <- GS_overlap$NRS[1:166]
overlap_mat <- t(overlap_mat)  # Transpose matrix before setting row names

# Generate the categories for each row according to your dataset
categories_vector <- c(rep("Breast cancer", 16), "Multi-cancer", rep("Immune system", 3))

# Transform the categories vector into a factor
categories_factor <- factor(categories_vector, levels = c("Breast cancer", "Multi-cancer", "Immune system"))

# Generate annotation dataframe
annotation_df_overlap <- data.frame(Annotations = categories_factor)
rownames(annotation_df_overlap) <- rownames(overlap_mat)  # Adjust row names of the annotation to match the matrix

# Generate annotation colors for heatmap
annotation_colors_overlap <- list(Annotations = c("Breast cancer" = rcartocolor::carto_pal(n = 7, name = "Burg")[6], 
                                                  "Multi-cancer" = rcartocolor::carto_pal(n = 7, name = "Teal")[7], 
                                                  "Immune system" = rcartocolor::carto_pal(n = 7, name = "Peach")[6]))

# Configure the color scheme and annotation
col_fun <- colorRamp2(c(0, 1), c("grey95", "dodgerblue4"))

# Heatmap annotation
ha = rowAnnotation(df = annotation_df_overlap, 
                   col = annotation_colors_overlap, 
                   show_annotation_name = FALSE)

# Generate the heatmap
heatmap <- Heatmap(overlap_mat,
                   name = "Overlaps",
                   show_row_names = TRUE,
                   show_column_names = TRUE,
                   cluster_rows = FALSE,
                   cluster_columns = FALSE,
                   row_title = "Signatures",
                   row_title_gp = gpar(fontface = "bold", fontsize = 15),
                   column_title = "NRS genes",
                   column_title_gp = gpar(fontface = "bold", fontsize = 15),
                   column_title_side = "bottom",
                   heatmap_legend_param = list(legend_direction = "horizontal",
                                               legend_width = unit(4, "cm"),
                                               color_bar = "discrete",
                                               at = c(0, 1),
                                               labels = c("Not overlapping", "Overlapping")),
                   col = col_fun,
                   border = "black",
                   row_names_side = "left",
                   row_names_gp = gpar(fontsize = 8),
                   column_names_gp = gpar(fontsize = 4),
                   left_annotation = ha)
png("Signatures/NRS_published_signatures_heatmap.png", res = 600, height = 3300, width = 6750)
plot.new()
heatmap
# Draw the heatmap
# draw(heatmap, heatmap_legend_side = "bottom")
dev.off()

# Venn diagrams #####
library(ggVennDiagram)
GS_Venn = c(list(NRS = NRS$Gene.Symbol), 
            lapply(Signatures, function(x) {x <- x %>% dplyr::select(Gene.Symbol)
            x <- as.character(unlist(x))}))
Entrez_Venn = c(list(NRS = NRS$EntrezGene.ID), 
                lapply(Signatures, function(x) {x <- x %>% dplyr::select(EntrezGene.ID)
                x <- as.character(unlist(x))}))

# Venn diagram list
GS_Venn_diagram_list = list()
Entrez_Venn_diagram_list = list()

for(i in 2:length(GS_Venn)){
  diagram = ggVennDiagram(
    GS_Venn[c(1,i)], label_alpha = 0,
    category.names = names(GS_Venn[c(1,i)]), set_size = 0.65, label_size = 1) +
    scale_color_manual(values = c("grey10", "grey15")) +
    scale_fill_gradient(low = "white", high = "grey30") + 
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 2.5),
          panel.background = element_rect(fill = "white", 
                                          colour = "white"),
          legend.key.size = unit(1.5, units = "mm"),
          legend.text = element_text(size = 2),
          legend.title = element_text(face = "bold", size = 3),
          legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.background = element_blank())+
    labs(title = paste0("Venn diagram of NRS and the ", names(GS_Venn)[i], " signature"))
  GS_Venn_diagram_list = c(GS_Venn_diagram_list, list(diagram)); rm(diagram)
  
  diagram = ggVennDiagram(
    Entrez_Venn[c(1,i)], label_alpha = 0,
    category.names = names(Entrez_Venn[c(1,i)]), set_size = 0.65, label_size = 1) +
    scale_color_manual(values = c("grey10", "grey15")) +
    scale_fill_gradient(low = "white", high = "grey30") + 
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 2.5),
          panel.background = element_rect(fill = "white", 
                                          colour = "white"),
          legend.key.size = unit(1.5, units = "mm"),
          legend.text = element_text(size = 2),
          legend.title = element_text(face = "bold", size = 3),
          legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.background = element_blank())+
    labs(title = paste0("Venn diagram of NRS and the ", names(Entrez_Venn)[i], " signature"))
  Entrez_Venn_diagram_list = c(Entrez_Venn_diagram_list, list(diagram)); rm(diagram)
}
rm(i); gc()
names(Entrez_Venn_diagram_list) = names(GS_Venn_diagram_list) = names(Signatures)

# Defining a multiplot function #####
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# Multiplot of Gene Symbol Venn diagrams
tiff("Signatures/Gene_Symbol_Venn_multiplot.tiff", 
     width = 3840, height = 2160, res = 700, compression = "lzw")
multiplot(GS_Venn_diagram_list[[1]], GS_Venn_diagram_list[[2]], GS_Venn_diagram_list[[3]],
          GS_Venn_diagram_list[[4]], GS_Venn_diagram_list[[5]], GS_Venn_diagram_list[[6]],
          GS_Venn_diagram_list[[7]], GS_Venn_diagram_list[[8]], GS_Venn_diagram_list[[9]],
          GS_Venn_diagram_list[[10]], GS_Venn_diagram_list[[11]], GS_Venn_diagram_list[[12]],
          GS_Venn_diagram_list[[13]], GS_Venn_diagram_list[[14]], GS_Venn_diagram_list[[15]],
          GS_Venn_diagram_list[[16]], GS_Venn_diagram_list[[17]], GS_Venn_diagram_list[[18]],
          GS_Venn_diagram_list[[19]], GS_Venn_diagram_list[[20]], cols = 4)
m = ggplot(multiplot(GS_Venn_diagram_list[[1]], GS_Venn_diagram_list[[2]], GS_Venn_diagram_list[[3]],
                     GS_Venn_diagram_list[[4]], GS_Venn_diagram_list[[5]], GS_Venn_diagram_list[[6]],
                     GS_Venn_diagram_list[[7]], GS_Venn_diagram_list[[8]], GS_Venn_diagram_list[[9]],
                     GS_Venn_diagram_list[[10]], GS_Venn_diagram_list[[11]], GS_Venn_diagram_list[[12]],
                     GS_Venn_diagram_list[[13]], GS_Venn_diagram_list[[14]], GS_Venn_diagram_list[[15]],
                     GS_Venn_diagram_list[[16]], GS_Venn_diagram_list[[17]], GS_Venn_diagram_list[[18]],
                     GS_Venn_diagram_list[[19]], GS_Venn_diagram_list[[20]], cols = 4))
dev.off(); rm(m)

# Multiplot of Entrez ID Venn diagrams
tiff("Signatures/Entrez_Venn_multiplot.tiff", 
     width = 3840, height = 2160, res = 700, compression = "lzw")
multiplot(Entrez_Venn_diagram_list[[1]], Entrez_Venn_diagram_list[[2]], Entrez_Venn_diagram_list[[3]],
          Entrez_Venn_diagram_list[[4]], Entrez_Venn_diagram_list[[5]], Entrez_Venn_diagram_list[[6]],
          Entrez_Venn_diagram_list[[7]], Entrez_Venn_diagram_list[[8]], Entrez_Venn_diagram_list[[9]],
          Entrez_Venn_diagram_list[[10]], Entrez_Venn_diagram_list[[11]], Entrez_Venn_diagram_list[[12]],
          Entrez_Venn_diagram_list[[13]], Entrez_Venn_diagram_list[[14]], Entrez_Venn_diagram_list[[15]],
          Entrez_Venn_diagram_list[[16]], Entrez_Venn_diagram_list[[17]], Entrez_Venn_diagram_list[[18]],
          Entrez_Venn_diagram_list[[19]], Entrez_Venn_diagram_list[[20]], cols = 4)
m = ggplot(multiplot(Entrez_Venn_diagram_list[[1]], Entrez_Venn_diagram_list[[2]], Entrez_Venn_diagram_list[[3]],
                     Entrez_Venn_diagram_list[[4]], Entrez_Venn_diagram_list[[5]], Entrez_Venn_diagram_list[[6]],
                     Entrez_Venn_diagram_list[[7]], Entrez_Venn_diagram_list[[8]], Entrez_Venn_diagram_list[[9]],
                     Entrez_Venn_diagram_list[[10]], Entrez_Venn_diagram_list[[11]], Entrez_Venn_diagram_list[[12]],
                     Entrez_Venn_diagram_list[[13]], Entrez_Venn_diagram_list[[14]], Entrez_Venn_diagram_list[[15]],
                     Entrez_Venn_diagram_list[[16]], Entrez_Venn_diagram_list[[17]], Entrez_Venn_diagram_list[[18]],
                     Entrez_Venn_diagram_list[[19]], Entrez_Venn_diagram_list[[20]], cols = 4))
dev.off(); rm(m)

# Venn diagrams of the intersections of signatures and our background genes (5,672) #####
background_genes = read.xlsx("data/Output sets/Downstream/DGEA.xlsx") %>%
  dplyr::select(Gene.Symbol, EntrezGene.ID)

GS_Venn_intersections = c(list(NRS = NRS$Gene.Symbol), 
                          lapply(Signatures, function(x) {x <- x$Gene.Symbol
                          x <- intersect(as.character(x), background_genes$Gene.Symbol)
                          x <- as.character(unlist(x))}))
Entrez_Venn_intersections = c(list(NRS = NRS$EntrezGene.ID), 
                              lapply(Signatures, function(x) {x <- x $EntrezGene.ID
                              x <- intersect(as.character(x), background_genes$EntrezGene.ID)
                              x <- as.character(unlist(x))}))

# Venn diagram list
GS_Venn_intersections_diagram_list = list()
Entrez_Venn_intersections_diagram_list = list()

for(i in 2:length(GS_Venn_intersections)){
  diagram = ggVennDiagram(
    GS_Venn_intersections[c(1,i)], label_alpha = 0,
    category.names = names(GS_Venn_intersections[c(1,i)]), set_size = 0.65, label_size = 1) +
    scale_color_manual(values = c("grey10", "grey15")) +
    scale_fill_gradient(low = "white", high = "grey30") + 
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 2.5),
          panel.background = element_rect(fill = "white", 
                                          colour = "white"),
          legend.key.size = unit(1.5, units = "mm"),
          legend.text = element_text(size = 2),
          legend.title = element_text(face = "bold", size = 3),
          legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.background = element_blank())+
    labs(title = paste0("Venn diagram of NRS and the ", names(GS_Venn_intersections)[i], " signature"))
  GS_Venn_intersections_diagram_list = c(GS_Venn_intersections_diagram_list, list(diagram)); rm(diagram)
  
  diagram = ggVennDiagram(
    Entrez_Venn_intersections[c(1,i)], label_alpha = 0,
    category.names = names(Entrez_Venn_intersections[c(1,i)]), set_size = 0.65, label_size = 1) +
    scale_color_manual(values = c("grey10", "grey15")) +
    scale_fill_gradient(low = "white", high = "grey30") + 
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 2.5),
          panel.background = element_rect(fill = "white", 
                                          colour = "white"),
          legend.key.size = unit(1.5, units = "mm"),
          legend.text = element_text(size = 2),
          legend.title = element_text(face = "bold", size = 3),
          legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "mm"),
          legend.spacing.y = unit(0.5, units = "mm"),
          legend.spacing.x = unit(0.5, units = "mm"),
          legend.background = element_blank())+
    labs(title = paste0("Venn diagram of NRS and the ", names(Entrez_Venn_intersections)[i], " signature"))
  Entrez_Venn_intersections_diagram_list = c(Entrez_Venn_intersections_diagram_list, list(diagram)); rm(diagram)
}
rm(i); gc()
names(Entrez_Venn_intersections_diagram_list) = names(GS_Venn_intersections_diagram_list) = names(Signatures)

# Multiplot of Gene Symbol Venn diagrams
tiff("Signatures/Gene_Symbol_Venn_intersections_multiplot.tiff", 
     width = 3840, height = 2160, res = 700, compression = "lzw")
multiplot(GS_Venn_intersections_diagram_list[[1]], GS_Venn_intersections_diagram_list[[2]], GS_Venn_intersections_diagram_list[[3]],
          GS_Venn_intersections_diagram_list[[4]], GS_Venn_intersections_diagram_list[[5]], GS_Venn_intersections_diagram_list[[6]],
          GS_Venn_intersections_diagram_list[[7]], GS_Venn_intersections_diagram_list[[8]], GS_Venn_intersections_diagram_list[[9]],
          GS_Venn_intersections_diagram_list[[10]], GS_Venn_intersections_diagram_list[[11]], GS_Venn_intersections_diagram_list[[12]],
          GS_Venn_intersections_diagram_list[[13]], GS_Venn_intersections_diagram_list[[14]], GS_Venn_intersections_diagram_list[[15]],
          GS_Venn_intersections_diagram_list[[16]], GS_Venn_intersections_diagram_list[[17]], GS_Venn_intersections_diagram_list[[18]],
          GS_Venn_intersections_diagram_list[[19]], GS_Venn_intersections_diagram_list[[20]], cols = 4)
m = ggplot(multiplot(GS_Venn_intersections_diagram_list[[1]], GS_Venn_intersections_diagram_list[[2]], GS_Venn_intersections_diagram_list[[3]],
                     GS_Venn_intersections_diagram_list[[4]], GS_Venn_intersections_diagram_list[[5]], GS_Venn_intersections_diagram_list[[6]],
                     GS_Venn_intersections_diagram_list[[7]], GS_Venn_intersections_diagram_list[[8]], GS_Venn_intersections_diagram_list[[9]],
                     GS_Venn_intersections_diagram_list[[10]], GS_Venn_intersections_diagram_list[[11]], GS_Venn_intersections_diagram_list[[12]],
                     GS_Venn_intersections_diagram_list[[13]], GS_Venn_intersections_diagram_list[[14]], GS_Venn_intersections_diagram_list[[15]],
                     GS_Venn_intersections_diagram_list[[16]], GS_Venn_intersections_diagram_list[[17]], GS_Venn_intersections_diagram_list[[18]],
                     GS_Venn_intersections_diagram_list[[19]], GS_Venn_intersections_diagram_list[[20]], cols = 4))
dev.off(); rm(m)

# Multiplot of Entrez ID Venn diagrams
tiff("Signatures/Entrez_Venn_intersections_multiplot.tiff", 
     width = 3840, height = 2160, res = 700, compression = "lzw")
multiplot(Entrez_Venn_intersections_diagram_list[[1]], Entrez_Venn_intersections_diagram_list[[2]], Entrez_Venn_intersections_diagram_list[[3]],
          Entrez_Venn_intersections_diagram_list[[4]], Entrez_Venn_intersections_diagram_list[[5]], Entrez_Venn_intersections_diagram_list[[6]],
          Entrez_Venn_intersections_diagram_list[[7]], Entrez_Venn_intersections_diagram_list[[8]], Entrez_Venn_intersections_diagram_list[[9]],
          Entrez_Venn_intersections_diagram_list[[10]], Entrez_Venn_intersections_diagram_list[[11]], Entrez_Venn_intersections_diagram_list[[12]],
          Entrez_Venn_intersections_diagram_list[[13]], Entrez_Venn_intersections_diagram_list[[14]], Entrez_Venn_intersections_diagram_list[[15]],
          Entrez_Venn_intersections_diagram_list[[16]], Entrez_Venn_intersections_diagram_list[[17]], Entrez_Venn_intersections_diagram_list[[18]],
          Entrez_Venn_intersections_diagram_list[[19]], Entrez_Venn_intersections_diagram_list[[20]], cols = 4)
m = ggplot(multiplot(Entrez_Venn_intersections_diagram_list[[1]], Entrez_Venn_intersections_diagram_list[[2]], Entrez_Venn_intersections_diagram_list[[3]],
                     Entrez_Venn_intersections_diagram_list[[4]], Entrez_Venn_intersections_diagram_list[[5]], Entrez_Venn_intersections_diagram_list[[6]],
                     Entrez_Venn_intersections_diagram_list[[7]], Entrez_Venn_intersections_diagram_list[[8]], Entrez_Venn_intersections_diagram_list[[9]],
                     Entrez_Venn_intersections_diagram_list[[10]], Entrez_Venn_intersections_diagram_list[[11]], Entrez_Venn_intersections_diagram_list[[12]],
                     Entrez_Venn_intersections_diagram_list[[13]], Entrez_Venn_intersections_diagram_list[[14]], Entrez_Venn_intersections_diagram_list[[15]],
                     Entrez_Venn_intersections_diagram_list[[16]], Entrez_Venn_intersections_diagram_list[[17]], Entrez_Venn_intersections_diagram_list[[18]],
                     Entrez_Venn_intersections_diagram_list[[19]], Entrez_Venn_intersections_diagram_list[[20]], cols = 4))
dev.off(); rm(m)

# Venn diagrams of selected lists only #####
# Full lists
# Venn diagrams of the intersections of signatures and our background genes (5,672) #####
GS_Venn_selected = c(list(NRS = NRS$Gene.Symbol), 
                                   lapply(Signatures[c("pam50", "Mammaprint", 
                                                       "iC10", "Immunophenoscore")], function(x) {x <- x$Gene.Symbol
                                                       x <- as.character(x)}))
Entrez_Venn_selected = c(list(NRS = NRS$EntrezGene.ID), 
                                       lapply(Signatures[c("pam50", "Mammaprint", 
                                                           "iC10", "Immunophenoscore")], function(x) {x <- x $EntrezGene.ID
                                                           x <- as.character(x)}))

# Venn diagram list
GS_Venn_selected_diagram_list = list()
Entrez_Venn_selected_diagram_list = list()

for(i in 2:length(GS_Venn_selected)){
  diagram = ggVennDiagram(
    GS_Venn_selected[c(1,i)], label_alpha = 0,
    category.names = names(GS_Venn_selected[c(1,i)]),
    label = "count", set_size = 2, label_size = 3) +
    scale_color_manual(values = c("grey10", "grey15")) +
    scale_fill_gradient(low = "white", high = "grey30") + 
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 6),
          panel.background = element_rect(fill = "white", 
                                          colour = "white"),
          legend.key.size = unit(3, units = "mm"),
          legend.text = element_text(size = 4),
          legend.title = element_text(face = "bold", size = 6),
          legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "mm"),
          legend.spacing.y = unit(2, units = "mm"),
          legend.spacing.x = unit(2, units = "mm"),
          legend.background = element_blank())+
    labs(title = paste0("Venn diagram of NRS and the ", names(GS_Venn_selected)[i], " signature"))
  GS_Venn_selected_diagram_list = c(GS_Venn_selected_diagram_list, list(diagram)); rm(diagram)
  
  diagram = ggVennDiagram(
    Entrez_Venn_selected[c(1,i)], label_alpha = 0,
    category.names = names(Entrez_Venn_selected[c(1,i)]),
    label = "count", set_size = 2, label_size = 3) +
    scale_color_manual(values = c("grey10", "grey15")) +
    scale_fill_gradient(low = "white", high = "grey30") + 
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 6),
          panel.background = element_rect(fill = "white", 
                                          colour = "white"),
          legend.key.size = unit(3, units = "mm"),
          legend.text = element_text(size = 4),
          legend.title = element_text(face = "bold", size = 6),
          legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "mm"),
          legend.spacing.y = unit(2, units = "mm"),
          legend.spacing.x = unit(2, units = "mm"),
          legend.background = element_blank())+
    labs(title = paste0("Venn diagram of NRS and the ", names(Entrez_Venn_selected)[i], " signature"))
  Entrez_Venn_selected_diagram_list = c(Entrez_Venn_selected_diagram_list, list(diagram)); rm(diagram)
}
rm(i); gc()
names(Entrez_Venn_selected_diagram_list) = names(GS_Venn_selected_diagram_list) =
  names(Signatures[c("pam50", "Mammaprint", 
                     "iC10", "Immunophenoscore")])

# Multiplot of Gene Symbol Venn diagrams
tiff("Signatures/Gene_Symbol_Venn_selected_multiplot.tiff", 
     width = 3840, height = 2160, res = 700, compression ="lzw")
multiplot(GS_Venn_selected_diagram_list[[1]], GS_Venn_selected_diagram_list[[2]], 
          GS_Venn_selected_diagram_list[[3]], GS_Venn_selected_diagram_list[[4]], 
          cols = 2)
m = ggplot(multiplot(GS_Venn_selected_diagram_list[[1]], 
                     GS_Venn_selected_diagram_list[[2]],
                     GS_Venn_selected_diagram_list[[3]], 
                     GS_Venn_selected_diagram_list[[4]], cols = 2))
dev.off(); rm(m)

# Multiplot of Entrez ID Venn diagrams
tiff("Signatures/Entrez_Venn_selected_multiplot.tiff", 
     width = 3840, height = 2160, res = 700, compression = "lzw")
multiplot(Entrez_Venn_selected_diagram_list[[1]], Entrez_Venn_selected_diagram_list[[2]],
          Entrez_Venn_selected_diagram_list[[3]], 
          Entrez_Venn_selected_diagram_list[[4]], cols = 2)
m = ggplot(multiplot(Entrez_Venn_selected_diagram_list[[1]], Entrez_Venn_selected_diagram_list[[2]],
                     Entrez_Venn_selected_diagram_list[[3]], 
                     Entrez_Venn_selected_diagram_list[[4]], cols = 2))
dev.off(); rm(m)

# Intersections
# Venn diagrams of the intersections of signatures and our background genes (5,672) #####
GS_Venn_selected_intersections = c(list(NRS = NRS$Gene.Symbol), 
                                   lapply(Signatures[c("pam50", "Mammaprint", 
                                                       "iC10", "Immunophenoscore")], function(x) {x <- x$Gene.Symbol
                                                       x <- intersect(as.character(x), background_genes$Gene.Symbol)
                                                       x <- as.character(unlist(x))}))
Entrez_Venn_selected_intersections = c(list(NRS = NRS$EntrezGene.ID), 
                                       lapply(Signatures[c("pam50", "Mammaprint", 
                                                           "iC10", "Immunophenoscore")], function(x) {x <- x $EntrezGene.ID
                                                           x <- intersect(as.character(x), background_genes$EntrezGene.ID)
                                                           x <- as.character(unlist(x))}))

# Venn diagram list
GS_Venn_selected_intersections_diagram_list = list()
Entrez_Venn_selected_intersections_diagram_list = list()

for(i in 2:length(GS_Venn_selected_intersections)){
  diagram = ggVennDiagram(
    GS_Venn_selected_intersections[c(1,i)], label_alpha = 0,
    category.names = names(GS_Venn_selected_intersections[c(1,i)]),
    label = "count", set_size = 2, label_size = 3) +
    scale_color_manual(values = c("grey10", "grey15")) +
    scale_fill_gradient(low = "white", high = "grey30") + 
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 6),
          panel.background = element_rect(fill = "white", 
                                          colour = "white"),
          legend.key.size = unit(3, units = "mm"),
          legend.text = element_text(size = 4),
          legend.title = element_text(face = "bold", size = 6),
          legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "mm"),
          legend.spacing.y = unit(2, units = "mm"),
          legend.spacing.x = unit(2, units = "mm"),
          legend.background = element_blank())+
    labs(title = paste0("Venn diagram of NRS and the ", names(GS_Venn_selected_intersections)[i], " signature"))
  GS_Venn_selected_intersections_diagram_list = c(GS_Venn_selected_intersections_diagram_list, list(diagram)); rm(diagram)
  
  diagram = ggVennDiagram(
    Entrez_Venn_selected_intersections[c(1,i)], label_alpha = 0,
    category.names = names(Entrez_Venn_selected_intersections[c(1,i)]),
    label = "count", set_size = 2, label_size = 3) +
    scale_color_manual(values = c("grey10", "grey15")) +
    scale_fill_gradient(low = "white", high = "grey30") + 
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 6),
          panel.background = element_rect(fill = "white", 
                                          colour = "white"),
          legend.key.size = unit(3, units = "mm"),
          legend.text = element_text(size = 4),
          legend.title = element_text(face = "bold", size = 6),
          legend.margin = ggplot2::margin(0, 0, 0, 0, unit = "mm"),
          legend.spacing.y = unit(2, units = "mm"),
          legend.spacing.x = unit(2, units = "mm"),
          legend.background = element_blank())+
    labs(title = paste0("Venn diagram of NRS and the ", names(Entrez_Venn_selected_intersections)[i], " signature"))
  Entrez_Venn_selected_intersections_diagram_list = c(Entrez_Venn_selected_intersections_diagram_list, list(diagram)); rm(diagram)
}
rm(i); gc()
names(Entrez_Venn_selected_intersections_diagram_list) = names(GS_Venn_selected_intersections_diagram_list) =
  names(Signatures[c("pam50", "Mammaprint", 
                     "iC10", "Immunophenoscore")])

# Multiplot of Gene Symbol Venn diagrams
tiff("Signatures/Gene_Symbol_Venn_selected_intersections_multiplot.tiff", 
     width = 3840, height = 2160, res = 700, compression = "lzw")
multiplot(GS_Venn_selected_intersections_diagram_list[[1]], 
          GS_Venn_selected_intersections_diagram_list[[2]], 
          GS_Venn_selected_intersections_diagram_list[[3]], 
          GS_Venn_selected_intersections_diagram_list[[4]], cols = 2)
m = ggplot(multiplot(GS_Venn_selected_intersections_diagram_list[[1]], 
                     GS_Venn_selected_intersections_diagram_list[[2]], 
                     GS_Venn_selected_intersections_diagram_list[[3]], 
                     GS_Venn_selected_intersections_diagram_list[[4]], cols = 2))
dev.off(); rm(m)

# Multiplot of Entrez ID Venn diagrams
tiff("Signatures/Entrez_Venn_selected_intersections_multiplot.tiff", 
     width = 3840, height = 2160, res = 700, compression = "lzw")
multiplot(Entrez_Venn_selected_intersections_diagram_list[[1]], 
          Entrez_Venn_selected_intersections_diagram_list[[2]], 
          Entrez_Venn_selected_intersections_diagram_list[[3]],
          Entrez_Venn_selected_intersections_diagram_list[[4]], cols = 2)
m = ggplot(multiplot(Entrez_Venn_selected_intersections_diagram_list[[1]],
                     Entrez_Venn_selected_intersections_diagram_list[[2]],
                     Entrez_Venn_selected_intersections_diagram_list[[3]],
                     Entrez_Venn_selected_intersections_diagram_list[[4]],  cols = 2))
dev.off(); rm(m)

# Immunophenoscore heatmap #####
library(pheatmap)
save_pheatmap_tiff = function(x, filename, width = 174, height = 120.8333, res = 650,
                              units = "mm") {
  tiff(filename, width = width, height = height, res = res, units = units, 
       compression = "lzw")
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# Import RDS objects
annotation_for_heatmap = readRDS("Signatures/annotation_for_heatmap.rds")
z_exprs = readRDS("Signatures/normalised_expression.rds")

# Samples ordered by response: Non-responders first
response_ordered = rownames(annotation_for_heatmap[order(annotation_for_heatmap$Response, 
                                                         decreasing = TRUE),])

# Genes ordered by Immunophenoscore group
Immunophenoscore = read.xlsx("Signatures/Immunophenoscore.xlsx")[-1,]
colnames(Immunophenoscore) = c("Gene.Symbol", "Category", "Immunity")
Immunophenoscore = Immunophenoscore[Immunophenoscore$Gene.Symbol %in% NRS$Gene.Symbol,] %>%
  dplyr::arrange(Category)
Immunophenoscore$Category = gsub("Effector memeory CD4 T cell", 
                                 "Effector memory CD4 T cell",
                                 Immunophenoscore$Category)

mock = NRS; rownames(mock) = mock$Gene.Symbol
genes_ordered = mock[Immunophenoscore$Gene.Symbol, ]
rm(mock)

plot_matrix = z_exprs[which(rownames(z_exprs) %in% genes_ordered$EntrezGene.ID), 
                      response_ordered]
rownames(plot_matrix) = Immunophenoscore$Gene.Symbol
plot_matrix = plot_matrix[genes_ordered$Gene.Symbol,]
plot_matrix = t(plot_matrix)

hmcol = rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdBu"))(255))
immune_categories = unlist(RColorBrewer::brewer.pal(10, "Paired"))
names(immune_categories) = unique(Immunophenoscore$Category)

annotation_for_heatmap = annotation_for_heatmap %>%
  dplyr::select(-scmod1, -IC10, -Mammaprint_risk, -rorS_risk, -Dataset)
ann_colors <- list(
  Response = c(Responder = "dodgerblue4", Non_responder = "deeppink4"),
  Timepoint = c(T1 = "goldenrod2", T2 = "purple4"),#, T2.5 = "green", T3 = "blue"),
  Treatment = c(Chemotherapy = "paleturquoise3", Endocrine_treatment = "sienna2"),
  pam50 = c(Basal = "red4", LumB = "skyblue",
            LumA = "darkblue", Her2 = "purple4", Normal = "lightgreen"),
  Cluster = c(`Cluster 1` = "orange", `Cluster 2` = "skyblue"),
  Category = immune_categories,
  Immunity = c(Innate = "cadetblue", Adaptive = "darkmagenta")
)

annotation_for_heatmap = annotation_for_heatmap %>%
  dplyr::select(Response, Cluster, Treatment, Timepoint, pam50)
rownames(Immunophenoscore) = Immunophenoscore$Gene.Symbol
Immunophenoscore = Immunophenoscore %>% dplyr::select(-Gene.Symbol)
heatmap = pheatmap(plot_matrix, col = (hmcol),
                   breaks = seq(-2, 2, 4/256),
                   annotation_col = Immunophenoscore,
                   annotation_row = annotation_for_heatmap,
                   annotation_colors = ann_colors,
                   cluster_cols = F,
                   cluster_rows = F,
                   gaps_row = 444,
                   treeheight_col =  0,
                   legend = TRUE,
                   show_rownames = F,
                   show_colnames = T,
                   fontsize = 5,
                   fontsize_col = 5,
                   legend_breaks = c(-2, # min(plot_matrix, na.rm = TRUE), 
                                     2), # max(plot_matrix, na.rm = TRUE)), 
                   legend_labels = (c("lower expression", "higher expression")),
                   main = "Immunophenoscore overlap with the NRS")
save_pheatmap_tiff(heatmap, "Signatures/Immunophenoscore_heatmap.tiff")

# Write out overlaps #####
write.xlsx(GS_overlap, "Signatures/Gene_Symbol_Signature_Overlaps.xlsx")
write.xlsx(Entrez_overlap, "Signatures/Entrez_Signature_Overlaps.xlsx")

# MCPcounter #####
library(MCPcounter)
library(openxlsx)
library(tidyverse)
library(broom)

trainval = read.xlsx("data/Output sets/Pheno.xlsx", sheet = 2)
responders_T1 = trainval$Sample.ID[trainval$Response == "Responder" &
                                     trainval$Timepoint_coded == "T1"]
responders_T2 = trainval$Sample.ID[trainval$Response == "Responder" &
                                     trainval$Timepoint_coded == "T2"]
non_responders_T1 = trainval$Sample.ID[trainval$Response == "Non_responder" &
                                         trainval$Timepoint_coded == "T1"]
non_responders_T2 = trainval$Sample.ID[trainval$Response == "Non_responder" &
                                         trainval$Timepoint_coded == "T2"]

MCP = list()
MCP[["Full"]] = as.data.frame(MCPcounter.estimate(z_exprs, featuresType = "ENTREZ_ID"))
# Found no markers for population(s): CD8 T cells, NK cells
MCP[["Responders"]] = MCP[["Full"]][, colnames(MCP[["Full"]]) %in% c(responders_T1, responders_T2)]
MCP[["Non_responders"]] = MCP[["Full"]][, colnames(MCP[["Full"]]) %in% c(non_responders_T1, non_responders_T2)]
MCP[["Responders_T1"]] = MCP[["Full"]][, colnames(MCP[["Full"]]) %in% responders_T1]
MCP[["Non_responders_T1"]] = MCP[["Full"]][, colnames(MCP[["Full"]]) %in% non_responders_T1]
MCP[["Responders_T2"]] = MCP[["Full"]][, colnames(MCP[["Full"]]) %in% responders_T2]
MCP[["Non_responders_T2"]] = MCP[["Full"]][, colnames(MCP[["Full"]]) %in% non_responders_T2]

# Calculate mean scores for each cell type in each dataset
MCP <- lapply(MCP, function(df) {
  df$mean_row <- rowMeans(df, na.rm = TRUE)
  return(df)
})

# Boxplots
# Converting the MCP dataframe into a long format
MCP_long <- lapply(MCP, function(df) {
  df$cell_type <- rownames(df)
  df_long <- df %>%
    tidyr::pivot_longer(-cell_type, names_to = "sample_id", values_to = "value")
  return(df_long)
})

# Assign response and timepoint groups
for (i in 1:length(MCP_long)) {
  MCP_long[[i]]$Response <- factor(ifelse(MCP_long[[i]]$sample_id %in% c(responders_T1, responders_T2),
                                                "Responder", "Non_responder"))
  MCP_long[[i]]$Timepoint <- factor(ifelse(MCP_long[[i]]$sample_id %in% c(responders_T1, non_responders_T1),
                                                 "T1", "T2"))
}

# Titles for the boxplots
titles = c("Responders vs. Non-responders",
            "Responders_T1 vs. Non_responders_T1",
            "Responders_T2 vs. Non_responders_T2",
            "Responders_T1 vs. Responders_T2",
            "Non_responders_T1 vs. Non_responders_T2")

# Create boxplots in a loop
boxplots = list(RespvsNonresp = list(), RespT1vsNonrespT1 = list(),
                RespT2vsNonrespT2 = list(), RespT1vsRespT2 = list(),
                NonrespT1vsNonrespT2 = list())
box_datasets = list(MCP_long[["Full"]], 
                    rbind(MCP_long[["Responders_T1"]], MCP_long[["Non_responders_T1"]]),
                    rbind(MCP_long[["Responders_T2"]], MCP_long[["Non_responders_T2"]]),
                    rbind(MCP_long[["Responders_T1"]], MCP_long[["Responders_T2"]]),
                    rbind(MCP_long[["Non_responders_T1"]], MCP_long[["Non_responders_T2"]]))
x_variables = c("Response", "Response","Response", "Timepoint", "Timepoint")
cell_types = sort(unique(MCP_long[["Full"]]$cell_type))

# Function to conduct t-tests between timepoints for each group and each cell type
perform_t_tests_timepoints <- function(data, filter) {
  result = data[filter, ] %>%
    dplyr::group_by(cell_type) %>%
    dplyr::do(tidy(t.test(value ~ Timepoint, data = .)))
  return(result)
}

perform_t_tests_response <- function(data, filter) {
  result = data[filter, ] %>%
    dplyr::group_by(cell_type) %>%
    dplyr::do(tidy(t.test(value ~ Response, data = .)))
  return(result)
}

# Create significance labels (*) for the plots
add_significance_labels <- function(t_test_result) {
  t_test_result %>%
    mutate(significance_label = case_when(
      p.value > 0.05 ~ "n.s.",
      p.value <= 0.05 & p.value > 0.01 ~ "*",
      p.value <= 0.01 & p.value > 0.001 ~ "**",
      p.value <= 0.001 & p.value > 0.0001 ~ "***",
      p.value <= 0.0001 ~ "****"
    )) %>%
    dplyr::arrange(cell_type)
}

# Generate t-test results
t_test_results = list()

t_test_results[["Resp.vs.Nonresp"]] = perform_t_tests_response(MCP_long[["Full"]], 
                                                               1:nrow(MCP_long[["Full"]]))
t_test_results[["Resp.vs.Nonresp"]] = add_significance_labels(as.data.frame(t_test_results[["Resp.vs.Nonresp"]]))

t_test_results[["RespT1.vs.Nonresp_T1"]] = perform_t_tests_response(MCP_long[["Full"]], 
                                                                    which(MCP_long[["Full"]]$Timepoint == "T1"))
t_test_results[["RespT1.vs.Nonresp_T1"]] = add_significance_labels(as.data.frame(t_test_results[["RespT1.vs.Nonresp_T1"]]))

t_test_results[["RespT2.vs.Nonresp_T2"]] = perform_t_tests_response(MCP_long[["Full"]], 
                                                                    which(MCP_long[["Full"]]$Timepoint == "T2"))
t_test_results[["RespT2.vs.Nonresp_T2"]] = add_significance_labels(as.data.frame(t_test_results[["RespT2.vs.Nonresp_T2"]]))

t_test_results[["RespT1.vs.RespT2"]] = perform_t_tests_timepoints(MCP_long[["Responders"]], 
                                                                  1:nrow(MCP_long[["Responders"]]))
t_test_results[["RespT1.vs.RespT2"]] = add_significance_labels(as.data.frame(t_test_results[["RespT1.vs.RespT2"]]))

t_test_results[["NonrespT1.vs.NonrespT2"]] = perform_t_tests_timepoints(MCP_long[["Non_responders"]], 
                                                                        1:nrow(MCP_long[["Non_responders"]]))
t_test_results[["NonrespT1.vs.NonrespT2"]] = add_significance_labels(as.data.frame(t_test_results[["NonrespT1.vs.NonrespT2"]]))


combined_results = rbind(t_test_results[[1]],
                         t_test_results[[2]],
                         t_test_results[[3]],
                         t_test_results[[4]],
                         t_test_results[[5]])

for (i in 1:length(titles)) {
  for (j in 1:nrow(MCP[["Full"]])) {
    box_dataset = box_datasets[[i]] %>% dplyr::filter(cell_type == cell_types[j])
    sig_label = t_test_results[[i]][["significance_label"]][t_test_results[[i]]$cell_type == cell_types[j]]
    box_dataset$sig_label = ifelse(nrow(box_dataset) > 0, sig_label, NA)
    boxplots[[i]][[j]] = ggplot(data = box_dataset, 
                                aes(x = !!sym(x_variables[i]), y = value, 
                                    fill = !!sym(x_variables[i]))) +
      geom_boxplot(alpha = 0.8, outlier.shape = NA, width = 0.3) +
      geom_text(aes(x = 1.5, y = max(value, na.rm = TRUE), 
                    label = paste(as.character(sig_label))), 
                size = 3, vjust = 1) +
      scale_fill_manual(values = c("Responder" = "dodgerblue4", "Non_responder" = "deeppink4",
                                   "T1" = "goldenrod2", "T2" = "purple4")) +
      labs(title = NULL, #paste0(cell_types[j], ": ", titles[i])
           x = NULL, y = cell_types[j]) +
      theme_classic() +
      theme(panel.background = element_rect(fill = "white", 
                                            colour = "white"),
            panel.grid = element_blank(),
            axis.line = element_line(),
            axis.title.y = element_text(size = 7, face = "bold"),
            axis.text.x = element_text(size = 5),
            legend.position = "none")
  }
  names(boxplots[[i]]) = cell_types
}
rm(box_dataset)

# Create and export ggarrange() objects
library(ggpubr)

# Create output directory
dir.create("Signatures/MCPcounter")
MCPcounter_full_plots = list()
for (i in 1:length(titles)) {
  MCPcounter_full_plots[[i]] = ggarrange(plotlist = boxplots[[i]], ncol = 2, 
                                         nrow = 4, labels = NULL)
  print(MCPcounter_full_plots[[i]])
  ggsave(filename = paste0("MCPcounter_", titles[i], ".tiff"),
         path = "Signatures/MCPcounter", 
         width = 3612, height = 6500, device = 'tiff', units = "px",
         dpi = 700, compression = "lzw")
  dev.off()
}

# Create a comprehensive plot of singificant and meaningful results
# i.e. results that are different between the groups of response/timepoint
# It is going to be a 3x1 grid of three 1x3 ggarrange() object (a 3x3 grid ultimately)

# A legend plot
p9 = ggplot() +
  geom_rect(aes(xmin = -1, xmax = 1, ymin = -1, ymax = 1), fill = NA) +
  annotate("richtext", x = -0.5, y = 0, 
           label = "<b>a, b, c:</b> All samples<br><br>
                    <b>d, e, f:</b> T2 samples<br><br>
                    <b>g:</b> Responders (T1 vs. T2)<br><br>
                    <b>h:</b> Non-responders (T1 vs. T2)",
           hjust = 0, vjust = 0.5, fill = NA, label.color = NA, label.size = 0, size = 2) +
  theme_void() +
  coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1))

ggarrange1 = ggarrange(boxplots$RespvsNonresp$`Endothelial cells`,
                       boxplots$RespvsNonresp$Fibroblasts,
                       boxplots$RespvsNonresp$Neutrophils,
                       ncol = 3, nrow = 1, labels = c("a", "b", "c"),
                       font.label = list(size = 8.5))
ggarrange2 = ggarrange(boxplots$RespT2vsNonrespT2$`Endothelial cells`,
                       boxplots$RespT2vsNonrespT2$Fibroblasts,
                       boxplots$RespT2vsNonrespT2$Neutrophils,
                       ncol = 3, nrow = 1, labels = c("d", "e", "f"),
                       font.label = list(size = 8.5))
ggarrange3 = ggarrange(boxplots$RespT1vsRespT2$`Endothelial cells`,
                       boxplots$NonrespT1vsNonrespT2$`Endothelial cells`,
                       p9,
                       ncol = 3, nrow = 1, labels = c("g", "h", "Legend"),
                       font.label = list(size = 8.5))
tiff("Signatures/MCPcounter/MCPcounter_results.tiff", width = 3500, height = 4800,
     res = 700, compression = "lzw")
ggarrange(ggarrange1, ggarrange2, ggarrange3,
          ncol = 1, nrow = 3, labels = NULL)
dev.off()

# MSigDB Oncogenic Collection Enrichment #####
library(msigdb)
library(ExperimentHub)
library(GSEABase)

# Query the database
eh = ExperimentHub()
query(eh , 'msigdb')

# Custom download
msigdb.hs = getMsigdb(org = 'hs', id = 'EZID', version = '7.2')
msigdb.hs

# Oncogenic signatures
onco = subsetCollection(msigdb.hs, 'c6')
onco_ids = geneIds(onco)

# Prepare for limma::fry or limma::camera
camera_indices = ids2indices(onco_ids, rownames(z_exprs))
z_design = readRDS("Signatures/design_matrix.rds") # design matrix
z_cm = readRDS("Signatures/contrast_matrix.rds") # contrast matrix

onco_enr = limma::camera(z_exprs, index = camera_indices, 
                  design = z_design, contrast = z_cm)
enriched_oncosets = onco_enr[onco_enr$FDR < 0.05,] # 5 sets

# Enriched oncogenic sets metadata: descriptions
onco_descriptions_long = list()
onco_descriptions_short = list()
for(i in rownames(onco_enr)){
  onco_descriptions_long[[i]] = msigdb.hs[[i]]@longDescription
  onco_descriptions_short[[i]] = msigdb.hs[[i]]@shortDescription
}
rm(i); gc()

onco_enr$Terms = rownames(onco_enr)
onco_enr$Long.Desc = onco_descriptions_long
onco_enr$Short.Desc = onco_descriptions_short

onco_descriptions_short[[1]]; onco_descriptions_long[[1]]
onco_descriptions_short[[2]]; onco_descriptions_long[[2]]
onco_descriptions_short[[3]]; onco_descriptions_long[[3]]
onco_descriptions_short[[4]]; onco_descriptions_long[[4]]
onco_descriptions_short[[5]]; onco_descriptions_long[[5]]

onco_enr = onco_enr %>%
  dplyr::select(Terms, NGenes, Direction, PValue, FDR, Short.Desc, Long.Desc)

# Write out enrichment results
write.xlsx(onco_enr, "Signatures/MSigDB_Oncogenic_Collection_Enrichment.xlsx",
           overwrite = TRUE)

# Session Info #####
sessionInfo()
