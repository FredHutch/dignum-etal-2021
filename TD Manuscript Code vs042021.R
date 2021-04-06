### Project: Analysis of scRNAseq data from E9 and E9.5 P-Sp/AGM V+E+61+ cells for T Dignum manuscript, Dev Cell submitted 4-2021
### Load packages
library(monocle3) 
library(leidenbase)
library(ggplot2)
library(dplyr)
library(tibble)
library(stringr)
library(pheatmap)
library(Rcpp)
library(reticulate)
library(VGAM)
library(viridis)
library(magrittr)
library(mygene)
library(ggpubr)
library(scales)
library(garnett)
library(org.Mm.eg.db)

###define file path
# DIR <- file.path("~/...")

simple_theme <-  theme(text = element_blank(),
                       panel.grid = element_blank(),
                       axis.title = element_blank(),
                       axis.text = element_blank(),
                       axis.ticks = element_blank(),
                       axis.line.x = element_blank(),
                       axis.line.y = element_blank(), legend.position = "none")  ### theme to remove legends and axis/font

simple_themeL <-  theme(text = element_text(size=25),
                        panel.grid = element_blank(),
                        axis.title = element_blank(),
                        axis.text = element_blank(),
                        axis.ticks = element_blank(),
                        axis.line.x = element_blank(),
                        axis.line.y = element_blank())   #### theme show the legend with large font

###load cds 
E9cds <- readRDS(file.path(DIR, "Sample1.RDS"))
E95cds <- readRDS(file.path(DIR, "Sample2.RDS"))

###combine E9 and E9.5 V+61+E+ only
cds <- combine_cds(list(E9cds,E95cds))  ### limited to v+61+E+ E9-9.5

### Rremove cells with low UMI and genes per cell
colData(cds)$n.umis <- Matrix::colSums(counts(cds))
cds <- cds[,Matrix::colSums(counts(cds)) > 5000]
cds <- detect_genes(cds, min_expr=0.1)    
cds <- cds[,colData(cds)$num_genes_expressed > 1000]

### plot UMI and genes per cell for each sample Fig S3B
simple_theme2 <-  theme(panel.background = element_rect(fill="white", color="white"),
                        legend.position = "none", axis.line = element_line("black"), aspect.ratio=0.5,
                        panel.grid = element_blank(), axis.title = element_blank(), axis.text = element_text(size=15))
df1 <- data.frame("Sample"=colData(cds)$Sample, "n.umis"=colData(cds)$n.umis)
df1$Sample <- factor(df1$Sample, levels = c("E9", "E95"))
ggplot(df1, aes(x=Sample, y=n.umis, fill=Sample)) +geom_boxplot() +  
  scale_y_continuous(trans='log10', limits=c(100,100000), labels=number) + simple_theme2

df2 <- data.frame("Sample"=colData(cds)$Sample, "num_genes_expressed"=colData(cds)$num_genes_expressed)
df2$Sample <- factor(df2$Sample, levels = c("E9", "E95"))
ggplot(df2, aes(x=Sample, y=num_genes_expressed, fill=Sample)) +geom_boxplot()  +
  scale_y_continuous(trans='log10', limits=c(100,10000), labels=number) + simple_theme2 

###pre-process cds
cds <- preprocess_cds(cds) 

###align cds to remove batch effects 
cds = align_cds(cds, alignment_group = "Sample")

### reduce dimensions
cds <- reduce_dimension(cds, umap.fast_sgd = FALSE, cores=1, n_sgd_threads=1) 

###plot in UMAP by sample Fig S3C
colData(cds)$Sample <- factor(colData(cds)$Sample, levels = c("E9", "E95"))
plot_cells(cds, color_cells_by = "Sample", label_cell_groups = F, cell_size = 1, show_trajectory_graph = FALSE) + simple_themeL

###cluster cells and plot in UMAP Fig S3D, 3A
set.seed(17)
cds = cluster_cells(cds, resolution =3E-2, random_seed = 17)
plot_cells(cds, label_cell_groups = T, cell_size = 1, group_label_size = 5, show_trajectory_graph = FALSE) +simple_themeL

###load cds saved following processing steps above
cds <- readRDS(file.path(DIR, "Sample1_2_P.RDS"))

###plot gene expression in UMAP: Fig 3B, S3E-F
plot_cells(cds, genes=c("Meox1"), cell_size = 1.5, label_cell_groups = F, show_trajectory_graph = FALSE) +simple_theme 
plot_cells(cds, genes=c("Pax3"), cell_size = 1.5, label_cell_groups = F, show_trajectory_graph = FALSE) +simple_theme 
plot_cells(cds, genes=c("Sox2"), cell_size = 1.5, label_cell_groups = F, show_trajectory_graph = FALSE) +simple_theme 
plot_cells(cds, genes=c("Cdh5"), cell_size = 1.5, label_cell_groups = F, show_trajectory_graph = FALSE) +simple_theme 
plot_cells(cds, genes=c("Procr"), cell_size = 1.5, label_cell_groups = F, show_trajectory_graph = FALSE) +simple_theme 
plot_cells(cds, genes=c("Etv2"), cell_size = 1.5, label_cell_groups = F, show_trajectory_graph = FALSE) +simple_theme 
plot_cells(cds, genes=c("Nr2f2"), cell_size = 1.5, label_cell_groups = F, show_trajectory_graph = FALSE) +simple_theme 
plot_cells(cds, genes=c("Dll4"), cell_size = 1.5, label_cell_groups = F, show_trajectory_graph = FALSE) +simple_theme 
plot_cells(cds, genes=c("Efnb2"), cell_size = 1.5, label_cell_groups = F, show_trajectory_graph = FALSE) +simple_theme 
plot_cells(cds, genes=c("Cxcr4"), cell_size = 1.5, label_cell_groups = F, show_trajectory_graph = FALSE) +simple_theme 
plot_cells(cds, genes=c("Cd44"), cell_size = 1.5, label_cell_groups = F, show_trajectory_graph = FALSE) +simple_theme 
plot_cells(cds, genes=c("Gja5"), cell_size = 1.5, label_cell_groups = F, show_trajectory_graph = FALSE) +simple_theme 
plot_cells(cds, genes=c("Runx1"), cell_size = 1.5, label_cell_groups = F, show_trajectory_graph = FALSE) +simple_theme 
plot_cells(cds, genes=c("Gfi1"), cell_size = 1.5, label_cell_groups = F, show_trajectory_graph = FALSE) +simple_theme 
plot_cells(cds, genes=c("Ptprc"), cell_size = 1.5, label_cell_groups = F, show_trajectory_graph = FALSE) +simple_theme 
plot_cells(cds, genes=c("Flt3"), cell_size = 1.5, label_cell_groups = F, show_trajectory_graph = FALSE) +simple_theme 
plot_cells(cds, genes=c("Il7r"), cell_size = 1.5, label_cell_groups = F, show_trajectory_graph = FALSE) +simple_theme 
plot_cells(cds, genes=c("Csf3r"), cell_size = 1.5, label_cell_groups = F, show_trajectory_graph = FALSE) +simple_theme 
plot_cells(cds, genes=c("Kit"), cell_size = 1.5, label_cell_groups = F, show_trajectory_graph = FALSE) +simple_theme 
plot_cells(cds, genes=c("Myb"), cell_size = 1.5, label_cell_groups = F, show_trajectory_graph = FALSE) +simple_theme 
plot_cells(cds, genes=c("Cd33"), cell_size = 1.5, label_cell_groups = F, show_trajectory_graph = FALSE) +simple_theme 


#### Psuedotime analysis
cds <- learn_graph(cds, learn_graph_control=list(ncenter=300, prune_graph=T))
cds <- order_cells(cds)   ### set root node in Etv2/Nr2f2+ EC cluster = least differentiated state

###select cells, remove outlying clusters (somitic)
cds<-choose_cells(cds) 

###assign cxcr4 expressing clusters based on cxcr4 expression in >5% of cells in a cluster 
###(accounts for low Cxcr4 transcript/dropout, Cxcr4+ by this appraoch ~%CXCR4+ by FACS)
### Fig 3C
pData(cds)$cxcr4 = Matrix::colSums(exprs(cds[fData(cds)$gene_short_name %in% c("Cxcr4")]))
df <- data.frame(monocle3::clusters(cds), (colData(cds)$cxcr4>0))
##renames the labels of the data frame 
names(df) <- c("cluster", "cxcr4")
df2 <- df %>% 
  group_by(cluster,cxcr4) %>% 
  summarise(count=n()) %>% 
  mutate(perc=count/sum(count))
df3<-df2 %>% filter(cxcr4 == "TRUE") 
df4<-df3 %>% filter(perc>0.05)
colData(cds)$cxcr4 = (monocle3::clusters(cds)  %in% df4$cluster)
cxcr4_color <- c("TRUE"= "darkmagenta", "FALSE"= "deepskyblue3")
plot_cells(cds, color_cells_by = "cxcr4", label_cell_groups = F, cell_size = 1, show_trajectory_graph = FALSE) +
  scale_color_manual(values = cxcr4_color) + simple_themeL

####cell type analysis using Garnette
set.seed(260)
marker_file_path <- "~/AGM_markers.txt"
hspc_classifier <- train_cell_classifier(cds = cds,
                                         marker_file = marker_file_path,
                                         db=org.Mm.eg.db,
                                         cds_gene_id_type = "ENSEMBL",
                                         num_unknown = 100,
                                         marker_file_gene_id_type = "SYMBOL")   ### trains the cell classifier

cds <- classify_cells(cds, hspc_classifier,
                      db = org.Mm.eg.db,
                      cds_gene_id_type = "ENSEMBL")              ### defines cell types 

### plot by cell type Fig 3D
tgray=alpha("gray90", 0.1)
cell_type_color2 <- c("HPC"= "green", "UEC" = "#F0E442", "Unknown"=tgray, "AEC1"= "#E69F00","AEC2"="#D55E00","HE"="#0072B2")
colData(cds)$cell_type <- factor(colData(cds)$cell_type, levels = c("UEC","AEC1","AEC2","HE", "HPC", "Unknown"))
plot_cells(cds,
           color_cells_by = "cell_type",
           show_trajectory_graph = FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1, cell_size = 1) + 
  scale_color_manual(values = cell_type_color2)  +simple_themeL

###determine relative cell type contributions by Cxcr4 cluster and plot by % in bar plots 
###Fig 3D, right upper panel
dft<-data.frame("cxcr4"=colData(cds)$cxcr4, "celltype"=colData(cds)$cell_type)
dft<-dft %>% filter(celltype !="Unknown")   
#dft$celltype<-factor(dft$celltype, levels = c("UEC","AEC1","AEC2","HE","HPC"))
df2 <- dft %>% 
  group_by(cxcr4,celltype) %>% 
  summarise(count=n()) %>% 
  mutate(perc=count/sum(count))
df2$celltype<-factor(df2$celltype, levels = c("HPC", "HE", "AEC2","AEC1","UEC"))
cell_type_color3 <- c("HPC"= "green", "UEC" = "#F0E442", "AEC1"= "#E69F00","AEC2"="#D55E00","HE"="#0072B2")
ggplot(df2, aes(x = factor(cxcr4), y = perc*100, fill = factor(celltype))) +
  geom_bar(stat="identity", width = 0.7) +
  labs(x = "cxcr4", y = "percent", fill = "celltype") +
  theme_minimal(base_size = 14) + scale_fill_manual(values = cell_type_color3) + theme(plot.margin = unit(c(0,3,0,3), "cm")) + ylim(0,101) 

###plot number of HE cell type in Cxcr4 positive vs negative clusters 
###Fig 3D right lower panel
df3<-df2 %>% filter(celltype =="HE")   ### keep HE only
ggplot(df3, aes(x = factor(cxcr4), y = count, fill = factor(cxcr4))) +
  geom_bar(stat="identity", width = 0.7) +
  labs(x = "cxcr4", y = "count", fill = "celltype") +
  theme_minimal(base_size = 14) + scale_fill_manual(values = cxcr4_color) + theme(plot.margin = unit(c(0,2,0,2), "cm")) + ylim(0,NA) 


#### Estimate gene set scores for signature gene lists
estimate_score <- function(cds, markers){
  cds_score = cds[fData(cds)$gene_short_name %in% markers,] 
  aggregate_score = exprs(cds_score)
  aggregate_score = Matrix::t(Matrix::t(aggregate_score) / pData(cds_score)$Size_Factor)
  aggregate_score = Matrix::colSums(aggregate_score)
  pData(cds)$score = log(aggregate_score +1) 
  return(cds)
}

estimate_score2 <- function(cds, markers){
  cds_score = cds[fData(cds)$gene_short_name %in% markers,] 
  aggregate_score = exprs(cds_score)
  aggregate_score = Matrix::t(Matrix::t(aggregate_score) / pData(cds_score)$Size_Factor)
  aggregate_score = Matrix::colSums(aggregate_score)
  pData(cds)$score2 = log(aggregate_score +1) 
  return(cds)
}


### signature gene sets:
### Arterial EC genes from Xu Nat Comm 2018
AEC_genes <- as.character(c(read.csv(file.path(DIR, "XU_AEC_GENES.csv"), header=F))$V1)
### AEC genes downregulated in Dll4 mutants (Luo, Benedito Nature 2020)
AEC2_genes <- as.character(c(read.csv(file.path(DIR, "LUO_DLL4_AEC_GENES.csv"), header=F))$V1)
#### TFs required for AEC fate in vitro (Aranguren 2013)
AEC3_genes <- as.character(c(read.csv(file.path(DIR, "ARANGUREN_AEC_TFS.csv"), header=F))$V1)
### AGM HSC-specific genes from Vink/Dzeirzak 2020:
AGM_HSC_sig_genes <- c("Hes1", "Hey1", "Notch4", "Sox17", "Il6st", "Nos3", "Flt1", "Kdr", "Ly6a", "Cdkn1c", "Pdzk1ip1", "Procr", "Ramp2", "Trim47", "Vwf", "Meis2", "Cdh5", "Gfi1", "Pbx1", "Ptprb", "Bmp2k", "Nrp1", "Ptpru", "Dll4")  
### HSC sig genes from Wilson 2015:
HSC_sig_genes <- c("Procr", "Pdzk1ip1", "Ltb", "Mllt3", "Ifitm1", "Gimap1", "Gimap6", "Limd2", "Trim47", "Neil2", "Vwf", "Pde1b", "Neo1", "Sqrdl", "Sult1a1", "Cd82", "Ramp2", "Ubl3", "Ly6a", "Cdkn1c", "Fgfr3", "Cldn10", "Ptpn14", "Mettl7a1", "Smtnl1", "Ctsf", "Gstm1", "Sox18", "Fads3")  
### Serially engrafting HSC signature genes from Rodriguez-Fraticelli et al 2020
HSC2_sig_genes <- as.character(c(read.csv(file.path(DIR, "Serial_HSC_genes_Camargo.csv"), header=F))$V1)
#HSC genes from Chanbers
HSC3_sig_genes  <- as.character(c(read.csv(file.path(DIR, "HSC_genes_Chambers.csv"), header=F))$V1)
##### dHSC genes from Cabezas
HSC4_sig_genes <- as.character(c(read.csv(file.path(DIR, "Cabezas_dHSC.csv"), header=F))$V1)
### HE sig genes from Hou 2020
HE_sig_genes <- c("Neurl3", "Phgdh", "Sfrp2", "Nupr1", "Mycn", "Hlf", "Gfi1", "Gck", "Ift57", "Eya2", "Lmo1")   

### Genes involved in HSC self renewal from published literature:
HSC_SR_genes <- c("Cdkn1c", "Mecom", "Pbx1", "Txnip", "Kmt2a", "Tcf15", "Erg", "Sirt6", "Kdm6b", "H19", "Id1", "Igf2", "Ndn", "Eng", "Cd81", "Prdm16")
### MYC pathways genes upregulated in Dll4 mutants (Luo, Benedito Nature 2020)
Myc_genes <- as.character(c(read.csv(file.path(DIR, "LUO_DLL4_REG_MYC_GENES.cs"), header=F))$V1)
####GSEA Hallmark Myc genes
Myc2_genes <- as.character(c(read.csv(file.path(DIR, "GSEA_Hallmark_Myc_V1.csv"), header=F))$V1)
#####GSEA Hallmark Myc genesV2
Myc3_genes <- as.character(c(read.csv(file.path(DIR, "GSEA_Hallmark_Myc_V2.csv"), header=F))$V1)
####Boroviak et al. 2015, diapaused Epiblast
Diapause_up_genes <- as.character(c(read.csv(file.path(DIR, "DIAPAUSE_UP_BOROVIAK.csv"), header=F))$V1)
####FRIDMAN_SENESCENCE_UP
Senesc_genes <- as.character(c(read.csv(file.path(DIR, "FRIDMAN_SENESCENCE_UP.csv"), header=F))$V1)

cds <- choose_cells(cds)      ## remove HPC cluster 23 not associated in pseudotime for GSS scores and module analysis below
cds <- estimate_score(cds, markers = AEC_genes)  ### repeat for each GSS

### plot scaled gene set scores in UMAP, Figs 3E, 4B, D, S3G, S4F
plot_cells(cds, color_cells_by = 'score', cell_size = 1, show_trajectory_graph = F) + 
  scale_color_gradient2(low="gray80",mid="gray80",high="red3",midpoint=(((max(pData(cds)$score)-min(pData(cds)$score))/2)+min(pData(cds)$score)), space = "Lab") +simple_theme

###makes a dataframe of gene set scores by cell type to generate a violin plot
###Figs 3E, 4B, D, S3G, S4F
dfx <- data.frame(colData(cds)$cell_type, colData(cds)$score)
names(dfx) <- c("celltype", "score")
dfx<-dfx %>% filter(celltype !="Unknown")     ### filter out Unknwon
dfx<-dfx %>% filter(celltype !="HPC")     ### filter out HPC
cell_type_color3 <- c("UEC" = "#F0E442", "AEC1"= "#E69F00","AEC2"="#D55E00","HE"="#0072B2")
p <- ggplot(dfx, aes(x= celltype, y=score, fill = celltype)) + geom_violin(trim=FALSE) + geom_boxplot(fill= "white", width=0.1, outlier.colour=NA) 
p <- p + theme_bw() + scale_fill_manual(values = cell_type_color3)
p <- p + xlab(NULL) + ylab("Gene-set Score") ##renames axis label
p <- p + theme(plot.margin = unit(c(0,4.5,0,4.5), "cm"))
p <- p + theme(panel.grid.major= element_blank(), panel.grid.minor = element_blank())
p <- p+theme(legend.position = "none") + theme(text = element_text(size=35))
p

###makes a dataframe of gene set score to of HE from cxcr4-expressing vs non-expressing clusters and generates a violin plot
###Figs 3E, 4B, D, S3G, S4F
df <- data.frame(colData(cds)$cxcr4, colData(cds)$cell_type, colData(cds)$score) 
names(df) <- c("cxcr4", "celltype", "score")
df<-df %>% filter(celltype =="HE")    
df$cxcr4<-factor(df$cxcr4, levels = c("TRUE", "FALSE"))
p <- ggplot(df, aes(x= cxcr4, y=score, fill = cxcr4)) + geom_violin(trim=FALSE) + geom_boxplot(fill= "white", width=0.1, outlier.colour=NA) 
p <- p + stat_compare_means() + stat_mean()   ###statistical comparison to calculate p value between 2 samples (Wilcoxon)
p <- p + theme_bw() + scale_fill_manual(values = cxcr4_color)
p <- p + xlab(NULL) + ylab("Gene-set Score") ##renames axis label
p <- p + theme(plot.margin = unit(c(0,5.5,0,5.5), "cm"))
p <- p + theme(panel.grid.major= element_blank(), panel.grid.minor = element_blank())
p <- p+theme(legend.position = "none") + theme(text = element_text(size=35))
p

#### generates a scatter plot of two GSS for cxcr4+ vs - HE with linear correlation (geom_smooth method lm)
#### Fig S4E, G
cds <- estimate_score(cds, markers = HSC4_sig_genes)   ### repeat for each GSS comparison
cds <- estimate_score2(cds, markers = Myc_genes)
df <- data.frame(colData(cds)$cxcr4, colData(cds)$cell_type, colData(cds)$score, colData(cds)$score2)
names(df) <- c("cxcr4", "celltype", "score1", "score2")
df<-df %>% filter(celltype =="HE")
df$cxcr4<-factor(df$cxcr4, levels = c("TRUE", "FALSE"))
p <- ggplot(df, aes(x=score1, y=score2)) +geom_point(aes(color=cxcr4, size=3)) +geom_smooth(method="lm", aes(x=score1, y=score2))
p <- p + theme_bw() + scale_color_manual(values = cxcr4_color)
p <- p + theme(plot.margin = unit(c(0.5,0,0.5,0), "cm"))
p <- p + theme(panel.grid.major= element_blank(), panel.grid.minor = element_blank())
p <- p+theme(legend.position = "none") + theme(text = element_text(size=20))
p

### Select genes that are significant DEG between Cxcr4+ vs Cxcr4- HE based on FDR<0.05
### Table S4
cds_sub <- cds[,(colData(cds)$cell_type == "HE"),]
cds_sub <- detect_genes(cds_sub, min_expr=10)
exp_genes <- row.names(subset(fData(cds_sub), fData(cds_sub)$num_cells_expressed>0))
cds_sub <- cds_sub[rowData(cds_sub)$id %in% exp_genes,]
gene_fits <- fit_models(cds_sub, model_formula_str = "~cxcr4")
fit_coefs <- coefficient_table(gene_fits)
identity_terms <- fit_coefs %>% filter(term != "(Intercept)")
identity_DEG <- identity_terms %>% filter (q_value < 0.05) %>%
  dplyr::select(gene_short_name, term, q_value, estimate)
write.csv(identity_DEG, file.path(DIR, "HE_DEG_CXCR4posVneg.csv"))

### generate violin plots for HSC self-renewal genes by cell type with separate Cxcr4+ vs Cxcr4- HE, Fig 4A
### create new cell type separating HE to CXCR4+ vs CXCR4-
colData(cds)$cell_type2 = as.character(paste(colData(cds)$cell_type, colData(cds)$cxcr4))
colData(cds)$cell_type2[colData(cds)$cell_type2 == "HE TRUE"] <- "HEcxcr4p"
colData(cds)$cell_type2[colData(cds)$cell_type2 == "HE FALSE"] <- "HEcxcr4n"
colData(cds)$cell_type2[colData(cds)$cell_type2 == "AEC1 TRUE"] <- "AEC1"
colData(cds)$cell_type2[colData(cds)$cell_type2 == "AEC1 FALSE"] <- "AEC1"
colData(cds)$cell_type2[colData(cds)$cell_type2 == "AEC2 TRUE"] <- "AEC2"
colData(cds)$cell_type2[colData(cds)$cell_type2 == "AEC2 FALSE"] <- "AEC2"
colData(cds)$cell_type2[colData(cds)$cell_type2 == "UEC TRUE"] <- "UEC"
colData(cds)$cell_type2[colData(cds)$cell_type2 == "UEC FALSE"] <- "UEC"
colData(cds)$cell_type2[colData(cds)$cell_type2 == "HPC TRUE"] <- "HPC"
colData(cds)$cell_type2[colData(cds)$cell_type2 == "HPC FALSE"] <- "HPC"
colData(cds)$cell_type2[colData(cds)$cell_type2 == "Unknown FALSE"] <- "Unknown"
colData(cds)$cell_type2[colData(cds)$cell_type2 == "Unknown TRUE"] <- "Unknown"
SR_genes <- c("H19", "Cdkn1c", "Ndn", "Mecom", "Cd81", "Eng", "Igf2")
cds_subset <- cds[,(colData(cds)$cell_type2 != "Unknown"),]  
cds_subset <- cds_subset[,(colData(cds_subset)$cell_type2 != "HPC"),] 
colData(cds_subset)$cell_type2<-factor(colData(cds_subset)$cell_type2, levels = c("UEC", "AEC1", "AEC2", "HEcxcr4p","HEcxcr4n"))
cell_type_color4 <- c("UEC" = "#F0E442", "AEC1"= "#E69F00","AEC2"="#D55E00","HEcxcr4p"= "darkmagenta", "HEcxcr4n"= "deepskyblue3")
cds_subset<-cds_subset[rowData(cds_subset)$gene_short_name %in% SR_genes,] 
plot_genes_violin(cds_subset, group_cells_by = "cell_type2",ncol=7) + 
  theme(axis.text.x=element_text(angle=0, hjust=1)) +scale_fill_manual(values = cell_type_color4) +
  theme(plot.margin = unit(c(2,0,2,0), "cm"))

### Estimate proliferation index, Fig 4C
estimate_cell_cycle <- function(cds, g1s_markers, g2m_markers){
  cds_g1s = cds[fData(cds)$gene_short_name %in% g1s_markers,]
  aggregate_g1s_expression = exprs(cds_g1s)
  aggregate_g1s_expression = Matrix::t(Matrix::t(aggregate_g1s_expression) / pData(cds_g1s)$Size_Factor)
  aggregate_g1s_expression = Matrix::colSums(aggregate_g1s_expression)
  
  cds_g2m = cds[fData(cds)$gene_short_name %in% g2m_markers,]
  aggregate_g2m_expression = exprs(cds_g2m)
  aggregate_g2m_expression = Matrix::t(Matrix::t(aggregate_g2m_expression) / pData(cds_g2m)$Size_Factor)
  aggregate_g2m_expression = Matrix::colSums(aggregate_g2m_expression)
  
  pData(cds)$g1s_score = log(aggregate_g1s_expression+1)
  pData(cds)$g2m_score = log(aggregate_g2m_expression+1)
  pData(cds)$proliferation_index = log(aggregate_g1s_expression + aggregate_g2m_expression + 1)
  return(cds)
}
s.genes <- c("Mcm5", "Pcna", "Tyms", "Fen1", "Mcm2", "Mcm4", "Rrm1", "Ung", "Gins2", "Mcm6", "Cdca7", "Dtl", "Prim1", "Uhrf1",
             "Mlf1ip", "Hells", "Rfc2", "Rpa2", "Nasp", "Rad51ap1", "Gmnn", "Wdr76", "Slbp", "Ccne2", "Ubr7", "Pold3", "Msh2", 
             "Atad2", "Rad51", "Rrm2", "Cdc45", "Cdc6", "Exo1", "Tipin", "Dscc1", "Blm", "Casp8ap2", "Usp1", "Clspn", "Pola1", 
             "Chaf1b", "Brip1", "E2f8")
g2m.genes <- c("Hmgb2", "Cdk1", "Nusap1", "Ube2c", "Birc5", "Tpx2", "Top2a", "Ndc80", "Cks2", "Nuf2", "Cks1b", "Mki67", "Tmpo",
               "Cenpf", "Tacc3", "Fam64a", "Smc4", "Ccnb2", "Ckap2l", "Ckap2", "Aurkb", "Bub1", "Kif11", "Anp32e", "Tubb4b", "Gtse1",
               "Kif20b", "Hjurp", "Cdca3", "Hn1", "Cdc20", "Ttk", "Cdc25c", "Kif2c", "Rangap1", "Ncapd2", "Dlgap5", "Cdca2", "Cdca8",
               "Ect2", "Kif23", "Hmmr", "Aurka", "Psrc1", "Anln", "Lbr", "Ckap5", "Cenpe", "Ctcf", "Nek2", "G2e3", "Gas2l3", "Cbx5", "Cenpa")

cds <- estimate_cell_cycle(cds, g1s_markers = s.genes, g2m_markers = g2m.genes)
mycol <- c("navy", "blue", "cyan", "lightcyan", "yellow", "red", "red4")
plot_cells(cds, color_cells_by = 'proliferation_index', cell_size = 1, show_trajectory_graph = F)  +
  scale_color_gradientn(colours = mycol) + simple_themeL

###generate a violin plot of proliferation index by cell type
dfx <- data.frame(colData(cds)$cell_type, colData(cds)$proliferation_index)
names(dfx) <- c("celltype", "score")
dfx<-dfx %>% filter(celltype !="Unknown")    
dfx<-dfx %>% filter(celltype !="HPC")     
cell_type_color3 <- c("UEC" = "#F0E442", "AEC1"= "#E69F00","AEC2"="#D55E00","HE"="#0072B2")
p <- ggplot(dfx, aes(x= celltype, y=score, fill = celltype)) + geom_violin(trim=FALSE) + geom_boxplot(fill= "white", width=0.1, outlier.colour=NA) 
p <- p + theme_bw() + scale_fill_manual(values = cell_type_color3)
p <- p + xlab(NULL) + ylab(NULL)
p <- p + theme(plot.margin = unit(c(0,4,0,4), "cm"))
p <- p + theme(panel.grid.major= element_blank(), panel.grid.minor = element_blank())
p <- p+theme(legend.position = "none") + theme(text = element_text(size=35))
p

###generate a violin plot of proliferation index for Cxcr4+ vs Cxcr4- HE 
df <- data.frame(colData(cds)$cxcr4, colData(cds)$cell_type, colData(cds)$proliferation_index) 
names(df) <- c("cxcr4", "celltype", "score")
df<-df %>% filter(celltype =="HE")
df$cxcr4<-factor(df$cxcr4, levels = c("TRUE", "FALSE"))
p <- ggplot(df, aes(x= cxcr4, y=score, fill = cxcr4)) + geom_violin(trim=FALSE) + geom_boxplot(fill= "white", width=0.1, outlier.colour=NA) 
p <- p + stat_compare_means() + stat_mean()   ### statistical comparison (Wilcoxon) to calculate p value
p <- p + theme_bw() + scale_fill_manual(values = cxcr4_color)
p <- p + xlab(NULL) + ylab(NULL)
p <- p + theme(plot.margin = unit(c(0,5.5,0,5.5), "cm"))
p <- p + theme(panel.grid.major= element_blank(), panel.grid.minor = element_blank())
p <- p+theme(legend.position = "none") + theme(text = element_text(size=35))
p


#### gene modules using principal_graph to identify modules as function of pseudotime
### plot cells in psuedotime
### Fig 3F
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE, trajectory_graph_segment_size = 1,
           label_branch_points=FALSE,
           graph_label_size=1,  cell_size = 1) +simple_theme

colData(cds)$Cluster = monocle3::clusters(cds) ## assign clusters in colData
subset_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=2)
pr_deg_ids <- row.names(subset(subset_pr_test_res, q_value < 0.05))
set.seed(7) 
gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution = 4E-3, random_seed = 7) 
### plot gene modules by cell type in pheatmap
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=colData(cds)$cell_type)  
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
colnames(agg_mat) <- stringr::str_c(" ", colnames(agg_mat))
pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                   scale="column", clustering_method="ward.D2",
                   fontsize=8)
### plot gene modules by cluster_ in pheatmap
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=clusters(cds)[colnames(cds)])
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
colnames(agg_mat) <- stringr::str_c(" ", colnames(agg_mat))
pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                   scale="column", clustering_method="ward.D2",
                   fontsize=8)

### plot gene modules in UMAP - adjust min_expr to set threshold for lower limit of expression
### Fig 3G, module 9 associated with Cxcr4+ HE, module 15 associated with Cxcr4- HE
plot_cells(cds,
           genes=gene_module_df %>% filter(module %in% c(9)),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE,
           cell_size = 1,
           scale_to_range = TRUE,
           min_expr = 35) +scale_color_gradientn(colors=c("gray80", "royalblue3", "royalblue3")) +simple_themeL

plot_cells(cds,
           genes=gene_module_df %>% filter(module %in% c(15)),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE,
           cell_size = 1,
           scale_to_range = TRUE,
           min_expr = 20) +simple_theme +scale_color_gradientn(colors=c("gray80", "royalblue3", "royalblue3")) +simple_themeL

### Save genes in gene modules with Moran's i and q-values
### Table S2
mod_genes<-gene_module_df%>%dplyr::filter(module==9)
names(mod_genes)[1]<-"gene_id"
mod_genes_names<-merge(mod_genes,subset_pr_test_res,by.x="gene_id",by.y="id")
mod_genes_names<-mod_genes_names[,c(1,9,10,11,12)]
write.csv(mod_genes_names, file.path(DIR, "mod9.csv"))

mod_genes<-gene_module_df%>%dplyr::filter(module==15)
names(mod_genes)[1]<-"gene_id"
mod_genes_names<-merge(mod_genes,subset_pr_test_res,by.x="gene_id",by.y="id")
mod_genes_names<-mod_genes_names[,c(1,9,10,11,12)]
write.csv(mod_genes_names, file.path(DIR, "mod15.csv"))






