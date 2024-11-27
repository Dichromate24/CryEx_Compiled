library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(biomaRt)

# Take the collapsed differential transcript analysis produced from de_analysis of SpliCeAT. 
# We take the "normal" set which did not have any modifications during prior augment transcriptome process

normal_collapse <- read_csv("/Users/darren/Desktop/NORMAL_collapsed_differential_transcript_analysis.csv")

# filter for qval_lrt and qval_wald <= 0.05
normal_collapse_filtered <- normal_collapse[(normal_collapse$qval_lrt <= 0.05),]
normal_collapse_filtered <- normal_collapse_filtered[(normal_collapse_filtered$qval_wald <= 0.05),]

# get upregulated and down regulated events

norm_up <- normal_collapse_filtered[(normal_collapse_filtered$b_wald >= 1),]
norm_up <- norm_up[order(norm_up$qval_wald),]
write.csv(norm_up, "/Users/darren/Desktop/norm_up.csv", row.names = F)
norm_up_genes <- norm_up$ext_gene

norm_down <- normal_collapse_filtered[(normal_collapse_filtered$b_wald <= -1),]
norm_down <- norm_down[order(norm_down$qval_wald),]
norm_down_genes <- norm_down$ext_gene

# get the ensembl ids, and remove version number for GO

norm_up_ensmb <- norm_up$ens_gene
norm_up_ensmb <- gsub("\\..*", "", norm_up_ensmb)

norm_down_ensmb <- norm_down$ens_gene
norm_down_ensmb <- gsub("\\..*", "", norm_down_ensmb)

# Perform GO analysis for upregulated genes and down regulated genes

go_up <- enrichGO(gene         = norm_up_ensmb,
                  OrgDb        = org.Hs.eg.db,
                  keyType      = "ENSEMBL",
                  ont          = "BP",  
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05)

View(as_tibble(go_up))
barplot(go_up,showCategory=12)
dotplot(go_up)
cnetplot(go_up,showCategory=10,cex_label_gene=0.5)

go_down <- enrichGO(gene       = norm_down_ensmb,
                  OrgDb        = org.Hs.eg.db,
                  keyType      = "ENSEMBL",
                  ont          = "BP",  
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05)

View(as_tibble(go_down))
barplot(go_down,showCategory=12)
dotplot(go_down)
cnetplot(go_down,showCategory=10,cex_label_gene=0.5)


# Perform GO analysis for all upregulated novel groups

norm_up_novel <- read_csv("/Users/darren/Desktop/norm_up_novel_ens.csv", col_names = F)
norm_up_novel_ensmb <- norm_up_novel |>
  mutate(X1 = gsub("\\..*", "", X1))
norm_up_novel_ensmb <- norm_up_novel_ensmb$X1

go_novel <- enrichGO(gene       = norm_up_novel_ensmb,
                    OrgDb        = org.Hs.eg.db,
                    keyType      = "ENSEMBL",
                    ont          = "BP",  
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05)

View(as_tibble(go_novel))                                  # There were no significant biological processes



# Zooming into top 20 upregulated / down regulated
norm_up_top20 <- norm_up[c(1:20),]
norm_up_top20 <- norm_up_top20 |> 
                  dplyr::select(ens_gene, ext_gene, target_id, qval_lrt, qval_wald, b_wald)

norm_down_top20 <- norm_down[c(1:20),]
norm_down_top20 <- norm_down_top20 |> 
  dplyr::select(ens_gene, ext_gene, target_id, qval_lrt, qval_wald, b_wald)



# Redo everything with the "union" set, which used the union of all genes during the augment transcriptome process

union_collapse <- read_csv("/Users/darren/Desktop/UNION_collapsed_differential_transcript_analysis.csv")

union_collapse_filtered <- union_collapse[(union_collapse$qval_lrt <= 0.05),]
union_collapse_filtered <- union_collapse_filtered[(union_collapse_filtered$qval_wald <= 0.05),]

union_up <- union_collapse_filtered[(union_collapse_filtered$b_wald >= 1),]
union_up <- union_up[order(union_up$qval_wald),]
union_up_genes <- union_up$ext_gene
write.csv(union_up_genes, "/Users/darren/Desktop/union_up_genes.csv", row.names = F)

union_down <- union_collapse_filtered[(union_collapse_filtered$b_wald <= -1),]
union_down <- union_down[order(union_down$qval_wald),]
union_down_genes <- union_down$ext_gene
write.csv(union_down_genes, "/Users/darren/Desktop/union_down_genes.csv", row.names = F)

union_up_ensmb <- union_up$ens_gene
union_up_ensmb <- gsub("\\..*", "", union_up_ensmb)

union_down_ensmb <- union_down$ens_gene
union_down_ensmb <- gsub("\\..*", "", union_down_ensmb)

union_go_up <- enrichGO(gene    = union_up_ensmb,
                  OrgDb        = org.Hs.eg.db,
                  keyType      = "ENSEMBL",
                  ont          = "BP",  
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05)

View(as_tibble(union_go_up))
barplot(union_go_up,showCategory=12)
dotplot(union_go_up)
cnetplot(union_go_up,showCategory=8,cex_label_gene=0.5)

union_go_down <- enrichGO(gene         = union_down_ensmb,
                        OrgDb        = org.Hs.eg.db,
                        keyType      = "ENSEMBL",
                        ont          = "BP",  
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05)

View(as_tibble(union_go_down))
barplot(union_go_down,showCategory=12)
dotplot(union_go_down)
cnetplot(union_go_down,showCategory=10,cex_label_gene=0.5)




# Redo everything with the "subset" set, which used a subset of union set during the augment transcriptome process

subset_collapse <- read_csv("/Users/darren/Desktop/SUBSET_collapsed_differential_transcript_analysis.csv")

subset_collapse_filtered <- subset_collapse[(subset_collapse$qval_lrt <= 0.05),]
subset_collapse_filtered <- subset_collapse_filtered[(subset_collapse_filtered$qval_wald <= 0.05),]

subset_up <- subset_collapse_filtered[(subset_collapse_filtered$b_wald >= 1),]
subset_up <- subset_up[order(subset_up$qval_wald),]
subset_up_genes <- subset_up$ext_gene
write.csv(subset_up_genes, "/Users/darren/Desktop/subset_up_genes.csv", row.names = F)

subset_down <- subset_collapse_filtered[(subset_collapse_filtered$b_wald <= -1),]
subset_down <- subset_down[order(subset_down$qval_wald),]
subset_down_genes <- subset_down$ext_gene
write.csv(subset_down_genes, "/Users/darren/Desktop/subset_down_genes.csv", row.names = F)

subset_up_ensmb <- subset_up$ens_gene
subset_up_ensmb <- gsub("\\..*", "", subset_up_ensmb)

subset_down_ensmb <- subset_down$ens_gene
subset_down_ensmb <- gsub("\\..*", "", subset_down_ensmb)

subset_go_up <- enrichGO(gene    = subset_up_ensmb,
                        OrgDb        = org.Hs.eg.db,
                        keyType      = "ENSEMBL",
                        ont          = "BP",  
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05)

View(as_tibble(subset_go_up))
barplot(union_go_up,showCategory=12)
dotplot(union_go_up)
cnetplot(subset_go_up,showCategory=5,cex_label_gene=0.5)

subset_go_down <- enrichGO(gene         = subset_down_ensmb,
                          OrgDb        = org.Hs.eg.db,
                          keyType      = "ENSEMBL",
                          ont          = "BP",  
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05)

View(as_tibble(subset_go_down))
barplot(union_go_down,showCategory=12)
dotplot(union_go_down)
cnetplot(union_go_down,showCategory=10,cex_label_gene=0.5)




# Peform GO analysis on the gene set that was published by CryEx 

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

CryEx_dpsi_up_genes <- read_csv("/Users/darren/Desktop/CryEx_up_genes.csv")

CryEx_gene_vector <- CryEx_dpsi_up_genes[[1]]

CryEx_gene_mapping <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                      filters = 'external_gene_name',
                      values = CryEx_gene_vector,
                      mart = ensembl)


CryEx_provided <- enrichGO(gene         = CryEx_gene_mapping$ensembl_gene_id,
                       OrgDb        = org.Hs.eg.db,
                       keyType      = "ENSEMBL",
                       ont          = "BP",  
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05)

View(as_tibble(CryEx_provided))
barplot(CryEx_provided,showCategory=12)





