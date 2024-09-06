library(rtracklayer)
library(tidyverse)
library(AnnotationDb)
library(org.Mm.eg.db)
library(tximport)
library(DESeq2)
library(clusterProfiler)

##############################################################################################################
#setup the directories
Basedirectory=file.path("~","Opto-CLIP","Figure_5_F_G")
list.files(Basedirectory)

Datadirectory=file.path(Basedirectory,"Data")
list.files(Datadirectory)

Outdirectory=file.path(Basedirectory,"Graphs")
dir.create(Outdirectory,showWarnings = TRUE)

##############################################################################################################
#reading in GTF for transcript information
mm10_gtf=rtracklayer::import(file.path(Datadirectory,"mm10.ensGene.gtf.gz"))
mm10_gtf_df=mm10_gtf %>% as_tibble
mm10_Tx=mm10_gtf_df %>% dplyr::filter(type=="transcript")
mm10_Tx_final=mm10_Tx %>% dplyr::select(seqnames,start,end,width,strand,gene_id,transcript_id)
mm10_Tx_final$gene_name<- mapIds(org.Mm.eg.db,
                                 keys = mm10_Tx_final$gene_id,
                                 column = "SYMBOL",
                                 keytype = "ENSEMBL",
                                 multiVals = "first")
mm10_Tx_final$entrez<- mapIds(org.Mm.eg.db,
                              keys = mm10_Tx_final$gene_id,
                              column = "ENTREZID",
                              keytype = "ENSEMBL",
                              multiVals = "first")

mm10_Tx_final_noNAs=subset(mm10_Tx_final, is.na(gene_name) == FALSE)

################################################################################
#Bring in Opto-RiboTag data
RiboTag_Metadata <- read.delim(file.path(Datadirectory,"RiboTag_Metadata_justIP.txt"))
RiboTag_files <- file.path(Datadirectory,RiboTag_Metadata$Filename, "quant.sf")
names(RiboTag_files) <- RiboTag_Metadata$Filename
all(file.exists(RiboTag_files))

RiboTag_txi.tx <- tximport(RiboTag_files, type = "salmon", txOut = TRUE)
rownames(RiboTag_Metadata) <- RiboTag_Metadata$Filename

RiboTag_Metadata <- mutate_at(RiboTag_Metadata, vars(Group), as.factor)

RT_dds <- DESeqDataSetFromTximport(RiboTag_txi.tx, RiboTag_Metadata, ~ Group)

keep <- rowSums(counts(RT_dds)) >= 1
RT_dds <- RT_dds[keep, ]
RT_dds <- DESeq(RT_dds)

res=results(RT_dds, contrast=c("Group", "ChR2","Control"))
res=res[order(res$padj), ]
resdata=as.data.frame(res)
resdata <- tibble::rownames_to_column(resdata, "Transcript_ID")

RiboTag_RSEM_TPM=as.data.frame(RiboTag_txi.tx$abundance)
RiboTag_RSEM_TPM <- tibble::rownames_to_column(RiboTag_RSEM_TPM, "trancript_ID")
colnames(RiboTag_RSEM_TPM) <- c("transcript_ID",paste(colnames(RiboTag_RSEM_TPM[2:ncol(RiboTag_RSEM_TPM)]), "TPM", sep = "."))

RiboTag_resdata_TPM=merge(resdata,RiboTag_RSEM_TPM, by.x="Transcript_ID", by.y="transcript_ID", all.x=TRUE, all.y=FALSE)

#add in gtf info
RiboTag_mm10_gtf_merged <- merge(RiboTag_resdata_TPM, mm10_Tx_final_noNAs, by.x="Transcript_ID", by.y="transcript_id")

#should be true
length(unique(RiboTag_mm10_gtf_merged$Transcript_ID)) == nrow(RiboTag_mm10_gtf_merged)

#should be false
length(unique(RiboTag_mm10_gtf_merged$gene_name)) == nrow(RiboTag_mm10_gtf_merged)

nrow(RiboTag_mm10_gtf_merged) #57064

RiboTag_mm10_gtf_merged_exp=RiboTag_mm10_gtf_merged %>% 
  dplyr::filter(RiboTag_ChR2_30min_IP_rep1_salmon.TPM > 0 & RiboTag_ChR2_30min_IP_rep2_salmon.TPM > 0 & RiboTag_Control_30min_IP_rep1_salmon.TPM > 0 & RiboTag_Control_30min_IP_rep2_salmon.TPM > 0)

nrow(RiboTag_mm10_gtf_merged_exp) #31892

RiboTag_mm10_gtf_merged_exp= RiboTag_mm10_gtf_merged_exp %>%
  rowwise() %>%
  mutate(avg_rest_RT_TPM    = mean(c_across(c("RiboTag_Control_30min_IP_rep1_salmon.TPM", "RiboTag_Control_30min_IP_rep2_salmon.TPM")), na.rm=TRUE),
         avg_opto_RT_TPM    = mean(c_across(c("RiboTag_ChR2_30min_IP_rep1_salmon.TPM", "RiboTag_ChR2_30min_IP_rep2_salmon.TPM")), na.rm=TRUE),
         avg_both_RT_TPM    = mean(c_across(c("RiboTag_Control_30min_IP_rep1_salmon.TPM","RiboTag_Control_30min_IP_rep2_salmon.TPM","RiboTag_ChR2_30min_IP_rep1_salmon.TPM", "RiboTag_ChR2_30min_IP_rep2_salmon.TPM")), na.rm=TRUE),
         avg_RT_TPM_FC      = avg_opto_RT_TPM/avg_rest_RT_TPM,
         avg_RT_TPM_FC_log2 = log2(avg_RT_TPM_FC))

#Figure 5F
RiboTag_mm10_gtf_merged_exp=RiboTag_mm10_gtf_merged_exp[order(RiboTag_mm10_gtf_merged_exp$padj), ]
voldf=data.frame(gene = RiboTag_mm10_gtf_merged_exp$gene_name,
                 padj = -log10(RiboTag_mm10_gtf_merged_exp$padj), 
                 lfc = RiboTag_mm10_gtf_merged_exp$log2FoldChange)
voldf = na.omit(voldf)
voldf = mutate(voldf, color = case_when(voldf$lfc > 0 & voldf$padj > 1.3 ~ "Increased",
                                        voldf$lfc < 0 & voldf$padj > 1.3 ~ "Decreased",
                                        voldf$padj < 1.3 ~ "NS"))
vol_plot=ggplot(voldf, aes(x = lfc, y = padj, color = color)) + 
  geom_point(size = 2, na.rm = T) + 
  xlim(-6,6) +
  scale_color_manual(name = "Directionality", values = c(Increased = "red3", Decreased = "dodgerblue", NS = "darkgray")) +
  scale_y_continuous(trans = "log1p") + 
  theme(axis.text=element_text(size=12), legend.position = "right") +  
  xlab(paste("Log2 Fold Change\n","Opto vs Rest", sep=" ")) + 
  ylab(expression(-log[10]("padj"))) + 
  geom_hline(yintercept = 1.3, linetype='dotted', colour = "darkgrey") +
  annotate("text", 
           x = 3, 
           y = 30, 
           label = paste0("Decreased (", nrow(voldf %>% filter(color == "Decreased")),")"),
           color = "dodgerblue") +
  annotate("text", 
           x = 3, 
           y = 25, 
           label = paste0("Increased (", nrow(voldf %>% filter(color == "Increased")),")"),
           color = "red3") +
  theme_classic()
ggsave(filename = paste0(Outdirectory,"/Figure5F.pdf"), plot =vol_plot)

################################################################################################################################
#Figure 5G- GSEA
RT_forGO_Up=subset(RiboTag_mm10_gtf_merged_exp, padj < 0.05 & log2FoldChange > 0 & is.na(entrez) == FALSE)
nrow(RT_forGO_Up) #139

df=RT_forGO_Up
go_enrich_all= enrichGO(gene = df$entrez,
                        OrgDb = 'org.Mm.eg.db', 
                        keyType = "ENTREZID",
                        readable = T,
                        ont = "all",
                        pvalueCutoff = 0.05, 
                        qvalueCutoff = 0.10)
godata_RT_forGO_Up=go_enrich_all@result
rownames(godata_RT_forGO_Up)=NULL

RT_forGO_Down=subset(RiboTag_mm10_gtf_merged_exp, padj < 0.05 & log2FoldChange < 0 & is.na(entrez) == FALSE)
nrow(RT_forGO_Down) #200

df=RT_forGO_Down
go_enrich_all= enrichGO(gene = df$entrez,
                        OrgDb = 'org.Mm.eg.db', 
                        keyType = "ENTREZID",
                        readable = T,
                        ont = "all",
                        pvalueCutoff = 0.05, 
                        qvalueCutoff = 0.10)
godata_RT_forGO_Down=go_enrich_all@result
rownames(godata_RT_forGO_Down)=NULL

RT_GOterms_list=list(godata_RT_forGO_Up,godata_RT_forGO_Down)
names(RT_GOterms_list)=c("RT_forGO_Up","RT_forGO_Down")

Allgoterms_RiboTag=do.call("rbind", RT_GOterms_list)

Allgoterms_RiboTag$sourcefile=rownames(Allgoterms_RiboTag)
row.names(Allgoterms_RiboTag) <- NULL

Allgoterms_RiboTag= Allgoterms_RiboTag %>% separate(sourcefile, into = c("Sample","two"), sep = "\\.", remove = TRUE) 
Allgoterms_RiboTag= Allgoterms_RiboTag %>% dplyr::select(!two)

RT_terms=c("negative regulation of response to external stimulus",
           "regulation of viral process",
           "response to interferon-beta",
           "ubiquitin-like protein ligase binding",
           "nuclear envelope",
           "GTPase activity",
           "regulation of RNA splicing",
           "regulation of mRNA processing",
           "neuron to neuron synapse",
           "postsynaptic density",
           "positive regulation of cell projection organization",
           "protein serine/threonine kinase activity",
           "learning or memory",
           "cognition",
           "peptidyl-serine phosphorylation",       
           "calmodulin binding")

Allgoterms_RiboTag_use=subset(Allgoterms_RiboTag, Description %in% RT_terms)
Allgoterms_RiboTag_use=Allgoterms_RiboTag_use[!duplicated(Allgoterms_RiboTag_use$Description), ]

RT_GO_plot=ggplot(Allgoterms_RiboTag_use, aes(Count, Description, fill = Sample)) +
  geom_bar(stat = "identity") +
  scale_y_discrete(limits=RT_terms) +
  theme_classic() +
  scale_fill_manual(name = "RiboTag", 
                    values = c(RT_forGO_Up ="red3",
                               RT_forGO_Down = "dodgerblue"))
ggsave(paste0(Outdirectory, "/Figure5G.pdf"), plot = RT_GO_plot, device = "pdf")
