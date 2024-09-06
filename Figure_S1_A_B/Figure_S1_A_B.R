library(tidyverse)
library(DESeq2)
library(corrplot)

##############################################################################################################
#setup the directory for all the files to be written to
Basedirectory=file.path("~","Opto-CLIP","Figure_S1_A_B")
list.files(Basedirectory)

Datadirectory=file.path(Basedirectory,"Data")
list.files(Datadirectory)

Outdirectory=file.path(Basedirectory,"Graphs")
dir.create(Outdirectory,showWarnings = TRUE)

##############################################################################################################
#Bring in FMRP CLIP files
FMRP_CLIP_STAR_files=grep("transcriptome.unique.bed$",list.files(Datadirectory),value=TRUE)

FMRP_CLIP_STAR_files

transcripts=NULL
for (i in 1:(length(FMRP_CLIP_STAR_files))) {
  transcripts[[i]]=read_delim(file.path(Datadirectory,FMRP_CLIP_STAR_files[[i]]), delim = "\t", col_names = FALSE)
  colnames(transcripts[[i]])= c("Transcript_id", "Txstart", "Txend", "read_ID", "Txscore", "Txstrand")
  names(transcripts)[[i]]=vapply(strsplit(FMRP_CLIP_STAR_files[[i]], ".", fixed = TRUE), "[", "", 1)
 }

All.FMRP.Tx=do.call("rbind", transcripts)

All.FMRP.Tx$sourcefile=rownames(All.FMRP.Tx)
row.names(All.FMRP.Tx) <- NULL

All.FMRP.Tx= All.FMRP.Tx %>% separate(sourcefile, into = c("Sample","two"), sep = "\\.", remove = TRUE) 
All.FMRP.Tx= All.FMRP.Tx %>% dplyr::select(!two)

counts.Tx <- All.FMRP.Tx %>%
  group_by(Transcript_id, Sample) %>%
  summarise(tags = dplyr::n()) %>% ungroup()

counts.Tx_wide=counts.Tx %>% 
  pivot_wider(names_from= Sample, values_from= tags) %>%
  replace(is.na(.), 0)

########################################################################################################################
#Figure S1A- PCA plot
FMRP_count_matrix=counts.Tx_wide

FMRP_count_matrix=column_to_rownames(FMRP_count_matrix, var = "Transcript_id")

FMRP_sampleTable <- data.frame(Sample= colnames(FMRP_count_matrix))

FMRP_sampleTable= FMRP_sampleTable %>% 
  separate(Sample, into = c("RBP","Treatment","Replicate"), sep = "_", remove = FALSE) 

FMRP_sampleTable=column_to_rownames(FMRP_sampleTable, var = "Sample")

all(rownames(FMRP_sampleTable) == colnames(FMRP_count_matrix))

FMRP_dds <- DESeqDataSetFromMatrix(countData = FMRP_count_matrix,
                                   colData = FMRP_sampleTable,
                                   design = ~ Treatment)

keep <- rowSums(counts(FMRP_dds)) >= 1
FMRP_dds <- FMRP_dds[keep,]
FMRP_dds <- DESeq(FMRP_dds)

rld <- rlog(FMRP_dds, blind=FALSE)

pcaData <- plotPCA(rld, intgroup="Treatment", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

PCA=ggplot(pcaData, aes(PC1, PC2, color=Treatment)) +
  geom_point(size=4) +
  scale_color_manual(values=c("red3","dodgerblue")) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_classic()
PCA
ggsave(paste0(Outdirectory,"/Figure_S1A.pdf"),plot= PCA,device="pdf")

########################################################################################################################
##Figure S1B- Correlation plot
my_corr_df=counts.Tx_wide %>% 
  ungroup() %>%
  select(FMRP_ChR2_rep1,
         FMRP_ChR2_rep2,
         FMRP_ChR2_rep3,
         FMRP_ChR2_rep4,
         FMRP_Control_rep1,
         FMRP_Control_rep2,
         FMRP_Control_rep3,
         FMRP_Control_rep4)

my_corr_matrix=cor(my_corr_df)

pdf(file = paste0(Outdirectory,"/Figure_S1B.pdf"))
corrplot.mixed(my_corr_matrix,tl.cex = 0.7,tl.col="black" , lower.col=COL2('PiYG', 100), upper.col=COL2('PiYG', 100))
dev.off()
