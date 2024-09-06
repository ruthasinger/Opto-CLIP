library(rtracklayer)
library(tximport)
library(tidyverse)
library(AnnotationDb)
library(org.Mm.eg.db)

##############################################################################################################
#setup the directories
Basedirectory=file.path("~","Opto-CLIP","Figure_5_E")
list.files(Basedirectory)

Datadirectory=file.path(Basedirectory,"Data")
list.files(Datadirectory)

Outdirectory=file.path(Basedirectory,"Graphs")
dir.create(Outdirectory,showWarnings = TRUE)

##############################################################################################################
#Bring in GTF for transcript information
mm10_gtf=rtracklayer::import(file.path(Datadirectory,"mm10.ensGene.gtf.gz"))
mm10_gtf_df=mm10_gtf %>% as_tibble
mm10_Tx=mm10_gtf_df %>% dplyr::filter(type=="transcript")
mm10_tx2geneid=mm10_Tx %>% dplyr::select(transcript_id,gene_id)

################################################################################
#Bring in Opto-RiboTag data
RiboTag_Metadata <- read.delim(file.path(Datadirectory,"RiboTag_Metadata_IPandInput_withnoCre.txt"))
RiboTag_files <- file.path(Datadirectory,RiboTag_Metadata$Filename, "quant.sf")
names(RiboTag_files) <- RiboTag_Metadata$Filename
RiboTag_files
all(file.exists(RiboTag_files))

RiboTag_txi <- tximport(RiboTag_files, type = "salmon", tx2gene = mm10_tx2geneid)

RiboTag_TPM=as.data.frame(RiboTag_txi$abundance)
RiboTag_TPM <- tibble::rownames_to_column(RiboTag_TPM, "gene_id")

colnames(RiboTag_TPM)

RiboTag_TPM=RiboTag_TPM %>%
  dplyr::mutate(Control_rep1_log2FC  = log2(RiboTag_Control_30min_IP_rep1_salmon/RiboTag_Control_30min_Input_rep1_salmon),
                Control_rep2_log2FC  = log2(RiboTag_Control_30min_IP_rep2_salmon/RiboTag_Control_30min_Input_rep2_salmon),
                ChR2_rep1_log2FC     = log2(RiboTag_ChR2_30min_IP_rep1_salmon/RiboTag_ChR2_30min_Input_rep1_salmon),
                ChR2_rep2_log2FC     = log2(RiboTag_ChR2_30min_IP_rep2_salmon/RiboTag_ChR2_30min_Input_rep2_salmon),
                noCre_rep1_log2FC    = log2(RiboTag_NoCre_30min_IP_rep1_salmon/RiboTag_NoCre_30min_input_rep1_salmon))

colnames(RiboTag_TPM)

RiboTag_TPM$gene_name<- mapIds(org.Mm.eg.db,
                               keys = RiboTag_TPM$gene_id,
                               column = "SYMBOL",
                               keytype = "ENSEMBL",
                               multiVals = "first")

GoI=c("Neurod6",
      "Camk2a",
      "Rbfox3",
      "Snap25",
      "Nrgn",
      "Hpca",
      "Crym",
      "Chn1",
      "Gad2",
      "Sst",
      "Calb2",
      "Pdgfra",
      "Ptprz1",
      "Mag",
      "Mal",
      "Mbp",
      "Mobp",
      "Plp1",
      "Gfap",
      "Glul",
      "Aqp4",
      "Aldh1l1",
      "Pla2g7",
      "Slc1a3",
      "Aldoc")

RiboTag_TPM_GOI=RiboTag_TPM[RiboTag_TPM$gene_name %in% GoI,]

RiboTag_TPM_GOI_log2FC=RiboTag_TPM_GOI %>% dplyr::select(c(gene_name,
                                                           Control_rep1_log2FC,
                                                           Control_rep2_log2FC,
                                                           ChR2_rep1_log2FC,
                                                           ChR2_rep2_log2FC,
                                                           noCre_rep1_log2FC))
RiboTag_TPM_GOI_long <- RiboTag_TPM_GOI_log2FC %>%
  pivot_longer(names_to = "Sample", values_to = "log2FC_TPM", cols = c(2:6)) 

RT_IPvsInput_plot=ggplot(RiboTag_TPM_GOI_long, aes(Sample, gene_name, fill= log2FC_TPM)) + 
  geom_tile() +
  scale_x_discrete(limits=c("noCre_rep1_log2FC","Control_rep1_log2FC","Control_rep2_log2FC","ChR2_rep1_log2FC", "ChR2_rep2_log2FC")) +
  scale_y_discrete(limits=rev(GoI)) +
  scale_fill_viridis_c(direction= 1) +
  coord_fixed(ratio=1/3) +
  theme_minimal()
RT_IPvsInput_plot
ggsave(paste0(Outdirectory, "/Figure5E.pdf"), plot = RT_IPvsInput_plot, device = "pdf")

