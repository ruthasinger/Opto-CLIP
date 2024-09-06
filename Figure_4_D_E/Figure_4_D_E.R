library(tidyverse)

##############################################################################################################
#setup the directories
Basedirectory=file.path("~","Opto-CLIP","Figure_4_D_E")
list.files(Basedirectory)

Datadirectory=file.path(Basedirectory,"Data")
list.files(Datadirectory)

Outdirectory=file.path(Basedirectory,"Graphs")
dir.create(Outdirectory,showWarnings = TRUE)

##############################################################################################################
#Bring in FMRP CLIP files
FMRP_CLIP_Tag_files=grep("tag.uniq.bed$",list.files(Datadirectory),value=TRUE)
FMRP_CLIP_Anno_files=grep("annot.txt$",list.files(Datadirectory),value=TRUE)
FMRP_CLIP_STAR_files=grep("transcriptome.unique.bed$",list.files(Datadirectory),value=TRUE)

FMRP_CLIP_Tag_files
FMRP_CLIP_Anno_files
FMRP_CLIP_STAR_files

tags=NULL
annos=NULL
tags_anno=NULL
tags_anno_use=NULL
transcripts=NULL
tags_anno_transcript=NULL
for (i in 1:(length(FMRP_CLIP_Tag_files))) {
  tags[[i]]=read_delim(file.path(Datadirectory,FMRP_CLIP_Tag_files[[i]]), delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
  colnames(tags[[i]])=c("chr","start","end","name","score","strand","thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts")
  names(tags)[[i]]=vapply(strsplit(FMRP_CLIP_Tag_files[[i]], ".", fixed = TRUE), "[", "", 1)
  annos[[i]]=read_delim(file.path(Datadirectory,FMRP_CLIP_Anno_files[[i]]), delim = "\t", escape_double = FALSE, col_names = TRUE, trim_ws = TRUE)
  colnames(annos[[i]])=c("name","gene_id","gene_name","region")
  names(annos)[[i]]=vapply(strsplit(FMRP_CLIP_Tag_files[[i]], ".", fixed = TRUE), "[", "", 1)
  tags_anno[[i]]=merge(tags[[i]],annos[[i]],by.x="name",by.y="name")
  names(tags_anno)[[i]]=vapply(strsplit(FMRP_CLIP_Tag_files[[i]], ".", fixed = TRUE), "[", "", 1)
  tags_anno_use[[i]]=tags_anno[[i]][, c("name", "gene_name","region")]
  names(tags_anno_use)[[i]]=vapply(strsplit(FMRP_CLIP_Tag_files[[i]], ".", fixed = TRUE), "[", "", 1)
  transcripts[[i]]=read_delim(file.path(Datadirectory,FMRP_CLIP_STAR_files[[i]]), delim = "\t", col_names = FALSE)
  colnames(transcripts[[i]])= c("Transcript_id", "Txstart", "Txend", "read_ID", "Txscore", "Txstrand")
  names(transcripts)[[i]]=vapply(strsplit(FMRP_CLIP_Tag_files[[i]], ".", fixed = TRUE), "[", "", 1)
  tags_anno_transcript[[i]]=merge(tags_anno_use[[i]],transcripts[[i]],by.x="name",by.y="read_ID")
  names(tags_anno_transcript)[[i]]=vapply(strsplit(FMRP_CLIP_Tag_files[[i]], ".", fixed = TRUE), "[", "", 1)
}

########################################################################################################################
#Figure 4D (actual graph made in prism)
No.UMTs <- sapply(tags, function(x) nrow(x))

No.UMTs.df=data.frame(Sample=names(No.UMTs),
                      UMTs=No.UMTs)

#export datafrane and make graph in Prism
write.csv(No.UMTs.df, file = paste0(Outdirectory, "/Figure4D_No.UMTs.df.csv"), row.names=FALSE)

########################################################################################################################
#Figure 4E
All.FMRP.tags=do.call("rbind", tags_anno_transcript)

All.FMRP.tags$sourcefile=rownames(All.FMRP.tags)
row.names(All.FMRP.tags) <- NULL

All.FMRP.tags= All.FMRP.tags %>% separate(sourcefile, into = c("Sample","two"), sep = "\\.", remove = TRUE) 
All.FMRP.tags= All.FMRP.tags %>% dplyr::select(!two)

All.FMRP.tags.Crepos=All.FMRP.tags %>% dplyr::filter(Sample %in% c("FMRP_ChR2_rep1", "FMRP_ChR2_rep2", "FMRP_ChR2_rep3", "FMRP_ChR2_rep4",
                                                                   "FMRP_Control_rep1","FMRP_Control_rep2","FMRP_Control_rep3","FMRP_Control_rep4"))

counts <- All.FMRP.tags.Crepos %>%
  group_by(Transcript_id, region, Sample) %>%
  summarise(tags = dplyr::n()) %>% ungroup()

counts2=counts %>% 
  separate(Sample, into = c("RBP","Treatment","rep"), sep = "\\_", remove = FALSE) %>%
  group_by(Transcript_id, region, Treatment) %>%                                                            
  summarise(BC = n())  %>%
  pivot_wider(names_from = Treatment, values_from = BC, names_prefix="BC_") %>% 
  replace(is.na(.), 0)

counts_wide=counts %>% 
  pivot_wider(names_from= Sample, values_from= tags) %>%
  replace(is.na(.), 0)

counts_BC <- cbind(counts_wide, counts2[,3:4])

BCs=colnames(counts_BC)[grepl("BC",colnames(counts_BC))]

toPlot_longer=counts_BC %>% 
  pivot_longer(names_to = "Sample", values_to = "BC", cols = BCs) 

toPlot_longer_use= toPlot_longer %>% 
  filter(region %in% c("3'UTR", "5'UTR", "CDS", "deep_intergenic", "intron")) %>%
  filter(BC > 0) %>% 
  mutate(BC= as.factor(BC))

CLIP_pipeline_annotations_withBC=ggplot(toPlot_longer_use, aes(x = BC, fill = region)) +
  geom_bar(stat = "count", position = "fill") +
  xlab("Number of Biological Replicates") + 
  facet_wrap(~Sample) +
  theme_classic() +
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07","#000066","#33CC99"))
CLIP_pipeline_annotations_withBC
ggsave(paste0(Outdirectory,"/","Figure_4E.pdf"),plot=CLIP_pipeline_annotations_withBC,device = "pdf")
