library("tidyr")
library("dplyr")
setwd("DEG/")
al<-read.delim("alldiff.txt", header = T)
al<-al[al$FDR < 0.05,]

al<-separate(al, R25C_Junc_Inclusive..Exclusive, 
             c("R25C_Junc_Inclusive","R25C_Junc_Exclusive"),
             sep = "::", remove = T)
al<-separate(al, R32C_Junc_Inclusive..Exclusive, 
             c("R32C_Junc_Inclusive","R32C_Junc_Exclusive"),
             sep = "::", remove = T)
al<-separate(al, R25C_Exp_Inclusive..Exclusive, 
             c("R25C_Exp_Inclusive","R25C_Exp_Exclusive"),
             sep = "::", remove = T)
al<-separate(al, R32C_Exp_Inclusive..Exclusive, 
             c("R32C_Exp_Inclusive","R32C_Exp_Exclusive"),
             sep = "::", remove = T)

hsplist<-read.delim("hspgene.txt", header = T)
al.hsp<-al[al$AccID %in% unique(hsplist$Genes),]

gene2name<-read.delim("D:/oyster/eggnog/gene2name.txt", header = TRUE)
names(gene2name)[1]="AccID"
DEG<-read.delim("DEG.allpairs.txt", header = TRUE)
DEG$AccID<-row.names(DEG)
DEG$logFC<-0-DEG$logFC

al.anno<-merge(al, DEG, by = "AccID", all.x = TRUE)

al.anno<-merge(al.anno, gene2name, by = "AccID", all.x = TRUE)

write.table(al.anno,
            "alternative_splicing_filter_anno.txt",
            sep = "\t", quote = F, row.names = F)

al.summary<-as.data.frame(xtabs(~SplicingType,al))
#pie chart
taxoncolor<-c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
              "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
              "#BC80BD", "#CCEBC5", "#FFED6F", "#D9D9D9",  "white")
piep<-ggplot(al.summary, aes(x="", y=Freq, fill=SplicingType))+
  geom_bar(stat = "identity", position = position_fill())+
  scale_fill_manual(values=taxoncolor)+
  geom_text(aes(label = Freq), position = position_fill(vjust = 0.5), color = "black")+
  labs(fill = "AS type")+
  coord_polar(theta = "y", start=0)+
  theme_void()

piep

f = "al.pie.pptx"
topptx(piep,f, width = 6, height = 5, units = "in")

##Volcano plot
a<-ggplot(data = al, aes(x = delta_PSI, y = -log10(FDR), 
                         colour = SplicingType))+
  geom_point(shape = 20, size = 3)+
  scale_colour_manual(values = taxoncolor)+
  scale_x_continuous(limits = c(-1,1), breaks = seq(-1,1,0.2))+
  scale_y_continuous(limits = c(0,12), breaks = seq(0,12,2))+
  labs(x = "delta_PSI", y = "-log10(FDR)", legend = "")+
  geom_hline(yintercept = -log10(0.05), color = "grey30", linetype = "dashed")+
  geom_vline(xintercept = 0.1, color = "grey30", linetype = "dashed")+
  geom_vline(xintercept = -0.1, color = "grey30", linetype = "dashed")+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 0),
        legend.position = "none",
        panel.grid.minor.x = element_line(linetype = "blank"),
        strip.text = element_text(size = 12))
a
ggsave("al.Volcano.png", width = 5, height = 5, units = "in", dpi = 300)
ggsave("al.Volcano.tiff", width = 5, height = 5, units = "in", dpi = 300)


a<-ggplot(data = al.anno, aes(x = delta_PSI, y = logFC, 
                          colour = SplicingType))+
  geom_point(shape = 20, size = 3)+
  scale_colour_manual(values = taxoncolor)+
  scale_x_continuous(limits = c(-1,1), breaks = seq(-1,1,0.2))+
  scale_y_continuous(limits = c(-4,4), breaks = seq(-4,4,1))+
  labs(x = "delta_PSI", y = "-log2(Fold change)", legend = "")+
  geom_hline(yintercept = -1, color = "grey30", linetype = "dashed")+
  geom_hline(yintercept = 1, color = "grey30", linetype = "dashed")+
  geom_vline(xintercept = 0.1, color = "grey30", linetype = "dashed")+
  geom_vline(xintercept = -0.1, color = "grey30", linetype = "dashed")+
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12, colour = "black"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 0),
        legend.position = "none",
        panel.grid.minor.x = element_line(linetype = "blank"),
        strip.text = element_text(size = 12))
a
ggsave("al.DEG.Volcano.png", width = 5, height = 5, units = "in", dpi = 300)
ggsave("al.DEG.Volcano.tiff", width = 5, height = 5, units = "in", dpi = 300)

##KOG enrich
GenesKOGpair.1v1<-read.delim("D:/oyster/eggnog/KOG.1v1.txt", header = TRUE)

kog2name<-read.delim("D:/3enrichment/kog2name.txt", header = TRUE)
kog2name$ONTOLOGY<-gsub("INFORMATION STORAGE AND PROCESSING", "Information storage\nand processing", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("CELLULAR PROCESSES AND SIGNALING", "Cellular processes\nand signaling", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("POORLY CHARACTERIZED", "Poor", kog2name$ONTOLOGY)
kog2name$ONTOLOGY<-gsub("METABOLISM", "Metabolism", kog2name$ONTOLOGY)

al<-read.delim("alldiff.txt", header = T)
al<-al[al$FDR < 0.05,]

DEG<-al
names(DEG)[names(DEG)=="AccID"]<-"Genes"
DEG$Groups<-"Alternative.Splicing"
ref="al"
{
  KOG.all.1<-compareCluster(Genes ~ Groups, 
                            data = DEG, 
                            fun = 'enricher',
                            TERM2GENE = GenesKOGpair.1v1,
                            TERM2NAME = kog2name,
                            pvalueCutoff = 1,
                            pAdjustMethod = "BH",
                            qvalueCutoff = 1,
                            minGSSize = 1,
                            maxGSSize = 200000)
  
  ##Plot
  plotin<-as.data.frame(KOG.all.1)
  View(plotin)
  plotin<-subset(plotin, is.na(Description)==FALSE)
  
  plotinsep<-separate(plotin, "GeneRatio", c("Genenumerator", "Genedenominator"),sep = "/", remove = FALSE, convert = TRUE)
  plotinsep<-separate(plotinsep, "BgRatio", c("BGnumerator", "BGdenominator"),sep = "/", remove = FALSE, convert = TRUE)
  
  colnames(plotinsep)[names(plotinsep)=="ID"]<-"kogClass"
  
  plotdata<-merge(plotinsep, kog2name, by = c("kogClass"), all.x = TRUE)
  
  plotdata$ratio1=plotdata$Genenumerator/plotdata$BGnumerator
  plotdata$ratio2=plotdata$Genenumerator/plotdata$Genedenominator
  plotdata$Group=paste(plotdata$Cluster, "(", plotdata$Genedenominator, ")")
  
  write.table(plotdata, file = paste0(ref,"KOGenrich.txt"),
              sep = '\t', row.names = FALSE,
              quote = FALSE)
  
  list(plotdata$Group)
  
  a<-plotdata %>%
    #mutate(Groups = factor(Groups, levels = group.order.sub)) %>% 
    ggplot(aes(y = paste0(Description,"(",Count,")"), x = Count, 
               fill = ONTOLOGY))+
    labs(title = "", 
         fill = "KOG ONTOLOGY",
         x = "", 
         y = "")+
    geom_bar(stat = "identity", position = "dodge")+
    #geom_vline(xintercept = 0, color = "black", linetype = "solid")+
    guides(size = guide_legend(order = 1))+
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
          axis.title.x = element_text(colour = "black", size = 12), 
          axis.title.y = element_text(colour = "white"), 
          axis.text = element_text(colour = "black", size = 12), 
          axis.text.x = element_text(size = 12
                                     #, angle = 90, hjust = 1, vjust = 0.5
                                     ), 
          axis.text.y = element_text(size = 12),
          legend.title = element_text(size = 12), 
          plot.title = element_text(size = 12),
          legend.position = "none")+
    facet_grid(ONTOLOGY~., scales = "free", space = "free")+
    theme(strip.text = element_text(size = 8))
  a
  ggsave(paste0(ref,".KOG.tiff"), width = 10, height = 6, units = "in", dpi = 300)
  ggsave(paste0(ref,".KOG.png"), width = 10, height = 6, units = "in", dpi = 300)
  
  f = paste0(ref,".KOG.pptx")
  topptx(a,f, width = 10, height = 6, units = "in")
}