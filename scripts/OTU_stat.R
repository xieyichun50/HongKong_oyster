library(tidyr)
library(dplyr)
library(ggplot2)
library(eoffice)


taxonomy<-read.delim("taxonomy.tsv", header = TRUE)
names(taxonomy)<-c("Feature.ID", "index", "Confidence")
taxonomy$index<-gsub(" ","",taxonomy$index)

OTU<-read.delim("OTU_L7.txt", header = FALSE)
OTU<-as.data.frame(t(OTU))
names(OTU)<-OTU[1,]
OTU<-OTU[-1,]
row.names(OTU)<-1:nrow(OTU)
OTU<-OTU[1:(nrow(OTU)-2),]

OTU<-separate(OTU, index, into = c("L1_Kingdom",
                                   "L2_Phylum",
                                   "L3_Class",
                                   "L4_Order",
                                   "L5_Family",
                                   "L6_Genus",
                                   "L7_Species"),
              sep = ";")
OTU$F1<-as.numeric(OTU$F1)
sum(OTU$F1)
OTU$F1[OTU$F1<10]<-0
sum(OTU$F1)

OTU$F2<-as.numeric(OTU$F2)
sum(OTU$F2)
OTU$F2[OTU$F2<10]<-0
sum(OTU$F2)

OTU$F3<-as.numeric(OTU$F3)
sum(OTU$F3)
OTU$F3[OTU$F3<10]<-0
sum(OTU$F3)

OTU.filter<-OTU[(OTU$F1+OTU$F2+OTU$F3)>0,]

OTU.filter$ratio.F1<-OTU.filter$F1/sum(OTU$F1)
OTU.filter$ratio.F2<-OTU.filter$F1/sum(OTU$F2)
OTU.filter$ratio.F3<-OTU.filter$F1/sum(OTU$F3)

OTU.genus<-OTU.filter %>% group_by(L1_Kingdom,
                                   L2_Phylum,
                                   L3_Class,
                                   L4_Order,
                                   L5_Family,
                                   L6_Genus) %>% summarise_if(is.numeric, sum)
OTU.genus<-as.data.frame(OTU.genus)
nrow(OTU.filter[OTU.filter$F1>0,])
nrow(OTU.filter[OTU.filter$F1>0 & OTU.filter$L7_Species != "s__" & OTU.filter$L7_Species != "__",])
nrow(OTU.filter[OTU.filter$F1>0 & OTU.filter$L6_Genus != "g__" & OTU.filter$L6_Genus != "__",])
sum(OTU.filter$F1[OTU.filter$F1>0 & OTU.filter$L6_Genus != "g__" & OTU.filter$L6_Genus != "__"])/sum(OTU.filter$F1)


nrow(OTU.filter[OTU.filter$F2>0,])
nrow(OTU.filter[OTU.filter$F2>0 & OTU.filter$L7_Species != "s__" & OTU.filter$L7_Species != "__",])
nrow(OTU.filter[OTU.filter$F2>0 & OTU.filter$L6_Genus != "g__" & OTU.filter$L6_Genus != "__",])
sum(OTU.filter$F2[OTU.filter$F2>0 & OTU.filter$L6_Genus != "g__" & OTU.filter$L6_Genus != "__"])/sum(OTU.filter$F2)

nrow(OTU.filter[OTU.filter$F3>0,])
nrow(OTU.filter[OTU.filter$F3>0 & OTU.filter$L7_Species != "s__" & OTU.filter$L7_Species != "__",])
nrow(OTU.filter[OTU.filter$F3>0 & OTU.filter$L6_Genus != "g__" & OTU.filter$L6_Genus != "__",])
sum(OTU.filter$F3[OTU.filter$F3>0 & OTU.filter$L6_Genus != "g__" & OTU.filter$L6_Genus != "__"])/sum(OTU.filter$F3)


##Venn diagram
library("VennDiagram")
venn.diagram(x = list(`Individual 1` = row.names(OTU.filter[OTU.filter$F1>0,]), 
                      `Individual 2` = row.names(OTU.filter[OTU.filter$F2>0,]),
                      `Individual 3` = row.names(OTU.filter[OTU.filter$F3>0,])),
             filename = "OTU.shared.L7.svg", imagetype = "svg", 
             cat.dist = 0.1, margin = 0.05,
             units = "in", height = 3, width = 3, resolution = 300,
             cat.fontface = 1,
             fill = c("burlywood2", "darkgray", "plum2"))

venn.diagram(x = list(`Individual 1` = row.names(OTU.genus[OTU.genus$F1>0,]), 
                      `Individual 2` = row.names(OTU.genus[OTU.genus$F2>0,]),
                      `Individual 3` = row.names(OTU.genus[OTU.genus$F3>0,])),
             filename = "OTU.shared.L6.png", imagetype = "png", 
             cat.dist = 0.1, margin = 0.05,
             units = "in", height = 3, width = 3, resolution = 300,
             fill = c("burlywood2", "darkgray", "plum2"))

##Phylum distribution
OTU.L2<-read.delim("OTU_L2.txt", header = FALSE)
OTU.L2<-as.data.frame(t(OTU.L2))
names(OTU.L2)<-OTU.L2[1,]
OTU.L2<-OTU.L2[-1,]
row.names(OTU.L2)<-1:nrow(OTU.L2)
OTU.L2<-OTU.L2[1:(nrow(OTU.L2)-2),]

OTU.L2<-separate(OTU.L2, index, into = c("L1_Kingdom","L2_Phylum"), sep = ";")

OTU.L2$F1<-as.numeric(OTU.L2$F1)
sum(OTU.L2$F1)
OTU.L2$F1[OTU.L2$F1<10]<-0
sum(OTU.L2$F1)

OTU.L2$F2<-as.numeric(OTU.L2$F2)
sum(OTU.L2$F2)
OTU.L2$F2[OTU.L2$F2<10]<-0
sum(OTU.L2$F2)

OTU.L2$F3<-as.numeric(OTU.L2$F3)
sum(OTU.L2$F3)
OTU.L2$F3[OTU.L2$F3<10]<-0
sum(OTU.L2$F3)

OTU.L2<-OTU.L2[(OTU.L2$F1+OTU.L2$F2+OTU.L2$F3)>0,]

OTU.L2$ratio.F1<-OTU.L2$F1/sum(OTU.L2$F1)
OTU.L2$ratio.F2<-OTU.L2$F2/sum(OTU.L2$F2)
OTU.L2$ratio.F3<-OTU.L2$F3/sum(OTU.L2$F3)

OTU.L2<-OTU.L2[,c("L2_Phylum","ratio.F1","ratio.F2","ratio.F3")]
OTU.L2.long<-as.data.frame(matrix(data = NA, nrow = 0, ncol = 3))
names(OTU.L2.long)<-c("taxon","ratio","sample")
for (i in 2:ncol(OTU.L2)) {
  OTU.L2.sub<-OTU.L2[,c(1,i)]
  names(OTU.L2.sub)<-c("taxon","ratio")
  
  OTU.L2.sub.minor<-OTU.L2.sub[OTU.L2.sub$ratio<=0.01 | OTU.L2.sub$taxon == "Bacteria[Unassigned]",]
  OTU.L2.sub<-OTU.L2.sub[-which(OTU.L2.sub$ratio<=0.01 | OTU.L2.sub$taxon == "Bacteria[Unassigned]"),]
  OTU.L2.sub<-rbind(OTU.L2.sub,c("Minor OTUs (<1%)", sum(OTU.L2.sub.minor$ratio)))
  OTU.L2.sub$sample<-names(OTU.L2)[i]
  OTU.L2.long<-rbind(OTU.L2.long, OTU.L2.sub)
}

OTU.L2.long$taxon<-gsub("p__","",OTU.L2.long$taxon)
OTU.L2.long$ratio<-as.numeric(OTU.L2.long$ratio)

taxonorder<-unique(OTU.L2.long[order(OTU.L2.long$ratio, 
                                     decreasing = T),"taxon"])
taxonorder<-taxonorder[-which(taxonorder == "Minor OTUs (<1%)")]
taxonorder<-c(taxonorder,"Minor OTUs (<1%)")
library(RColorBrewer)
taxoncolor<-c("Red", brewer.pal(12, "Set3"), "white")
taxoncolor<-c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072",
              "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
              "#BC80BD", "#CCEBC5", "#FFED6F", "#D9D9D9",  "white")
a<-ggplot(data = OTU.L2.long, aes(x = sample, y = ratio, fill = factor(taxon, levels = taxonorder)))+
  geom_bar(stat = "identity", position = "stack")+
  geom_text(stat = "identity", aes(label= scales::percent(ratio, accuracy = 0.01)), 
            position = position_stack(0.5), size = 3.5)+
  scale_fill_manual(values = taxoncolor)+
  scale_y_continuous(labels = scales::percent, breaks = seq(0,1,0.1))+
  labs(x="",y="", title = "", fill = "Phylum")+
  theme(panel.background = element_rect(fill = NA),
        axis.line = element_line(linetype = "solid", size = 0.5),
        axis.ticks = element_line(colour = "black", size = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10, colour = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 0),
        legend.position = "right")

a
f="OTU.phylum.pptx"
topptx(a,f, width = 4, height = 6, units = "in")

##Family distribution
OTU.L2<-read.delim("OTU_L5.txt", header = FALSE)
OTU.L2<-as.data.frame(t(OTU.L2))

names(OTU.L2)<-OTU.L2[1,]
OTU.L2<-OTU.L2[-1,]

row.names(OTU.L2)<-1:nrow(OTU.L2)
OTU.L2<-OTU.L2[1:(nrow(OTU.L2)-2),]

OTU.L2<-separate(OTU.L2, index, into = c("L1_Kingdom",
                                         "L2_Phylum",
                                         "L3_Class",
                                         "L4_Order",
                                         "L5_Family"), sep = ";")

OTU.L2$F1<-as.numeric(OTU.L2$F1)
sum(OTU.L2$F1)
OTU.L2$F1[OTU.L2$F1<10]<-0
sum(OTU.L2$F1)

OTU.L2$F2<-as.numeric(OTU.L2$F2)
sum(OTU.L2$F2)
OTU.L2$F2[OTU.L2$F2<10]<-0
sum(OTU.L2$F2)

OTU.L2$F3<-as.numeric(OTU.L2$F3)
sum(OTU.L2$F3)
OTU.L2$F3[OTU.L2$F3<10]<-0
sum(OTU.L2$F3)

OTU.L2<-OTU.L2[(OTU.L2$F1+OTU.L2$F2+OTU.L2$F3)>0,]

OTU.L2$ratio.F1<-OTU.L2$F1/sum(OTU.L2$F1)
OTU.L2$ratio.F2<-OTU.L2$F2/sum(OTU.L2$F2)
OTU.L2$ratio.F3<-OTU.L2$F3/sum(OTU.L2$F3)

#OTU.L2<-OTU.L2[,c("L6_Genus","ratio.F1","ratio.F2","ratio.F3")]
OTU.L2.long<-as.data.frame(matrix(data = NA, nrow = 0, ncol = 7))
names(OTU.L2.long)<-c("L1_Kingdom",
                      "L2_Phylum",
                      "L3_Class",
                      "L4_Order",
                      "L5_Family",
                      "ratio","sample")
for (i in 9:ncol(OTU.L2)) {
  OTU.L2.sub<-OTU.L2[,c(1:5,i)]
  names(OTU.L2.sub)[6]<-"ratio"
  
  OTU.L2.sub<-OTU.L2.sub[OTU.L2.sub$ratio>=0.01 & OTU.L2.sub$L4_Order != "o__" & OTU.L2.sub$L4_Order != "__",]

  OTU.L2.sub<-rbind(OTU.L2.sub,c(rep("Others",5), 1-sum(OTU.L2.sub$ratio)))
  OTU.L2.sub$sample<-names(OTU.L2)[i]
  OTU.L2.long<-rbind(OTU.L2.long, OTU.L2.sub)
}

OTU.L2.long$ratio<-as.numeric(OTU.L2.long$ratio)
OTU.L2.long$L5_Family[OTU.L2.long$L5_Family=="f__"]<-paste0("Unclassified_",OTU.L2.long$L4_Order[OTU.L2.long$L5_Family=="f__"])
OTU.L2.long$L5_Family<-gsub("\\[","",OTU.L2.long$L5_Family)
OTU.L2.long$L5_Family<-gsub("\\]","",OTU.L2.long$L5_Family)

taxonorder<-unique(OTU.L2.long[order(OTU.L2.long$L5_Family, 
                                     decreasing = F),"L5_Family"])
taxonorder<-taxonorder[-which(taxonorder == "Others")]
taxonorder<-c(taxonorder,"Others")
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set3"))

a<-ggplot(data = OTU.L2.long, aes(x = sample, y = ratio, fill = factor(L5_Family, levels = taxonorder)))+
  geom_bar(stat = "identity", position = "stack")+
  geom_text(stat = "identity", aes(label= scales::percent(ratio, accuracy = 0.01)), 
            position = position_stack(0.5), size = 3.5)+
  scale_fill_manual(values = getPalette(length(taxonorder)))+
  scale_y_continuous(labels = scales::percent, breaks = seq(0,1,0.1))+
  labs(x="",y="", title = "", fill = "Phylum")+
  theme(panel.background = element_rect(fill = NA),
        axis.line = element_line(linetype = "solid", size = 0.5),
        axis.ticks = element_line(colour = "black", size = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10, colour = "black"),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 0),
        legend.position = "right")

a
f="OTU.family.pptx"
topptx(a,f, width = 5, height = 6, units = "in")

##Genus heatmap
OTU.L2<-read.delim("OTU_L6.txt", header = T)
row.names(OTU.L2)<-1:nrow(OTU.L2)
OTU.L2<-OTU.L2[1:(nrow(OTU.L2)-2),]
OTU.L2<-separate(OTU.L2, index, into = c("L1_Kingdom",
                                      "L2_Phylum",
                                      "L3_Class",
                                      "L4_Order",
                                      "L5_Family",
                                      "L6_Genus"), sep = ";")


OTU.L2$F1<-as.numeric(OTU.L2$F1)
sum(OTU.L2$F1)
OTU.L2$F1[OTU.L2$F1<10]<-0
sum(OTU.L2$F1)

OTU.L2$F2<-as.numeric(OTU.L2$F2)
sum(OTU.L2$F2)
OTU.L2$F2[OTU.L2$F2<10]<-0
sum(OTU.L2$F2)

OTU.L2$F3<-as.numeric(OTU.L2$F3)
sum(OTU.L2$F3)
OTU.L2$F3[OTU.L2$F3<10]<-0
sum(OTU.L2$F3)

OTU.L2<-OTU.L2[(OTU.L2$F1+OTU.L2$F2+OTU.L2$F3)>0,]

OTU.L2$ratio.F1<-OTU.L2$F1/sum(OTU.L2$F1)
OTU.L2$ratio.F2<-OTU.L2$F2/sum(OTU.L2$F2)
OTU.L2$ratio.F3<-OTU.L2$F3/sum(OTU.L2$F3)

#OTU.L2<-OTU.L2[OTU.L2$L6_Genus!="__",]
#OTU.L2<-OTU.L2[OTU.L2$L6_Genus!="g__",]

OTU.L2.heat<-OTU.L2[OTU.L2$L5_Family != "__" & OTU.L2$L5_Family != "f__",]

OTU.L2.heat$L6_Genus[OTU.L2.heat$L6_Genus=="g__"]<-"__Unclassified"
OTU.L2.heat$L6_Genus[OTU.L2.heat$L6_Genus=="__"]<-"__Unclassified"

OTU.L2.heat.u<-OTU.L2.heat[OTU.L2.heat$L6_Genus=="__Unclassified",] %>% group_by(L1_Kingdom,
                                                                                L2_Phylum,
                                                                                L3_Class,
                                                                                L4_Order,
                                                                                L5_Family,
                                                                                L6_Genus) %>% summarise_if(is.numeric, sum)
OTU.L2.heat.u<-OTU.L2.heat.u[order(OTU.L2.heat.u$ratio.F1,
                                   OTU.L2.heat.u$ratio.F2,
                                   OTU.L2.heat.u$ratio.F3, 
                                   decreasing = T),]

OTU.L2.heat.c<-OTU.L2.heat[OTU.L2.heat$L6_Genus!="__Unclassified",]
OTU.L2.heat.c<-OTU.L2.heat.c[order(OTU.L2.heat.c$ratio.F1,
                                   OTU.L2.heat.c$ratio.F2,
                                   OTU.L2.heat.c$ratio.F3, 
                                   decreasing = T),]

OTU.L2.heat.c<-rbind(OTU.L2.heat.c,OTU.L2.heat.u)

OTU.L2.heat.c$L5_Family<-gsub("f__","",OTU.L2.heat.c$L5_Family)
OTU.L2.heat.c$L6_Genus<-gsub("g__","__",OTU.L2.heat.c$L6_Genus)
OTU.L2.heat.c$label<-paste0(OTU.L2.heat.c$L5_Family, "__",OTU.L2.heat.c$L6_Genus)

OTU.L2.heat.c$label<-gsub("\\[","",OTU.L2.heat.c$label)
OTU.L2.heat.c$label<-gsub("\\]","",OTU.L2.heat.c$label)

co<-0.01
OTU.L2.heat.p<-OTU.L2.heat.c[OTU.L2.heat.c$ratio.F1 > co | OTU.L2.heat.c$ratio.F2 > co | OTU.L2.heat.c$ratio.F3 > co,]

row.names(OTU.L2.heat.p)<-OTU.L2.heat.p$label
OTU.L2.heat.p<-OTU.L2.heat.p[,c("ratio.F1","ratio.F2","ratio.F3")]
names(OTU.L2.heat.p)<-c("Individual 1", "Individual 2", "Individual 3")

#Genus heatmap
library(pheatmap)
OTU.L2.heat.p<-log2(OTU.L2.heat.p)
bk.limit<-ceiling(max(abs(OTU.L2.heat.p)))
bk<-c(seq(-bk.limit,-3,by=0.1))
heatp<-pheatmap(OTU.L2.heat.p,
                color = c(colorRampPalette(colors = c("white","red"))(length(bk))),
                #legend_breaks=seq(-bk.limit,bk.limit,1), breaks=bk,
                legend = FALSE,
                scale = "none",
                cluster_rows = FALSE, cluster_cols = FALSE,
                treeheight_row = 0,
                border_color = NA,
                cellwidth = 16, cellheight = 12, fontsize = 10,
                angle_col = c("45"),
                width = 6, height = 6,
                show_colnames = TRUE,
                show_rownames = TRUE,
                filename = "OTU.genus.m.png")
